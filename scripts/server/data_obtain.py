# -*- coding: utf-8 -*-
"""
data_obtain.py —— 因果考古 (RLCD) 项目的 TNG50 数据提取脚本
基于 001_orbit_new.py (最终版), 在 TNG JupyterLab 服务器上运行。

与 001_orbit_new.py 的区别
------------------------------------------------------------
  沿用: 出生半径 (galpos + TNG50_disk_trajectory.hdf5)、丰度公式 [I][J]、
        风粒子剔除、tform 上边界夹取 [B]、对齐断言 [A]、年轻星自检 [L]、
        BirthPos 必须在 face_on 之前提取 [F]。
  新增:
    [N] 保存 z=0 盘坐标系下恒星的 pos + vel (算作用量 Jz/Jphi、偏心率 e、Rg 用)
    [O] 保存暗物质 / 气体 / 黑洞 / 风粒子的 pos + mass (本地用 AGAMA 拟合引力势用)
    [P] 保留"早于轨道起点"的恒星 (birth_R = NaN, has_birth = False):
        它们没有出生半径, 但它们的质量是引力势的一部分, 不能丢
    [Q] 另存原始质量分数 X_H / X_Mg / X_Fe, 以后换太阳丰度标度不用重跑
    [R] 输出改为每星系一个 HDF5 (float32 + gzip), 比 txt 小且快;
        先写临时文件再改名, 中断后重跑会自动跳过已完成的星系

输出: ./data_obtain/<ID>.hdf5
  /header (attrs)  ID, snap, a, h, time_Gyr, face_on 参数, 单位说明
  /frame           center_pos, center_vel, rotation   (z=0 face_on 返回值)
  /stars           iord, pos(N,3), vel(N,3), mass, init_mass, age, aform,
                   FeH, MgFe, X_H, X_Mg, X_Fe,
                   birth_pos(N,3), birth_R_2D, birth_R_3D, has_birth
  /wind            pos, mass                 (风粒子, 只用于引力势)
  /dm              pos, mass
  /gas             pos, mass
  /bh              pos, mass
单位: pos kpc (物理), vel km/s, mass Msol, age Gyr。
坐标系: z=0 时以星系中心为原点、恒星角动量 (r<8 kpc) 为 z 轴 —— 与
001_orbit_new.py 的 Current_R_2D / Current_z 完全相同。
birth_pos 是 galpos 给出的"出生时刻盘坐标系"坐标, 与 z=0 坐标系不同。
"""

from AnastrisTNG import TNGsimulation
import numpy as np
import os
import h5py
from galpos import GalaxyPoseTrajectory
from galpos import pynbody_decorate

# ==================== 设置 ====================
run = 'TNG50'
BasePath = '/home/tnguser/sims.TNG/' + run + '-1/output'
snap99 = TNGsimulation.Snapshot(BasePath, 99)

# 默认已加载 dm,star,gas,bh 四类粒子的 Coordinates/Velocities/Masses/ParticleIDs;
# 这里只给恒星追加额外字段
snap99.load_particle_para['star_fields'] = [
    'GFM_StellarFormationTime', 'BirthPos', 'BirthVel', 'GFM_InitialMass',
    'GFM_Metals'
]

orbit_file_path = 'TNG50_disk_trajectory.hdf5'   # 相对当前运行目录; 可改成绝对路径
if not os.path.exists(orbit_file_path):
    raise FileNotFoundError(
        f"找不到轨道文件 {os.path.abspath(orbit_file_path)}。请把 TNG50_disk_trajectory.hdf5 "
        f"软链接到当前目录 (ln -s <原路径> .), 或修改 orbit_file_path。")

# 607654 不在 423 样本 (轨道文件) 中, 已删除
halo_ids = np.array([487742, 530330, 538905, 545437, 546870,
                     547293, 552581, 557721, 560421, 562029,
                     567607, 569251])

FACE_ON_RMAX = 8          # kpc, 与 001_orbit_new.py 一致
FACE_ON_ALIGN = 'star'

output_dir = './data_obtain/'
os.makedirs(output_dir, exist_ok=True)

print(f"将处理 {len(halo_ids)} 个星系: {halo_ids}")

# ==================== [I] 丰度计算常数 ====================
A_H, A_Mg, A_Fe = 1.008, 24.305, 55.845
LOGEPS_MG_SUN = 7.60      # Asplund et al. 2009
LOGEPS_FE_SUN = 7.50
SUN_FEH_TERM  = LOGEPS_FE_SUN - 12.0            # = -4.50
SUN_MGFE_TERM = LOGEPS_MG_SUN - LOGEPS_FE_SUN   # = +0.10

T_EDGE_TOL = 1e-3         # [B] Gyr

processing_errors = []
ledger_all = {}


def load_galaxy_trajectory(orbit_file, galaxy_id):
    with h5py.File(orbit_file, "r") as file:
        if str(galaxy_id) not in file:
            return None
        g = file[str(galaxy_id)]
        snap_all = np.sort(np.int16(list(g.keys())))
        return {
            'snap_all': snap_all,
            'ang_mom': np.array([g[f"{i}"]['ang_mom'][...] for i in snap_all]),
            'cpos':    np.array([g[f"{i}"]['cpos'][...]    for i in snap_all]),
            'cvel':    np.array([g[f"{i}"]['cvel'][...]    for i in snap_all]),
            't':       np.array([g[f"{i}"]['t'][...]       for i in snap_all]),
        }


def arr(x, unit=None, dtype=np.float32):
    """SimArray -> 普通 numpy 数组 (可选单位换算)"""
    if unit is not None:
        x = x.in_units(unit)
    return np.asarray(x.view(np.ndarray), dtype=dtype).copy()


def write_group(h5, name, data):
    grp = h5.create_group(name)
    for k, v in data.items():
        v = np.asarray(v)
        if v.ndim >= 1 and v.size > 1000:
            grp.create_dataset(k, data=v, compression='gzip',
                               compression_opts=4, shuffle=True)
        else:
            grp.create_dataset(k, data=v)


# ==================== 主循环 ====================
for ID in halo_ids:
    ID = int(ID)
    out_file = os.path.join(output_dir, f'{ID}.hdf5')
    if os.path.exists(out_file):
        print(f"Subhalo {ID}: 输出已存在, 跳过")
        continue
    print(f"Processing Subhalo ID: {ID}")

    try:
        # ---------- 轨道数据 ----------
        trajectory_data = load_galaxy_trajectory(orbit_file_path, ID)
        if trajectory_data is None:
            raise RuntimeError(f"Galaxy {ID}: 轨道文件中无此星系")

        t_min_orbit = trajectory_data['t'].min()
        t_max_orbit = trajectory_data['t'].max()
        print(f"    轨道时间范围: {t_min_orbit:.3f} - {t_max_orbit:.3f} Gyr")

        orbits = GalaxyPoseTrajectory(
            trajectory_data['t'],
            trajectory_data['cpos'],
            trajectory_data['cvel'],
            angular_momentum=trajectory_data['ang_mom'],
            box_size=35000/0.6774
        )

        # ---------- 粒子数据 (dm, star, gas, bh 全部加载) ----------
        testsub = snap99.load_particle(ID)
        testsub.physical_units()

        n_total = len(testsub.s)
        if n_total == 0:
            print(f"    No stellar particles for {ID}, skip")
            continue

        # ==========================================================
        # [F] 以下恒星数组必须在 face_on 之前提取
        # ==========================================================
        birth_pos  = testsub.s['BirthPos'].in_units('a kpc').view(np.ndarray).copy()
        birth_vel  = testsub.s['BirthVel'].in_units('kpc Gyr**-1 a**1/2').view(np.ndarray).copy()
        mass       = testsub.s['mass'].in_units('Msol').view(np.ndarray).copy()
        init_mass  = testsub.s['GFM_InitialMass'].in_units('Msol').view(np.ndarray).copy()
        age   = testsub.s['age'].in_units('Gyr').view(np.ndarray).copy()
        aform = testsub.s['aform'].view(np.ndarray).copy()
        tform = testsub.s['tform'].in_units('Gyr').view(np.ndarray).copy()
        iord  = testsub.s['iord'].view(np.ndarray).copy()
        metals_all = testsub.s['GFM_Metals'].view(np.ndarray).copy()

        X_H_all  = metals_all[:, 0]
        X_Mg_all = metals_all[:, 6]
        X_Fe_all = metals_all[:, 8]

        # ---------- 剔除风粒子 ----------
        valid_mask = (
            (aform > 0) & (aform <= 1) &
            (tform > 0) & np.isfinite(tform) & np.isfinite(aform)
        )
        n_valid = int(np.sum(valid_mask))
        n_wind  = n_total - n_valid
        print(f"    恒星族粒子: {n_total}, 真实恒星: {n_valid}, 风粒子: {n_wind}")
        if n_valid == 0:
            print(f"    警告: {ID} 没有有效恒星, 跳过")
            continue

        mass_valid      = mass[valid_mask]
        age_valid       = age[valid_mask]
        aform_valid     = aform[valid_mask]
        tform_valid     = tform[valid_mask]
        X_H_v  = X_H_all[valid_mask]
        X_Mg_v = X_Mg_all[valid_mask]
        X_Fe_v = X_Fe_all[valid_mask]

        # ---------- [I][J] 丰度 ----------
        feh_ok  = (X_H_v > 0)  & (X_Fe_v > 0)
        mgfe_ok = (X_Fe_v > 0) & (X_Mg_v > 0)
        FeH  = np.full(n_valid, np.nan)
        MgFe = np.full(n_valid, np.nan)
        FeH[feh_ok]   = np.log10((X_Fe_v[feh_ok] / A_Fe) / (X_H_v[feh_ok] / A_H)) - SUN_FEH_TERM
        MgFe[mgfe_ok] = np.log10((X_Mg_v[mgfe_ok] / A_Mg) / (X_Fe_v[mgfe_ok] / A_Fe)) - SUN_MGFE_TERM

        # ---------- [B] tform 上边界夹取 ----------
        edge_mask = (tform_valid > t_max_orbit) & (tform_valid < t_max_orbit + T_EDGE_TOL)
        n_edge = int(np.sum(edge_mask))
        tform_for_orbit = tform_valid.copy()
        tform_for_orbit[edge_mask] = t_max_orbit
        n_early = int(np.sum(tform_valid < t_min_orbit))
        print(f"    上边界夹取: {n_edge}, 早于轨道起点: {n_early}")

        # ---------- 出生坐标变换 (galpos) ----------
        birth = pynbody_decorate.make_star_birth(
            orbits,
            birth_pos=birth_pos[valid_mask],
            mass=mass_valid,
            birth_velocity=birth_vel[valid_mask],
            birth_time=tform_for_orbit,
            birth_pos_units="a kpc",
            birth_velocity_units='kpc Gyr**-1 a**1/2',
            scale_factor=aform_valid
        )
        birth.align_with_galaxy(orientation_align=True)

        # [A] 对齐断言
        assert len(birth.s) == n_valid, f"{ID}: make_star_birth 改变了粒子数"
        assert np.allclose(np.asarray(mass_valid, dtype=np.float64),
                           np.asarray(birth.s['mass'], dtype=np.float64),
                           rtol=1e-6), f"{ID}: make_star_birth 可能重排了粒子顺序"

        bpos = np.asarray(birth.s['pos'], dtype=np.float64)
        birth_r_2d = np.sqrt(bpos[:, 0]**2 + bpos[:, 1]**2)
        birth_r_3d = np.sqrt(np.sum(bpos**2, axis=1))
        has_birth = np.isfinite(birth_r_2d) & np.isfinite(birth_r_3d)   # [P]
        n_nan = int(np.sum(~has_birth))
        if n_nan != n_early:
            print(f"    *** 警告 *** NaN 数({n_nan}) != 早于轨道起点数({n_early})")

        # ---------- z=0 face_on (整个子晕一起旋转: dm/star/gas/bh) ----------
        center_pos, center_vel, rot = testsub.face_on(rmax=FACE_ON_RMAX,
                                                      alignwith=FACE_ON_ALIGN)

        s_pos = arr(testsub.s['pos'], 'kpc', np.float64)
        s_vel = arr(testsub.s['vel'], 'km s**-1', np.float64)

        # [L] 年轻星自检: 出生半径 vs 当前半径
        cur_r_2d = np.sqrt(s_pos[valid_mask, 0]**2 + s_pos[valid_mask, 1]**2)
        for lo, hi in [(0, 0.02), (0.02, 0.05), (0.05, 0.1)]:
            m = (age_valid >= lo) & (age_valid < hi) & has_birth
            if np.sum(m) > 10:
                d = birth_r_2d[m] - cur_r_2d[m]
                print(f"    [自检] age {lo:.2f}-{hi:.2f} Gyr ({np.sum(m):>5d}个): "
                      f"中位数 {np.median(d):+.3f} kpc, |dr|95% = "
                      f"{np.percentile(np.abs(d), 95):.3f} kpc")

        # ---------- 收集各类粒子 ----------
        stars = {
            'iord':       iord[valid_mask].astype(np.uint64),
            'pos':        s_pos[valid_mask].astype(np.float32),
            'vel':        s_vel[valid_mask].astype(np.float32),
            'mass':       mass_valid.astype(np.float32),
            'init_mass':  init_mass[valid_mask].astype(np.float32),
            'age':        age_valid.astype(np.float32),
            'aform':      aform_valid.astype(np.float32),
            'FeH':        FeH.astype(np.float32),
            'MgFe':       MgFe.astype(np.float32),
            'X_H':        X_H_v.astype(np.float32),
            'X_Mg':       X_Mg_v.astype(np.float32),
            'X_Fe':       X_Fe_v.astype(np.float32),
            'birth_pos':  bpos.astype(np.float32),
            'birth_R_2D': birth_r_2d.astype(np.float32),
            'birth_R_3D': birth_r_3d.astype(np.float32),
            'has_birth':  has_birth,
        }
        wind = {
            'pos':  s_pos[~valid_mask].astype(np.float32).reshape(-1, 3),
            'mass': mass[~valid_mask].astype(np.float32),
        }

        others = {}
        for fam_name, fam_attr in [('dm', 'dm'), ('gas', 'g'), ('bh', 'bh')]:
            try:
                sub = getattr(testsub, fam_attr)
                n = len(sub)
            except Exception:
                n = 0
            if n > 0:
                others[fam_name] = {
                    'pos':  arr(sub['pos'], 'kpc'),
                    'mass': arr(sub['mass'], 'Msol'),
                }
            else:
                others[fam_name] = {
                    'pos':  np.zeros((0, 3), np.float32),
                    'mass': np.zeros(0, np.float32),
                }

        # ---------- 写 HDF5 (先写临时文件再改名) ----------
        tmp_file = out_file + '.tmp'
        with h5py.File(tmp_file, 'w') as h5:
            hdr = h5.create_group('header')
            hdr.attrs['ID'] = ID
            hdr.attrs['snap'] = 99
            for key, getter in [('a', lambda: float(testsub.properties['a'])),
                                ('h', lambda: float(testsub.properties['h'])),
                                ('time_Gyr', lambda: float(testsub.properties['t'].in_units('Gyr')))]:
                try:
                    hdr.attrs[key] = getter()
                except Exception:
                    hdr.attrs[key] = np.nan
            hdr.attrs['face_on_rmax_kpc'] = FACE_ON_RMAX
            hdr.attrs['face_on_alignwith'] = FACE_ON_ALIGN
            hdr.attrs['units'] = 'pos: kpc (physical); vel: km/s; mass: Msol; age: Gyr'
            hdr.attrs['solar_scale'] = 'Asplund 2009: logeps_Fe=7.50, logeps_Mg=7.60'
            hdr.attrs['script'] = 'data_obtain.py'

            write_group(h5, 'frame', {
                'center_pos': arr(center_pos, None, np.float64),
                'center_vel': arr(center_vel, None, np.float64),
                'rotation':   np.asarray(rot, dtype=np.float64),
            })
            h5['frame'].attrs['center_pos_units'] = str(center_pos.units)
            h5['frame'].attrs['center_vel_units'] = str(center_vel.units)

            write_group(h5, 'stars', stars)
            write_group(h5, 'wind', wind)
            for fam_name, d in others.items():
                write_group(h5, fam_name, d)
        os.replace(tmp_file, out_file)

        size_mb = os.path.getsize(out_file) / 1e6
        ledger_all[ID] = {
            'star_family': n_total, 'wind': n_wind, 'valid': n_valid,
            'edge_clipped': n_edge, 'no_birth': n_nan,
            'dm': len(others['dm']['mass']), 'gas': len(others['gas']['mass']),
            'bh': len(others['bh']['mass']), 'MB': size_mb,
        }
        print(f"    已保存 {out_file} ({size_mb:.1f} MB)")
        print(f"    dm {ledger_all[ID]['dm']}, gas {ledger_all[ID]['gas']}, "
              f"bh {ledger_all[ID]['bh']}, stars {n_valid} (其中无出生半径 {n_nan})")
        print()

        del testsub, birth

    except Exception as e:
        print(f"*** ERROR *** Subhalo {ID}: {e}")
        processing_errors.append((ID, str(e)))
        tmp_file = out_file + '.tmp'
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        continue

# =================== 最终报告 ===================
print("\n" + "=" * 90)
if processing_errors:
    print(f"*** {len(processing_errors)} 个星系出错:")
    for gid, err in processing_errors:
        print(f"  - {gid}: {err}")
print(f"{'ID':>8} {'恒星族':>9} {'风':>6} {'恒星':>9} {'无出生R':>7} "
      f"{'DM':>9} {'gas':>9} {'bh':>4} {'MB':>7}")
for gid, L in ledger_all.items():
    print(f"{gid:>8} {L['star_family']:>9} {L['wind']:>6} {L['valid']:>9} "
          f"{L['no_birth']:>7} {L['dm']:>9} {L['gas']:>9} {L['bh']:>4} {L['MB']:>7.1f}")
print(f"\n本次成功: {len(ledger_all)} / {len(halo_ids)}, 输出目录: {output_dir}")
print("=" * 90)
