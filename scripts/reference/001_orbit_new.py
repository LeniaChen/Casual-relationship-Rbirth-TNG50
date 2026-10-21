# -*- coding: utf-8 -*-
"""
001_orbit 数据提取脚本 —— 最终版 (第 1 类 + 第 2 类修改 + Current_z 列)

第 1 类修改 (逻辑/对齐, 已验证):
  [A] make_star_birth 后对齐断言 (487742 验证: 保长保序)
  [B] tform 上边界浮点夹取 (z=0 出生恒星舍入误判)
  [C] 缺轨道数据直接 raise, 删除回退分支 (13 星系轨道已确认齐全)
  [D] 输出真实粒子 ID (iord)
  [E] 空星系检查移到 load_particle 之后
  [F] BirthPos/BirthVel 必须在 face_on 之前提取 (锁死顺序)
  [G] 删除未使用的 gas_fields / 轨道 align,a,re 读取
  [H] 每星系粒子账本

第 2 类修改 (数值正确性):
  [I] 丰度公式: 数密度比 + Asplund 2009 logε 太阳参考
      [Fe/H]  = log10( (X_Fe/A_Fe) / (X_H /A_H ) ) + 4.50
      [Mg/Fe] = log10( (X_Mg/A_Mg) / (X_Fe/A_Fe) ) - 0.10
  [J] FeH / MgFe 掩码解耦
  [K] 输出增加 InitMass 列
  [L] 年轻星自检 (年龄分箱): 验证坐标链路端到端正确性

  [M] 输出增加 Current_z 列 (z=0 face_on 盘坐标系垂向坐标, 物理 kpc)。
      输出全部有效恒星, 不做 |z| 筛选 —— 盘星选择 (如 |z|<0.3 kpc)
      留给下游分析脚本, 避免改阈值时重跑昂贵的提取流程。

已核实 (基于轨道文件生成代码):
  cpos = 无h共动kpc (physical_units 后 in_units("a kpc")), 与 box_size=35000/0.6774 一致
  t 与主脚本 tform 同源 (pynbody 同一套宇宙学换算)
  ang_mom 为盒子坐标系方向矢量, 与 BirthPos 同系
备忘 (不影响本脚本输出):
  cvel 存为 "a kpc Gyr**-1", 与 BirthVel 的 "kpc Gyr**-1 a**1/2" 约定不同;
  当前只输出半径/坐标不用速度。将来若使用出生速度, 必须先核对
  GalaxyPoseTrajectory 对 cvel 单位的期望。
  轨道早期 snapshot 的 alignwith 会切换到 baryon/all, 最老恒星的"盘面"
  定义与后期不同, 早期分析建议用 Birth_R_3D。
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
# [G] 本脚本不使用气体, 不加载 gas_fields
snap99.load_particle_para['star_fields'] = [
    'GFM_StellarFormationTime', 'BirthPos', 'BirthVel', 'GFM_InitialMass',
    'GFM_Metals'
]

orbit_file_path = 'TNG50_disk_trajectory.hdf5'

halo_ids = np.array([487742, 530330, 538905, 545437, 546870,
                     547293, 552581, 557721, 560421, 562029,
                     567607, 569251, 607654])

print(f"将处理 {len(halo_ids)} 个星系: {halo_ids}")

output_dir = './001_orbit_new/'
os.makedirs(output_dir, exist_ok=True)

# ==================== [I] 丰度计算常数 ====================
# 原子量
A_H, A_Mg, A_Fe = 1.008, 24.305, 55.845
# Asplund et al. 2009 光球层丰度 logε (logε_H = 12.00 标度)
LOGEPS_MG_SUN = 7.60
LOGEPS_FE_SUN = 7.50
SUN_FEH_TERM  = LOGEPS_FE_SUN - 12.0            # = -4.50 → 公式中为 +4.50
SUN_MGFE_TERM = LOGEPS_MG_SUN - LOGEPS_FE_SUN   # = +0.10 → 公式中为 -0.10

# [B] tform 上边界夹取容差 (Gyr): 远大于浮点舍入(~1e-6), 远小于物理时标
T_EDGE_TOL = 1e-3

processing_errors = []
ledger_all = {}   # [H] 粒子账本

def load_galaxy_trajectory(orbit_file, galaxy_id):
    """[G] 只读实际用到的字段: t, cpos, cvel, ang_mom"""
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

# ==================== 主循环 ====================
for ID in halo_ids:
    ID = int(ID)
    print(f"Processing Subhalo ID: {ID}")

    try:
        # ---------- 轨道数据 ----------
        trajectory_data = load_galaxy_trajectory(orbit_file_path, ID)
        # [C] 缺轨道 = 硬错误
        if trajectory_data is None:
            raise RuntimeError(
                f"Galaxy {ID}: 轨道文件中无此星系。已确认 13 个目标星系轨道齐全, "
                f"出现此错误说明轨道文件或 ID 列表被改动, 请检查。")

        t_min_orbit = trajectory_data['t'].min()
        t_max_orbit = trajectory_data['t'].max()
        print(f"    轨道时间范围: {t_min_orbit:.3f} - {t_max_orbit:.3f} Gyr")

        orbits = GalaxyPoseTrajectory(
            trajectory_data['t'],
            trajectory_data['cpos'],
            trajectory_data['cvel'],
            angular_momentum=trajectory_data['ang_mom'],
            box_size=35000/0.6774   # 已核实: 与轨道文件 cpos 同为无h共动kpc
        )

        # ---------- 粒子数据 ----------
        testsub = snap99.load_particle(ID)
        testsub.physical_units()

        # [E] 立刻检查空星系
        n_total = len(testsub.s['mass'])
        if n_total == 0:
            print(f"    No stellar particles found for Subhalo ID: {ID}, skip")
            continue

        # ==========================================================
        # [F] !!! 提取顺序锁死: 以下所有数组必须在 face_on 之前提取 !!!
        # pynbody 的 face_on 会就地旋转 SimSnap 中所有 3-vector 数组
        # (包括 BirthPos/BirthVel)。BirthPos/BirthVel 的坐标变换由
        # galpos 负责 (用各恒星出生时刻的中心和盘面), 绝不能先被
        # face_on 用 z=0 的中心/盘面污染。
        # ==========================================================
        birth_pos  = testsub.s['BirthPos'].in_units('a kpc').view(np.ndarray).copy()
        birth_vel  = testsub.s['BirthVel'].in_units('kpc Gyr**-1 a**1/2').view(np.ndarray).copy()
        mass       = testsub.s['mass'].in_units('Msol').view(np.ndarray).copy()
        init_mass  = testsub.s['GFM_InitialMass'].in_units('Msol').view(np.ndarray).copy()  # [K]
        age   = testsub.s['age'].in_units('Gyr').view(np.ndarray).copy()
        aform = testsub.s['aform'].view(np.ndarray).copy()
        tform = testsub.s['tform'].in_units('Gyr').view(np.ndarray).copy()
        iord  = testsub.s['iord'].view(np.ndarray).copy()      # [D] 真实粒子 ID
        metals_all = testsub.s['GFM_Metals'].view(np.ndarray).copy()  # (N, 10)

        # TNG GFM_Metals 顺序: [H, He, C, N, O, Ne, Mg, Si, Fe, other]
        X_H_all  = metals_all[:, 0]
        X_Mg_all = metals_all[:, 6]
        X_Fe_all = metals_all[:, 8]

        # ---------- 筛选: 剔除风粒子 (aform<0 为风粒子, TNG 官方约定) ----------
        valid_mask = (
            (aform > 0) & (aform <= 1) &
            (tform > 0) & np.isfinite(tform) & np.isfinite(aform)
        )
        n_valid = int(np.sum(valid_mask))
        n_wind  = n_total - n_valid
        print(f"    过滤前恒星粒子数: {n_total}")
        print(f"    过滤后恒星粒子数: {n_valid} (剔除风粒子等 {n_wind})")

        if n_valid == 0:
            print(f"    警告: Galaxy ID {ID} 没有有效的恒星粒子, 跳过")
            continue

        birth_pos_valid = birth_pos[valid_mask]
        birth_vel_valid = birth_vel[valid_mask]
        mass_valid      = mass[valid_mask]
        init_mass_valid = init_mass[valid_mask]     # [K]
        age_valid   = age[valid_mask]
        aform_valid = aform[valid_mask]
        tform_valid = tform[valid_mask]
        iord_valid  = iord[valid_mask]

        X_H_v  = X_H_all[valid_mask]
        X_Mg_v = X_Mg_all[valid_mask]
        X_Fe_v = X_Fe_all[valid_mask]

        # ---------- [I][J] 丰度计算 (数密度比 + Asplund 2009) ----------
        feh_ok  = (X_H_v > 0)  & (X_Fe_v > 0)
        mgfe_ok = (X_Fe_v > 0) & (X_Mg_v > 0)

        FeH  = np.full(n_valid, np.nan)
        MgFe = np.full(n_valid, np.nan)
        FeH[feh_ok]   = np.log10((X_Fe_v[feh_ok]  / A_Fe) / (X_H_v[feh_ok]   / A_H)) \
                        - SUN_FEH_TERM                     # 即 + 4.50
        MgFe[mgfe_ok] = np.log10((X_Mg_v[mgfe_ok] / A_Mg) / (X_Fe_v[mgfe_ok] / A_Fe)) \
                        - SUN_MGFE_TERM                    # 即 - 0.10

        # ---------- [B] tform 上边界浮点夹取 ----------
        # z=0 出生恒星 (aform=1) 的 tform 与轨道终点物理上同一时刻,
        # 但两套换算的浮点舍入可差 ~1e-6 Gyr → 插值器误判越界输出 NaN。
        # 夹回 t_max 是消除舍入差, 不是外推。
        # 下边界 (早于轨道起点) 不夹取: 真实超出覆盖范围, 保持 NaN 剔除。
        edge_mask = (tform_valid > t_max_orbit) & (tform_valid < t_max_orbit + T_EDGE_TOL)
        n_edge = int(np.sum(edge_mask))
        tform_for_orbit = tform_valid.copy()
        tform_for_orbit[edge_mask] = t_max_orbit

        n_early = int(np.sum(tform_valid < t_min_orbit))
        print(f"    上边界夹取 (z=0 出生, 浮点舍入): {n_edge}")
        print(f"    早于轨道起点 (将被 NaN 剔除):   {n_early}")

        # ---------- 出生坐标变换 (galpos) ----------
        birth = pynbody_decorate.make_star_birth(
            orbits,
            birth_pos=birth_pos_valid,
            mass=mass_valid,
            birth_velocity=birth_vel_valid,
            birth_time=tform_for_orbit,          # [B] 夹取后的时间
            birth_pos_units="a kpc",
            birth_velocity_units='kpc Gyr**-1 a**1/2',
            scale_factor=aform_valid
        )
        birth.align_with_galaxy(orientation_align=True)

        # [A] 对齐保险 (已验证; galpos 若更新导致行为改变会立刻在此报错)
        assert len(birth.s) == n_valid, (
            f"Galaxy {ID}: make_star_birth 改变了粒子数 "
            f"{n_valid} -> {len(birth.s)}, 逐行配对失效!")
        assert np.allclose(np.asarray(mass_valid, dtype=np.float64),
                           np.asarray(birth.s['mass'], dtype=np.float64),
                           rtol=1e-6), (
            f"Galaxy {ID}: make_star_birth 输出质量与输入不逐元素相等, "
            f"可能重排了粒子顺序, 逐行配对失效!")

        # ---------- 当前坐标 (face_on 之后提取) ----------
        testsub.face_on(rmax=8, alignwith="star")   # 就地旋转到 z=0 盘坐标系
        _pos = np.array(testsub.s['pos'])
        current_r_2d = np.sqrt(_pos[valid_mask, 0]**2 + _pos[valid_mask, 1]**2)
        current_r_3d = np.sqrt(np.sum(_pos[valid_mask]**2, axis=1))
        current_z    = _pos[valid_mask, 2]          # [M] z=0 盘坐标系垂向坐标 (物理 kpc)

        # ---------- 出生半径 ----------
        _bpos = np.asarray(birth.s['pos'], dtype=np.float64)
        birth_r_2d = np.sqrt(_bpos[:, 0]**2 + _bpos[:, 1]**2)
        birth_r_3d = np.sqrt(np.sum(_bpos**2, axis=1))

        # ---------- 剔除 NaN (应仅为早于轨道起点的粒子) ----------
        finite_mask = np.isfinite(birth_r_2d) & np.isfinite(birth_r_3d)
        n_nan = int(np.sum(~finite_mask))
        n_final = int(np.sum(finite_mask))
        print(f"    NaN 剔除: {n_nan} 个, 剩余: {n_final} 个")

        # [B] 一致性检查: NaN 只能来自早于轨道起点的粒子
        if n_nan != n_early:
            print(f"    *** 警告 *** NaN 数({n_nan}) != 早于轨道起点数({n_early}), "
                  f"存在未解释的 NaN 来源, 需要检查!")

        # [L] 自检: dr 随年龄分箱 (真实运动→随龄增大; 坐标误差→与龄无关)
        for lo, hi in [(0, 0.02), (0.02, 0.05), (0.05, 0.1)]:
            m = (age_valid >= lo) & (age_valid < hi) & finite_mask
            if np.sum(m) > 10:
                d = birth_r_2d[m] - current_r_2d[m]
                print(f"    [自检] age {lo:.2f}-{hi:.2f} Gyr ({np.sum(m):>5d}个): "
                      f"中位数 {np.median(d):+.3f} kpc, "
                      f"|dr|95% = {np.percentile(np.abs(d), 95):.3f} kpc")

        # ---------- [H] 粒子账本 ----------
        ledger_all[ID] = {
            'total': n_total, 'wind_removed': n_wind, 'valid': n_valid,
            'edge_clipped': n_edge, 'early_nan_removed': n_nan, 'final': n_final,
        }

        # ---------- 保存 (全部有效恒星, 盘星筛选留给下游) ----------
        output_file = os.path.join(output_dir, f'{ID}.txt')
        with open(output_file, 'w') as f:
            # [D] ParticleID = TNG 真实 ParticleIDs (iord)
            # [K] Mass = z=0 当前质量; InitMass = GFM_InitialMass (形成时质量)
            # [M] Current_z = z=0 face_on 盘坐标系垂向坐标 (物理 kpc, 带符号)
            #     盘星选择示例 (下游): |Current_z| < 0.3
            f.write(f"{'ParticleID':>15} {'Age_Gyr':>10} {'Current_R_2D':>14} {'Birth_R_2D':>12} {'Current_R_3D':>14} {'Birth_R_3D':>12} {'Current_z':>11} {'Mass':>12} {'InitMass':>12} {'FeH':>10} {'MgFe':>10}\n")
            for i in np.where(finite_mask)[0]:
                feh_str  = f"{FeH[i]:>10.6f}"  if np.isfinite(FeH[i])  else f"{'nan':>10}"
                mgfe_str = f"{MgFe[i]:>10.6f}" if np.isfinite(MgFe[i]) else f"{'nan':>10}"
                f.write(f"{iord_valid[i]:>15d} {age_valid[i]:>10.6f} {current_r_2d[i]:>14.6f} {birth_r_2d[i]:>12.6f} {current_r_3d[i]:>14.6f} {birth_r_3d[i]:>12.6f} {current_z[i]:>11.6f} {mass_valid[i]:>12.6f} {init_mass_valid[i]:>12.6f} {feh_str} {mgfe_str}\n")

        print(f"Subhalo {ID}: {n_final} stellar particles saved")
        print(f"  账本: 总数 {n_total} → 剔风粒子 {n_wind} → 有效 {n_valid} "
              f"→ 边界夹回 {n_edge} → NaN 剔除 {n_nan} → 最终 {n_final}")
        print(f"  Age range:          {np.nanmin(age_valid[finite_mask]):.2f} - {np.nanmax(age_valid[finite_mask]):.2f} Gyr")
        print(f"  Current_R_2D range: {np.nanmin(current_r_2d[finite_mask]):.2f} - {np.nanmax(current_r_2d[finite_mask]):.2f} kpc")
        print(f"  Birth_R_2D range:   {np.nanmin(birth_r_2d[finite_mask]):.2f} - {np.nanmax(birth_r_2d[finite_mask]):.2f} kpc")
        print(f"  Current_z range:    {np.nanmin(current_z[finite_mask]):+.3f} ~ {np.nanmax(current_z[finite_mask]):+.3f} kpc")
        n_disk_est = int(np.sum(finite_mask & (np.abs(current_z) < 0.3)))
        print(f"  |z|<0.3 kpc 粒子数 (供参考): {n_disk_est} / {n_final}")
        print(f"  Mass range:         {np.nanmin(mass_valid[finite_mask]):.3f} - {np.nanmax(mass_valid[finite_mask]):.3f} Msol")
        print(f"  InitMass range:     {np.nanmin(init_mass_valid[finite_mask]):.3f} - {np.nanmax(init_mass_valid[finite_mask]):.3f} Msol")
        feh_fin  = FeH[finite_mask]
        mgfe_fin = MgFe[finite_mask]
        print(f"  [Fe/H]:  有效 {np.sum(np.isfinite(feh_fin))} / {n_final}, "
              f"中位数 = {np.nanmedian(feh_fin):+.3f}, "
              f"范围 {np.nanmin(feh_fin):+.2f} ~ {np.nanmax(feh_fin):+.2f}")
        print(f"  [Mg/Fe]: 有效 {np.sum(np.isfinite(mgfe_fin))} / {n_final}, "
              f"中位数 = {np.nanmedian(mgfe_fin):+.3f}, "
              f"范围 {np.nanmin(mgfe_fin):+.2f} ~ {np.nanmax(mgfe_fin):+.2f}")
        print()

    except Exception as e:
        error_info = f"Error processing Subhalo ID {ID}: {str(e)}"
        print(f"*** ERROR *** {error_info}")
        processing_errors.append((ID, str(e)))
        continue

# =================== 最终报告 ===================
print("\n" + "="*80)
print("PROCESSING COMPLETED - FINAL REPORT")
print("="*80)

if processing_errors:
    print(f"*** ERROR *** {len(processing_errors)} galaxies had PROCESSING ERRORS:")
    for gal_id, error in processing_errors:
        print(f"  - Galaxy ID {gal_id}: {error}")
    print()

print(f"{'ID':>8} {'总数':>9} {'剔风':>6} {'有效':>9} {'夹回':>5} {'NaN剔':>6} {'最终':>9}")
for gid, L in ledger_all.items():
    print(f"{gid:>8} {L['total']:>9} {L['wind_removed']:>6} {L['valid']:>9} "
          f"{L['edge_clipped']:>5} {L['early_nan_removed']:>6} {L['final']:>9}")

print(f"\nSUMMARY:")
print(f"  Total galaxies: {len(halo_ids)}")
print(f"  Successfully processed: {len(ledger_all)}")
print(f"\nData saved in directory: {output_dir}")
print("="*80)