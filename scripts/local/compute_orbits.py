# -*- coding: utf-8 -*-
"""
compute_orbits.py —— 用 AGAMA 为每个星系拟合引力势, 并对所有恒星计算
    作用量 (J_R, J_z, J_phi)、能量 E、导引半径 Rg、近心距 / 远心距 / 偏心率 (积分轨道)、z_max

运行 (在 casuality 文件夹里, 用 venv 的 Python):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py --test     # 测试: 1 个星系, 随机 5000 颗星
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py            # 全部 12 个星系, 全部恒星
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py --ids 487742 530330   # 只算指定星系

输入:  data/tng50_extract/data_obtain_all.hdf5
输出:  data/tng50_derived/orbits_all.hdf5          (测试模式: orbits_test.hdf5)
         /<ID>/Jr, Jz, Jphi, E, Rg, rperi, rapo, ecc, zmax, n_period  —— 每个数组长度 = 该星系恒星数,
                                                                        顺序与 data_obtain_all.hdf5 的 /<ID>/stars 完全相同
         /<ID>/rotation_curve/R, vc                  —— 拟合势的圆周速度曲线 (检查用)
       data/tng50_derived/potentials/<ID>_halo.coef, <ID>_disk.coef  —— 势的展开系数, 可用 agama.Potential(file=...) 重新读入
已经算完的星系会自动跳过; 中途中断后重跑即可继续。

方法 (与 Jin et al. 论文一致: AGAMA + 轴对称多极展开, 包含暗物质、气体、恒星):
  - 坐标系: data_obtain 保存的 z=0 盘坐标系 (中心为原点, 恒星角动量为 +z)
  - 引力势 = 两部分之和, 都强制轴对称; 网格内边界 rmin = 0.2 kpc (约为 TNG50 的引力软化长度),
    更小半径处由 AGAMA 平滑外推 (若 ActionFinder 仍报错, 自动改用 0.5、1.0 kpc):
      "halo" 部分: Multipole 展开 (lmax=8), 包含暗物质、黑洞、风粒子、盘外的气体和恒星
      "disk" 部分: CylSpline 展开 (柱坐标网格), 包含 r<50 kpc 的恒星和 R<30 kpc 且 |z|<5 kpc 的气体
    (这是 AGAMA 处理模拟星系的标准做法: 扁平的盘用柱坐标网格描述更准确)
  - 作用量: AGAMA ActionFinder (Staeckel fudge); J_phi 就是角动量 z 分量 L_z
  - Rg: 满足 L_z = Rg * v_c(Rg) 的半径, AGAMA pot.Rcirc(L=|L_z|)
  - 偏心率: 在该势中积分每颗星的轨道 N_PERIOD 个圆轨道周期, 取三维距离的最小/最大值
        e = (r_apo - r_peri) / (r_apo + r_peri)
  - 未束缚 (E >= 0) 的恒星: 所有轨道量记为 NaN
单位: 长度 kpc, 速度 km/s, 质量 Msun, 作用量 kpc*km/s, 能量 (km/s)^2
"""
import os
import sys
import time
import argparse
import numpy as np
import h5py
import agama

# ==================== 参数 ====================
IN_FILE = 'data/tng50_extract/data_obtain_all.hdf5'
OUT_DIR = 'data/tng50_derived'
POT_DIR = os.path.join(OUT_DIR, 'potentials')

N_PERIOD = 15        # 积分时长 = N_PERIOD * 圆轨道周期 Tcirc(E)
TRAJ_SIZE = 1501     # 每条轨道记录的点数 (每个圆周期约 100 个点)
CHUNK = 5000         # 每批积分的轨道数 (控制内存: 5000*1501*6*4 字节 ≈ 180 MB)
ACCURACY = 1e-7      # 轨道积分相对精度

# 分配给 "disk" (CylSpline) 部分的粒子
DISK_STAR_RMAX = 50.0    # 恒星: 三维距离 r < 50 kpc
DISK_GAS_RMAX = 30.0     # 气体: 柱半径 R < 30 kpc
DISK_GAS_ZMAX = 5.0      #        且 |z| < 5 kpc

N_TEST = 5000

# 引力势网格的内边界 (kpc)。TNG50 在 z=0 的引力软化长度约 0.3 kpc, 更小的尺度上质量分布不可靠;
# 自动网格会延伸到 ~0.01 kpc, 在那里粒子噪声会让圆轨道角动量随半径不单调, 导致 ActionFinder 报错。
# 依次尝试这些内边界, 直到 ActionFinder 能正常初始化。
RMIN_TRY = [0.2, 0.5, 1.0]
RMAX_HALO = 500.0        # Multipole 外边界 (kpc)
RMAX_DISK = 60.0         # CylSpline R 和 z 的外边界 (kpc)
ZMIN_DISK = 0.05         # CylSpline z 方向第一个格点 (kpc)
POT_VERSION = 'v2-rmin-grid'   # 势的设置改变时更新; 结果文件里版本不同的星系会被重算

agama.setUnits(mass=1, length=1, velocity=1)   # Msun, kpc, km/s (时间单位 = 0.978 Gyr)
TIME_UNIT_GYR = 0.9778


def build_potential(g, gid, rmin):
    """用 dm / gas / stars / wind / bh 粒子构建轴对称引力势 (网格内边界 rmin kpc)"""
    rd = lambda grp: (g[grp]['pos'][:].astype(np.float64), g[grp]['mass'][:].astype(np.float64))
    s_pos, s_m = rd('stars')
    g_pos, g_m = rd('gas')
    d_pos, d_m = rd('dm')
    w_pos, w_m = rd('wind')
    b_pos, b_m = rd('bh')

    s_disk = np.linalg.norm(s_pos, axis=1) < DISK_STAR_RMAX
    g_disk = (np.hypot(g_pos[:, 0], g_pos[:, 1]) < DISK_GAS_RMAX) & (np.abs(g_pos[:, 2]) < DISK_GAS_ZMAX)

    disk_pos = np.vstack([s_pos[s_disk], g_pos[g_disk]])
    disk_m = np.concatenate([s_m[s_disk], g_m[g_disk]])
    halo_pos = np.vstack([d_pos, b_pos, w_pos, s_pos[~s_disk], g_pos[~g_disk]])
    halo_m = np.concatenate([d_m, b_m, w_m, s_m[~s_disk], g_m[~g_disk]])

    print(f"    势 (rmin = {rmin} kpc): disk 部分 {len(disk_m):,} 个粒子 ({disk_m.sum():.3e} Msun), "
          f"halo 部分 {len(halo_m):,} 个粒子 ({halo_m.sum():.3e} Msun)")

    pot_halo = agama.Potential(type='Multipole', particles=(halo_pos, halo_m),
                               symmetry='axisymmetric', lmax=8, gridSizeR=40,
                               rmin=rmin, rmax=RMAX_HALO)
    pot_disk = agama.Potential(type='CylSpline', particles=(disk_pos, disk_m),
                               symmetry='axisymmetric', mmax=0, gridSizeR=30, gridSizez=30,
                               rmin=rmin, rmax=RMAX_DISK, zmin=ZMIN_DISK, zmax=RMAX_DISK)
    os.makedirs(POT_DIR, exist_ok=True)
    pot_halo.export(os.path.join(POT_DIR, f'{gid}_halo.coef'))
    pot_disk.export(os.path.join(POT_DIR, f'{gid}_disk.coef'))
    return agama.Potential(pot_halo, pot_disk)


def peri_apo(pot, xv, T):
    """分批积分轨道, 返回 r_peri, r_apo, z_max"""
    n = len(xv)
    rperi = np.full(n, np.nan)
    rapo = np.full(n, np.nan)
    zmax = np.full(n, np.nan)
    t0 = time.time()
    for i0 in range(0, n, CHUNK):
        i1 = min(i0 + CHUNK, n)
        _, traj = agama.orbit(potential=pot, ic=xv[i0:i1], time=T[i0:i1],
                              trajsize=TRAJ_SIZE, accuracy=ACCURACY,
                              separateTime=True, verbose=False)
        traj = np.asarray(traj, dtype=np.float32)          # (m, TRAJ_SIZE, 6)
        r = np.sqrt(np.sum(traj[:, :, :3].astype(np.float64)**2, axis=2))
        rperi[i0:i1] = r.min(axis=1)
        rapo[i0:i1] = r.max(axis=1)
        zmax[i0:i1] = np.abs(traj[:, :, 2]).max(axis=1)
        done = i1
        el = time.time() - t0
        print(f"\r    轨道积分 {done:,}/{n:,}  已用 {el/60:.1f} 分钟, "
              f"预计还需 {el/done*(n-done)/60:.1f} 分钟   ", end='', flush=True)
    print()
    return rperi, rapo, zmax


def process_galaxy(fin, fout, gid, test=False, rng=None):
    g = fin[gid]
    t_start = time.time()
    print(f"\n==== 星系 {gid} ====")

    pot = af = None
    for rmin in RMIN_TRY:
        pot = build_potential(g, gid, rmin)
        try:
            af = agama.ActionFinder(pot, interp=False)
            break
        except RuntimeError as err:
            print(f"    rmin = {rmin} kpc 时 ActionFinder 初始化失败: {err}")
            print("    换更大的 rmin 重试")
    if af is None:
        raise RuntimeError(f"星系 {gid}: 所有 rmin 都失败 {RMIN_TRY}")

    pos = g['stars']['pos'][:].astype(np.float64)
    vel = g['stars']['vel'][:].astype(np.float64)
    n_all = len(pos)
    if test:
        idx = np.sort(rng.choice(n_all, size=min(N_TEST, n_all), replace=False))
    else:
        idx = np.arange(n_all)
    xv = np.hstack([pos[idx], vel[idx]])
    n = len(idx)
    print(f"    恒星数: {n:,}" + (f" (测试: 从 {n_all:,} 颗中随机抽取)" if test else ""))

    # ---- 能量、束缚与否 ----
    E = pot.potential(xv[:, :3]) + 0.5 * np.sum(xv[:, 3:]**2, axis=1)
    bound = np.isfinite(E) & (E < 0)
    print(f"    束缚恒星: {bound.sum():,} / {n:,}")

    out = {k: np.full(n, np.nan) for k in ['Jr', 'Jz', 'Jphi', 'Rg', 'rperi', 'rapo', 'ecc', 'zmax']}
    out['E'] = E
    out['Jphi'] = xv[:, 0] * xv[:, 4] - xv[:, 1] * xv[:, 3]     # L_z, 对所有恒星都有定义

    # ---- 作用量 ----
    t0 = time.time()
    J = af(xv[bound])
    out['Jr'][bound] = J[:, 0]
    out['Jz'][bound] = J[:, 1]
    print(f"    作用量完成 ({time.time()-t0:.1f} 秒)")

    # ---- 导引半径 ----
    out['Rg'][bound] = pot.Rcirc(L=np.abs(out['Jphi'][bound]))

    # ---- 轨道积分: 近心距、远心距、偏心率 ----
    T = N_PERIOD * pot.Tcirc(xv[bound])
    ok = np.isfinite(T) & (T > 0)
    b_idx = np.where(bound)[0][ok]
    rp, ra, zm = peri_apo(pot, xv[b_idx], T[ok])
    out['rperi'][b_idx] = rp
    out['rapo'][b_idx] = ra
    out['zmax'][b_idx] = zm
    out['ecc'][b_idx] = (ra - rp) / (ra + rp)

    # ---- 圆周速度曲线 ----
    Rgrid = np.linspace(0.5, 30, 60)
    pts = np.column_stack([Rgrid, 0 * Rgrid, 0 * Rgrid])
    vc = np.sqrt(-Rgrid * pot.force(pts)[:, 0])

    # ---- 写入 ----
    grp = fout.create_group(gid)
    grp.attrs['n_stars'] = n
    grp.attrs['test_mode'] = bool(test)
    grp.attrs['pot_version'] = POT_VERSION
    grp.attrs['pot_rmin_kpc'] = rmin
    grp.attrs['N_PERIOD'] = N_PERIOD
    grp.attrs['TRAJ_SIZE'] = TRAJ_SIZE
    grp.attrs['units'] = 'kpc, km/s, Msun; actions kpc*km/s; E (km/s)^2'
    if test:
        grp.create_dataset('star_index', data=idx)
    for k, v in out.items():
        grp.create_dataset(k, data=v.astype(np.float32), compression='gzip',
                           compression_opts=4, shuffle=True)
    rc = grp.create_group('rotation_curve')
    rc.create_dataset('R', data=Rgrid)
    rc.create_dataset('vc', data=vc)

    # ---- 简单诊断 ----
    R = np.hypot(xv[:, 0], xv[:, 1])
    disk = bound & (R > 4) & (R < 12) & (np.abs(xv[:, 2]) < 1)
    def med(a): return np.nanmedian(a[disk]) if disk.sum() else np.nan
    print(f"    v_c(8 kpc) = {np.interp(8, Rgrid, vc):.1f} km/s,  v_c 最大值 = {vc.max():.1f} km/s")
    print(f"    盘内恒星 (4<R<12, |z|<1): {disk.sum():,} 颗; 中位数: "
          f"J_z = {med(out['Jz']):.1f}, J_phi = {med(out['Jphi']):.0f}, e = {med(out['ecc']):.3f}, "
          f"Rg - R = {med(out['Rg'] - R):+.2f} kpc")
    nan_e = np.sum(bound & ~np.isfinite(out['ecc']))
    if nan_e:
        print(f"    注意: {nan_e} 颗束缚恒星没有算出偏心率")
    el = time.time() - t_start
    print(f"    星系 {gid} 完成, 用时 {el/60:.1f} 分钟")
    return el, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--test', action='store_true', help='只算第 1 个 (或 --ids 指定的) 星系的 5000 颗随机恒星')
    ap.add_argument('--ids', nargs='*', help='只算这些星系')
    args = ap.parse_args()

    if not os.path.exists(IN_FILE):
        sys.exit(f"找不到 {IN_FILE}, 请在 casuality 文件夹里运行")
    os.makedirs(OUT_DIR, exist_ok=True)
    out_file = os.path.join(OUT_DIR, 'orbits_test.hdf5' if args.test else 'orbits_all.hdf5')
    rng = np.random.default_rng(42)

    with h5py.File(IN_FILE, 'r') as fin:
        ids = sorted(fin.keys())
        if args.ids:
            ids = [i for i in ids if i in args.ids]
        if args.test:
            ids = ids[:1]
            if os.path.exists(out_file):
                os.remove(out_file)
        print(f"输出: {out_file}\n星系: {ids}")

        for gid in ids:
            with h5py.File(out_file, 'a') as fout:
                if gid in fout:
                    if fout[gid].attrs.get('pot_version', '') == POT_VERSION:
                        print(f"\n星系 {gid} 已算过, 跳过")
                        continue
                    print(f"\n星系 {gid} 是用旧版势设置算的, 删除后重算")
                    del fout[gid]
            # 计算全部完成后才写入该星系的组; 出错时删除, 不留半成品
            with h5py.File(out_file, 'a') as fout:
                try:
                    el, n = process_galaxy(fin, fout, gid, test=args.test, rng=rng)
                except BaseException:
                    if gid in fout:
                        del fout[gid]
                    raise
            if args.test:
                print(f"\n[测试] {n} 颗星用时 {el:.0f} 秒。按此速度, 全部恒星 (约 60–100 万颗/星系) "
                      f"每个星系约需 {el/n*8e5/60:.0f} 分钟 (粗略估计, 势的构建只算一次)。")
    print("\n全部完成。")


if __name__ == '__main__':
    main()
