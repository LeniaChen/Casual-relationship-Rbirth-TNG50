# -*- coding: utf-8 -*-
"""
check_vs_001_orbit.py —— 核对 data_obtain_all.hdf5 与以前的 001_orbit 结果是否一致

在 casuality 文件夹里运行 (用你平时的 anaconda Python, 需要 numpy / pandas / h5py):
    cd /Volumes/001HardDisk/MW/casuality
    python3 scripts/local/check_vs_001_orbit.py

报告同时打印在屏幕上, 并写入 results/check_vs_001_orbit.txt

一、结构检查: 每个星系组是否齐全, 数组长度是否一致, 数值范围是否合理
二、与 423_sample/001_orbit/<ID>.txt (001_orbit.py, 2025-09) 逐行比较
    该文件每行 = 一颗真实恒星 (已去风粒子), ParticleID 是行号, 顺序与 /stars 相同
    Age_Gyr   <-> stars/age
    Current_R <-> |stars/pos|        (注意: 旧脚本的 Current_R 是三维距离)
    Birth_R   <-> stars/birth_R_3D   (核对发现: 旧脚本的 Birth_R 也是三维距离, 与 Current_R 一样)
    Mass      <-> stars/mass
三、若存在, 再与 Downloads/001_orbit_s99_3D_z/<ID>.txt 比较 (多了 Current_z, Birth_R_3D)

预期差异 (不算错误):
  - 新脚本把 "z=0 附近出生" 的恒星的出生时间夹回轨道终点 (修改 [B]),
    旧脚本里这些星的 Birth_R 是 NaN, 新脚本里有数值
  - float32 存储带来 ~1e-6 相对误差
"""
import os
import sys
import numpy as np
import pandas as pd
try:
    import h5py
except ImportError:          # 没有 h5py 时用纯 Python 的 pyfive 读取 (接口相同)
    import pyfive as h5py

H5 = 'data/tng50_extract/data_obtain_all.hdf5'
def first_existing(*cands):
    for c in cands:
        if os.path.isdir(c):
            return c
    return cands[0]

OLD_DIRS = {
    '001_orbit': first_existing('/Volumes/001HardDisk_6/TNG_server/423_sample/001_orbit',
                                os.path.expanduser('~/mnt/423_sample/001_orbit')),
    '001_orbit_s99_3D_z': first_existing(os.path.expanduser('~/Downloads/001_orbit_s99_3D_z'),
                                         os.path.expanduser('~/mnt/Downloads/001_orbit_s99_3D_z')),
}
REPORT = 'results/check_vs_001_orbit.txt'

EXPECTED_GROUPS = ['header', 'frame', 'stars', 'wind', 'dm', 'gas', 'bh']
TOL = {'age': 1e-4, 'R': 1e-3, 'Rbirth': 1e-2, 'mass_rel': 1e-5}   # Gyr, kpc, kpc, 相对
# 出生半径容差放宽到 10 pc: 不同时间运行的结果之间有 ~几 pc 的数值差 (相对差 ~2e-4)

lines = []
def log(msg=''):
    print(msg)
    lines.append(msg)


def compare(name, new, old, tol, rel=False):
    """返回一行比较结果; 只在两边都有限的位置比较"""
    both = np.isfinite(new) & np.isfinite(old)
    only_new = np.sum(np.isfinite(new) & ~np.isfinite(old))
    only_old = np.sum(~np.isfinite(new) & np.isfinite(old))
    if both.sum() == 0:
        return f"    {name:<22} 没有可比较的数据", False
    d = np.abs(new[both] - old[both])
    if rel:
        d = d / np.maximum(np.abs(old[both]), 1e-30)
    n_bad = int(np.sum(d > tol))
    ok = (n_bad == 0) and (only_old == 0)
    flag = 'OK ' if ok else '!! '
    return (f"  {flag}{name:<22} 比较 {both.sum():>8d} 个, 最大差 {d.max():.2e}, "
            f"超差 {n_bad:>6d}, 仅新有值 {only_new:>5d}, 仅旧有值 {only_old:>5d}"), ok


def read_old(path):
    return pd.read_csv(path, sep=r'\s+')


if not os.path.exists(H5):
    sys.exit(f"找不到 {H5}, 请在 casuality 文件夹里运行, 并确认下载已完成")

all_ok = True
with h5py.File(H5, 'r') as f:
    ids = sorted(f.keys())
    log(f"文件: {H5}  ({os.path.getsize(H5)/1e9:.2f} GB), 星系数: {len(ids)}")
    log(f"星系: {', '.join(ids)}")

    for gid in ids:
        g = f[gid]
        log('')
        log('=' * 100)
        log(f"星系 {gid}")

        # ---------- 一、结构 ----------
        missing = [k for k in EXPECTED_GROUPS if k not in g]
        if missing:
            log(f"  !! 缺少组: {missing}")
            all_ok = False
            continue
        s = g['stars']
        n = s['age'].shape[0]
        lens = {k: s[k].shape[0] for k in s.keys()}
        if len(set(lens.values())) != 1:
            log(f"  !! stars 各数组长度不一致: {lens}")
            all_ok = False
        hdr = dict(g['header'].attrs)
        log(f"  header: a={hdr.get('a')}, h={hdr.get('h')}, t={hdr.get('time_Gyr')} Gyr")
        log(f"  粒子数: stars {n}, wind {g['wind/mass'].shape[0]}, dm {g['dm/mass'].shape[0]}, "
            f"gas {g['gas/mass'].shape[0]}, bh {g['bh/mass'].shape[0]}")
        if g['bh/mass'].shape[0] == 0:
            log("  !! 没有黑洞粒子 (引力势里少了黑洞, 影响很小, 但说明 bh 没读到)")

        age = s['age'][:].astype(np.float64)
        pos = s['pos'][:].astype(np.float64)
        vel = s['vel'][:].astype(np.float64)
        mass = s['mass'][:].astype(np.float64)
        feh = s['FeH'][:].astype(np.float64)
        mgfe = s['MgFe'][:].astype(np.float64)
        bR2 = s['birth_R_2D'][:].astype(np.float64)
        bR3 = s['birth_R_3D'][:].astype(np.float64)
        has = s['has_birth'][:].astype(bool)

        R_cyl = np.hypot(pos[:, 0], pos[:, 1])
        r_3d = np.linalg.norm(pos, axis=1)
        vphi = (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]) / np.maximum(R_cyl, 1e-9)
        disk = (R_cyl > 2) & (R_cyl < 15) & (np.abs(pos[:, 2]) < 1)

        log(f"  age: {np.nanmin(age):.3f} – {np.nanmax(age):.3f} Gyr")
        log(f"  [Fe/H]: 中位数 {np.nanmedian(feh):+.3f}, NaN {np.sum(~np.isfinite(feh))};  "
            f"[Mg/Fe]: 中位数 {np.nanmedian(mgfe):+.3f}, NaN {np.sum(~np.isfinite(mgfe))}")
        log(f"  has_birth: {has.sum()} / {n}  (NaN 出生半径 {np.sum(~np.isfinite(bR2))})")
        if np.sum(~has) != np.sum(~np.isfinite(bR2)):
            log("  !! has_birth 与 birth_R_2D 的 NaN 个数不一致")
            all_ok = False
        log(f"  盘内恒星 (2<R<15, |z|<1 kpc): {disk.sum()}, 其 v_phi 中位数 "
            f"{np.median(vphi[disk]):.1f} km/s  (应为正且约 150–250)")
        young = (age < 0.1) & has
        if young.sum() > 10:
            dR = bR2[young] - R_cyl[young]
            log(f"  年轻星 (<0.1 Gyr) 出生半径 - 当前半径: 中位数 {np.median(dR):+.3f} kpc, "
                f"|差|95% {np.percentile(np.abs(dR), 95):.3f} kpc  (应接近 0)")

        # ---------- 二、三、与旧结果比较 ----------
        for label, d in OLD_DIRS.items():
            path = os.path.join(d, f'{gid}.txt')
            if not os.path.exists(path):
                log(f"  [{label}] 没有 {gid}.txt, 跳过")
                continue
            old = read_old(path)
            log(f"  [{label}] 旧文件 {len(old)} 行, 新文件 {n} 行")
            if len(old) != n:
                log(f"  !! 行数不同, 无法逐行比较 (旧脚本可能剔除了 NaN 行)")
                all_ok = False
                continue
            checks = [
                ('Age', age, old['Age_Gyr'].values, TOL['age'], False),
                ('Mass', mass, old['Mass'].values, TOL['mass_rel'], True),
            ]
            if label == '001_orbit':
                checks += [
                    ('Current_R(3D)', r_3d, old['Current_R'].values, TOL['R'], False),
                    ('Birth_R(=3D)', bR3, old['Birth_R'].values, TOL['Rbirth'], False),
                ]
            else:
                checks += [
                    ('Current_R(2D)', R_cyl, old['Current_R'].values, TOL['R'], False),
                    ('Current_R_3D', r_3d, old['Current_R_3D'].values, TOL['R'], False),
                    ('Current_z', pos[:, 2], old['Current_z'].values, TOL['R'], False),
                    ('Birth_R(2D)', bR2, old['Birth_R'].values, TOL['Rbirth'], False),
                    ('Birth_R_3D', bR3, old['Birth_R_3D'].values, TOL['Rbirth'], False),
                ]
            for name, a, b, tol, rel in checks:
                msg, ok = compare(name, a, b.astype(np.float64), tol, rel)
                log(msg)
                all_ok &= ok

log('')
log('=' * 100)
log('结论: 全部一致 (OK)' if all_ok else '结论: 有标 !! 的项目需要查看')

os.makedirs(os.path.dirname(REPORT), exist_ok=True)
with open(REPORT, 'w') as fo:
    fo.write('\n'.join(lines) + '\n')
print(f"\n报告已写入 {REPORT}")
