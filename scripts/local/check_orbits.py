# -*- coding: utf-8 -*-
"""
check_orbits.py —— 检查 compute_orbits.py 对 12 个星系的计算结果 (第 1 阶段 1.5 步)

运行 (在 casuality 文件夹里; 有 h5py 用 h5py, 没有就用 pyfive):
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/check_orbits.py

输入:  data/tng50_extract/data_obtain_all.hdf5      (恒星位置、速度、年龄等)
       data/tng50_derived/orbits_all.hdf5            (作用量、Rg、偏心率等)
输出:  results/orbit_check/summary.txt               每个星系的检查结果 (文字)
       results/orbit_check/summary.csv               同上 (表格)
       results/orbit_check/rotation_curves.png       12 个星系的圆周速度曲线
       results/orbit_check/ecc_disk.png              盘内恒星的偏心率分布
       results/orbit_check/Rg_vs_R.png               导引半径 Rg 与当前半径 R
       results/orbit_check/logJz_disk.png            盘内恒星的 log10(J_z) 分布

检查项 (每一项都对应一个"算错了就会露馅"的物理关系):
  1. 束缚恒星比例、各量的 NaN 比例 (束缚恒星中应为 0)
  2. J_phi 是否等于直接算的 L_z = x*vy - y*vx
  3. 当前三维距离 r 是否落在 [r_peri, r_apo] 之间;  当前 |z| 是否 <= z_max
  4. Rg * v_c(Rg) / |L_z| 是否 ≈ 1 (导引半径的定义)
  5. 偏心率是否在 [0, 1]
  6. 盘内 (4<R<12 kpc, |z|<1 kpc) 恒星的中位数: J_z, J_phi, J_R, e, Rg-R, 逆向转动比例
  7. 引力势的设置: pot_version、pot_rmin_kpc
"""
import os
import numpy as np
try:
    import h5py
except ImportError:
    import pyfive as h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_PART = 'data/tng50_extract/data_obtain_all.hdf5'
IN_ORB = 'data/tng50_derived/orbits_all.hdf5'
OUT = 'results/orbit_check'
os.makedirs(OUT, exist_ok=True)

INK = '#1f2937'
MUTED = '#9ca3af'
BLUE = '#2563eb'

lines, rows, curves, eccs, rgr, ljz = [], [], {}, {}, {}, {}


def log(s=''):
    print(s)
    lines.append(s)


fp = h5py.File(IN_PART, 'r')
fo = h5py.File(IN_ORB, 'r')
ids = sorted(fo.keys())
log(f"orbits_all.hdf5 中的星系: {len(ids)} 个: {', '.join(ids)}")
all_ok = True

for gid in ids:
    o, s = fo[gid], fp[gid]['stars']
    log('')
    log('=' * 90)
    at = dict(o.attrs)
    ver = at.get('pot_version', b'?')
    ver = ver.decode() if isinstance(ver, bytes) else ver
    log(f"星系 {gid}   pot_version={ver}   pot_rmin={at.get('pot_rmin_kpc')} kpc   "
        f"test_mode={at.get('test_mode')}")

    pos = s['pos'][:].astype(float)
    vel = s['vel'][:].astype(float)
    age = s['age'][:].astype(float)
    q = {k: o[k][:].astype(float) for k in ['Jr', 'Jz', 'Jphi', 'E', 'Rg', 'rperi', 'rapo', 'ecc', 'zmax']}
    n = len(pos)
    if len(q['E']) != n:
        log(f"  !! 长度不一致: 恒星 {n}, 轨道结果 {len(q['E'])}")
        all_ok = False
        continue

    R = np.hypot(pos[:, 0], pos[:, 1])
    r = np.linalg.norm(pos, axis=1)
    z = pos[:, 2]
    Lz = pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]
    bound = np.isfinite(q['E']) & (q['E'] < 0)
    nb = bound.sum()

    # 1. NaN
    nan_counts = {k: int(np.sum(bound & ~np.isfinite(v))) for k, v in q.items()}
    bad_nan = {k: v for k, v in nan_counts.items() if v > 0}
    log(f"  恒星 {n:,}, 束缚 {nb:,} ({nb/n*100:.2f}%);  束缚恒星中 NaN: "
        f"{bad_nan if bad_nan else '无'}")
    ok = bound & np.all([np.isfinite(v) for v in q.values()], axis=0)

    # 2. Jphi == Lz
    dL = np.nanmax(np.abs(q['Jphi'] - Lz))
    # 3. r in [rperi, rapo], |z| <= zmax
    f_r = np.mean((q['rperi'][ok] <= r[ok] * 1.001 + 1e-3) & (r[ok] <= q['rapo'][ok] * 1.001 + 1e-3))
    f_z = np.mean(q['zmax'][ok] >= np.abs(z[ok]) * 0.999 - 1e-3)
    # 4. Rg * vc(Rg) / |Lz|
    Rc, vc = o['rotation_curve']['R'][:], o['rotation_curve']['vc'][:]
    m4 = ok & (q['Rg'] > 1) & (q['Rg'] < 29) & (np.abs(Lz) > 10)
    ratio = q['Rg'][m4] * np.interp(q['Rg'][m4], Rc, vc) / np.abs(Lz[m4])
    # 5. ecc range
    emin, emax = np.nanmin(q['ecc'][ok]), np.nanmax(q['ecc'][ok])
    log(f"  |J_phi - L_z| 最大 {dL:.3g};  r∈[r_peri,r_apo] {f_r*100:.3f}%;  |z|≤z_max {f_z*100:.3f}%;  "
        f"Rg·v_c/|L_z| 中位数 {np.median(ratio):.4f} (5–95%: {np.percentile(ratio,5):.3f}–{np.percentile(ratio,95):.3f});  "
        f"e ∈ [{emin:.4f}, {emax:.4f}]")
    checks_ok = (not bad_nan) and dL < 0.05 and f_r > 0.999 and f_z > 0.999 \
        and abs(np.median(ratio) - 1) < 0.01 and emin >= 0 and emax <= 1
    if not checks_ok:
        log("  !! 有检查项没通过")
        all_ok = False

    # 6. disk medians
    disk = ok & (R > 4) & (R < 12) & (np.abs(z) < 1)
    med = lambda a: float(np.median(a[disk]))
    v8 = float(np.interp(8, Rc, vc))
    row = dict(ID=gid, n_stars=n, bound_frac=nb / n, pot_rmin=at.get('pot_rmin_kpc'),
               vc8=v8, vcmax=float(vc.max()), n_disk=int(disk.sum()),
               Jz_med=med(q['Jz']), Jphi_med=med(q['Jphi']), Jr_med=med(q['Jr']),
               ecc_med=med(q['ecc']), RgR_med=med(q['Rg'] - R),
               counterrot_frac=float(np.mean(Lz[disk] < 0)),
               ecc_young_med=float(np.median(q['ecc'][disk & (age < 3)])) if np.any(disk & (age < 3)) else np.nan,
               checks_ok=bool(checks_ok))
    rows.append(row)
    log(f"  v_c(8 kpc) = {v8:.1f} km/s, v_c 最大 {vc.max():.1f};  盘内恒星 {disk.sum():,}: 中位数 "
        f"J_z {row['Jz_med']:.1f}, J_phi {row['Jphi_med']:.0f}, J_R {row['Jr_med']:.1f}, e {row['ecc_med']:.3f} "
        f"(年轻星<3 Gyr: {row['ecc_young_med']:.3f}), Rg−R {row['RgR_med']:+.2f} kpc, 逆向转动 {row['counterrot_frac']*100:.1f}%")

    curves[gid] = (Rc, vc)
    eccs[gid] = q['ecc'][disk]
    sub = ok & (R < 25) & (np.abs(z) < 3)
    rgr[gid] = (R[sub], q['Rg'][sub])
    ljz[gid] = np.log10(np.clip(q['Jz'][disk], 1e-3, None))

fp.close()
fo.close()

log('')
log('=' * 90)
log('结论: 全部检查通过' if all_ok else '结论: 有星系没通过检查, 见上面标 !! 的行')

with open(os.path.join(OUT, 'summary.txt'), 'w') as f:
    f.write('\n'.join(lines) + '\n')
keys = list(rows[0].keys())
with open(os.path.join(OUT, 'summary.csv'), 'w') as f:
    f.write(','.join(keys) + '\n')
    for rw in rows:
        f.write(','.join(f"{rw[k]:.6g}" if isinstance(rw[k], float) else str(rw[k]) for k in keys) + '\n')


# ---------------- 图 (图中文字用英文, 避免系统缺中文字体时显示成方块) ----------------
def style(ax):
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=8)
    ax.grid(color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)


def grid_fig(title):
    fig, axes = plt.subplots(3, 4, figsize=(13, 8.5), sharex=True, sharey=True)
    fig.suptitle(title, color=INK, fontsize=13)
    return fig, axes.ravel()


# 1. rotation curves: one panel, every galaxy in muted gray, median in blue
fig, ax = plt.subplots(figsize=(7, 4.5))
style(ax)
allv = []
for gid, (Rc, vc) in curves.items():
    ax.plot(Rc, vc, color=MUTED, lw=1)
    allv.append(vc)
ax.plot(Rc, np.median(allv, axis=0), color=BLUE, lw=2, label='median of 12 galaxies')
ax.set_xlabel('R [kpc]', color=INK)
ax.set_ylabel('v_c [km/s]', color=INK)
ax.set_title('Circular velocity curves (gray: each galaxy)', color=INK, fontsize=11)
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'rotation_curves.png'), dpi=130)
plt.close(fig)

# 2. ecc distributions (small multiples)
fig, axs = grid_fig('Eccentricity of disk stars (4<R<12 kpc, |z|<1 kpc)')
bins = np.linspace(0, 1, 41)
for ax, gid in zip(axs, eccs):
    style(ax)
    ax.hist(eccs[gid], bins=bins, color=BLUE, alpha=0.85, edgecolor='white', linewidth=0.4)
    ax.axvline(np.median(eccs[gid]), color=INK, lw=1, ls='--')
    ax.set_title(f"{gid}  median {np.median(eccs[gid]):.2f}", fontsize=9, color=INK)
for ax in axs[8:]:
    ax.set_xlabel('e', color=INK)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'ecc_disk.png'), dpi=120)
plt.close(fig)

# 3. Rg vs R
fig, axs = grid_fig('Guiding radius Rg vs current R (R<25, |z|<3 kpc; color = log star count)')
for ax, gid in zip(axs, rgr):
    style(ax)
    Rr, Rg = rgr[gid]
    ax.hexbin(Rr, Rg, gridsize=50, extent=(0, 25, 0, 25), bins='log', cmap='Blues', mincnt=1)
    ax.plot([0, 25], [0, 25], color=INK, lw=0.8)
    ax.set_title(gid, fontsize=9, color=INK)
for ax in axs[8:]:
    ax.set_xlabel('R [kpc]', color=INK)
for ax in axs[::4]:
    ax.set_ylabel('Rg [kpc]', color=INK)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'Rg_vs_R.png'), dpi=120)
plt.close(fig)

# 4. log10 Jz
fig, axs = grid_fig('log10 J_z [kpc km/s] of disk stars (4<R<12 kpc, |z|<1 kpc)')
bins = np.linspace(-2, 3.5, 45)
for ax, gid in zip(axs, ljz):
    style(ax)
    ax.hist(ljz[gid], bins=bins, color=BLUE, alpha=0.85, edgecolor='white', linewidth=0.4)
    ax.set_title(f"{gid}  median {np.median(ljz[gid]):.2f}", fontsize=9, color=INK)
for ax in axs[8:]:
    ax.set_xlabel('log10 J_z', color=INK)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'logJz_disk.png'), dpi=120)
plt.close(fig)
print(f"\n输出在 {OUT}/")
