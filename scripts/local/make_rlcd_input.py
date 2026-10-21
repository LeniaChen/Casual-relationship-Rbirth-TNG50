# -*- coding: utf-8 -*-
"""
make_rlcd_input.py —— 选星, 并生成 RLCD 的输入 csv 和对应的真值 csv (第 1 阶段 1.6 步)

运行 (在 casuality 文件夹里; 有 h5py 用 h5py, 没有就用 pyfive; 需要 numpy、pandas):
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/make_rlcd_input.py

输入:  data/tng50_extract/data_obtain_all.hdf5   (位置、年龄、丰度、出生半径)
       data/tng50_derived/orbits_all.hdf5        (作用量、偏心率、Rg)
输出:  data/rlcd_input/nonoise/<ID>_e.csv        RLCD 输入, 6 列, 格式与作者仓库 data/g2.79e12_e.csv 相同:
                                                  logJZ, JPhi, ecc, age, FeH, MgH
       data/rlcd_input/nonoise/<ID>_truth.csv    与 _e.csv 逐行对应的真值和辅助量
       data/rlcd_input/noise/<ID>_e.csv          加了观测误差的版本 (年龄 10%, [Fe/H] 0.02 dex, [Mg/Fe] 0.06 dex)
       data/rlcd_input/noise/<ID>_truth.csv      真值 (同时保存加误差前后的 age / FeH / MgFe / MgH)
       results/rlcd_input_selection.txt          每个星系每一步选星后剩多少颗星

选星条件 (2026-09-26 与用户确定):
  1. 束缚 (E < 0), 且 J_R、J_z、J_phi、Rg、偏心率全部有值
  2. 有出生半径 (has_birth), 且出生半径 R_bir < 20 kpc (同作者仓库数据)
  3. 当前柱半径 7 < R < 10 kpc              (论文正文的范围)
  4. 离盘面的高度 |z| < 0.3 kpc
  5. [Fe/H] > -0.5                          (同作者仓库数据; 不切 [Mg/Fe], 因为多数 TNG50 星系没有清楚的高/低 α 两条序列)
  6. J_z > 0                                (要取对数)
  逆向转动的星 (J_phi < 0) 保留, 与作者仓库数据一致。
  加误差的版本: 在同一批星上加误差, 选星用的是不加误差的值 (所以两个版本逐行对应同一批星)。

各列的含义和单位:
  logJZ = log10(J_z / (kpc km/s));  JPhi = J_phi = L_z [kpc km/s];  ecc = 偏心率 (三维距离的近心距/远心距);
  age [Gyr];  FeH = [Fe/H];  MgH = [Mg/H] = [Mg/Fe] + [Fe/H]
  truth: iord (TNG 粒子 ID), R, z [kpc, z=0 盘坐标系], Rg [kpc], R_bir [kpc, 出生时刻盘坐标系的柱半径],
         R_bir_3D, z_bir [kpc], JR, JZ, JPhi, rperi, rapo, zmax, E, star_index (在 data_obtain_all.hdf5 /stars 里的行号)
"""
import os
import numpy as np
import pandas as pd
try:
    import h5py
except ImportError:
    import pyfive as h5py

IN_PART = 'data/tng50_extract/data_obtain_all.hdf5'
IN_ORB = 'data/tng50_derived/orbits_all.hdf5'
OUT_DIR = 'data/rlcd_input'
REPORT = 'results/rlcd_input_selection.txt'

# ---- 选星条件 ----
R_MIN, R_MAX = 7.0, 10.0
Z_MAX = 0.3
FEH_MIN = -0.5
RBIR_MAX = 20.0

# ---- 观测误差 (论文 2.1 节) ----
AGE_FRAC_ERR = 0.10     # 年龄: 相对误差 10%
FEH_ERR = 0.02          # [Fe/H]: 0.02 dex
MGFE_ERR = 0.06         # [Mg/Fe]: 0.06 dex
NOISE_SEED = 20260926   # 固定随机种子, 每个星系再加上它的 ID, 保证可重复

lines = []


def log(s=''):
    print(s)
    lines.append(s)


os.makedirs(os.path.join(OUT_DIR, 'nonoise'), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'noise'), exist_ok=True)
os.makedirs(os.path.dirname(REPORT), exist_ok=True)

fp = h5py.File(IN_PART, 'r')
fo = h5py.File(IN_ORB, 'r')
ids = sorted(fo.keys())

log(f"选星条件: 束缚且轨道量齐全; has_birth 且 R_bir < {RBIR_MAX}; {R_MIN} < R < {R_MAX} kpc; |z| < {Z_MAX} kpc; "
    f"[Fe/H] > {FEH_MIN}; J_z > 0; 不切 [Mg/Fe]")
log(f"观测误差 (noise 版本): 年龄 {AGE_FRAC_ERR*100:.0f}%, [Fe/H] {FEH_ERR} dex, [Mg/Fe] {MGFE_ERR} dex, 种子 {NOISE_SEED}+ID")
log('')
hdr = f"{'ID':>7} {'全部':>9} {'轨道齐全':>9} {'有出生R':>9} {'R范围':>8} {'|z|':>7} {'[Fe/H]':>7} {'Jz>0':>7}"
log(hdr)

totals = []
for gid in ids:
    s, o = fp[gid]['stars'], fo[gid]
    pos = s['pos'][:].astype(np.float64)
    R = np.hypot(pos[:, 0], pos[:, 1])
    z = pos[:, 2]
    age = s['age'][:].astype(np.float64)
    feh = s['FeH'][:].astype(np.float64)
    mgfe = s['MgFe'][:].astype(np.float64)
    rbir = s['birth_R_2D'][:].astype(np.float64)
    rbir3 = s['birth_R_3D'][:].astype(np.float64)
    zbir = s['birth_pos'][:][:, 2].astype(np.float64)
    has = s['has_birth'][:].astype(bool)
    iord = s['iord'][:]
    q = {k: o[k][:].astype(np.float64) for k in ['Jr', 'Jz', 'Jphi', 'E', 'Rg', 'rperi', 'rapo', 'ecc', 'zmax']}
    n = len(R)

    steps = []
    m = np.isfinite(q['E']) & (q['E'] < 0)
    for k in ['Jr', 'Jz', 'Jphi', 'Rg', 'ecc']:
        m &= np.isfinite(q[k])
    steps.append(m.sum())
    m &= has & np.isfinite(rbir) & (rbir < RBIR_MAX)
    steps.append(m.sum())
    m &= (R > R_MIN) & (R < R_MAX)
    steps.append(m.sum())
    m &= np.abs(z) < Z_MAX
    steps.append(m.sum())
    m &= np.isfinite(feh) & (feh > FEH_MIN) & np.isfinite(mgfe)
    steps.append(m.sum())
    m &= q['Jz'] > 0
    steps.append(m.sum())
    log(f"{gid:>7} {n:>9,} " + ' '.join(f"{v:>9,}" if i < 2 else f"{v:>8,}" if i == 2 else f"{v:>7,}" for i, v in enumerate(steps)))
    totals.append(steps[-1])

    idx = np.where(m)[0]
    mgh = mgfe[idx] + feh[idx]
    e_df = pd.DataFrame({
        'logJZ': np.log10(q['Jz'][idx]),
        'JPhi': q['Jphi'][idx],
        'ecc': q['ecc'][idx],
        'age': age[idx],
        'FeH': feh[idx],
        'MgH': mgh,
    })
    truth = pd.DataFrame({
        'iord': iord[idx].astype(np.uint64),
        'star_index': idx,
        'R': R[idx], 'z': z[idx],
        'Rg': q['Rg'][idx],
        'R_bir': rbir[idx], 'R_bir_3D': rbir3[idx], 'z_bir': zbir[idx],
        'age': age[idx], 'FeH': feh[idx], 'MgFe': mgfe[idx], 'MgH': mgh,
        'JR': q['Jr'][idx], 'JZ': q['Jz'][idx], 'JPhi': q['Jphi'][idx],
        'ecc': q['ecc'][idx], 'rperi': q['rperi'][idx], 'rapo': q['rapo'][idx],
        'zmax': q['zmax'][idx], 'E': q['E'][idx],
    })
    e_df.to_csv(os.path.join(OUT_DIR, 'nonoise', f'{gid}_e.csv'), index=False)
    truth.to_csv(os.path.join(OUT_DIR, 'nonoise', f'{gid}_truth.csv'), index=False)

    # ---- 加观测误差的版本 ----
    rng = np.random.default_rng(NOISE_SEED + int(gid))
    age_n = age[idx] * (1 + AGE_FRAC_ERR * rng.standard_normal(len(idx)))
    feh_n = feh[idx] + FEH_ERR * rng.standard_normal(len(idx))
    mgfe_n = mgfe[idx] + MGFE_ERR * rng.standard_normal(len(idx))
    e_n = e_df.copy()
    e_n['age'] = age_n
    e_n['FeH'] = feh_n
    e_n['MgH'] = mgfe_n + feh_n
    truth_n = truth.copy()
    truth_n['age_obs'] = age_n
    truth_n['FeH_obs'] = feh_n
    truth_n['MgFe_obs'] = mgfe_n
    truth_n['MgH_obs'] = mgfe_n + feh_n
    e_n.to_csv(os.path.join(OUT_DIR, 'noise', f'{gid}_e.csv'), index=False)
    truth_n.to_csv(os.path.join(OUT_DIR, 'noise', f'{gid}_truth.csv'), index=False)

fp.close()
fo.close()

log('')
log(f"最终样本: 每个星系 {min(totals):,}–{max(totals):,} 颗, 合计 {sum(totals):,} 颗")
log(f"输出: {OUT_DIR}/nonoise/<ID>_e.csv, <ID>_truth.csv;  {OUT_DIR}/noise/<ID>_e.csv, <ID>_truth.csv")
with open(REPORT, 'w') as f:
    f.write('\n'.join(lines) + '\n')
