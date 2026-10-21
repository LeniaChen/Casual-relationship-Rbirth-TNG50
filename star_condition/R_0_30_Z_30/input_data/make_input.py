# -*- coding: utf-8 -*-
"""
make_input.py —— 选星条件 R_0_30_Z_30 (0 < R < 30 kpc, |z| < 30 kpc) 下, 生成 RLCD 的输入和真值 (2026-09-27, 用户指定)

与 scripts/local/make_rlcd_input.py 完全相同, 只改了两处:
  1. 选星范围: 0 < R < 30 kpc, |z| < 30 kpc (原来是 7 < R < 10、|z| < 0.3)。其余条件不变:
     束缚且轨道量齐全; 有出生半径且 R_bir < 20 kpc; [Fe/H] > -0.5; J_z > 0; 不切 [Mg/Fe]; 逆向转动的星保留。
  2. 输出位置 (都在本条件的文件夹 star_condition/R_0_30_Z_30/ 里):
       input_data/nonoise/<ID>_e.csv, input_data/noise/<ID>_e.csv      RLCD 输入, 6 列 logJZ, JPhi, ecc, age, FeH, MgH
       true_value/nonoise/<ID>_truth.csv, true_value/noise/<ID>_truth.csv   与 _e.csv 逐行对应的真值 (R_bir, Rg 等)
       input_data/selection.txt                                          每个星系每一步选星后剩多少颗星

运行 (有 h5py 用 h5py, 没有就用 pyfive; 需要 numpy、pandas):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_30/input_data/make_input.py            # 12 个星系
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_30/input_data/make_input.py 538905     # 只做指定星系

输入: data/tng50_extract/data_obtain_all.hdf5, data/tng50_derived/orbits_all.hdf5 (与原来相同)
观测误差 (noise 版本) 与原来相同: 年龄 10%, [Fe/H] 0.02 dex, [Mg/Fe] 0.06 dex, 随机种子 20260926 + 星系 ID。
各列的含义和单位见 scripts/local/make_rlcd_input.py。
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
try:
    import h5py
except ImportError:
    import pyfive as h5py

HERE = Path(__file__).resolve()
COND = HERE.parents[1]                         # star_condition/R_0_30_Z_30/
ROOT = HERE.parents[3]                         # casuality/
IN_PART = str(ROOT / 'data' / 'tng50_extract' / 'data_obtain_all.hdf5')
IN_ORB = str(ROOT / 'data' / 'tng50_derived' / 'orbits_all.hdf5')
OUT_E = COND / 'input_data'                    # <版本>/<ID>_e.csv
OUT_T = COND / 'true_value'                    # <版本>/<ID>_truth.csv
REPORT = str(COND / 'input_data' / 'selection.txt')

# ---- 选星条件 ----
R_MIN, R_MAX = 0.0, 30.0
Z_MAX = 30.0
FEH_MIN = -0.5
RBIR_MAX = 20.0

# ---- 观测误差 ----
AGE_FRAC_ERR = 0.10     # 年龄: 相对误差 10%
FEH_ERR = 0.02          # [Fe/H]: 0.02 dex
MGFE_ERR = 0.06         # [Mg/Fe]: 0.06 dex
NOISE_SEED = 20260926   # 固定随机种子, 每个星系再加上它的 ID, 保证可重复

lines = []


def log(s=''):
    print(s)
    lines.append(s)


for d in [OUT_E / 'nonoise', OUT_E / 'noise', OUT_T / 'nonoise', OUT_T / 'noise']:
    os.makedirs(d, exist_ok=True)

fp = h5py.File(IN_PART, 'r')
fo = h5py.File(IN_ORB, 'r')
ids = sys.argv[1:] or sorted(fo.keys())
append = len(sys.argv) > 1 and os.path.exists(REPORT)   # 分批运行时把结果追加到 selection.txt

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
    e_df.to_csv(OUT_E / 'nonoise' / f'{gid}_e.csv', index=False)
    truth.to_csv(OUT_T / 'nonoise' / f'{gid}_truth.csv', index=False)

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
    e_n.to_csv(OUT_E / 'noise' / f'{gid}_e.csv', index=False)
    truth_n.to_csv(OUT_T / 'noise' / f'{gid}_truth.csv', index=False)

fp.close()
fo.close()

log('')
log(f"最终样本: 每个星系 {min(totals):,}–{max(totals):,} 颗, 合计 {sum(totals):,} 颗")
log(f"输出: input_data/<版本>/<ID>_e.csv;  true_value/<版本>/<ID>_truth.csv  (版本 = nonoise, noise)")
with open(REPORT, 'a' if append else 'w') as f:
    f.write('\n'.join(lines if not append else [l for l in lines if l[:7].strip().isdigit()]) + '\n')
