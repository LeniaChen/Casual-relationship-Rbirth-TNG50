# -*- coding: utf-8 -*-
"""
make_input.py — RLCD inputs and truth for the selection R_0_30_Z_30 (0 < R < 30 kpc, |z| < 30 kpc) (2026-09-27)

Identical to scripts/local/make_rlcd_input.py except for two changes:
  1. selection range: 0 < R < 30 kpc, |z| < 30 kpc (instead of 7 < R < 10, |z| < 0.3). Other cuts unchanged:
     bound with all orbit quantities; birth radius with R_bir < 20 kpc; [Fe/H] > -0.5; J_z > 0; no [Mg/Fe] cut; counter-rotating stars kept.
  2. output location (all inside this selection's folder star_condition/R_0_30_Z_30/):
       input_data/nonoise/<ID>_e.csv, input_data/noise/<ID>_e.csv      RLCD input, 6 columns logJZ, JPhi, ecc, age, FeH, MgH
       true_value/nonoise/<ID>_truth.csv, true_value/noise/<ID>_truth.csv   truth row-matched to _e.csv (R_bir, Rg, ...)
       input_data/selection.txt                                          number of stars left after each selection step, per galaxy

Run (uses h5py, or pyfive if h5py is missing; needs numpy, pandas):
    cd <repo root>
    $PY star_condition/R_0_30_Z_30/input_data/make_input.py            # 12 galaxies
    $PY star_condition/R_0_30_Z_30/input_data/make_input.py 538905     # selected galaxies only

Inputs: data/tng50_extract/data_obtain_all.hdf5, data/tng50_derived/orbits_all.hdf5 (as before)
Observational errors (noise version) as before: age 10%, [Fe/H] 0.02 dex, [Mg/Fe] 0.06 dex, random seed 20260926 + galaxy ID.
Column meanings and units: see scripts/local/make_rlcd_input.py.
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
ROOT = HERE.parents[3]                         # repo root
IN_PART = str(ROOT / 'data' / 'tng50_extract' / 'data_obtain_all.hdf5')
IN_ORB = str(ROOT / 'data' / 'tng50_derived' / 'orbits_all.hdf5')
OUT_E = COND / 'input_data'                    # <version>/<ID>_e.csv
OUT_T = COND / 'true_value'                    # <version>/<ID>_truth.csv
REPORT = str(COND / 'input_data' / 'selection.txt')

# ---- selection ----
R_MIN, R_MAX = 0.0, 30.0
Z_MAX = 30.0
FEH_MIN = -0.5
RBIR_MAX = 20.0

# ---- observational errors ----
AGE_FRAC_ERR = 0.10     # age: 10% relative error
FEH_ERR = 0.02          # [Fe/H]: 0.02 dex
MGFE_ERR = 0.06         # [Mg/Fe]: 0.06 dex
NOISE_SEED = 20260926   # fixed random seed; the galaxy ID is added to it, for reproducibility

lines = []


def log(s=''):
    print(s)
    lines.append(s)


for d in [OUT_E / 'nonoise', OUT_E / 'noise', OUT_T / 'nonoise', OUT_T / 'noise']:
    os.makedirs(d, exist_ok=True)

fp = h5py.File(IN_PART, 'r')
fo = h5py.File(IN_ORB, 'r')
ids = sys.argv[1:] or sorted(fo.keys())
append = len(sys.argv) > 1 and os.path.exists(REPORT)   # when run in batches, append to selection.txt

log(f"selection: bound with all orbit quantities; has_birth and R_bir < {RBIR_MAX}; {R_MIN} < R < {R_MAX} kpc; |z| < {Z_MAX} kpc; "
    f"[Fe/H] > {FEH_MIN}; J_z > 0; no [Mg/Fe] cut")
log(f"observational errors (noise version): age {AGE_FRAC_ERR*100:.0f}%, [Fe/H] {FEH_ERR} dex, [Mg/Fe] {MGFE_ERR} dex, seed {NOISE_SEED}+ID")
log('')
hdr = f"{'ID':>7} {'all':>9} {'orbits ok':>9} {'has Rbir':>9} {'R range':>8} {'|z|':>7} {'[Fe/H]':>7} {'Jz>0':>7}"
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

    # ---- version with observational errors ----
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
log(f"final sample: {min(totals):,}–{max(totals):,} stars per galaxy, {sum(totals):,} in total")
log(f"outputs: input_data/<version>/<ID>_e.csv;  true_value/<version>/<ID>_truth.csv  (version = nonoise, noise)")
with open(REPORT, 'a' if append else 'w') as f:
    f.write('\n'.join(lines if not append else [l for l in lines if l[:7].strip().isdigit()]) + '\n')
