# -*- coding: utf-8 -*-
"""
make_rlcd_input.py — select stars and write the RLCD input csv and the matching truth csv (stage 1, step 1.6)

Run (from the repo root; uses h5py, or pyfive if h5py is missing; needs numpy, pandas):
    $PY scripts/local/make_rlcd_input.py

Inputs:  data/tng50_extract/data_obtain_all.hdf5   (positions, ages, abundances, birth radii)
         data/tng50_derived/orbits_all.hdf5        (actions, eccentricity, Rg)
Outputs: data/rlcd_input/nonoise/<ID>_e.csv        RLCD input, 6 columns, same format as data/g2.79e12_e.csv of the authors' repo:
                                                  logJZ, JPhi, ecc, age, FeH, MgH
         data/rlcd_input/nonoise/<ID>_truth.csv    truth and auxiliary quantities, row-matched to _e.csv
         data/rlcd_input/noise/<ID>_e.csv          version with observational errors (age 10%, [Fe/H] 0.02 dex, [Mg/Fe] 0.06 dex)
         data/rlcd_input/noise/<ID>_truth.csv      truth (also stores age / FeH / MgFe / MgH before and after adding errors)
         results/rlcd_input_selection.txt          number of stars left after each selection step, per galaxy

Selection (fixed 2026-09-26):
  1. bound (E < 0), with J_R, J_z, J_phi, Rg and eccentricity all defined
  2. has a birth radius (has_birth) and R_bir < 20 kpc (as in the authors' data)
  3. current cylindrical radius 7 < R < 10 kpc   (range in the paper)
  4. height above the disk plane |z| < 0.3 kpc
  5. [Fe/H] > -0.5                               (as in the authors' data; no [Mg/Fe] cut, since most TNG50 galaxies show no clear high-/low-alpha sequences)
  6. J_z > 0                                     (logarithm is taken)
  Counter-rotating stars (J_phi < 0) are kept, as in the authors' data.
  Noisy version: errors are added to the same stars; selection uses the error-free values (both versions are row-matched).

Column meanings and units:
  logJZ = log10(J_z / (kpc km/s));  JPhi = J_phi = L_z [kpc km/s];  ecc = eccentricity (from 3D peri-/apocentre distances);
  age [Gyr];  FeH = [Fe/H];  MgH = [Mg/H] = [Mg/Fe] + [Fe/H]
  truth: iord (TNG particle ID), R, z [kpc, z=0 disk frame], Rg [kpc], R_bir [kpc, cylindrical radius in the disk frame at birth],
         R_bir_3D, z_bir [kpc], JR, JZ, JPhi, rperi, rapo, zmax, E, star_index (row in data_obtain_all.hdf5 /stars)
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

# ---- selection ----
R_MIN, R_MAX = 7.0, 10.0
Z_MAX = 0.3
FEH_MIN = -0.5
RBIR_MAX = 20.0

# ---- observational errors (paper Sec. 2.1) ----
AGE_FRAC_ERR = 0.10     # age: 10% relative error
FEH_ERR = 0.02          # [Fe/H]: 0.02 dex
MGFE_ERR = 0.06         # [Mg/Fe]: 0.06 dex
NOISE_SEED = 20260926   # fixed random seed; the galaxy ID is added to it, for reproducibility

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
    e_df.to_csv(os.path.join(OUT_DIR, 'nonoise', f'{gid}_e.csv'), index=False)
    truth.to_csv(os.path.join(OUT_DIR, 'nonoise', f'{gid}_truth.csv'), index=False)

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
    e_n.to_csv(os.path.join(OUT_DIR, 'noise', f'{gid}_e.csv'), index=False)
    truth_n.to_csv(os.path.join(OUT_DIR, 'noise', f'{gid}_truth.csv'), index=False)

fp.close()
fo.close()

log('')
log(f"final sample: {min(totals):,}–{max(totals):,} stars per galaxy, {sum(totals):,} in total")
log(f"outputs: {OUT_DIR}/nonoise/<ID>_e.csv, <ID>_truth.csv;  {OUT_DIR}/noise/<ID>_e.csv, <ID>_truth.csv")
with open(REPORT, 'w') as f:
    f.write('\n'.join(lines) + '\n')
