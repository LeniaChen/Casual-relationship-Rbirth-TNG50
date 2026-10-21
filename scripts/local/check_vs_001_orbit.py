# -*- coding: utf-8 -*-
"""
check_vs_001_orbit.py — check that data_obtain_all.hdf5 agrees with the earlier 001_orbit results

Run from the repo root (needs numpy / pandas / h5py):
    cd <repo root>
    python3 scripts/local/check_vs_001_orbit.py

The report is printed and also written to results/check_vs_001_orbit.txt

1. Structure: every galaxy group complete, array lengths consistent, value ranges sensible
2. Row-by-row comparison with 423_sample/001_orbit/<ID>.txt (001_orbit.py, 2025-09)
    each row of that file = one real star (wind particles removed); ParticleID is the row index, same order as /stars
    Age_Gyr   <-> stars/age
    Current_R <-> |stars/pos|        (note: Current_R of the old script is a 3D distance)
    Birth_R   <-> stars/birth_R_3D   (found during the check: Birth_R of the old script is also a 3D distance)
    Mass      <-> stars/mass
3. If present, comparison with 001_orbit_s99_3D_z/<ID>.txt (adds Current_z, Birth_R_3D)

Expected differences (not errors):
  - the new script clamps the formation time of stars born "near z=0" back to the end of the trajectory (change [B]);
    their Birth_R is NaN in the old script and finite in the new one
  - float32 storage gives ~1e-6 relative differences
"""
import os
import sys
import numpy as np
import pandas as pd
try:
    import h5py
except ImportError:          # fall back to the pure-Python pyfive reader (same interface) if h5py is missing
    import pyfive as h5py

H5 = 'data/tng50_extract/data_obtain_all.hdf5'
def first_existing(*cands):
    for c in cands:
        if os.path.isdir(c):
            return c
    return cands[0]

# Earlier birth-radius results (local files of the author's earlier project, not distributed). Set their location with:
#   LEGACY_001_ORBIT_DIR      → 423_sample/001_orbit/<ID>.txt
#   LEGACY_001_ORBIT_3DZ_DIR  → 001_orbit_s99_3D_z/<ID>.txt
# This is only an internal consistency check (result in results/check_vs_001_orbit.txt); skip it if you do not have these files.
OLD_DIRS = {
    '001_orbit': os.environ.get('LEGACY_001_ORBIT_DIR', 'legacy/001_orbit'),
    '001_orbit_s99_3D_z': os.environ.get('LEGACY_001_ORBIT_3DZ_DIR', 'legacy/001_orbit_s99_3D_z'),
}
REPORT = 'results/check_vs_001_orbit.txt'

EXPECTED_GROUPS = ['header', 'frame', 'stars', 'wind', 'dm', 'gas', 'bh']
TOL = {'age': 1e-4, 'R': 1e-3, 'Rbirth': 1e-2, 'mass_rel': 1e-5}   # Gyr, kpc, kpc, relative
# birth-radius tolerance relaxed to 10 pc: runs at different times differ by a few pc (relative ~2e-4)

lines = []
def log(msg=''):
    print(msg)
    lines.append(msg)


def compare(name, new, old, tol, rel=False):
    """Return one comparison line; compare only where both sides are finite"""
    both = np.isfinite(new) & np.isfinite(old)
    only_new = np.sum(np.isfinite(new) & ~np.isfinite(old))
    only_old = np.sum(~np.isfinite(new) & np.isfinite(old))
    if both.sum() == 0:
        return f"    {name:<22} nothing to compare", False
    d = np.abs(new[both] - old[both])
    if rel:
        d = d / np.maximum(np.abs(old[both]), 1e-30)
    n_bad = int(np.sum(d > tol))
    ok = (n_bad == 0) and (only_old == 0)
    flag = 'OK ' if ok else '!! '
    return (f"  {flag}{name:<22} compared {both.sum():>8d}, max diff {d.max():.2e}, "
            f"out of tol {n_bad:>6d}, only new finite {only_new:>5d}, only old finite {only_old:>5d}"), ok


def read_old(path):
    return pd.read_csv(path, sep=r'\s+')


if not os.path.exists(H5):
    sys.exit(f"{H5} not found; run from the repo root and make sure the download is complete")

all_ok = True
with h5py.File(H5, 'r') as f:
    ids = sorted(f.keys())
    log(f"file: {H5}  ({os.path.getsize(H5)/1e9:.2f} GB), galaxies: {len(ids)}")
    log(f"galaxies: {', '.join(ids)}")

    for gid in ids:
        g = f[gid]
        log('')
        log('=' * 100)
        log(f"galaxy {gid}")

        # ---------- 1. structure ----------
        missing = [k for k in EXPECTED_GROUPS if k not in g]
        if missing:
            log(f"  !! missing groups: {missing}")
            all_ok = False
            continue
        s = g['stars']
        n = s['age'].shape[0]
        lens = {k: s[k].shape[0] for k in s.keys()}
        if len(set(lens.values())) != 1:
            log(f"  !! inconsistent array lengths in stars: {lens}")
            all_ok = False
        hdr = dict(g['header'].attrs)
        log(f"  header: a={hdr.get('a')}, h={hdr.get('h')}, t={hdr.get('time_Gyr')} Gyr")
        log(f"  particles: stars {n}, wind {g['wind/mass'].shape[0]}, dm {g['dm/mass'].shape[0]}, "
            f"gas {g['gas/mass'].shape[0]}, bh {g['bh/mass'].shape[0]}")
        if g['bh/mass'].shape[0] == 0:
            log("  !! no black-hole particle (tiny effect on the potential, but means bh was not read)")

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
        log(f"  [Fe/H]: median {np.nanmedian(feh):+.3f}, NaN {np.sum(~np.isfinite(feh))};  "
            f"[Mg/Fe]: median {np.nanmedian(mgfe):+.3f}, NaN {np.sum(~np.isfinite(mgfe))}")
        log(f"  has_birth: {has.sum()} / {n}  (NaN birth radii {np.sum(~np.isfinite(bR2))})")
        if np.sum(~has) != np.sum(~np.isfinite(bR2)):
            log("  !! number of has_birth does not match the NaNs of birth_R_2D")
            all_ok = False
        log(f"  disk stars (2<R<15, |z|<1 kpc): {disk.sum()}, median v_phi "
            f"{np.median(vphi[disk]):.1f} km/s  (should be positive, ~150–250)")
        young = (age < 0.1) & has
        if young.sum() > 10:
            dR = bR2[young] - R_cyl[young]
            log(f"  young stars (<0.1 Gyr) birth radius - current radius: median {np.median(dR):+.3f} kpc, "
                f"95% of |diff| {np.percentile(np.abs(dR), 95):.3f} kpc  (should be ~0)")

        # ---------- 2./3. comparison with the old results ----------
        for label, d in OLD_DIRS.items():
            path = os.path.join(d, f'{gid}.txt')
            if not os.path.exists(path):
                log(f"  [{label}] no {gid}.txt, skipped")
                continue
            old = read_old(path)
            log(f"  [{label}] old file {len(old)} rows, new file {n} rows")
            if len(old) != n:
                log(f"  !! different number of rows, cannot compare row by row (the old script may have dropped NaN rows)")
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
log('Verdict: all consistent (OK)' if all_ok else 'Verdict: items marked !! need attention')

os.makedirs(os.path.dirname(REPORT), exist_ok=True)
with open(REPORT, 'w') as fo:
    fo.write('\n'.join(lines) + '\n')
print(f"\nreport written to {REPORT}")
