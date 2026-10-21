# -*- coding: utf-8 -*-
"""
make_input.py — RLCD inputs and truth for the selection R_0_30_Z_0.3 (0 < R < 30 kpc, |z| < 0.3 kpc) (2026-09-27)

Instead of selecting from the simulation again, the stars with |z| < 0.3 kpc are taken from the results of star_condition/R_0_30_Z_30/:
  - all other cuts of R_0_30_Z_30 are identical (0 < R < 30 kpc; bound with all orbit quantities; birth radius with R_bir < 20 kpc;
    [Fe/H] > -0.5; J_z > 0; no [Mg/Fe] cut); the only difference is the |z| limit, 30 → 0.3 kpc;
  - true_value/<version>/<ID>_truth.csv contains the current height z of every star and is row-matched to input_data/<version>/<ID>_e.csv;
  - the noisy version reuses the errors already drawn in R_0_30_Z_30, so a star has identical values in both selections.
  Rerunning make_rlcd_input.py from scratch gives the same stars, but the noisy version would draw different random errors (different random-number order).

Run (numpy, pandas):
    cd <repo root>
    $PY star_condition/R_0_30_Z_0.3/input_data/make_input.py            # 12 galaxies
    $PY star_condition/R_0_30_Z_0.3/input_data/make_input.py 538905     # selected galaxies only
Outputs: input_data/<version>/<ID>_e.csv, true_value/<version>/<ID>_truth.csv (version = nonoise, noise), input_data/selection.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
COND = HERE.parents[1]                                   # star_condition/R_0_30_Z_0.3/
SRC = COND.parent / 'R_0_30_Z_30'
Z_MAX = 0.3
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    ids = sys.argv[1:] or IDS
    rows = {}
    for gid in ids:
        for v in ['nonoise', 'noise']:
            e = pd.read_csv(SRC / 'input_data' / v / f'{gid}_e.csv')
            t = pd.read_csv(SRC / 'true_value' / v / f'{gid}_truth.csv')
            assert len(e) == len(t)
            m = (np.abs(t['z'].to_numpy()) < Z_MAX)
            for sub, df in [('input_data', e), ('true_value', t)]:
                d = COND / sub / v
                d.mkdir(parents=True, exist_ok=True)
                df[m].to_csv(d / (f'{gid}_e.csv' if sub == 'input_data' else f'{gid}_truth.csv'), index=False)
            rows[gid] = (len(t), int(m.sum()))
        print(f'{gid}: R_0_30_Z_30 {rows[gid][0]:,} stars → {rows[gid][1]:,} after |z|<{Z_MAX}', flush=True)
    rep = COND / 'input_data' / 'selection.txt'
    old = {}
    if rep.exists() and len(sys.argv) > 1:                       # merge when run in batches
        for l in rep.read_text().splitlines():
            p = l.split()
            if p and p[0].isdigit():
                old[p[0]] = (int(p[1].replace(',', '')), int(p[2].replace(',', '')))
    old.update(rows)
    lines = [f'selection: 0 < R < 30 kpc, |z| < {Z_MAX} kpc; otherwise as R_0_30_Z_30 (bound with all orbit quantities; has_birth and R_bir < 20; '
             '[Fe/H] > -0.5; J_z > 0; no [Mg/Fe] cut); selected from the R_0_30_Z_30 files',
             '', f"{'ID':>7} {'R_0_30_Z_30':>12} {f'|z|<{Z_MAX}':>10}"]
    for g in sorted(old):
        lines.append(f'{g:>7} {old[g][0]:>12,} {old[g][1]:>10,}')
    tot = [v[1] for v in old.values()]
    lines += ['', f'final sample: {min(tot):,}–{max(tot):,} stars per galaxy, {sum(tot):,} in total ({len(old)} galaxies)']
    rep.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
