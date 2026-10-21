# -*- coding: utf-8 -*-
"""
make_input.py — Approach B, step 0: prepare the input data and the truth (2026-09-27)

Approach B first uses the same main selection as the paper: 7 < R < 10 kpc, |z| < 0.3 kpc (plus: bound with all orbit quantities; birth radius with R_bir < 20 kpc;
[Fe/H] > -0.5; J_z > 0), i.e. the same stars as the main result of Approach A. These stars were selected by scripts/local/make_rlcd_input.py
and stored in data/rlcd_input/<version>/; this script just copies them unchanged into PlanB/, so that all Approach-B files live in PlanB:
    data/rlcd_input/<version>/<ID>_e.csv      → PlanB/input_data/<version>/<ID>_e.csv     (6-column RLCD input)
    data/rlcd_input/<version>/<ID>_truth.csv  → PlanB/true_value/<version>/<ID>_truth.csv (true birth radius R_bir, true guiding radius Rg, ...)
The two files are row-matched (row i is the same star). Versions: nonoise (no observational errors), noise (with observational errors).

Run (numpy, pandas):
    cd <repo root>
    $PY PlanB/input_data/make_input.py
Outputs: the csv files in the two folders above, and PlanB/input_data/selection.txt (number of stars per galaxy)
"""
import shutil
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
SRC = ROOT / 'data' / 'rlcd_input'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    lines = ['Approach B selection (as in the paper): 7 < R < 10 kpc, |z| < 0.3 kpc; bound with all orbit quantities; has_birth and R_bir < 20 kpc; '
             '[Fe/H] > -0.5; J_z > 0', f'source: {SRC.relative_to(ROOT)}/<version>/ (copied unchanged)', '',
             f"{'ID':>7} {'nonoise':>9} {'noise':>9}"]
    tot = 0
    for gid in IDS:
        n = {}
        for v in ['nonoise', 'noise']:
            for sub, name in [('input_data', f'{gid}_e.csv'), ('true_value', f'{gid}_truth.csv')]:
                d = PLANB / sub / v
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(SRC / v / name, d / name)
            e = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')
            t = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir'])
            assert len(e) == len(t), gid
            n[v] = len(e)
        tot += n['nonoise']
        lines.append(f"{gid:>7} {n['nonoise']:>9,} {n['noise']:>9,}")
        print(lines[-1], flush=True)
    lines += ['', f'total {tot:,} stars (12 galaxies)']
    (PLANB / 'input_data' / 'selection.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
