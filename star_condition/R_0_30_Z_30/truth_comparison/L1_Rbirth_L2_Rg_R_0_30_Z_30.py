# -*- coding: utf-8 -*-
"""
L1_Rbirth_L2_Rg_R_0_30_Z_30.py — 12 TNG50 galaxies for the selection R_0_30_Z_30 (0<R<30 kpc, |z|<30 kpc):
f(L1) vs true Rbirth and f(L2) vs true Rg (same as results/TNG/truth_comparison/tng_L1_Rbirth_L2_Rg.py with different data)

Run (numpy, pandas, matplotlib; the analysis environment is enough):
    cd <repo root>
    $PY "star_condition/R_0_30_Z_30/truth_comparison/L1_Rbirth_L2_Rg_R_0_30_Z_30.py"                  # no errors
    $PY "star_condition/R_0_30_Z_30/truth_comparison/L1_Rbirth_L2_Rg_R_0_30_Z_30.py" --version noise  # → ..._noise.png
Output: L1_Rbirth_L2_Rg_R_0_30_Z_30.png next to this script

Procedure (identical to reproduce_Jin/truth_comparison/jin_L1_Rbirth_L2_Rg.py, with TNG50 data and coefficients):
  1. data: input_data/<version>/<ID>_e.csv (6 columns logJZ, JPhi, ecc, age, FeH, MgH), every column standardized;
     truth: R_bir (birth radius) and Rg (guiding radius) from true_value/<version>/<ID>_truth.csv, row-matched to _e.csv.
  2. edge coefficients: the authors' graph frozen (as in paper Secs. 4.2, 4.3), coefficients from
     the smallest-NLL solution in output_data/<version>_<ID>_ci.csv (maximum likelihood; same as drawn in causal_graph/).
  3. L2, L1 of every star (L2 first) by hierarchical least squares:
         r_c = c − Σ_{p≠L} a_pc · p          (child c minus the contributions of its other parents p)
         L   = Σ_c a_Lc · r_c / Σ_c a_Lc²    (the L minimizing Σ_c (r_c − a_Lc·L)²)
       c: a child of L (standardized value); p: another parent of c; a_pc: coefficient p→c; a_Lc: coefficient L→c
  4. linear conversion: Rbirth ≈ α·L1 + β, Rg ≈ α'·L2 + β' (least squares; uses the truth only to set the units), Pearson r annotated.
  5. each panel also shows (REPORT.md Sec. 5.6):
       range: |r| range over all solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84) — the coefficients are not identifiable, only one is drawn;
       age only / [Fe/H] only / J_phi only: |r| of age, [Fe/H] or J_phi alone with the true Rbirth (baselines);
         ⚠ this selection covers the whole disk, where J_phi (≈ current radius) is already strongly correlated with the birth radius; L1 must at least beat it to be meaningful;
       L1 resembles: |r| of L1 with [Fe/H] and age, i.e. which observable L1 actually "looks like".
  How to read it: a narrow vertical band in an L1 panel is not a bug but the signature of a weak correlation: the slope of the linear conversion f = α·L + β is
  α = r·σ_R/σ_L (r: correlation, σ_R: std of the truth, σ_L: std of L), so the spread of f is only |r| times the spread of the truth;
  for r ≈ 0.1 all values of f pile up near the mean.
  The L1 computation was checked on the results/TNG/truth_comparison/ version (differs by < 1e-7 from the authors' per-star Nelder–Mead).
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
COND = HERE.parents[1]                                   # star_condition/R_0_30_Z_30/
ROOT = HERE.parents[3]                                   # repo root
IDENT = COND / 'output_data'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LR95 = 3.84
INK = '#1f2937'


def coefficients(row):
    """one row of <tag>_ci.csv → {(parent, child): coefficient}, X_ prefix removed"""
    out = {}
    for k, v in row.items():
        if '->' in str(k):
            a, b = k.split('->')
            out[(a.replace('X_', ''), b.replace('X_', ''))] = float(v)
    return out


def solve_latent(L, coef, val):
    num, den = 0.0, 0.0
    for (p, c), a in coef.items():
        if p != L:
            continue
        r = val[c].copy()
        for (q, c2), b in coef.items():
            if c2 == c and q != L:
                r = r - b * val[q]
        num, den = num + a * r, den + a * a
    return num / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    out = HERE.with_suffix('.png') if args.version == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')

    fig, axs = plt.subplots(6, 4, figsize=(17, 25))
    lines = []
    for k, gid in enumerate(IDS):
        x = pd.read_csv(COND / 'input_data' / args.version / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(COND / 'true_value' / args.version / f'{gid}_truth.csv', usecols=['R_bir', 'Rg'])
        z = (x - x.mean()) / x.std()
        ci = pd.read_csv(IDENT / f'{args.version}_{gid}_ci.csv').sort_values('nll')
        coef = coefficients(ci.iloc[0])
        inn = ci[ci.lr < LR95]
        val = {c: z[c].to_numpy() for c in COLS}
        val['L2'] = solve_latent('L2', coef, val)
        val['L1'] = solve_latent('L1', coef, val)
        r_age = abs(np.corrcoef(x['age'], t['R_bir'])[0, 1])
        r_feh = abs(np.corrcoef(x['FeH'], t['R_bir'])[0, 1])
        r_jphi = abs(np.corrcoef(x['JPhi'], t['R_bir'])[0, 1])
        row, c0 = divmod(k, 2)
        for j, (L, truth, lab, rng) in enumerate([('L1', 'R_bir', 'R_b', ('abs_r_L1_Rbir')),
                                                  ('L2', 'Rg', 'R_g', ('abs_r_L2_Rg'))]):
            ax = axs[row, 2 * c0 + j]
            y = t[truth].to_numpy()
            slope, icpt = np.polyfit(val[L], y, 1)
            f = slope * val[L] + icpt
            r = np.corrcoef(f, y)[0, 1]
            lim = 20 if L == 'L1' else 30          # birth radius < 20 kpc; guiding radius up to 30 kpc
            ax.hexbin(f, y, bins='log', extent=(0, lim, 0, lim), gridsize=60, cmap='Blues', mincnt=1)
            ax.plot([0, lim], [0, lim], color=INK, lw=1)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_aspect('equal')
            ax.set_xlabel(f'$f(L_{L[1]})$ [kpc]', fontsize=9)
            ax.set_ylabel(f'${lab}$ (true) [kpc]', fontsize=9)
            ax.tick_params(labelsize=8)
            txt = (f'$r$ = {r:.2f}   ($f$ = {slope:+.2f} $L_{L[1]}$ {icpt:+.2f})\n'
                   f'range (indistinguishable fits): {inn[rng].min():.2f}–{inn[rng].max():.2f}')
            if L == 'L1':
                txt += (f'\nage only: {r_age:.2f}   [Fe/H] only: {r_feh:.2f}   ' + r'$J_\phi$' + f' only: {r_jphi:.2f}'
                        f'\nL1 resembles: |r(L1,[Fe/H])| {abs(np.corrcoef(val["L1"], x["FeH"])[0, 1]):.2f},'
                        f' |r(L1,age)| {abs(np.corrcoef(val["L1"], x["age"])[0, 1]):.2f}')
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=7.8, color=INK,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.85))
            ax.set_title(f'{gid}  ({len(x):,} stars)  —  {"birth radius" if L == "L1" else "guiding radius"}',
                         fontsize=9.5, color=INK)
            lines.append(f'{gid} {truth} = {slope:+.3f}*{L} {icpt:+.3f}, r = {r:.3f}')
    fig.suptitle(f'TNG50, stars 0<R<30 kpc, |z|<30 kpc ({args.version}): author\'s causal graph frozen, maximum-likelihood edge coefficients per galaxy\n'
                 'f(L1) vs true birth radius, f(L2) vs true guiding radius.  Coefficients are not unique: '
                 '"range" = |r| over all statistically indistinguishable fits.\n'
                 'A narrow vertical band = weak correlation: the linear map f = a·L + b shrinks toward the mean '
                 '(spread of f = |r| × spread of the truth).', fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.968))
    fig.savefig(out, dpi=100)
    print('\n'.join(lines))
    print(f'saved {out}')


if __name__ == '__main__':
    main()
