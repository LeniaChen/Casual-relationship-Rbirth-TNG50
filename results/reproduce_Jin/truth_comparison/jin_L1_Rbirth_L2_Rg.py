# -*- coding: utf-8 -*-
"""
jin_L1_Rbirth_L2_Rg.py — reproduce the two panels of Fig. 3 of Jin et al.: f(L1) vs true Rbirth, f(L2) vs true Rg (simulated galaxy g2.79e12)

Run (numpy, pandas, matplotlib; the analysis environment is enough):
    cd <repo root>
    $PY "results/reproduce_Jin/truth_comparison/jin_L1_Rbirth_L2_Rg.py"
Output: jin_L1_Rbirth_L2_Rg.png next to this script

Identical to the authors' notebook (Causal-Galactic-Archaeology/L1L2calculation/GetL1_L2_279sim.ipynb):
  1. data: data/g2.79e12_multi_ele_action_Jin_SN.csv (119,236 stars); MgH = MgFe + FeH; logjz = ln(JZ);
     every column standardized (subtract mean, divide by standard deviation).
  2. edge coefficients: the authors' published parameter-learning result
     rlcd/results_causal_vertical_heating_g279_e/result_..._param_by_trek_nll1.5861642734160428.dot
     (identical to the numbers hard-coded in the notebook; these are the numbers of paper Fig. 1).
  3. L2, L1 of every star (L2 first): the notebook minimizes χ² per star with Nelder-Mead;
     χ² is quadratic, so the minimum has an exact solution, used here (same result, thousands of times faster):
         r_c = c − Σ_{p≠L} a_pc · p          (child c minus the contributions of its other parents p)
         L   = Σ_c a_Lc · r_c / Σ_c a_Lc²    (the L minimizing Σ_c (r_c − a_Lc·L)²)
       c: a child of L (standardized); p: another parent of c; a_pc: coefficient p→c; a_Lc: coefficient L→c
  4. linear conversion (fitpoints in the notebook): Rbirth ≈ α·L1 + β, Rg ≈ α'·L2 + β' (least squares),
     f(L1) = α·L1 + β, f(L2) = α'·L2 + β'.
  5. plot as in the notebook: hexbin (log colour scale), range 0–20 kpc, black 1:1 line, Pearson r annotated.
Paper: Rb = −1.97 L1 + 9.31, r = 0.64;  Rg = −1.29 L2 + 9.52 (r ≈ 0.98).

Note (REPORT.md Secs. 4.2, 4.3): the authors' coefficients are not the maximum-likelihood solution (NLL 1.593, best ~ −1.272);
refitting to convergence with the maximum-likelihood method described in the paper gives different results. This script only "reproduces the authors' figure with the authors' coefficients".
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # repo root
REPO = ROOT / 'Causal-Galactic-Archaeology'
DATA = REPO / 'data' / 'g2.79e12_multi_ele_action_Jin_SN.csv'
DOT = REPO / 'rlcd' / 'results_causal_vertical_heating_g279_e' / \
    'result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1_param_by_trek_nll1.5861642734160428.dot'
OUT = HERE.with_suffix('.png')

# coefficients hard-coded in the authors' notebook (used to check the dot file)
NOTEBOOK = {('L2', 'JPhi'): -0.57, ('MgH', 'JPhi'): -0.082, ('age', 'JPhi'): -0.18,
            ('L2', 'ecc'): 0.085, ('MgH', 'ecc'): 0.098, ('age', 'ecc'): 0.4,
            ('L1', 'logJZ'): -0.099, ('age', 'logJZ'): 0.42, ('MgH', 'logJZ'): -0.13,
            ('L1', 'FeH'): 0.84, ('age', 'FeH'): 0.0032, ('MgH', 'FeH'): 0.13,
            ('L1', 'age'): -0.27, ('L1', 'MgH'): 0.83, ('L1', 'L2'): 0.43}


def read_coefficients(dot):
    """Read the coefficient of every edge in a .dot; return {(parent, child): coefficient}, X_ prefix removed"""
    coef = {}
    for a, b, lab in re.findall(r'(\w+)\s*->\s*(\w+)\s*\[[^\]]*label="?([-+\deE.]+)', open(dot).read()):
        coef[(a.replace('X_', ''), b.replace('X_', ''))] = float(lab)
    return coef


def solve_latent(L, coef, val):
    """Hierarchical least squares: solve L from all its children (other parents of the children must be known)"""
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
    coef = read_coefficients(DOT)
    assert coef == NOTEBOOK, f'coefficients in the dot file differ from the notebook: {coef}'

    df = pd.read_csv(DATA)
    df['MgH'] = df['MgFe'] + df['FeH']
    df['logJZ'] = np.log(df['JZ'])
    z = (df - df.mean()) / df.std()
    val = {c: z[c].to_numpy() for c in ['JPhi', 'ecc', 'age', 'FeH', 'MgH', 'logJZ']}
    val['L2'] = solve_latent('L2', coef, val)          # children: JPhi, ecc
    val['L1'] = solve_latent('L1', coef, val)          # children: L2, age, MgH, FeH, logJZ

    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = [f'g2.79e12: {len(df):,} stars']
    for ax, L, truth, lab in [(axs[0], 'L1', 'R_bir', 'R_b'), (axs[1], 'L2', 'Rg', 'R_g')]:
        y = df[truth].to_numpy()
        slope, icpt = np.polyfit(val[L], y, 1)
        f = slope * val[L] + icpt
        r = np.corrcoef(f, y)[0, 1]
        hb = ax.hexbin(f, y, bins='log', extent=(0, 20, 0, 20), gridsize=80, cmap='Blues', mincnt=1)
        ax.plot([0, 20], [0, 20], color='#1f2937', lw=1)
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 20)
        ax.set_aspect('equal')
        ax.set_xlabel(f'$f(L_{L[1]})$ [kpc]')
        ax.set_ylabel(f'${lab}$ (true) [kpc]')
        ax.text(0.04, 0.95, f'$r$ = {r:.3f}\n$f(L_{L[1]})$ = {slope:+.3f} $L_{L[1]}$ {icpt:+.3f}',
                transform=ax.transAxes, va='top', fontsize=11, color='#1f2937')
        fig.colorbar(hb, ax=ax, label='# of stars')
        lines.append(f'{truth} = {slope:+.3f}*{L} {icpt:+.3f}, Pearson r = {r:.3f}')
    axs[0].set_title('birth radius vs $L_1$ (paper Fig. 3 left: r = 0.64)', fontsize=11)
    axs[1].set_title('guiding radius vs $L_2$ (paper Fig. 3 right)', fontsize=11)
    fig.suptitle("Reproduction of Jin et al. Fig. 3 — g2.79e12, author's released edge coefficients (NLL 1.586)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT, dpi=130)
    print('\n'.join(lines))
    print(f'saved {OUT}')


if __name__ == '__main__':
    main()
