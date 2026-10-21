# -*- coding: utf-8 -*-
"""
planB_latent_vs_truth.py — Approach B: latents of the causal graph each TNG50 galaxy "learns by itself", compared with the true birth and guiding radii

Run (numpy, pandas, matplotlib; the analysis environment is enough) after PlanB/output_data/rlcd_B.py and rlcd_B_eval.py:
    cd <repo root>
    $PY "PlanB/truth_comparison/planB_latent_vs_truth.py"                  # no errors
    $PY "PlanB/truth_comparison/planB_latent_vs_truth.py" --version noise  # → planB_latent_vs_truth_noise.png
Output: planB_latent_vs_truth.png next to this script

Procedure (each galaxy uses its most frequent structure; if that structure has undirected edges or no latent, the panel says so and stays empty):
  1. read the 6 inputs of PlanB/input_data/<version>/<ID>_e.csv and standardize each column (subtract mean, divide by std) to get z_i;
  2. each latent L = Σ_i w_i·z_i. w comes from eval_<version>_best.json: the weights computed by rlcd_B_eval.py with the same hierarchical least squares as the authors' notebook
     (from the lowest latent upwards, inferring each latent from its children) at the maximum-likelihood edge coefficients;
  3. ⚠ graphs learned in Approach B carry no label saying "this latent is the birth radius" (the numbering is arbitrary). The latent with the largest |r| with the true
     birth radius is taken as the "birth-radius candidate" L_Rb, and the one with the largest |r| with the true guiding radius as L_Rg. This "peeks at the answer": a best case that is impossible with real observations;
  4. linear conversion (as in the authors' notebook): f(L) = α·L + β, where α, β define the line minimizing Σ (truth − α·L − β)²
     (α: slope, β: intercept); x axis f(L_Rb), y axis true birth radius; black line 1:1;
  5. baselines in each panel: |r| of age, [Fe/H], J_φ alone with the true birth radius, and the linear-regression ceiling of the 6 inputs (no "linear combination of the 6 inputs"
     can exceed it in |r| with the birth radius).
  The horizontal width of f(L) = |r| × width of the truth (r: Pearson correlation of L with the truth), so for small |r| the points form a narrow vertical band.
"""
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    v = args.version
    st = pd.read_csv(PLANB / 'output_data' / f'structures_{v}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    best = json.load(open(PLANB / 'output_data' / f'eval_{v}_best.json'))
    fig, axs = plt.subplots(6, 4, figsize=(17, 25))
    lines = []
    for k, gid in enumerate(IDS):
        axR, axG = axs.flat[2 * k], axs.flat[2 * k + 1]
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        sid, cnt = freq.index[0], freq.iloc[0]
        info = best.get(gid, {}).get(str(sid))
        e = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir', 'Rg'])
        z = ((e - e.mean()) / e.std()).to_numpy()
        rb = lambda c: abs(np.corrcoef(e[c], t.R_bir)[0, 1])
        base = f'age only {rb("age"):.2f}   [Fe/H] only {rb("FeH"):.2f}   $J_\\phi$ only {rb("JPhi"):.2f}'
        beta, *_ = np.linalg.lstsq(np.column_stack([z, np.ones(len(z))]), t.R_bir.to_numpy(), rcond=None)
        ceil = np.corrcoef(np.column_stack([z, np.ones(len(z))]) @ beta, t.R_bir)[0, 1]
        head = f'{gid} ({len(e):,} stars) — structure {sid} ({cnt}/{len(g)} seeds)'
        if info is None or not info['weights']:
            why = 'no latent variable' if info is not None else 'structure has undirected edges / not evaluated'
            for ax, nm in [(axR, 'birth radius'), (axG, 'guiding radius')]:
                ax.text(0.5, 0.5, why, ha='center', va='center', transform=ax.transAxes, fontsize=10, color='#6b7280')
                ax.set_title(f'{head}\n{nm}', fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            axR.text(0.04, 0.95, base + f'\nlinear ceiling {ceil:.2f}', transform=axR.transAxes, va='top', fontsize=7.5)
            lines.append(f'{gid}: structure {sid} — {why}')
            continue
        Lval = {L: z @ np.array(w) for L, w in info['weights'].items()}
        for ax, col, nm, lim in [(axR, 'R_bir', 'birth radius', 20), (axG, 'Rg', 'guiding radius', 20)]:
            y = t[col].to_numpy()
            rs = {L: abs(np.corrcoef(x, y)[0, 1]) for L, x in Lval.items()}
            L = max(rs, key=rs.get)
            a, b = np.polyfit(Lval[L], y, 1)
            f = a * Lval[L] + b
            r = np.corrcoef(f, y)[0, 1]
            ax.hexbin(f, y, bins='log', extent=(0, lim, 0, lim), gridsize=60, cmap='Blues', mincnt=1)
            ax.plot([0, lim], [0, lim], color='#1f2937', lw=1)
            ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal')
            ax.set_xlabel(f'$f({L})$ [kpc]', fontsize=9)
            ax.set_ylabel(('$R_b$' if col == 'R_bir' else '$R_g$') + ' (true) [kpc]', fontsize=9)
            txt = f'$r$ = {r:.2f}   ($f$ = {a:+.2f} {L} {b:+.2f})\nall latents |r|: ' + ', '.join(f'{k2} {v2:.2f}' for k2, v2 in sorted(rs.items()))
            if col == 'R_bir':
                txt += f'\n{base}\nlinear ceiling {ceil:.2f}'
            ax.text(0.04, 0.96, txt, transform=ax.transAxes, va='top', fontsize=7.3)
            ax.set_title(f'{head}\n{nm}: best-matching latent = {L}', fontsize=9)
            lines.append(f'{gid} {col}: best-matching latent {L}, f = {a:+.3f}*{L} {b:+.3f}, r = {r:.3f}')
    fig.suptitle(f'Plan B — TNG50, 7<R<10 kpc, |z|<0.3 kpc ({v}): latents of the graph RLCD learned in each galaxy (most frequent of 30 seeds) vs truth\n'
                 'Latent names are arbitrary, so the latent with the highest |r| is chosen using the TRUE values (best case, not possible for real data).\n'
                 'f(L) = linear least-squares map of L onto the truth; spread of f = |r| × spread of the truth.', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = HERE.with_suffix('.png') if v == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')
    fig.savefig(out, dpi=100)
    print('\n'.join(lines))
    print(f'saved {out}')


if __name__ == '__main__':
    main()
