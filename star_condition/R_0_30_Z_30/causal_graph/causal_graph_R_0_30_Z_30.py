# -*- coding: utf-8 -*-
"""
causal_graph_R_0_30_Z_30.py — causal graphs of the 12 TNG50 galaxies for the selection R_0_30_Z_30 (0<R<30 kpc, |z|<30 kpc)
(Approach A: the authors' graph frozen, edge coefficients re-estimated per galaxy; after paper Fig. 5; same as results/TNG/causal_graph/tng_causal_graph.py with different data)

Run (numpy, pandas, matplotlib; the analysis environment is enough):
    cd <repo root>
    $PY "star_condition/R_0_30_Z_30/causal_graph/causal_graph_R_0_30_Z_30.py"                  # no errors
    $PY "star_condition/R_0_30_Z_30/causal_graph/causal_graph_R_0_30_Z_30.py" --version noise  # → ..._noise.png
Output: causal_graph_R_0_30_Z_30.png next to this script

Inputs (all in this repo):
  - graph structure: the authors' final graph results/rlcd_g279/author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot
    (same as paper Fig. 1, 15 edges; paper Secs. 4.2, 4.3 freeze this graph and transfer it to other galaxies)
  - edge coefficients: the row with the smallest NLL in star_condition/R_0_30_Z_30/output_data/<version>_<ID>_ci.csv, i.e.
    the maximum-likelihood solution of output_data/identifiability_moments.py (230 starts, strict convergence, noise variances >= 0, coefficients within ±0.999).
    This is what "the maximum-likelihood method of the paper, run to convergence" gives.

⚠ Important (REPORT.md Sec. 4.3):
  - with this graph frozen the coefficients are not identifiable: many solutions are statistically indistinguishable from the best (n_eq in each panel title)
    and their coefficients can differ a lot; only the one with the smallest NLL is drawn;
  - observables marked "ω≈0": the optimum pushes their noise variance to 0 (the model treats them as fully determined by their parents, which is unphysical);
  - a coefficient of ±0.999 means it hit the bound.
"""
import re
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
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDENT = COND / 'output_data'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
LR95 = 3.84

# ---------------------------------------------------------------- drawing (same style as the reproduce_Jin/causal_graph script)
NODES = {  # name: (x, y, label, is latent)
    'L1': (0.55, 0.93, 'L1\nBirth Radius', True),
    'L2': (0.08, 0.62, 'L2\nGuiding Radius', True),
    'MgH': (0.36, 0.60, '[Mg/H]', False),
    'age': (0.74, 0.60, 'Stellar Age', False),
    'JPhi': (0.06, 0.12, r'$J_\phi$' + '\nAzimuthal Action', False),
    'ecc': (0.36, 0.12, r'$e$' + '\nEccentricity', False),
    'FeH': (0.70, 0.12, '[Fe/H]', False),
    'logJZ': (1.10, 0.60, r'$J_z$' + '\nVertical Action', False),
}
EW, EH = 0.105, 0.065          # ellipse half width, half height
BEND = {('MgH', 'logJZ'): -0.22, ('L1', 'logJZ'): 0.04}   # curvature: offset of the control point along the left normal of the edge (negative bends right/down)
LABEL_T = {('age', 'JPhi'): 0.22, ('MgH', 'FeH'): 0.55, ('age', 'ecc'): 0.30, ('MgH', 'ecc'): 0.55,
           ('L2', 'ecc'): 0.35, ('age', 'FeH'): 0.55, ('MgH', 'JPhi'): 0.62, ('MgH', 'logJZ'): 0.80,
           ('L1', 'FeH'): 0.30, ('L2', 'JPhi'): 0.5}
LATENT_FC, LATENT_EC = '#fde2e1', '#c2362f'
OBS_FC, OBS_EC = '#dbe8fb', '#2563eb'
INK, MUTED = '#1f2937', '#6b7280'


def _edge_point(center, toward):
    import numpy as np
    d = np.array(toward, float) - np.array(center, float)
    s = 1.0 / np.sqrt((d[0] / EW) ** 2 + (d[1] / EH) ** 2)
    return np.array(center, float) + d * s


def draw_graph(ax, coef, title=None, fontsize=8.5):
    """coef: {(parent, child): coefficient or None}; variable names without the X_ prefix"""
    import numpy as np
    from matplotlib.patches import Ellipse, FancyArrowPatch
    from matplotlib.path import Path as MPath
    ax.set_xlim(-0.07, 1.24)
    ax.set_ylim(0.03, 1.02)
    ax.set_aspect('equal')
    ax.axis('off')
    for name, (x, y, lab, lat) in NODES.items():
        ax.add_patch(Ellipse((x, y), 2 * EW, 2 * EH, fc=LATENT_FC if lat else OBS_FC,
                             ec=LATENT_EC if lat else OBS_EC, lw=1.4, zorder=3))
        ax.text(x, y, lab, ha='center', va='center', fontsize=fontsize - 0.5, color=INK, zorder=4)
    vals = [abs(v) for v in coef.values() if v is not None]
    vmax = max(vals) if vals else 1
    for (p, c), v in coef.items():
        P, Q = np.array(NODES[p][:2]), np.array(NODES[c][:2])
        d = Q - P
        C = (P + Q) / 2 + BEND.get((p, c), 0.0) * np.array([-d[1], d[0]]) / np.hypot(*d)
        P1, Q1 = _edge_point(P, C), _edge_point(Q, C)
        path = MPath([P1, C, Q1], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
        lw = 0.7 + 3.3 * (abs(v) / vmax if v is not None else 0.3)
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=10 + 2 * lw, lw=lw,
                                     color=INK, alpha=0.8, zorder=2))
        if v is not None:
            t = LABEL_T.get((p, c), 0.5)
            B = (1 - t) ** 2 * P1 + 2 * (1 - t) * t * C + t ** 2 * Q1
            ax.text(B[0], B[1], f'{v:.3f}' if abs(v) >= 0.995 else f'{v:.2g}', ha='center', va='center', fontsize=fontsize - 1, color=INK, zorder=5,
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.92))
    if title:
        ax.set_title(title, fontsize=fontsize + 1.5, color=INK)


def author_edges():
    txt = open(AUTHOR_DOT).read()
    return [(a.replace('X_', ''), b.replace('X_', '')) for a, b in re.findall(r'^\s*(\w+)\s*->\s*(\w+)', txt, re.M)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    out = HERE.with_suffix('.png') if args.version == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')
    edges = author_edges()
    col = {e: ('X_' if e[0] not in ('L1', 'L2') else '') + e[0] + '->' + ('X_' if e[1] not in ('L1', 'L2') else '') + e[1]
           for e in edges}

    fig, axs = plt.subplots(3, 4, figsize=(22, 15.5))
    for ax, gid in zip(axs.flat, IDS):
        ci = pd.read_csv(IDENT / f'{args.version}_{gid}_ci.csv').sort_values('nll')
        best = ci.iloc[0]
        coef = {e: float(best[col[e]]) for e in edges}
        n_eq = int((ci.lr < LR95).sum())
        draw_graph(ax, coef, fontsize=8.5)
        node = str(best.omega_min_node).replace('X_', '')
        if best.omega_min < 0.005 and node in NODES:
            x, y = NODES[node][:2]
            ax.text(x, y - EH - 0.035, 'ω≈0', ha='center', va='center', fontsize=8, color='#c2362f', fontweight='bold')
        ax.set_title(f'{gid}   (NLL {best.nll:.3f};  statistically indistinguishable fits n_eq = {n_eq})',
                     fontsize=10.5, color=INK)
    fig.suptitle(f'TNG50, stars 0<R<30 kpc, |z|<30 kpc ({args.version}): author\'s causal graph frozen, edge coefficients = maximum-likelihood fit per galaxy\n'
                 'ω≈0: noise variance of that observable driven to 0 at the optimum;  coefficients are NOT unique '
                 '(see REPORT.md §4.3)', fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=110)
    print(f'saved {out}')


if __name__ == '__main__':
    main()
