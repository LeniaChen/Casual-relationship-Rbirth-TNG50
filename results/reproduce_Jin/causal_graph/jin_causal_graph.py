# -*- coding: utf-8 -*-
"""
jin_causal_graph.py — reproduce Fig. 1 of Jin et al.: the causal graph learned by RLCD on g2.79e12 (with edge coefficients)

Run (must use the RLCD environment, which has causal-learn, statsmodels, ...):
    cd <repo root>
    $PY_RLCD "results/reproduce_Jin/causal_graph/jin_causal_graph.py"
Output: jin_causal_graph.png next to this script (about 1 minute)

What it does:
  1. Structure learning (actually rerun): the repo's RLCD on data/g2.79e12_e.csv (6 columns logJZ, JPhi, ecc, age, FeH, MgH),
     with the parameters in the authors' result file name: alpha0..3 = 1e-5, 1e-15, 1e-5, 1e-15, stage-1 GES (sparsity 2),
     exactly the defaults of rlcd/run_rlcd.py.
     ⚠ The RLCD result depends on Python's string-hash randomization (issue 1): the authors' graph is obtained only with PYTHONHASHSEED=1.
       This script restarts itself with PYTHONHASHSEED=1; use --seed 3 to see another seed.
       With a fixed seed the result is deterministic (same seed, same graph); "random" means it changes with the seed.
       Scan of PYTHONHASHSEED = 0–29 (30 seeds, 2026-09-27): the authors' graph appears 7 times (1, 7, 8, 11, 14, 22, 25);
       the most common result is another 15-edge graph (12 times, e.g. seeds 0, 2, 3); the other 11 runs give 18-edge graphs (4 kinds).
     Structure learning runs in memory and writes nothing into the repo's results_* folders (the authors' files are not overwritten).
  2. Edge-by-edge comparison with the authors' committed graph (results/rlcd_g279/author_result_..._ranktestscale1.dot)
     (latent labels aligned automatically) and check of the prior of paper Sec. 3.2 (age has no observed parent).
  3. Numbers on the edges: the authors' published parameter-learning result (…_param_by_trek_nll1.5861642734160428.dot, i.e. the numbers of paper Fig. 1).
     ⚠ These numbers are NOT computed here; they are read from the authors' file, so of course they match the paper.
       Parameter learning uses unseeded random starts + SGD, so the authors' run cannot be reproduced bit by bit (the 13 parameter-learning
       results for the same graph in the authors' repo all differ); and this set is not the maximum-likelihood solution (issue 2). Our own fits are in results/rlcd_A/identifiability/g279_ci.*.
  4. Drawing: nodes as in paper Fig. 1 (layout slightly adjusted to reduce crossings); red = latents, blue = observables; line width ∝ |coefficient|;
     the text below the graph states whether the structure is identical to the authors'.
"""
import os
import re
import sys
import argparse
import itertools
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # repo root
REPO = ROOT / 'Causal-Galactic-Archaeology'
DATA = REPO / 'data' / 'g2.79e12_e.csv'
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
COEF_DOT = REPO / 'rlcd' / 'results_causal_vertical_heating_g279_e' / \
    'result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1_param_by_trek_nll1.5861642734160428.dot'
OUT = HERE.with_suffix('.png')
ALPHA = (1e-5, 1e-15, 1e-5, 1e-15)

# ---------------------------------------------------------------- drawing (same style as the TNG/causal_graph script)
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


def read_dot_edges(path, with_coef=False):
    """Read the directed edges of a .dot → {(parent, child): coefficient or None}, X_ prefix removed"""
    txt = open(path).read()
    out = {}
    for m in re.finditer(r'^\s*(\w+)\s*->\s*(\w+)\s*(\[[^\]]*\])?', txt, re.M):
        a, b, attr = m.group(1).replace('X_', ''), m.group(2).replace('X_', ''), m.group(3) or ''
        lab = re.search(r'label="?([-+\deE.]+)', attr)
        out[(a, b)] = float(lab.group(1)) if (with_coef and lab) else None
    return out


def align_latents(ours, ref):
    """Latent labels are arbitrary: try all relabelings and keep the one sharing most edges with the reference"""
    lat = sorted({v for e in ours for v in e if v.startswith('L')})
    best = None
    for perm in itertools.permutations(lat):
        mp = dict(zip(lat, perm))
        ren = {(mp.get(a, a), mp.get(b, b)) for a, b in ours}
        score = len(ren & set(ref))
        if best is None or score > best[0]:
            best = (score, ren, mp)
    return best[1], best[2]


# ---------------------------------------------------------------- structure learning
def run_rlcd(seed):
    """Run RLCD in this process (same settings as run_rlcd.py); return (set of directed edges, set of undirected edges)"""
    os.chdir(REPO / 'rlcd')                    # the repo code finds its files by relative paths
    sys.path.insert(0, str(REPO / 'rlcd'))
    import pandas as pd
    from DGM.DataModel import DataModel
    from utils.Chi2RankTest import Chi2RankTest
    from StructureLearning.RLCD.RLCD_alg import RLCD
    import logging
    logging.getLogger().setLevel(logging.WARNING)   # RLCD logs a lot; keep warnings only

    df = pd.read_csv(DATA).dropna()                 # same as scenarios.causal_vertical_heating_g279_e()
    df = df.rename(columns={k: 'X_' + k.replace('-', '_') for k in df.columns})
    dgm = DataModel(df, df)
    df_x, _ = dgm.generate_data(N=-1, normalized=True)
    params = {'ranktest_method': Chi2RankTest(df_x.to_numpy(), 1), 'citest_method': None,
              'stage1_method': 'ges', 'alpha_dict': dict(enumerate(ALPHA)),
              'stage1_ges_sparsity': 2, 'stage1_partition_thres': 3}
    graph, _, _, _ = RLCD(1, dgm.xvars, df_x, params)
    strip = lambda v: v.replace('X_', '')
    directed = {(strip(a), strip(b)) for a, b in graph.dirEdges}
    undirected = {tuple(sorted(map(strip, e))) for e in graph.undirEdges}
    return directed, undirected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', default='1', help='PYTHONHASHSEED (the authors\' graph corresponds to 1)')
    args = ap.parse_args()
    if os.environ.get('PYTHONHASHSEED') != args.seed:          # the hash seed can only be set before Python starts
        env = dict(os.environ, PYTHONHASHSEED=args.seed)
        os.execve(sys.executable, [sys.executable, str(HERE)] + sys.argv[1:], env)

    ours, undirected = run_rlcd(args.seed)
    author = read_dot_edges(AUTHOR_DOT)
    ours_aligned, mapping = align_latents(ours, author)
    same = ours_aligned == set(author) and not undirected
    extra, missing = sorted(ours_aligned - set(author)), sorted(set(author) - ours_aligned)
    obs = {'JPhi', 'ecc', 'age', 'FeH', 'MgH', 'logJZ'}
    prior_ok = not any(c == 'age' and p in obs for p, c in ours_aligned)

    coef = read_dot_edges(COEF_DOT, with_coef=True)
    assert set(coef) == set(author), 'edges in the coefficient file differ from the authors\' graph'

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 7.6))
    if same:
        draw_graph(ax, coef, fontsize=10)
        note = (f'Structure re-learned with RLCD (PYTHONHASHSEED={args.seed}, alpha={ALPHA}): '
                f'IDENTICAL to the author\'s graph ({len(ours)} edges);\nage-has-no-observed-parent prior: '
                f'{"passed" if prior_ok else "FAILED"}.\n'
                'Edge numbers: READ from the author\'s released parameter-learning file (NLL 1.586), not re-computed '
                '(the paper\'s Fig. 1 numbers).')
    else:
        for L in ('L1', 'L2'):                     # if the structure differs, the latents no longer mean "birth radius / guiding radius"
            x, y, _, lat = NODES[L]
            NODES[L] = (x, y, L, lat)
        draw_graph(ax, {e: None for e in ours_aligned if e[0] in NODES and e[1] in NODES}, fontsize=10)
        note = (f'Structure re-learned with RLCD (PYTHONHASHSEED={args.seed}): DIFFERENT from the author\'s graph.\n'
                f'extra edges: {extra}   missing edges: {missing}   undirected: {sorted(undirected)}   '
                f'prior {"passed" if prior_ok else "FAILED"}')
    fig.text(0.5, 0.015, note, ha='center', va='bottom', fontsize=7.5, color=MUTED, wrap=True)
    fig.suptitle('Reproduction of Jin et al. Fig. 1 — causal graph of g2.79e12', fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(OUT, dpi=150)

    print(f'PYTHONHASHSEED = {args.seed}; latent alignment {mapping}')
    print(f'our graph: {len(ours)} directed edges, {len(undirected)} undirected edges; authors\' graph: {len(author)}')
    print('identical to the authors\' graph' if same else f'different from the authors\' graph: extra {extra}, missing {missing}')
    print(f'prior (age has no observed parent): {"pass" if prior_ok else "FAIL"}')
    print(f'saved {OUT}')


if __name__ == '__main__':
    main()
