# -*- coding: utf-8 -*-
"""
planB_causal_graph.py — Approach B: the causal graph each TNG50 galaxy "learns by itself" (its most frequent one) with maximum-likelihood edge coefficients

Run (numpy, pandas, matplotlib; the analysis environment is enough) after PlanB/output_data/rlcd_B.py and rlcd_B_eval.py:
    cd <repo root>
    $PY "PlanB/causal_graph/planB_causal_graph.py"                  # no errors
    $PY "PlanB/causal_graph/planB_causal_graph.py" --version noise  # with errors → planB_causal_graph_noise.png
    option: --rank 2 draws the second most frequent structure of each galaxy (→ planB_causal_graph_rank2.png)
Output: planB_causal_graph.png next to this script

Files read (PlanB/output_data/):
  structures_<version>.csv/.txt   graphs learned with 30 random seeds (PYTHONHASHSEED 0–29) per galaxy; structure ids are shared by all galaxies
  eval_<version>_best.json        maximum-likelihood solution of every (galaxy, structure): edge coefficients, noise variance of each observable
  eval_<version>.csv              |r| of every latent with the true birth and guiding radii at the maximum-likelihood solution

What the figure shows:
  - red ellipses = latents (unobserved quantities found by RLCD); blue = the 6 observables;
  - red dashed lines = undirected edges (RLCD could not orient them); graphs with undirected edges cannot be fitted, only the structure is drawn;
  - latent labels L1, L2, ... are arbitrary (chosen as the lexicographically smallest labeling, only to unify results across seeds) and do NOT mean "birth radius" or "guiding radius";
    each latent ellipse shows its correlations with the truth at the maximum-likelihood solution: Rb = |r(L, true birth radius)|, Rg = |r(L, true guiding radius)|;
  - numbers on the arrows = edge coefficients (standardized: how many standard deviations the child changes on average per standard deviation of the parent);
  - panel title: galaxy, structure id, number of seeds (out of 30) giving it, number of edges, Jacobian rank (coefficients are determined by the data only if it equals the number of edges),
    "prior ✗" marks a violation of the prior of paper Sec. 3.2 (age has an observed parent);
  - layout: layered in causal order (nodes without parents on top), unrelated to the positions in the authors' figure.
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
OUT_DATA = PLANB / 'output_data'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LABEL = {'logJZ': r'$J_z$', 'JPhi': r'$J_\phi$', 'ecc': r'$e$', 'age': 'Age', 'FeH': '[Fe/H]', 'MgH': '[Mg/H]'}
EW, EH = 0.085, 0.052
LATENT_FC, LATENT_EC = '#fde2e1', '#c2362f'
OBS_FC, OBS_EC = '#dbe8fb', '#2563eb'
INK, MUTED = '#1f2937', '#6b7280'


def layout(nodes, edges):
    """Layered layout: layer = length of the longest path from a parentless node; within a layer, sorted by the mean x of the parents (fewer crossings)"""
    par = {v: [a for a, b in edges if b == v] for v in nodes}
    depth = {}

    def d(v):
        if v not in depth:
            depth[v] = 0 if not par[v] else 1 + max(d(p) for p in par[v])
        return depth[v]
    for v in nodes:
        d(v)
    nl = max(depth.values()) + 1
    layers = [[v for v in nodes if depth[v] == k] for k in range(nl)]
    pos = {}
    for k, lay in enumerate(layers):
        if k > 0:
            lay.sort(key=lambda v: np.mean([pos[p][0] for p in par[v]]))
        y = 0.9 - 0.8 * k / max(nl - 1, 1)
        for i, v in enumerate(lay):
            pos[v] = ((i + 1) / (len(lay) + 1), y)
    return pos


def boundary(p, q):
    """Point where the ray from the centre of ellipse p towards q crosses the ellipse boundary"""
    dx, dy = q[0] - p[0], q[1] - p[1]
    t = 1 / np.sqrt((dx / EW) ** 2 + (dy / EH) ** 2)
    return p[0] + t * dx, p[1] + t * dy


def draw(ax, edges, coef, undirected, omega, ev, title):
    """edges: directed edges; coef: {edge: coefficient} (empty if not fitted); undirected: undirected edges (dashed, not oriented); ev: |r| of every latent"""
    lat = sorted({v for e in list(edges) + list(undirected) for v in e if v not in COLS})
    nodes = lat + COLS
    pos = layout(nodes, edges)
    vmax = max([abs(c) for c in coef.values()] + [1e-9]) if coef else 1
    nout = {}
    for (a, b) in edges:
        c = coef.get((a, b))
        j = nout.get(a, 0)
        nout[a] = j + 1
        pa, pb = pos[a], pos[b]
        same_layer = abs(pa[1] - pb[1]) < 1e-6
        rad = 0.35 if same_layer else 0.08
        s, e = boundary(pa, pb), boundary(pb, pa)
        lw = 0.6 + 2.4 * abs(c) / vmax if c is not None else 1.0
        ax.annotate('', xy=e, xytext=s, arrowprops=dict(arrowstyle='-|>', color=INK, lw=lw, shrinkA=0, shrinkB=0,
                    mutation_scale=9, connectionstyle=f'arc3,rad={rad}'))
        if c is None:
            continue
        # label position: control point of arc3 = midpoint + rad·(dy, −dx); take the point at parameter t on the quadratic Bezier curve;
        # for several edges from the same parent, t cycles through 0.3 / 0.47 / 0.64 so the labels do not overlap
        dx, dy = e[0] - s[0], e[1] - s[1]
        cpt = ((s[0] + e[0]) / 2 + rad * dy, (s[1] + e[1]) / 2 - rad * dx)
        tt = 0.3 + 0.17 * (j % 3)
        lx = (1 - tt) ** 2 * s[0] + 2 * tt * (1 - tt) * cpt[0] + tt ** 2 * e[0]
        ly = (1 - tt) ** 2 * s[1] + 2 * tt * (1 - tt) * cpt[1] + tt ** 2 * e[1]
        ax.text(lx, ly, f'{c:.2f}', fontsize=6.3, ha='center', va='center', color=INK, zorder=5,
                bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85))
    for (a, b) in undirected:
        s, e = boundary(pos[a], pos[b]), boundary(pos[b], pos[a])
        ax.annotate('', xy=e, xytext=s, arrowprops=dict(arrowstyle='-', color=LATENT_EC, lw=1.3, ls='--',
                    shrinkA=0, shrinkB=0, connectionstyle='arc3,rad=0.25'))
    for v in nodes:
        x, y = pos[v]
        is_lat = v not in COLS
        h = EH * (1.45 if is_lat and v in ev.index else 1)
        ax.add_patch(matplotlib.patches.Ellipse((x, y), 2 * EW * (1.12 if is_lat else 1), 2 * h,
                                                fc=LATENT_FC if is_lat else OBS_FC,
                                                ec=LATENT_EC if is_lat else OBS_EC, lw=1.1, zorder=6))
        if is_lat and v in ev.index:
            ax.text(x, y + 0.022, v, ha='center', va='center', fontsize=8.5, zorder=7, color=INK, fontweight='bold')
            ax.text(x, y - 0.02, f'Rb {ev.loc[v, "r_Rbir_best"]:.2f}\nRg {ev.loc[v, "r_Rg_best"]:.2f}',
                    ha='center', va='center', fontsize=6.3, zorder=7, color=LATENT_EC, linespacing=1.0)
        else:
            ax.text(x, y, v if is_lat else LABEL[v], ha='center', va='center', fontsize=8, zorder=7, color=INK)
        if not is_lat and omega.get('X_' + v, omega.get(v, 1)) < 1e-3:
            ax.text(x, y - EH - 0.012, 'ω≈0', ha='center', va='top', fontsize=6.5, color=LATENT_EC, fontweight='bold', zorder=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.03, 1.0)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=8.2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--rank', type=int, default=1, help='1 = most frequent structure of each galaxy, 2 = second most frequent, ...')
    args = ap.parse_args()
    st = pd.read_csv(OUT_DATA / f'structures_{args.version}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    best = json.load(open(OUT_DATA / f'eval_{args.version}_best.json'))
    ev = pd.read_csv(OUT_DATA / f'eval_{args.version}.csv')
    ev['galaxy'] = ev.galaxy.astype(str)
    fig, axs = plt.subplots(3, 4, figsize=(20, 16))
    for ax, gid in zip(axs.flat, IDS):
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        if len(freq) < args.rank:
            ax.axis('off')
            ax.set_title(f'{gid}: only {len(freq)} structure(s) found', fontsize=8.2)
            continue
        sid, cnt = freq.index[args.rank - 1], freq.iloc[args.rank - 1]
        row = g[g.structure_id == sid].iloc[0]
        info = best.get(gid, {}).get(str(sid))
        head = f'{gid}   structure {sid}: {cnt}/{len(g)} seeds, {len(freq)} distinct structures'
        prior = '' if row.prior_ok else '   prior ✗'
        edges = [tuple(x) for x in json.loads(row.structure)]
        und = [tuple(x) for x in json.loads(row.undirected)]
        if info is None:          # undirected edges (not oriented) or not evaluated: structure only, no coefficients
            why = (f'{len(und)} undirected edge(s) (dashed red): direction not determined → no coefficients'
                   if row.n_undirected > 0 else 'not evaluated (freq < 20%)')
            draw(ax, edges, {}, und, {}, pd.DataFrame(), head + f'\n{why}{prior}')
            continue
        e = ev[(ev.galaxy == gid) & (ev.structure_id == sid)].set_index('latent')
        title = (head + f'\n{info["n_edges"]} edges, Jacobian rank {info["jacobian_rank"]}/{info["n_edges"]}'
                 f'{" (identifiable)" if info["jacobian_rank"] == info["n_edges"] else " (NOT identifiable)"}{prior}')
        coef = {tuple(k): c for k, c in zip(info['edges'], info['coef'])}
        draw(ax, edges, coef, [], info['omega'], e, title)
    rk = '' if args.rank == 1 else f' (rank-{args.rank} most frequent structure)'
    fig.suptitle(f'Plan B — TNG50, 7<R<10 kpc, |z|<0.3 kpc ({args.version}): causal graph learned from scratch by RLCD in each galaxy'
                 f'{rk}\nmost frequent structure over 30 hash seeds; edge numbers = maximum-likelihood coefficients; '
                 'inside each latent: |r| with the TRUE birth radius (Rb) and guiding radius (Rg).  Latent names are arbitrary.',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    stem = HERE.stem + ('' if args.version == 'nonoise' else '_noise') + ('' if args.rank == 1 else f'_rank{args.rank}')
    out = HERE.with_name(stem + '.png')
    fig.savefig(out, dpi=110)
    print(f'saved {out}')


if __name__ == '__main__':
    main()
