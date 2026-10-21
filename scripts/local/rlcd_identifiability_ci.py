# -*- coding: utf-8 -*-
"""
rlcd_identifiability_ci.py — strict version of the "possible range": 95% profile-likelihood region (complements rlcd_identifiability.py)

Why this script is needed:
  rlcd_identifiability.py defines "equally good solutions" as NLL within 0.01 of the best. But the NLL is averaged per star,
  and for 7 thousand to 120 thousand stars a difference of 0.01 is statistically distinguishable. The proper criterion is a likelihood-ratio test:
      2·N·ΔNLL < 3.84
  meaning: the difference in fit quality of two coefficient sets is too small to say, at 95% confidence, that one is better.
      N: number of stars;  ΔNLL: difference from the best NLL (NLL averaged per star);  3.84: 95% quantile of the χ² distribution with 1 degree of freedom
  Also, Jacobian rank 14 means there is a direction of exactly equal fit; but the optimum usually sits on the boundary "some noise variance = 0"
  (a Heywood case in factor analysis), and moving along that direction quickly makes a noise variance negative, so how far one can move depends on how fast the fit degrades.

Procedure:
  1. re-optimize every valid solution in <tag>_multistart.csv (the 30 random starts) with very strict convergence (polish),
     removing NLL differences caused by unconverged optimizers; --extra adds more random starts (g2.79e12 used --extra 200,
     because its likelihood surface has many close but statistically distinguishable local optima);
  2. keep the solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84);
  3. report the range of |r(L1, true Rbirth)| and |r(L2, true Rg)| among them.
  This is the range "found so far": starts are finite; if a better optimum exists that was not found, the best solution and this range can change (g2.79e12 is an example).
  Constraints: coefficients within ±0.999; every node's noise variance >= 0 (physical requirement; rlcd_identifiability.py used a soft >= 0.01,
  which hardly matters for TNG50, but the optimum of g2.79e12 needs a noise variance of [Fe/H] ≈ 0.0005, hence >= 0 here).

Run (after rlcd_identifiability.py):
    $PY scripts/local/rlcd_identifiability_ci.py --version g279 --extra 200
    $PY scripts/local/rlcd_identifiability_ci.py --version nonoise --ids 538905
    $PY scripts/local/rlcd_identifiability_ci.py --version nonoise      # 12 galaxies
Outputs: results/rlcd_A/identifiability/<tag>_ci.csv (every polished solution), <tag>_ci.txt; summary summary_ci_<version>.csv/.txt/.png
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlcd_identifiability import COLS, DOT, OUT, load, setup, implied, nll, latents  # noqa: E402

IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
LR95 = 3.84      # 95% quantile of the χ² distribution with 1 degree of freedom
PEN = 1e7        # penalty coefficient for negative noise variances
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['g279', 'nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--dot', default=DOT)
    ap.add_argument('--summary-only', action='store_true', help='only regenerate the summary tables')
    ap.add_argument('--extra', type=int, default=0, help='number of additional random starts (optimized directly with the strict settings)')
    args = ap.parse_args()
    ids = ['g279'] if args.version == 'g279' else (args.ids or IDS)
    if args.summary_only:
        write_summary(args.version)
        return
    names, idx, edges, order = setup(args.dot)
    n, ne = len(names), len(edges)
    xi = [idx['X_' + c] for c in COLS]
    lab = [f'{names[a]}->{names[b]}' for a, b in edges]

    for gid in ids:
        tag = gid if args.version == 'g279' else f'{args.version}_{gid}'
        X, T = load(args.version, gid)
        N = len(X)
        Semp = X[COLS].cov().to_numpy()
        ms = pd.read_csv(os.path.join(OUT, f'{tag}_multistart.csv'))
        ms = ms[ms.valid]

        # 1. polish (physical constraint: every node noise variance >= 0, a large penalty approximates the hard constraint)
        def nll_phys(th):
            _, _, om = implied(th, n, edges, order)
            return nll(th, n, edges, order, xi, Semp, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

        starts = [row[lab].to_numpy(float) for _, row in ms.iterrows()]
        starts += [np.tanh(np.random.default_rng(5000 + i).normal(0, 0.5, ne)) for i in range(args.extra)]
        pol = []
        for i_start, th_start in enumerate(starts):
            r = minimize(nll_phys, th_start, method='L-BFGS-B',
                         bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
            F, S, om = implied(r.x, n, edges, order)
            if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
                continue
            pol.append((nll(r.x, n, edges, order, xi, Semp, pen=0), r.x, i_start))
        pol.sort(key=lambda z: z[0])
        best = pol[0][0]

        rows = []
        for v, th, i_start in pol:
            F, S, om = implied(th, n, edges, order)
            lat = latents(X, names, idx, F)
            rows.append(dict(start=i_start, nll=v, lr=2 * N * (v - best), omega_min=om.min(), omega_min_node=names[int(np.argmin(om))],
                             abs_r_L1_Rbir=abs(np.corrcoef(lat['L1'], T['R_bir'])[0, 1]),
                             abs_r_L2_Rg=abs(np.corrcoef(lat['L2'], T['Rg'])[0, 1]),
                             **{l: val for l, val in zip(lab, th)}))

        ci = pd.DataFrame(rows).sort_values('nll')
        ci.to_csv(os.path.join(OUT, f'{tag}_ci.csv'), index=False)
        inn = ci[ci.lr < LR95]
        lines = [f'{tag}: {N} stars; {len(ci)} valid solutions after polishing, best NLL {best:.6f}',
                 f'{len(inn)} solutions statistically indistinguishable from the best (2·N·ΔNLL < {LR95}): '
                 f'L1->age from {inn["L1->X_age"].min():+.2f} to {inn["L1->X_age"].max():+.2f}',
                 f'their |r(L1, Rbirth)| range {inn.abs_r_L1_Rbir.min():.3f} – {inn.abs_r_L1_Rbir.max():.3f}; '
                 f'|r(L2, Rg)| range {inn.abs_r_L2_Rg.min():.3f} – {inn.abs_r_L2_Rg.max():.3f}',
                 (f'best solution using only the {len(ms)} valid first-round multistart solutions: NLL {ci[ci.start < len(ms)].nll.iloc[0]:.6f}, '
                  f'|r(L1, Rbirth)| {ci[ci.start < len(ms)].abs_r_L1_Rbir.iloc[0]:.3f}, |r(L2, Rg)| {ci[ci.start < len(ms)].abs_r_L2_Rg.iloc[0]:.3f}'
                  if args.extra > 0 else '(no extra starts)'),
                 f'smallest noise variance at the best solution {ci.omega_min.iloc[0]:.4f} ({ci.omega_min_node.iloc[0]}; ≈0 means the optimum sits on the boundary "noise variance ≥ 0", i.e. that observable is almost fully determined by its parents)']
        txt = '\n'.join(lines)
        print(txt + '\n')
        open(os.path.join(OUT, f'{tag}_ci.txt'), 'w').write(txt + '\n')

    if args.version != 'g279':
        write_summary(args.version)


def write_summary(version):
    """Read existing <tag>_ci.csv files and summarize them into summary_ci_<version>.csv/.txt (g2.79e12 first)"""
    rows = []
    for tag in ['g279'] + [f'{version}_{g}' for g in IDS]:
        f = os.path.join(OUT, f'{tag}_ci.csv')
        if not os.path.exists(f):
            continue
        c = pd.read_csv(f).sort_values('nll')
        inn = c[c.lr < LR95]
        rows.append(dict(tag=tag, n_polished=len(c), n_indistinguishable=len(inn),
                         r1_min=inn.abs_r_L1_Rbir.min(), r1_max=inn.abs_r_L1_Rbir.max(),
                         r2_min=inn.abs_r_L2_Rg.min(), r2_max=inn.abs_r_L2_Rg.max(),
                         best_nll=c.nll.iloc[0], omega_min_at_best=c.omega_min.iloc[0],
                         omega_min_node=c.omega_min_node.iloc[0]))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, f'summary_ci_{version}.csv'), index=False)
    head = (f'Strict criterion: after polishing, solutions statistically indistinguishable from the best (2·N·ΔNLL < {LR95}, 95% likelihood-ratio test; noise variances >= 0); '
            f'range of |r(L1, true Rbirth)| and |r(L2, true Rg)| (finite starts: if a better optimum exists that was not found, the best solution and the range can change)\n'
            f'omega_min_node: node with the smallest noise variance at the best solution (≈0 means it is almost fully determined by its parents)\n\n')
    txt = head + tab.round(3).to_string(index=False) + '\n'
    open(os.path.join(OUT, f'summary_ci_{version}.txt'), 'w').write(txt)
    print(txt)

    # figure: one row per galaxy; thick blue line = strict range of |r(L1,Rb)|, thin grey line = loose range (ΔNLL<0.01),
    #     black triangle = age alone |r(age,Rb)|, cross = [Fe/H] alone, open circle = linear-regression ceiling of the 6 inputs
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
    loose = os.path.join(OUT, f'summary_{version}.csv')
    loose = pd.read_csv(loose).set_index('tag') if os.path.exists(loose) else None
    ceil = os.path.join(OUT, 'linear_ceiling.csv')
    ceil = pd.read_csv(ceil).set_index('tag') if os.path.exists(ceil) else None
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(tab) + 1.8))
    y = np.arange(len(tab))[::-1]
    for yi, (_, r) in zip(y, tab.iterrows()):
        if loose is not None and r.tag in loose.index:
            ax.plot([loose.loc[r.tag, 'r1_min'], loose.loc[r.tag, 'r1_max']], [yi] * 2, color=MUTED, lw=1.5)
        x0, x1 = r.r1_min, max(r.r1_max, r.r1_min + 0.006)
        ax.plot([x0, x1], [yi] * 2, color=BLUE, lw=6, solid_capstyle='round')
        if ceil is not None and r.tag in ceil.index:
            ax.plot(ceil.loc[r.tag, 'r_age_Rbir'], yi, marker='v', color=INK, ms=6, ls='none')
            ax.plot(ceil.loc[r.tag, 'r_FeH_Rbir'], yi, marker='x', color=INK, ms=6, ls='none')
            ax.plot(ceil.loc[r.tag, 'r_ceiling_Rbir'], yi, marker='o', mfc='none', mec=INK, ms=6, ls='none')
    ax.axvline(0.64, color=MUTED, ls='--', lw=1)
    ax.text(0.645, y.max() + 0.55, 'paper 0.64', color=MUTED, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace(f'{version}_', '') for t in tab.tag], fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel('|Pearson r| with the true birth radius')
    ax.plot([], [], color=BLUE, lw=6, label='L1, statistically indistinguishable fits')
    ax.plot([], [], color=MUTED, lw=1.5, label='L1, loose criterion (dNLL<0.01)')
    ax.plot([], [], marker='v', color=INK, ls='none', label='age alone')
    ax.plot([], [], marker='x', color=INK, ls='none', label='[Fe/H] alone')
    ax.plot([], [], marker='o', mfc='none', mec=INK, ls='none', label='linear regression on all 6 inputs')
    ax.set_ylim(-2.6, len(tab) - 0.2)
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f'How well can L1 track the birth radius? ({version}, author graph frozen)', fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f'summary_ci_{version}.png'), dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
