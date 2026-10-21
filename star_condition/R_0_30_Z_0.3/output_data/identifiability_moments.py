# -*- coding: utf-8 -*-
"""
identifiability_moments.py — Approach A (authors' causal graph frozen) + strict identifiability analysis for the selection R_0_30_Z_0.3 (step 3)

Method identical to scripts/local/rlcd_identifiability.py + rlcd_identifiability_ci.py (REPORT.md Sec. 3.9),
but the input is the covariance matrix computed by compute_moments.py (number of stars per galaxy in input_data/selection.txt; star-by-star computation is unnecessary):
  1. Jacobian rank: rank of the Jacobian of "15 edge coefficients → 15 observed correlations" at 3 random coefficient sets (14 expected for the authors' graph).
  2. first round: 30 random starts (seeds 1000+s), L-BFGS-B maximum likelihood, soft constraint noise variance >= 0.01 (as in the original script);
  3. strict convergence: these 30 solutions + 200 new random starts (seeds 5000+i), all strictly converged (polished) with noise variances >= 0 (penalty 1e7);
  4. statistically indistinguishable from the best: 2·N·ΔNLL < 3.84 (95% likelihood-ratio test);
       N: number of stars;  ΔNLL: NLL difference from the best solution (NLL averaged per star);  3.84: 95% quantile of the χ² distribution with 1 degree of freedom
  5. |r(L1, true Rbirth)| and |r(L2, true Rg)| of each solution: L1, L2 from the hierarchical least squares of the authors' notebook are
     linear combinations L = w·z of the standardized inputs z, so
         r(L, R) = (w·c) / sqrt(w·C·w · var(R))
     meaning: correlation = covariance of L and R ÷ (std of L × std of R), all computable from the covariance matrix.
       w: weights of L on the 6 standardized inputs; c: covariances of the standardized inputs with R; C: correlation matrix of the standardized inputs; var(R): variance of R
     (--selftest checks this on the earlier 7<R<10 data: identical to the star-by-star computation)
  6. baselines: |r| of age, [Fe/H], [Mg/H], J_phi alone with Rbirth; linear-regression ceiling of the 6 inputs sqrt(c·C⁻¹·c / var(R))
     (meaning: the correlation reached by the best linear combination of the 6 inputs predicting the birth radius; c, C, var(R) as above).
     ⚠ This selection covers the whole disk (0–30 kpc, |z|<0.3 kpc), where the current angular momentum J_phi (≈ current radius) is already strongly correlated with the birth radius (mostly 0.8–0.9),
       so L1 must be compared with "J_phi alone"; the paper uses only stars at 7–10 kpc precisely to remove this "trivial" information.

Run (numpy, scipy, pandas, matplotlib; analysis environment):
    cd <repo root>
    $PY star_condition/R_0_30_Z_0.3/output_data/identifiability_moments.py --version nonoise
    $PY star_condition/R_0_30_Z_0.3/output_data/identifiability_moments.py --version noise
    (about 5–10 minutes per galaxy; --ids 538905 runs a single galaxy; --extra 200 is the number of new starts)
Outputs (this folder):
  <version>_<ID>_ci.csv      all valid strictly converged solutions: start (start index, <30 = first round), nll, lr = 2·N·ΔNLL, omega_min,
                             omega_min_node (node with the smallest noise variance), abs_r_L1_Rbir, abs_r_L2_Rg, the 15 edge coefficients
  <version>_<ID>_summary.txt Jacobian rank, non-identifiable direction, range among indistinguishable solutions
  summary_ci_<version>.csv/.txt/.png   summary of the 12 galaxies (strict ranges + baselines)
"""
import os
import sys
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve()
COND = HERE.parents[1]
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'scripts' / 'local'))
from rlcd_identifiability import setup, implied, nll  # noqa: E402

DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LR95 = 3.84
PEN = 1e7
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def standardized_moments(m):
    """From the 8×8 covariance: C (correlation matrix of the 6 inputs), c_R / c_G (covariances of the standardized inputs with R_bir / Rg), var_R, var_G"""
    cov = np.array(m['cov'])
    sd = np.sqrt(np.diag(cov)[:6])
    C = cov[:6, :6] / np.outer(sd, sd)
    c_R, c_G = cov[:6, 6] / sd, cov[:6, 7] / sd
    return C, c_R, c_G, cov[6, 6], cov[7, 7]


def latent_weights(coef):
    """Hierarchical least squares of the authors' notebook, written as weights on the 6 standardized inputs. coef: {(parent, child): coefficient}, names without X_"""
    w = {c: np.eye(6)[i] for i, c in enumerate(COLS)}

    def solve(L):
        num, den = np.zeros(6), 0.0
        for (p, c), a in coef.items():
            if p != L:
                continue
            r = w[c].copy()
            for (q, c2), b in coef.items():
                if c2 == c and q != L:
                    r = r - b * w[q]
            num, den = num + a * r, den + a * a
        return num / den
    w['L2'] = solve('L2')
    w['L1'] = solve('L1')
    return w['L1'], w['L2']


def corr_from_weights(wv, C, c, var):
    return abs(wv @ c) / np.sqrt((wv @ C @ wv) * var)


def analyse(tag, m, names, idx, edges, order, extra):
    n, ne = len(names), len(edges)
    xi = [idx['X_' + c] for c in COLS]
    iu = np.triu_indices(6, 1)
    lab = [f'{names[a]}->{names[b]}' for a, b in edges]
    N = m['N']
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    Semp = C
    lines = [f'{tag}: {N} stars; frozen graph {DOT.name}; {ne} edges, {len(iu[0])} observed correlations',
             f'theoretical minimum NLL = 0.5*logdet(S)+3 = {0.5*np.linalg.slogdet(Semp)[1]+3:.4f}']

    # 1. Jacobian rank
    cv = lambda th: implied(th, n, edges, order)[1][np.ix_(xi, xi)][iu]
    rng = np.random.default_rng(0)
    for trial in range(3):
        th = np.tanh(rng.normal(0, 0.5, ne))
        J = np.array([(cv(th + 1e-6 * np.eye(ne)[k]) - cv(th - 1e-6 * np.eye(ne)[k])) / 2e-6 for k in range(ne)]).T
        sv = np.linalg.svd(J, compute_uv=False)
        lines.append(f'Jacobian rank (random point {trial+1}): {int(np.sum(sv > 1e-8 * sv[0]))} / {ne}')
    null = np.linalg.svd(J)[2][-1]
    lines.append('non-identifiable direction (changing the coefficients along it leaves the observed correlation matrix unchanged), largest components first:')
    for l, x in sorted(zip(lab, null), key=lambda z: -abs(z[1]))[:6]:
        lines.append(f'    {l:18s} {x:+.3f}')

    # 2. first-round multistart (soft constraint, as in rlcd_identifiability.py)
    starts = []
    for s in range(30):
        th0 = np.tanh(np.random.default_rng(1000 + s).normal(0, np.sqrt(2 / n), ne))
        r = minimize(nll, th0, args=(n, edges, order, xi, Semp), method='L-BFGS-B',
                     bounds=[(-0.999, 0.999)] * ne, options={'maxiter': 5000})
        F, S, om = implied(r.x, n, edges, order)
        if om.min() > -1e-6 and np.linalg.slogdet(S[np.ix_(xi, xi)])[0] > 0:
            starts.append(r.x)
    n_first = len(starts)
    starts += [np.tanh(np.random.default_rng(5000 + i).normal(0, 0.5, ne)) for i in range(extra)]

    # 3. strict convergence, noise variances >= 0
    def nll_phys(th):
        _, _, om = implied(th, n, edges, order)
        return nll(th, n, edges, order, xi, Semp, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

    rows = []
    for i, th0 in enumerate(starts):
        r = minimize(nll_phys, th0, method='L-BFGS-B', bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
        F, S, om = implied(r.x, n, edges, order)
        if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
            continue
        coef = {(names[a].replace('X_', ''), names[b].replace('X_', '')): v for (a, b), v in zip(edges, r.x)}
        w1, w2 = latent_weights(coef)
        rows.append(dict(start=i, nll=nll(r.x, n, edges, order, xi, Semp, pen=0), omega_min=om.min(),
                         omega_min_node=names[int(np.argmin(om))],
                         abs_r_L1_Rbir=corr_from_weights(w1, C, c_R, var_R),
                         abs_r_L2_Rg=corr_from_weights(w2, C, c_G, var_G),
                         **{l: v for l, v in zip(lab, r.x)}))
    ci = pd.DataFrame(rows).sort_values('nll').reset_index(drop=True)
    ci.insert(2, 'lr', 2 * N * (ci.nll - ci.nll.iloc[0]))
    ci.to_csv(HERE.with_name(f'{tag}_ci.csv'), index=False)
    inn = ci[ci.lr < LR95]
    first = ci[ci.start < n_first]
    lines += [f'{len(ci)} valid solutions after strict convergence ({n_first} from the first round + {extra} new starts), best NLL {ci.nll.iloc[0]:.6f}',
              f'{len(inn)} solutions statistically indistinguishable from the best (2·N·ΔNLL < {LR95}): '
              f'L1->age from {inn["L1->X_age"].min():+.2f} to {inn["L1->X_age"].max():+.2f}',
              f'their |r(L1, Rbirth)| range {inn.abs_r_L1_Rbir.min():.3f} – {inn.abs_r_L1_Rbir.max():.3f}; '
              f'|r(L2, Rg)| range {inn.abs_r_L2_Rg.min():.3f} – {inn.abs_r_L2_Rg.max():.3f}',
              f'best solution of the first round only: NLL {first.nll.iloc[0]:.6f}, |r(L1, Rbirth)| {first.abs_r_L1_Rbir.iloc[0]:.3f}, '
              f'|r(L2, Rg)| {first.abs_r_L2_Rg.iloc[0]:.3f}',
              f'smallest noise variance at the best solution {ci.omega_min.iloc[0]:.4f} ({ci.omega_min_node.iloc[0]}; ≈0 means this observable is treated as an exact linear combination of its parents)']
    txt = '\n'.join(lines)
    print(txt + '\n', flush=True)
    open(HERE.with_name(f'{tag}_summary.txt'), 'w').write(txt + '\n')


def baselines(m):
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    r = lambda i: abs(c_R[i]) / np.sqrt(var_R)
    ceil = np.sqrt(c_R @ np.linalg.solve(C, c_R) / var_R)
    ceil_g = np.sqrt(c_G @ np.linalg.solve(C, c_G) / var_G)
    return dict(r_age_Rbir=r(3), r_FeH_Rbir=r(4), r_MgH_Rbir=r(5), r_JPhi_Rbir=r(1),
                r_best_single_Rbir=max(r(i) for i in range(6)), best_single=COLS[int(np.argmax([r(i) for i in range(6)]))],
                r_ceiling_Rbir=ceil, r_ceiling_Rg=ceil_g,
                r_JPhi_Rg=abs(c_G[1]) / np.sqrt(var_G), std_Rbir=np.sqrt(var_R))


def write_summary(version, moments):
    rows = []
    for gid in IDS:
        f = HERE.with_name(f'{version}_{gid}_ci.csv')
        if not f.exists():
            continue
        c = pd.read_csv(f)
        inn = c[c.lr < LR95]
        rows.append(dict(galaxy=gid, n_stars=moments[gid]['N'], n_valid=len(c), n_indistinguishable=len(inn),
                         r1_best=c.abs_r_L1_Rbir.iloc[0], r1_min=inn.abs_r_L1_Rbir.min(), r1_max=inn.abs_r_L1_Rbir.max(),
                         r2_best=c.abs_r_L2_Rg.iloc[0], r2_min=inn.abs_r_L2_Rg.min(), r2_max=inn.abs_r_L2_Rg.max(),
                         omega_min_at_best=c.omega_min.iloc[0], omega_min_node=c.omega_min_node.iloc[0],
                         **baselines(moments[gid])))
    tab = pd.DataFrame(rows)
    tab.to_csv(HERE.with_name(f'summary_ci_{version}.csv'), index=False)
    head = (f'selection R_0_30_Z_0.3 (0<R<30 kpc, |z|<0.3 kpc), {version}; authors\' causal graph frozen\n'
            f'r1/r2_best: |r(L1,Rbirth)| / |r(L2,Rg)| of the maximum-likelihood solution;  r1/r2_min–max: range among solutions statistically indistinguishable from the best (2·N·ΔNLL<{LR95})\n'
            'r_age / r_FeH / r_MgH / r_JPhi _Rbir: |r| of a single input with the birth radius;  r_best_single: the best single input (best_single says which)\n'
            'r_ceiling: linear-regression ceiling of the 6 inputs;  std_Rbir: standard deviation of the birth radius [kpc]\n'
            '⚠ this selection covers the whole disk (0–30 kpc, |z|<0.3 kpc); the current angular momentum J_phi (≈ current radius) is already strongly correlated with the birth radius, '
            'so compare with "J_phi alone", not only with age and [Fe/H]\n\n')
    txt = head + tab.round(3).to_string(index=False) + '\n'
    open(HERE.with_name(f'summary_ci_{version}.txt'), 'w').write(txt)
    print(txt)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(tab) + 2.2))
    y = np.arange(len(tab))[::-1]
    for yi, (_, r) in zip(y, tab.iterrows()):
        ax.plot([r.r1_min, max(r.r1_max, r.r1_min + 0.006)], [yi] * 2, color=BLUE, lw=6, solid_capstyle='round')
        ax.plot(r.r1_best, yi, marker='|', color='white', ms=9, mew=2)
        ax.plot(r.r_age_Rbir, yi, marker='v', color=INK, ms=6, ls='none')
        ax.plot(r.r_FeH_Rbir, yi, marker='x', color=INK, ms=6, ls='none')
        ax.plot(r.r_JPhi_Rbir, yi, marker='s', mfc='#f59e0b', mec=INK, ms=6, ls='none')
        ax.plot(r.r_ceiling_Rbir, yi, marker='o', mfc='none', mec=INK, ms=6, ls='none')
    ax.axvline(0.64, color=MUTED, ls='--', lw=1)
    ax.text(0.645, y.max() + 0.55, 'paper 0.64', color=MUTED, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(tab.galaxy, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(-2.6, len(tab) - 0.2)
    ax.set_xlabel('|Pearson r| with the true birth radius')
    ax.plot([], [], color=BLUE, lw=6, label='L1, statistically indistinguishable fits (white tick: best fit)')
    ax.plot([], [], marker='v', color=INK, ls='none', label='age alone')
    ax.plot([], [], marker='x', color=INK, ls='none', label='[Fe/H] alone')
    ax.plot([], [], marker='s', mfc='#f59e0b', mec=INK, ls='none', label=r'$J_\phi$ alone (≈ present radius)')
    ax.plot([], [], marker='o', mfc='none', mec=INK, ls='none', label='linear regression on all 6 inputs')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f'R_0_30_Z_0.3 ({version}): how well can L1 track the birth radius?\n(author graph frozen)',
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(HERE.with_name(f'summary_ci_{version}.png'), dpi=130)
    plt.close(fig)


def selftest():
    """Check on the earlier 7<R<10 data: |r| from the covariance method equals the star-by-star computation (<tag>_ci.csv)"""
    ci = pd.read_csv(ROOT / 'results' / 'rlcd_A' / 'identifiability' / 'nonoise_538905_ci.csv').sort_values('nll')
    e = pd.read_csv(ROOT / 'data' / 'rlcd_input' / 'nonoise' / '538905_e.csv')[COLS]
    t = pd.read_csv(ROOT / 'data' / 'rlcd_input' / 'nonoise' / '538905_truth.csv')
    X = np.column_stack([e.to_numpy(float), t[['R_bir', 'Rg']].to_numpy(float)])
    m = dict(N=len(X), cov=np.cov(X.T, ddof=1).tolist())
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    for k in range(3):
        row = ci.iloc[k]
        coef = {tuple(x.replace('X_', '') for x in col.split('->')): row[col] for col in ci.columns if '->' in col}
        w1, w2 = latent_weights(coef)
        print(f'solution {k}: star-by-star |r(L1,Rb)| {row.abs_r_L1_Rbir:.6f}  covariance method {corr_from_weights(w1, C, c_R, var_R):.6f};  '
              f'star-by-star |r(L2,Rg)| {row.abs_r_L2_Rg:.6f}  covariance method {corr_from_weights(w2, C, c_G, var_G):.6f}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--extra', type=int, default=200)
    ap.add_argument('--summary-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    moments = json.load(open(HERE.with_name(f'moments_{args.version}.json')))
    if not args.summary_only:
        names, idx, edges, order = setup(str(DOT))
        for gid in args.ids or IDS:
            analyse(f'{args.version}_{gid}', moments[gid], names, idx, edges, order, args.extra)
    write_summary(args.version, moments)


if __name__ == '__main__':
    main()
