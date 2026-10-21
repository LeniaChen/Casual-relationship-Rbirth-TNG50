# -*- coding: utf-8 -*-
"""
rlcd_identifiability.py — can the edge coefficients be determined uniquely from the data when the authors' causal graph is frozen? (identifiability)

Needs only numpy, scipy, pandas, matplotlib (either environment works):
    cd <repo root>
    $PY scripts/local/rlcd_identifiability.py --version g279
    $PY scripts/local/rlcd_identifiability.py --version nonoise --ids 538905

Three things are done:
  1. direct identifiability test: at random coefficient values, the rank of the Jacobian of the map "15 edge coefficients -> 15 observed correlations".
     rank = 15: locally identifiable (coefficients fixed by the data); rank < 15: there is at least one direction along which the coefficients can change without changing the observed correlation matrix at all
  2. maximum-likelihood fits with L-BFGS from many random starts (truly converged), to see whether the optimum is unique
  3. profile: fix L1->age at a series of values, refit the other 14 coefficients, plot the NLL and |r(L1, true Rbirth)| against it
The model is exactly that of parameter_by_trek in the repo: all variables standardized to unit variance, latent variance 1, covariance from the trek rule (equivalent to layer-by-layer recursion);
in addition every node's implied noise variance is required to be >= 0 (the repo code does not check this).

Outputs: results/rlcd_A/identifiability/<ID>_{profile.csv, multistart.csv, profile.png, summary.txt}
"""
import os
import re
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize

COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
REPO = 'Causal-Galactic-Archaeology'
DOT = 'results/rlcd_g279/author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
OUT = 'results/rlcd_A/identifiability'


def load(version, gid):
    if version == 'g279':
        e = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_e.csv'))
        t = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_multi_ele_action_Jin_SN.csv'))
    else:
        e = pd.read_csv(f'data/rlcd_input/{version}/{gid}_e.csv')
        t = pd.read_csv(f'data/rlcd_input/{version}/{gid}_truth.csv')
    ok = np.isfinite(e[COLS].to_numpy()).all(axis=1)
    e, t = e.loc[ok, COLS].reset_index(drop=True), t.loc[ok].reset_index(drop=True)
    return (e - e.mean()) / e.std(), t


def setup(dotfile):
    E = re.findall(r'^(\w+)\s*->\s*(\w+)', open(dotfile).read(), re.M)
    names = sorted({a for a, _ in E} | {b for _, b in E})
    idx = {v: i for i, v in enumerate(names)}
    edges = [(idx[a], idx[b]) for a, b in E]
    order, rem = [], list(range(len(names)))
    while rem:
        for k in rem:
            if all((p, k) not in edges or p in order for p in range(len(names))):
                order.append(k)
                rem.remove(k)
                break
    return names, idx, edges, order


def implied(theta, n, edges, order):
    F = np.zeros((n, n))
    for (a, b), w in zip(edges, theta):
        F[a, b] = w
    S = np.zeros((n, n))
    om = np.ones(n)
    for pos, k in enumerate(order):
        prev = order[:pos]
        par = [p for p in prev if F[p, k] != 0]
        if par:
            Fp = F[par, k]
            S[prev, k] = S[np.ix_(prev, par)] @ Fp
            S[k, prev] = S[prev, k]
            om[k] = 1 - Fp @ S[np.ix_(par, par)] @ Fp
        S[k, k] = 1
    return F, S, om


def nll(theta, n, edges, order, xi, Semp, pen=1e3):
    F, S, om = implied(theta, n, edges, order)
    Sx = S[np.ix_(xi, xi)]
    sg, ld = np.linalg.slogdet(Sx)
    bad = np.sum(np.clip(0.01 - om, 0, None) ** 2)
    if sg <= 0:
        return 1e6 + pen * bad
    return 0.5 * ld + 0.5 * np.trace(np.linalg.solve(Sx, Semp)) + pen * bad


def latents(X, names, idx, F):
    val = {f'X_{c}': X[c].to_numpy() for c in COLS}
    lat = [v for v in names if not v.startswith('X_')]
    done = {}
    while len(done) < len(lat):
        for L in lat:
            if L in done:
                continue
            ch = [c for c in names if F[idx[L], idx[c]] != 0]
            if any(c in lat and c not in done for c in ch):
                continue
            num, den = 0.0, 0.0
            for c in ch:
                r = (val[c] if c in val else done[c]).copy()
                for p in names:
                    if p != L and F[idx[p], idx[c]] != 0:
                        r = r - F[idx[p], idx[c]] * (val[p] if p in val else done[p])
                a = F[idx[L], idx[c]]
                num, den = num + a * r, den + a * a
            done[L] = num / den
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['g279', 'nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--nstart', type=int, default=40)
    ap.add_argument('--dot', default=DOT)
    args = ap.parse_args()
    ids = ['g279'] if args.version == 'g279' else args.ids
    os.makedirs(OUT, exist_ok=True)
    names, idx, edges, order = setup(args.dot)
    n = len(names)
    xi = [idx['X_' + c] for c in COLS]
    iu = np.triu_indices(6, 1)
    lab = [f'{names[a]}->{names[b]}' for a, b in edges]

    for gid in ids:
        X, T = load(args.version, gid)
        Semp = X[COLS].cov().to_numpy()
        tag = gid if args.version == 'g279' else f'{args.version}_{gid}'
        lines = [f'{tag}: {len(X)} stars; frozen graph {args.dot}; {len(edges)} edges, {len(iu[0])} observed correlations']
        lines.append(f'theoretical minimum NLL = 0.5*logdet(S)+3 = {0.5*np.linalg.slogdet(Semp)[1]+3:.4f}')

        # 1. Jacobian rank
        cv = lambda th: implied(th, n, edges, order)[1][np.ix_(xi, xi)][iu]
        rng = np.random.default_rng(0)
        for trial in range(3):
            th = np.tanh(rng.normal(0, 0.5, len(edges)))
            J = np.array([(cv(th + 1e-6 * np.eye(len(edges))[k]) - cv(th - 1e-6 * np.eye(len(edges))[k])) / 2e-6
                          for k in range(len(edges))]).T
            sv = np.linalg.svd(J, compute_uv=False)
            rank = int(np.sum(sv > 1e-8 * sv[0]))
            lines.append(f'Jacobian rank (random point {trial+1}): {rank} / {len(edges)}')
        null = np.linalg.svd(J)[2][-1]
        lines.append('non-identifiable direction (changing the coefficients along it leaves the observed correlation matrix unchanged), largest components first:')
        for l, x in sorted(zip(lab, null), key=lambda z: -abs(z[1]))[:8]:
            lines.append(f'    {l:18s} {x:+.3f}')

        # 2. multistart ML
        rows = []
        for s in range(args.nstart):
            r0 = np.random.default_rng(1000 + s)
            th0 = np.tanh(r0.normal(0, np.sqrt(2 / n), len(edges)))
            r = minimize(nll, th0, args=(n, edges, order, xi, Semp), method='L-BFGS-B',
                         bounds=[(-0.999, 0.999)] * len(edges), options={'maxiter': 5000})
            F, S, om = implied(r.x, n, edges, order)
            Sx = S[np.ix_(xi, xi)]
            valid = bool(om.min() > -1e-6 and np.linalg.slogdet(Sx)[0] > 0)
            lat = latents(X, names, idx, F)
            rows.append(dict(start=s, valid=valid, nll=nll(r.x, n, edges, order, xi, Semp, pen=0),
                             rmse=float(np.sqrt(np.mean((Sx - Semp)[iu] ** 2))),
                             abs_r_L1_Rbir=abs(np.corrcoef(lat['L1'], T['R_bir'])[0, 1]),
                             abs_r_L2_Rg=abs(np.corrcoef(lat['L2'], T['Rg'])[0, 1]),
                             **{l: v for l, v in zip(lab, r.x)}))
        ms = pd.DataFrame(rows).sort_values('nll')
        ms.to_csv(os.path.join(OUT, f'{tag}_multistart.csv'), index=False)
        good = ms[ms.valid & (ms.nll < ms[ms.valid].nll.min() + 0.01)]
        lines.append(f'multistart maximum likelihood ({args.nstart} starts): {int(ms.valid.sum())} valid solutions; '
                     f'{len(good)} within 0.01 of the best NLL, their |r(L1,Rb)| range '
                     f'{good.abs_r_L1_Rbir.min():.2f}–{good.abs_r_L1_Rbir.max():.2f}, '
                     f'|r(L2,Rg)| range {good.abs_r_L2_Rg.min():.2f}–{good.abs_r_L2_Rg.max():.2f}')

        # 3. profile along L1->age
        k = lab.index('L1->X_age')
        prow = []
        for g in np.round(np.linspace(-0.95, 0.95, 39), 3):
            best = None
            for s in range(8):
                r0 = np.random.default_rng(2000 + s)
                th0 = np.tanh(r0.normal(0, 0.5, len(edges)))
                f = lambda t_: nll(np.insert(t_, k, g), n, edges, order, xi, Semp)
                r = minimize(f, np.delete(th0, k), method='L-BFGS-B',
                             bounds=[(-0.999, 0.999)] * (len(edges) - 1), options={'maxiter': 5000})
                if best is None or r.fun < best.fun:
                    best = r
            th = np.insert(best.x, k, g)
            F, S, om = implied(th, n, edges, order)
            lat = latents(X, names, idx, F)
            prow.append(dict(L1_age=g, valid=bool(om.min() > -1e-6), nll=nll(th, n, edges, order, xi, Semp, pen=0),
                             abs_r_L1_Rbir=abs(np.corrcoef(lat['L1'], T['R_bir'])[0, 1]),
                             abs_r_L2_Rg=abs(np.corrcoef(lat['L2'], T['Rg'])[0, 1]),
                             L1_MgH=F[idx['L1'], idx['X_MgH']]))
        pf = pd.DataFrame(prow)
        pf.to_csv(os.path.join(OUT, f'{tag}_profile.csv'), index=False)
        v = pf[pf.valid]
        flat = v[v.nll < v.nll.min() + 0.01]
        lines.append(f'profile (L1->age fixed, rest refitted): points within 0.01 of the best NLL have L1->age from {flat.L1_age.min():.2f} to '
                     f'{flat.L1_age.max():.2f}; their |r(L1,Rb)| from {flat.abs_r_L1_Rbir.min():.2f} to {flat.abs_r_L1_Rbir.max():.2f}')

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
        mn = 0.5 * np.linalg.slogdet(Semp)[1] + 3
        axs[0].plot(v.L1_age, v.nll, 'o-', color=BLUE, ms=4, lw=1.5)
        axs[0].axhline(mn, color=MUTED, ls='--', lw=1)
        axs[0].set_ylim(mn - 0.1, mn + 0.8)
        axs[0].set_xlabel('coefficient L1 -> age (held fixed)')
        axs[0].set_ylabel('best NLL (other 14 refit)')
        axs[0].set_title('fit quality (dashed: theoretical minimum)', fontsize=10, color=INK)
        axs[1].plot(v.L1_age, v.abs_r_L1_Rbir, 'o-', color=BLUE, ms=4, lw=1.5, label='|r(L1, true Rbirth)|')
        axs[1].plot(v.L1_age, v.abs_r_L2_Rg, 's-', color=INK, ms=3, lw=1, label='|r(L2, true Rg)|')
        axs[1].set_ylim(0, 1)
        axs[1].set_xlabel('coefficient L1 -> age (held fixed)')
        axs[1].set_title('meaning of the latents', fontsize=10, color=INK)
        axs[1].legend(frameon=False, fontsize=8, loc='lower center')
        for ax in axs:
            for sp in ['top', 'right']:
                ax.spines[sp].set_visible(False)
            ax.grid(color='#e5e7eb', lw=0.6)
            ax.set_axisbelow(True)
        fig.suptitle(f'{tag}: author\'s causal graph frozen', fontsize=11, color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f'{tag}_profile.png'), dpi=130)
        plt.close(fig)

        txt = '\n'.join(lines)
        print(txt + '\n')
        with open(os.path.join(OUT, f'{tag}_summary.txt'), 'w') as f:
            f.write(txt + '\n')


if __name__ == '__main__':
    main()
