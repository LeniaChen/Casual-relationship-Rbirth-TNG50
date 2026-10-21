# -*- coding: utf-8 -*-
"""
rlcd_A.py — Approach A: freeze the authors' causal graph, re-estimate the edge coefficients on each galaxy, compute L1, L2 for every star and compare with the truth

Must be run with the RLCD environment (needs torch, pydot, pandas, matplotlib), from the repo root:
    cd <repo root>
    $PY_RLCD scripts/local/rlcd_A.py --version g279            # first validate the whole chain on the authors' data
    $PY_RLCD scripts/local/rlcd_A.py --version nonoise --ids 538905   # trial run on one galaxy
    $PY_RLCD scripts/local/rlcd_A.py --version nonoise          # 12 galaxies, no errors
    $PY_RLCD scripts/local/rlcd_A.py --version noise            # 12 galaxies, with errors

Procedure (as in Sec. 3.3 of the paper and the authors' notebook):
  1. read the 6 observables and standardize them (subtract mean, divide by std), as DataModel.generate_data(normalized=True) in the repo
  2. causal graph: the authors' final g2.79e12 graph (results/rlcd_g279/author_...dot), structure kept fixed
  3. parameter learning: call parameter_by_trek from ParameterLearning/linear_parameter_identification.py of the repo
     (as in its README: lr=0.005, epochs=600, lr_decay=plateau). It starts from random initial values and the result depends on them,
     so it is repeated N_REPEAT times (torch seeds 0..N-1, default 50) and one solution is picked:
       (a) first drop physically impossible solutions: an implied node noise variance < 0, or a non-positive-definite covariance.
           The repo code does not check this, and it computes the NLL with the ABSOLUTE value of the determinant, so impossible
           solutions can get a "very good" NLL (NLL = -7.03 occurred on g2.79e12, below the theoretical minimum -1.27, which is impossible)
       (b) among the remaining ones, pick by default the smallest RMSE between model and data correlation matrices (--select rmse);
           or the smallest NLL (--select nll).
       Reason: the coefficients used in the authors' notebook and paper are the repo's set with the SECOND smallest NLL (1.586), which has the smallest RMSE (0.117);
       the set with the smallest NLL (1.498) has RMSE 0.260, its L1 is almost disconnected from the chemistry and correlates with the birth radius at only 0.01.
       The NLL is dominated by the 0.997 correlation between [Fe/H] and [Mg/H] and is not a good selection criterion here.
     The correlations of L1, L2 with the truth for every repeat are stored in nll_repeats.csv for diagnostics only; they are never used for selection
  4. latent values for every star (paper Sec. 3.3, from the lowest latent upwards):
       every child c of a latent L satisfies  c = a_Lc * L + Σ(other parents p) a_pc * p + noise
       let r_c = c - Σ(other parents) a_pc * p  (what is left for L to explain after removing the other parents);
       the L minimizing Σ_c (r_c - a_Lc * L)^2 is   L = Σ_c a_Lc * r_c / Σ_c a_Lc^2
     This is identical to the per-star Nelder-Mead minimization of chisq in the authors' notebook (a quadratic has a closed-form minimum)
  5. linear conversion to physical units: f(L1) = α1*L1 + β1 fitted to the true Rbirth, f(L2) = α2*L2 + β2 fitted to the true Rg (paper pp. 5-6);
     reports the Pearson correlation r, the residual standard deviation and the median absolute deviation

Outputs: results/rlcd_A/<version>/<ID>/
        nll_repeats.csv      NLL and all edge coefficients of every repeat
        best_params.dot/.png the selected coefficients drawn on the causal graph
        latents.csv          per star: iord (if available), L1, L2, f(L1), f(L2), R_bir, Rg
        fit.txt              conversion formulas and metrics
      results/rlcd_A/<version>/summary.csv, summary.txt
      results/rlcd_A/<version>/L1_vs_Rbir.png, L2_vs_Rg.png  (one panel per galaxy)
"""
import os
import io
import sys
import time
import argparse
import contextlib
import numpy as np
import pandas as pd

REPO = 'Causal-Galactic-Archaeology'
sys.path.insert(0, os.path.join(REPO, 'rlcd', 'ParameterLearning'))
sys.path.insert(0, os.path.join(REPO, 'rlcd'))
import torch                                              # noqa: E402
import pydot                                              # noqa: E402
from linear_parameter_identification import (             # noqa: E402
    parameter_by_trek, get_F_bin_and_varnames_from_dotgraph, add_edge_coefficients_to_graph)

AUTHOR_DOT = 'results/rlcd_g279/author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def load(version, gid):
    """Return (observables DataFrame, truth DataFrame)"""
    if version == 'g279':
        e = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_e.csv'))
        t = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_multi_ele_action_Jin_SN.csv'))
    else:
        e = pd.read_csv(f'data/rlcd_input/{version}/{gid}_e.csv')
        t = pd.read_csv(f'data/rlcd_input/{version}/{gid}_truth.csv')
    assert len(e) == len(t)
    ok = np.isfinite(e[COLS].to_numpy()).all(axis=1)
    return e.loc[ok, COLS].reset_index(drop=True), t.loc[ok].reset_index(drop=True)


def estimate_latents(Xn, F, names):
    """Solve the latents of every star layer by layer as in paper Sec. 3.3 (closed-form least squares)"""
    val = {f'X_{c}': Xn[c].to_numpy() for c in COLS}
    lat = [v for v in names if not v.startswith('X_')]
    idx = {v: i for i, v in enumerate(names)}
    done = {}
    while len(done) < len(lat):
        progress = False
        for L in lat:
            if L in done:
                continue
            children = [c for c in names if F[idx[L], idx[c]] != 0]
            if any(c in lat and c not in done for c in children):
                continue                       # solve the lower latents first
            num = 0.0
            den = 0.0
            for c in children:
                xc = val[c] if c in val else done[c]
                r = xc.copy()
                for p in names:
                    if p != L and F[idx[p], idx[c]] != 0:
                        xp = val[p] if p in val else done.get(p)
                        if xp is None:
                            raise RuntimeError(f"parent {p} of {c} not solved yet; cannot solve {L}")
                        r = r - F[idx[p], idx[c]] * xp
                a = F[idx[L], idx[c]]
                num = num + a * r
                den = den + a * a
            done[L] = num / den
            progress = True
        if not progress:
            raise RuntimeError('cycle among the latents; cannot solve layer by layer')
    return done


def implied_cov(F, names):
    """Model-implied covariance from the edge coefficients (all variables standardized to unit variance; equivalent to
    the trek rule of parameter_by_trek), plus the noise variance of each node ω_k = 1 - Var(contribution of its parents).
    ω_k < 0 means the coefficients are physically impossible (the parents "explain" more than the node's total variance)."""
    n = len(names)
    order, rem = [], list(range(n))
    while rem:
        for k in rem:
            if all(F[p, k] == 0 or p in order for p in range(n)):
                order.append(k)
                rem.remove(k)
                break
        else:
            raise RuntimeError('the causal graph has a cycle')
    S = np.zeros((n, n))
    omega = np.zeros(n)
    for pos, k in enumerate(order):
        par = [p for p in order[:pos] if F[p, k] != 0]
        for j in order[:pos]:
            S[j, k] = S[k, j] = sum(F[p, k] * S[j, p] for p in par)
        omega[k] = 1 - sum(F[p, k] * F[q, k] * S[p, q] for p in par for q in par)
        S[k, k] = 1.0
    return S, omega


def evaluate_solution(F, names, S_emp):
    """Quality of a coefficient set: physically valid or not, the correct NLL, RMSE against the data correlation matrix"""
    S, omega = implied_cov(F, names)
    xi = [i for i, v in enumerate(names) if v.startswith('X_')]
    Sx = S[np.ix_(xi, xi)]
    sign, logdet = np.linalg.slogdet(Sx)
    valid = bool(sign > 0 and omega.min() > -1e-6)
    nll = 0.5 * logdet + 0.5 * np.trace(np.linalg.solve(Sx, S_emp)) if sign > 0 else np.inf
    iu = np.triu_indices(len(xi), 1)
    rmse = float(np.sqrt(np.mean((Sx - S_emp)[iu] ** 2)))
    return dict(valid=valid, nll_true=float(nll), rmse=rmse, min_noise_var=float(omega.min()))


def fit_linear(L, truth):
    ok = np.isfinite(L) & np.isfinite(truth)
    a, b = np.polyfit(L[ok], truth[ok], 1)
    f = a * L + b
    res = truth[ok] - f[ok]
    r = np.corrcoef(f[ok], truth[ok])[0, 1]
    return a, b, f, r, float(np.std(res)), float(np.median(np.abs(res)))


def run_one(version, gid, dot_path, n_repeat, epochs, lr, lr_decay, outdir, select='rmse'):
    t0 = time.time()
    X, T = load(version, gid)
    Xn = (X - X.mean()) / X.std()
    df = Xn.rename(columns={c: f'X_{c}' for c in COLS})

    graph = pydot.graph_from_dot_file(dot_path)[0]
    F_bin, names = get_F_bin_and_varnames_from_dotgraph(graph)
    df = df[[v for v in names if v.startswith('X_')]]

    S_emp = df.cov().to_numpy()
    rows, sols = [], []
    for seed in range(n_repeat):
        torch.manual_seed(seed)
        np.random.seed(seed)
        with contextlib.redirect_stdout(io.StringIO()):
            F, nll = parameter_by_trek(df, F_bin, names, lr=lr, epochs=epochs, lr_decay=lr_decay)
        ev = evaluate_solution(F, names, S_emp)
        # diagnostics only (NOT used for selection!): correlation of L1, L2 from these coefficients with the truth
        lat_i = estimate_latents(Xn, F, names)
        diag_r1 = abs(np.corrcoef(lat_i['L1'], T['R_bir'])[0, 1])
        diag_r2 = abs(np.corrcoef(lat_i['L2'], T['Rg'])[0, 1])
        row = {'seed': seed, 'nll_code': float(nll), **ev,
               'diag_abs_r_L1_Rbir': diag_r1, 'diag_abs_r_L2_Rg': diag_r2}
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if F_bin[i, j]:
                    row[f'{a}->{b}'] = F[i, j]
        rows.append(row)
        sols.append((F, row))
        print(f"    [{gid}] repeat {seed+1:>2}/{n_repeat}: {'valid' if ev['valid'] else 'invalid'}  "
              f"NLL {ev['nll_true']:8.4f}  RMSE {ev['rmse']:.4f}   (diagnostic |r(L1,Rb)| {diag_r1:.2f}, |r(L2,Rg)| {diag_r2:.2f})",
              flush=True)

    os.makedirs(outdir, exist_ok=True)
    tab = pd.DataFrame(rows)
    valid = [x for x in sols if x[1]['valid']]
    if not valid:
        raise RuntimeError(f"{gid}: none of the {n_repeat} parameter-learning runs is physically valid")
    key = 'rmse' if select == 'rmse' else 'nll_true'
    F, brow = min(valid, key=lambda x: x[1][key])
    nll, seed = brow['nll_true'], brow['seed']
    tab['selected'] = tab['seed'] == seed
    tab.sort_values(['valid', key], ascending=[False, True]).to_csv(os.path.join(outdir, 'nll_repeats.csv'), index=False)
    n_valid = len(valid)
    top = sorted(valid, key=lambda x: x[1][key])[:5]
    top_txt = ', '.join(f"{x[1]['diag_abs_r_L1_Rbir']:.2f}" for x in top)
    add_edge_coefficients_to_graph(graph, F_bin, F, names)
    graph.write_raw(os.path.join(outdir, 'best_params.dot'))
    try:
        graph.write_png(os.path.join(outdir, 'best_params.png'))
    except Exception as err:
        print(f"    (png not drawn, graphviz may be missing: {err})")

    lat = estimate_latents(Xn, F, names)
    L1, L2 = lat['L1'], lat['L2']
    a1, b1, f1, r1, s1, m1 = fit_linear(L1, T['R_bir'].to_numpy())
    a2, b2, f2, r2, s2, m2 = fit_linear(L2, T['Rg'].to_numpy())

    out = pd.DataFrame({'L1': L1, 'L2': L2, 'f_L1': f1, 'f_L2': f2,
                        'R_bir': T['R_bir'], 'Rg': T['Rg'], 'age': X['age']})
    if 'iord' in T:
        out.insert(0, 'iord', T['iord'])
    out.to_csv(os.path.join(outdir, 'latents.csv'), index=False)

    txt = (f"galaxy {gid} ({version}), {len(X)} stars; {n_valid} of {n_repeat} parameter-learning runs physically valid; "
           f"selected by {key}: seed {seed}, NLL {nll:.4f}, RMSE {brow['rmse']:.4f}\n"
           f"  Rbirth = {a1:+.3f} * L1 {b1:+.3f}   Pearson r = {r1:.3f}   residual std {s1:.3f} kpc   median abs dev {m1:.3f} kpc\n"
           f"  Rg     = {a2:+.3f} * L2 {b2:+.3f}   Pearson r = {r2:.3f}   residual std {s2:.3f} kpc   median abs dev {m2:.3f} kpc\n"
           f"  (stability diagnostic: top-5 valid solutions by {key}, |r(L1,Rb)| = {top_txt})\n"
           f"  time {time.time()-t0:.0f} s\n")
    with open(os.path.join(outdir, 'fit.txt'), 'w') as f:
        f.write(txt)
    print(txt, flush=True)
    return dict(ID=gid, n=len(X), n_valid=n_valid, best_nll=nll, best_rmse=brow['rmse'], best_seed=seed,
                a1=a1, b1=b1, r_L1_Rbir=r1, std_Rbir=s1, mad_Rbir=m1,
                a2=a2, b2=b2, r_L2_Rg=r2, std_Rg=s2, mad_Rg=m2), out


def plots(results, outroot, version):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(results)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    for key, truth, lab, rkey in [('f_L1', 'R_bir', 'Rbirth', 'r_L1_Rbir'), ('f_L2', 'Rg', 'Rg', 'r_L2_Rg')]:
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.2 * nrow), squeeze=False)
        axs = axs.ravel()
        for ax, (summ, out) in zip(axs, results):
            lim = (0, 16)
            ax.hexbin(out[key], out[truth], gridsize=45, extent=lim + lim, bins='log', cmap='Blues', mincnt=1)
            ax.plot(lim, lim, color='#1f2937', lw=0.8)
            ax.set_title(f"{summ['ID']}  r = {summ[rkey]:.2f}", fontsize=9)
            ax.set_xlabel(f'f({key[2:]}) [kpc]', fontsize=8)
            ax.set_ylabel(f'true {lab} [kpc]', fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in axs[n:]:
            ax.axis('off')
        fig.suptitle(f'{version}: f({key[2:]}) vs true {lab}', fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(outroot, f'{key[2:]}_vs_{truth}.png'), dpi=120)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['g279', 'nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', help='only these galaxies (default: all 12)')
    ap.add_argument('--n-repeat', type=int, default=50, help='number of parameter-learning repeats')
    ap.add_argument('--select', default='rmse', choices=['rmse', 'nll'],
                    help='selection among physically valid solutions: rmse = smallest difference between model and data correlation matrices (default); nll = smallest negative log-likelihood')
    ap.add_argument('--epochs', type=int, default=600)
    ap.add_argument('--lr', type=float, default=0.005)
    ap.add_argument('--lr-decay', default='plateau')
    ap.add_argument('--dot', default=AUTHOR_DOT, help='frozen causal graph (default: the authors\' g2.79e12 graph)')
    args = ap.parse_args()

    if not os.path.isdir(REPO):
        sys.exit('run from the repo root')
    ids = ['g279'] if args.version == 'g279' else (args.ids or IDS)
    outroot = os.path.join('results', 'rlcd_A', args.version)
    os.makedirs(outroot, exist_ok=True)
    print(f"Approach A: frozen causal graph {args.dot}\nversion {args.version}, galaxies {ids}, {args.n_repeat} parameter-learning repeats per galaxy\n")

    results = []
    for gid in ids:
        summ, out = run_one(args.version, gid, args.dot, args.n_repeat, args.epochs, args.lr,
                            args.lr_decay, os.path.join(outroot, gid), select=args.select)
        results.append((summ, out))

    summ = pd.DataFrame([s for s, _ in results])
    tag = '' if (args.version == 'g279' or not args.ids) else '_' + '_'.join(ids)
    summ.to_csv(os.path.join(outroot, f'summary{tag}.csv'), index=False)
    with open(os.path.join(outroot, f'summary{tag}.txt'), 'w') as f:
        f.write(summ.to_string(index=False) + '\n')
    print(summ[['ID', 'n', 'n_valid', 'best_nll', 'best_rmse', 'r_L1_Rbir', 'std_Rbir', 'r_L2_Rg', 'std_Rg']].to_string(index=False))
    if args.version == 'g279':
        print("\npaper (g2.79e12): Rb = -1.97 L1 + 9.31, r(L1, Rb) = 0.64;  Rg = -1.29 L2 + 9.52, r(L2, Rg) = 0.98")
        print("the authors actually used the coefficient set with the second smallest NLL (repo file nll1.586), which has the smallest RMSE (0.117)")
        print("note: the sign of a latent is arbitrary, so the slope may have the opposite sign to the paper; this does not matter")
    plots(results, outroot, args.version)
    print(f"\noutputs in {outroot}/")


if __name__ == '__main__':
    main()
