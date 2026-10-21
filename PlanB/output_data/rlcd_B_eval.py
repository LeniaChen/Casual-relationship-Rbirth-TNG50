# -*- coding: utf-8 -*-
"""
rlcd_B_eval.py — Approach B, step 2: for the graph(s) learned on each galaxy, fit the edge coefficients and relate every latent to the true birth and guiding radii

Run rlcd_B.py in the same folder first to get PlanB/output_data/structures_<version>.csv, then (numpy, scipy, pandas; the analysis environment is enough):
    cd <repo root>
    $PY PlanB/output_data/rlcd_B_eval.py --version nonoise
    $PY PlanB/output_data/rlcd_B_eval.py --version noise
    (--min-frac 0.2: evaluate only structures with frequency >= 20% for that galaxy, and at least the most common one; --nstart 50: number of random starts)

For each (galaxy, structure):
  1. Jacobian rank: number of edges vs number of independent observed correlations (15); can the coefficients be determined by the data (as in rlcd_identifiability.py);
  2. maximum likelihood: nstart random starts, noise variances >= 0, strict convergence (as in rlcd_identifiability_ci.py);
     all solutions statistically indistinguishable from the best (2·N·ΔNLL < 3.84) are kept;
  3. latent values: the same hierarchical least squares as the authors' notebook (from the lowest latent upwards); each latent is a linear combination of the 6 standardized inputs,
     so its correlation with the truth follows directly from the covariance matrix (see star_condition/R_0_30_Z_30/output_data/identifiability_moments.py);
  4. for every latent, the range of |r(L, true birth radius)| and |r(L, true guiding radius)| among the indistinguishable solutions;
     baselines: |r| of age, [Fe/H], J_phi alone with the birth radius, and the linear-regression ceiling of the 6 inputs.
Outputs (PlanB/output_data/; other --cond values in subfolders):
  eval_<version>.csv   one row per galaxy, structure and latent (|r| ranges, Jacobian rank, baselines)
  eval_<version>.txt   same, as text
  eval_<version>_best.json  maximum-likelihood solution of every galaxy and structure: edge coefficients and the weights w of every latent on the 6 standardized inputs
                    (latent L = Σ_i w_i·z_i, z_i = input i minus its mean divided by its std); the plotting scripts use it to compute the latents of every star
"""
import sys
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
sys.path.insert(0, str(ROOT / 'scripts' / 'local'))      # rlcd_identifiability.py (implied, nll)
sys.path.insert(0, str(HERE.parent))                     # rlcd_B.py of the same folder (first on the path, so the older copy in scripts/local is not used)
from rlcd_identifiability import implied, nll  # noqa: E402
from rlcd_B import cond_paths, COLS, OBS, IDS  # noqa: E402

LR95 = 3.84
PEN = 1e7
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def moments(cond, version, gid):
    """8×8 covariance (6 inputs + R_bir + Rg) and number of stars"""
    src, _ = cond_paths(cond)
    if cond == 'main':
        e = pd.read_csv(PLANB / 'input_data' / version / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(PLANB / 'true_value' / version / f'{gid}_truth.csv', usecols=['R_bir', 'Rg'])
        X = np.column_stack([e.to_numpy(float), t[['R_bir', 'Rg']].to_numpy(float)])
        return np.cov(X.T, ddof=1), len(X)
    m = json.load(open(src / 'output_data' / f'moments_{version}.json'))[gid]
    return np.array(m['cov']), m['N']


def graph_setup(edges):
    names = COLS + sorted({v for e in edges for v in e if v not in OBS})
    idx = {v: i for i, v in enumerate(names)}
    E = [(idx[a], idx[b]) for a, b in edges]
    order, rem = [], list(range(len(names)))
    while rem:
        for k in rem:
            if all((p, k) not in E or p in order for p in range(len(names))):
                order.append(k)
                rem.remove(k)
                break
        else:
            raise ValueError('the graph has a cycle')
    return names, idx, E, order


def latent_weights(names, edges_named, coef):
    """Weights of all latents on the 6 standardized inputs.
    Same hierarchical least squares as the authors' notebook: starting from the bottom, take the group S of latents whose latent children are all solved,
    subtract from every child c of S the contribution of its known parents outside S, r_c = c − Σ a_qc·q, then find the S minimizing Σ_c (r_c − Σ_{L∈S} a_Lc·L)².
    With a single latent in S this is the authors' formula L = Σ a_Lc·r_c / Σ a_Lc²; two latents sharing children in the same layer (not present in the authors' graph)
    are solved jointly. Contributions of higher (not yet solved) latents to the children are set to 0, as in the authors' layer-by-layer procedure."""
    w = {c: np.eye(6)[i] for i, c in enumerate(COLS)}
    lat = [v for v in names if v not in OBS]
    todo = set(lat)
    while todo:
        S = [L for L in sorted(todo) if not any(c in todo for (p, c) in edges_named if p == L)]
        if not S:
            raise ValueError('cycle among the latents')
        ch = sorted({c for (p, c) in edges_named if p in S})
        A = np.array([[coef.get((L, c), 0.0) for L in S] for c in ch])
        R = []
        for c in ch:
            r = w[c].copy()
            for (q, c2) in edges_named:
                if c2 == c and q not in S and q in w:
                    r = r - coef[(q, c2)] * w[q]
            R.append(r)
        Wm = np.linalg.pinv(A) @ np.array(R)            # |S| × 6
        for L, row in zip(S, Wm):
            w[L] = row
        todo -= set(S)
    return {L: w[L] for L in lat}


def evaluate(edges_named, cov8, N, nstart):
    names, idx, E, order = graph_setup(edges_named)
    n, ne = len(names), len(E)
    xi = [idx[c] for c in COLS]
    iu = np.triu_indices(6, 1)
    sd = np.sqrt(np.diag(cov8))
    C = cov8[:6, :6] / np.outer(sd[:6], sd[:6])
    cR, cG = cov8[:6, 6] / sd[:6], cov8[:6, 7] / sd[:6]
    vR, vG = cov8[6, 6], cov8[7, 7]
    cv = lambda th: implied(th, n, E, order)[1][np.ix_(xi, xi)][iu]
    th = np.tanh(np.random.default_rng(0).normal(0, 0.5, ne))
    J = np.array([(cv(th + 1e-6 * np.eye(ne)[k]) - cv(th - 1e-6 * np.eye(ne)[k])) / 2e-6 for k in range(ne)]).T
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(sv > 1e-8 * sv[0]))

    def f(t):
        _, _, om = implied(t, n, E, order)
        return nll(t, n, E, order, xi, C, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

    sols = []
    for s in range(nstart):
        t0 = np.tanh(np.random.default_rng(7000 + s).normal(0, 0.5, ne))
        r = minimize(f, t0, method='L-BFGS-B', bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
        F, S, om = implied(r.x, n, E, order)
        if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
            continue
        sols.append((nll(r.x, n, E, order, xi, C, pen=0), r.x, om.min(), names[int(np.argmin(om))]))
    if not sols:
        return rank, ne, 0, np.nan, pd.DataFrame(columns=['latent', 'nll', 'omega_min', 'omega_min_node', 'r_Rbir', 'r_Rg']), None
    sols.sort(key=lambda z: z[0])
    best = sols[0][0]
    out, best_info = [], None
    for v, t, ommin, omnode in sols:
        if 2 * N * (v - best) >= LR95:
            continue
        coef = {(names[a], names[b]): x for (a, b), x in zip(E, t)}
        W = latent_weights(names, edges_named, coef)
        if best_info is None:                       # sols are sorted by NLL, the first one is the maximum-likelihood solution
            _, _, om = implied(t, n, E, order)
            best_info = dict(nll=v, edges=[list(e) for e in edges_named], coef=[coef[e] for e in edges_named],
                             omega={names[k]: float(om[k]) for k in range(n)},
                             weights={L: [float(x) for x in wv] for L, wv in W.items()})
        for L, wv in W.items():
            q = wv @ C @ wv
            out.append(dict(latent=L, nll=v, omega_min=ommin, omega_min_node=omnode,
                            r_Rbir=abs(wv @ cR) / np.sqrt(q * vR) if q > 0 else np.nan,
                            r_Rg=abs(wv @ cG) / np.sqrt(q * vG) if q > 0 else np.nan))
    cols = ['latent', 'nll', 'omega_min', 'omega_min_node', 'r_Rbir', 'r_Rg']
    return rank, ne, len(sols), best, pd.DataFrame(out, columns=cols), best_info


def baselines(cov8):
    sd = np.sqrt(np.diag(cov8))
    R = cov8 / np.outer(sd, sd)
    C, cR = R[:6, :6], R[:6, 6]
    return dict(age=abs(R[3, 6]), FeH=abs(R[4, 6]), JPhi=abs(R[1, 6]), ceiling=np.sqrt(cR @ np.linalg.solve(C, cR)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cond', default='main', choices=['main', 'R_0_30_Z_30', 'R_0_30_Z_3', 'R_0_30_Z_0.3'])
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--min-frac', type=float, default=0.2)
    ap.add_argument('--nstart', type=int, default=50)
    ap.add_argument('--ids', nargs='*', default=None)
    args = ap.parse_args()
    _, out = cond_paths(args.cond)
    st = pd.read_csv(out / f'structures_{args.version}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    rows, lines = [], [f'Approach B evaluation ({args.cond}, {args.version}): per galaxy, structures with frequency >= {args.min_frac:.0%} (at least the most common one)',
                        'for each latent: |r(L, birth radius)| and |r(L, guiding radius)| at the maximum-likelihood solution; in brackets the range among solutions statistically indistinguishable from the best (2·N·ΔNLL<3.84)', '']
    bestjson = {}
    for gid in args.ids or IDS:
        g = st[st.galaxy == gid]
        if not len(g):
            continue
        cov8, N = moments(args.cond, args.version, gid)
        b = baselines(cov8)
        freq = g.structure_id.value_counts()
        pick = [sid for sid, c in freq.items() if c / len(g) >= args.min_frac] or [freq.index[0]]
        lines.append(f'== {gid} (N={N:,}); baselines: age alone {b["age"]:.2f}, [Fe/H] {b["FeH"]:.2f}, J_phi {b["JPhi"]:.2f}, linear ceiling {b["ceiling"]:.2f}')
        for sid in pick:
            gg = g[g.structure_id == sid]
            if gg.n_undirected.iloc[0] > 0:
                lines.append(f'  structure {sid} ({len(gg)}/{len(g)} times): has undirected edges, skipped')
                continue
            edges = [tuple(e) for e in json.loads(gg.structure.iloc[0])]
            rank, ne, nvalid, best, sol, info = evaluate(edges, cov8, N, args.nstart)
            lines.append(f'  structure {sid} ({len(gg)}/{len(g)} times, {gg.n_latent.iloc[0]} latents, {ne} edges): Jacobian rank {rank}/{ne}'
                         f'{" (identifiable)" if rank == ne else " (NOT identifiable)"}; valid solutions {nvalid}, best NLL {best:.4f}')
            if info is not None:
                bestjson.setdefault(gid, {})[str(sid)] = dict(freq=len(gg) / len(g), n_seeds=len(g), jacobian_rank=rank,
                                                              n_edges=ne, N=N, **info)
            for L, s in sol.groupby('latent', sort=True):
                lines.append(f'    {L}: maximum likelihood |r(L,birth radius)| {s.r_Rbir.iloc[0]:.2f} (indistinguishable {s.r_Rbir.min():.2f}–{s.r_Rbir.max():.2f})'
                             f'   |r(L,guiding radius)| {s.r_Rg.iloc[0]:.2f} ({s.r_Rg.min():.2f}–{s.r_Rg.max():.2f})')
                rows.append(dict(galaxy=gid, structure_id=sid, freq=len(gg) / len(g), n_latent=gg.n_latent.iloc[0], n_edges=ne,
                                 jacobian_rank=rank, latent=L, n_indistinguishable=len(s),
                                 r_Rbir_best=s.r_Rbir.iloc[0], r_Rbir_min=s.r_Rbir.min(), r_Rbir_max=s.r_Rbir.max(),
                                 r_Rg_best=s.r_Rg.iloc[0], r_Rg_min=s.r_Rg.min(), r_Rg_max=s.r_Rg.max(),
                                 omega_min_best=s.omega_min.iloc[0], omega_min_node=s.omega_min_node.iloc[0],
                                 base_age=b['age'], base_FeH=b['FeH'], base_JPhi=b['JPhi'], ceiling=b['ceiling']))
            print('\n'.join(lines[-(1 + sol.latent.nunique()):]), flush=True)
    pd.DataFrame(rows).to_csv(out / f'eval_{args.version}.csv', index=False)
    (out / f'eval_{args.version}.txt').write_text('\n'.join(lines) + '\n')
    json.dump(bestjson, open(out / f'eval_{args.version}_best.json', 'w'), indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
