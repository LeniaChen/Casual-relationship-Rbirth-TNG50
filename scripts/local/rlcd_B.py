# -*- coding: utf-8 -*-
"""
rlcd_B.py — Approach B: on each TNG50 galaxy, run RLCD structure learning from scratch with many random seeds and see which causal graphs the data produce
(older copy; the version used for the results is PlanB/output_data/rlcd_B.py)

Run with the RLCD environment (needs causal-learn etc.):
    cd <repo root>
    $PY_RLCD scripts/local/rlcd_B.py --cond main        --version nonoise
    $PY_RLCD scripts/local/rlcd_B.py --cond R_0_30_Z_30 --version nonoise
    $PY_RLCD scripts/local/rlcd_B.py --cond R_0_30_Z_3  --version nonoise
    options: --seeds 20 (PYTHONHASHSEED = 0..19), --ids 538905 (selected galaxies), --alpha 1e-5 1e-15 1e-5 1e-15

Star selection (--cond):
  main         7 < R < 10 kpc, |z| < 0.3 kpc  (data data/rlcd_input/<version>/)       → output results/TNG/planB/
  R_0_30_Z_30  0 < R < 30 kpc, |z| < 30 kpc   (star_condition/R_0_30_Z_30/)          → output star_condition/R_0_30_Z_30/planB/
  R_0_30_Z_3   0 < R < 30 kpc, |z| < 3 kpc    (star_condition/R_0_30_Z_3/)           → output star_condition/R_0_30_Z_3/planB/

Procedure:
  1. RLCD settings identical to the authors' (defaults of rlcd/run_rlcd.py: stage-1 GES, sparsity 2; alpha defaults to the values in the authors'
     result file name, 1e-5, 1e-15, 1e-5, 1e-15). Each seed runs in its own Python process with PYTHONHASHSEED set (the hash seed can only be set at start-up);
     the learned graph (directed edges, undirected edges, number of latents) is recorded.
  2. Why the covariance matrix suffices: RLCD uses the data in only two places — the rank tests (canonical correlation analysis, which uses only the covariance matrix and N)
     and the BIC score of the stage-1 GES (np.cov and N). So for the large data sets (0.37–0.72 million stars) this script replaces the raw data by a Gaussian data set
     with exactly the same sample covariance and the same N; the result is identical to using the raw data (--selftest checks this on 538905 with 7<R<10).
     The main-selection data are small and are used directly.
  3. Each graph is compared with the authors' graph (latent labels aligned automatically) and checked against the prior of paper Sec. 3.2 (age has no observed parent).
Outputs (<output dir>/):
  structures_<version>.csv   one row per galaxy and seed: number of latents, directed edges, undirected edges, same as the authors' graph?, prior passed?, structure id
  structures_<version>.txt   the structures found for each galaxy and their frequencies
"""
import os
import sys
import json
import argparse
import itertools
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
REPO = ROOT / 'Causal-Galactic-Archaeology'
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
OBS = set(COLS)


def cond_paths(cond):
    if cond == 'main':
        return ROOT / 'data' / 'rlcd_input', ROOT / 'results' / 'TNG' / 'planB'
    base = ROOT / 'star_condition' / cond
    return base, base / 'planB'


def load_data(cond, version, gid):
    """Return (6-column DataFrame, whether it is a covariance-matched surrogate)"""
    src, _ = cond_paths(cond)
    if cond == 'main':
        return pd.read_csv(src / version / f'{gid}_e.csv')[COLS], False
    m = json.load(open(src / 'output_data' / f'moments_{version}.json'))[gid]
    return surrogate(np.array(m['cov'])[:6, :6], np.array(m['mean'])[:6], m['N']), True


def surrogate(cov, mean, N, seed=0):
    """Build an N-row Gaussian data set whose sample mean and sample covariance (denominator N-1) equal the given ones exactly"""
    Z = np.random.default_rng(seed).standard_normal((N, len(mean)))
    Z -= Z.mean(0)
    Z = Z @ np.linalg.inv(np.linalg.cholesky(np.cov(Z.T, ddof=1))).T      # sample covariance = identity
    X = Z @ np.linalg.cholesky(cov).T + mean
    return pd.DataFrame(X, columns=COLS)


def run_rlcd(df, alpha):
    """Run RLCD in this process (same settings as run_rlcd.py); return (directed edges, undirected edges) with the X_ prefix removed"""
    os.chdir(REPO / 'rlcd')
    sys.path.insert(0, str(REPO / 'rlcd'))
    import logging
    from DGM.DataModel import DataModel
    from utils.Chi2RankTest import Chi2RankTest
    from StructureLearning.RLCD.RLCD_alg import RLCD
    logging.getLogger().setLevel(logging.WARNING)
    df = df.rename(columns={k: 'X_' + k for k in df.columns})
    dgm = DataModel(df, df)
    df_x, _ = dgm.generate_data(N=-1, normalized=True)
    params = {'ranktest_method': Chi2RankTest(df_x.to_numpy(), 1), 'citest_method': None, 'stage1_method': 'ges',
              'alpha_dict': dict(enumerate(alpha)), 'stage1_ges_sparsity': 2, 'stage1_partition_thres': 3}
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):          # RLCD prints a lot of intermediate output
        graph, _, _, _ = RLCD(1, dgm.xvars, df_x, params)
    s = lambda v: v.replace('X_', '')
    return sorted((s(a), s(b)) for a, b in graph.dirEdges), sorted(tuple(sorted(map(s, e))) for e in graph.undirEdges)


def author_edges():
    import re
    txt = open(AUTHOR_DOT).read()
    return {(a.replace('X_', ''), b.replace('X_', '')) for a, b in re.findall(r'^\s*(\w+)\s*->\s*(\w+)', txt, re.M)}


def canonical(edges):
    """Latent labels are arbitrary: use the lexicographically smallest edge list over all relabelings as the structure's fingerprint"""
    lat = sorted({v for e in edges for v in e if v not in OBS})
    best = None
    for perm in itertools.permutations(range(1, len(lat) + 1)):
        mp = {L: f'L{p}' for L, p in zip(lat, perm)}
        key = tuple(sorted((mp.get(a, a), mp.get(b, b)) for a, b in edges))
        if best is None or key < best:
            best = key
    return best


def worker(args):
    df, _ = load_data(args.cond, args.version, args.ids[0])
    d, u = run_rlcd(df, args.alpha)
    print('RESULT ' + json.dumps({'directed': d, 'undirected': u}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cond', default='main', choices=['main', 'R_0_30_Z_30', 'R_0_30_Z_3'])
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--seeds', type=int, default=20)
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--alpha', nargs=4, type=float, default=[1e-5, 1e-15, 1e-5, 1e-15])
    ap.add_argument('--worker', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.worker:
        return worker(args)
    if args.selftest:
        return selftest(args)

    _, out = cond_paths(args.cond)
    out.mkdir(parents=True, exist_ok=True)
    csv = out / f'structures_{args.version}.csv'
    done = pd.read_csv(csv) if csv.exists() else pd.DataFrame()
    author = canonical(author_edges())
    rows = done.to_dict('records')
    for gid in args.ids or IDS:
        for s in range(args.seeds):
            if len(done) and ((done.galaxy.astype(str) == gid) & (done.seed == s)).any():
                continue
            env = dict(os.environ, PYTHONHASHSEED=str(s))
            cmd = [sys.executable, str(HERE), '--worker', '--cond', args.cond, '--version', args.version,
                   '--ids', gid, '--alpha', *map(str, args.alpha)]
            p = subprocess.run(cmd, env=env, capture_output=True, text=True)
            line = [l for l in p.stdout.splitlines() if l.startswith('RESULT ')]
            if not line:
                print(p.stderr[-2000:])
                raise RuntimeError(f'{gid} seed {s} failed')
            res = json.loads(line[0][7:])
            d = [tuple(e) for e in res['directed']]
            key = canonical(d)
            lat = sorted({v for e in d for v in e if v not in OBS})
            rows.append(dict(galaxy=gid, seed=s, n_latent=len(lat), n_directed=len(d), n_undirected=len(res['undirected']),
                             same_as_author=(key == author and not res['undirected']),
                             prior_ok=not any(c == 'age' and p in OBS for p, c in d),
                             structure=json.dumps(key), undirected=json.dumps(res['undirected'])))
            print(f'{args.cond} {args.version} {gid} seed {s}: {len(lat)} latents, {len(d)} directed edges'
                  f'{", same as the authors graph" if rows[-1]["same_as_author"] else ""}', flush=True)
            pd.DataFrame(rows).to_csv(csv, index=False)
    summarize(out, args.version)


def summarize(out, version):
    df = pd.read_csv(out / f'structures_{version}.csv')
    codes = {k: i + 1 for i, k in enumerate(df.structure.value_counts().index)}   # structure id: by total count over all galaxies, most frequent first
    df['structure_id'] = df.structure.map(codes)
    df.to_csv(out / f'structures_{version}.csv', index=False)
    lines = [f'Approach B ({version}): {df.groupby("galaxy").size().max()} seeds per galaxy; structure ids sorted by total count over all galaxies',
             f'identical to the authors\' graph: {int(df.same_as_author.sum())} / {len(df)}', '']
    for gid, g in df.groupby('galaxy'):
        lines.append(f'== {gid}')
        for sid, gg in g.groupby('structure_id'):
            r = gg.iloc[0]
            lines.append(f'  structure {sid:>2}: {len(gg):>2} times  latents {r.n_latent}  directed {r.n_directed}  undirected {r.n_undirected}  '
                         f'prior {"pass" if r.prior_ok else "FAIL"}{"  <- authors graph" if r.same_as_author else ""}  seeds {sorted(gg.seed.tolist())}')
    lines += ['', 'edges of each structure (latents relabeled consistently):']
    for k, i in sorted(codes.items(), key=lambda x: x[1]):
        lines.append(f'  structure {i}: ' + ', '.join(f'{a}->{b}' for a, b in json.loads(k)))
    (out / f'structures_{version}.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:3]))


def selftest(args):
    """Check: RLCD learns exactly the same graph from the covariance-matched surrogate as from the raw data (538905 with 7<R<10, seed 1)"""
    gid = '538905'
    x = pd.read_csv(ROOT / 'data' / 'rlcd_input' / 'nonoise' / f'{gid}_e.csv')[COLS]
    sur = surrogate(np.cov(x.to_numpy().T, ddof=1), x.mean().to_numpy(), len(x))
    print('max covariance difference', np.abs(np.cov(sur.to_numpy().T) - np.cov(x.to_numpy().T)).max())
    a = run_rlcd(x, args.alpha)
    b = run_rlcd(sur, args.alpha)
    print('raw data:   ', a)
    print('surrogate:  ', b)
    print('identical' if a == b else 'DIFFERENT!')


if __name__ == '__main__':
    main()
