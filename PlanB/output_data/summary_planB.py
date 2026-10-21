# -*- coding: utf-8 -*-
"""
summary_planB.py — Approach B, step 3: one-row-per-galaxy summary of the 12 galaxies, compared with Approach A (frozen authors' graph)

Run rlcd_B.py and rlcd_B_eval.py in the same folder first, then (numpy, pandas):
    cd <repo root>
    $PY PlanB/output_data/summary_planB.py --version nonoise
    $PY PlanB/output_data/summary_planB.py --version noise
Outputs: PlanB/output_data/summary_planB_<version>.csv and .txt

For each galaxy only the most frequent structure is used (main selection: 7<R<10 kpc, |z|<0.3 kpc, as in the paper and the main result of Approach A):
  n_structures   number of different graphs learned with the 30 seeds;  top_freq  how often the most common graph appears
  n_latent / n_edges / n_undirected   number of latents / directed edges / undirected edges (no coefficients can be fitted with undirected edges)
  same_as_author is it the authors' graph;  prior_ok does it satisfy the prior of paper Sec. 3.2 (age has no observed parent)
  jac_rank       Jacobian rank (= number of edges when the coefficients are determined by the data)
  B_best_L_Rb    the latent with the largest |r| with the true birth radius (chosen using the truth: best case) and its |r| at the maximum-likelihood solution
  B_best_L_Rg    same, for the guiding radius
  A_r1_range     range of |r(L1, birth radius)| in Approach A (authors' graph frozen; results/rlcd_A/identifiability/summary_ci_<version>.csv)
  NLL_B / NLL_A  maximum-likelihood NLL of the two graphs (mean negative log-likelihood per star; smaller = closer to the data; both use the same observed correlation matrix, so they are directly comparable)
  dBIC_per_star  (NLL_B − NLL_A) + (edges_B − 15)·ln N / (2N): the difference after the BIC penalty for extra edges (N = number of stars, 15 = edges of the authors' graph);
                 negative = learned graph better, positive = authors' graph better
  age_only / JPhi_only / ceiling   |r| of age alone / J_phi alone with the birth radius; linear-regression ceiling of the 6 inputs
"""
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    v = ap.parse_args().version
    st = pd.read_csv(HERE.with_name(f'structures_{v}.csv'))
    st['galaxy'] = st.galaxy.astype(str)
    ev = pd.read_csv(HERE.with_name(f'eval_{v}.csv'))
    ev['galaxy'] = ev.galaxy.astype(str)
    best = json.load(open(HERE.with_name(f'eval_{v}_best.json')))
    A = pd.read_csv(ROOT / 'results' / 'rlcd_A' / 'identifiability' / f'summary_ci_{v}.csv').set_index('tag')
    rows = []
    for gid in IDS:
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        sid = freq.index[0]
        r = g[g.structure_id == sid].iloc[0]
        e = ev[(ev.galaxy == gid) & (ev.structure_id == sid)]
        info = best.get(gid, {}).get(str(sid))
        a = A.loc[f'{v}_{gid}']
        row = dict(galaxy=gid, n_structures=len(freq), top_structure=sid, top_freq=f'{freq.iloc[0]}/{len(g)}',
                   n_latent=r.n_latent, n_edges=r.n_directed, n_undirected=r.n_undirected,
                   same_as_author=bool(r.same_as_author), prior_ok=bool(r.prior_ok),
                   jac_rank=info['jacobian_rank'] if info else None)
        if len(e):
            iR, iG = e.r_Rbir_best.idxmax(), e.r_Rg_best.idxmax()
            row.update(B_best_L_Rb=f'{e.latent[iR]} {e.r_Rbir_best[iR]:.2f}', B_r_Rb=e.r_Rbir_best[iR],
                       B_best_L_Rg=f'{e.latent[iG]} {e.r_Rg_best[iG]:.2f}', B_r_Rg=e.r_Rg_best[iG])
        else:
            row.update(B_best_L_Rb='—', B_r_Rb=None, B_best_L_Rg='—', B_r_Rg=None)
        x = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')[COLS]
        y = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir']).R_bir.to_numpy()
        z = ((x - x.mean()) / x.std()).to_numpy()
        N = len(z)
        rb = lambda c: abs(np.corrcoef(x[c], y)[0, 1])
        X1 = np.column_stack([z, np.ones(N)])
        ceil = np.corrcoef(X1 @ np.linalg.lstsq(X1, y, rcond=None)[0], y)[0, 1]
        # BIC difference per star: (NLL_B − NLL_A) + (edges_B − 15)·ln N / (2N); negative means the learned graph is better
        dbic = (info['nll'] - a.best_nll + (info['n_edges'] - 15) * np.log(N) / (2 * N)) if info else None
        row.update(A_r1_range=f'{a.r1_min:.2f}–{a.r1_max:.2f}', A_r2_range=f'{a.r2_min:.2f}–{a.r2_max:.2f}',
                   NLL_B=info['nll'] if info else None, NLL_A=a.best_nll, dBIC_per_star=dbic,
                   age_only=rb('age'), JPhi_only=rb('JPhi'), ceiling=ceil)
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(HERE.with_name(f'summary_planB_{v}.csv'), index=False)
    show = tab.drop(columns=['B_r_Rb', 'B_r_Rg'])
    txt = (f'Approach B summary (main selection 7<R<10 kpc, |z|<0.3 kpc; {v}); per galaxy only the most frequent of the 30 seeds\' structures\n'
           'B_best_L_Rb / _Rg: latent with the largest |r| with the true birth / guiding radius (best case, chosen with the truth) and its |r| (maximum likelihood)\n'
           'A_r1_range / A_r2_range: Approach A (authors\' graph frozen), range of |r(L1,birth radius)| / |r(L2,guiding radius)| among indistinguishable solutions\n'
           'NLL_B / NLL_A: maximum-likelihood NLL of the learned / authors\' graph (smaller = closer to the data); dBIC_per_star > 0 means the authors\' graph still fits better after accounting for the number of edges\n\n'
           + show.round(3).to_string(index=False) + '\n')
    if tab.B_r_Rb.notna().any():
        ok = tab.dropna(subset=['B_r_Rb'])
        txt += (f'\ngalaxies where coefficients can be fitted: {len(ok)}; median |r(L,birth radius)| of the best-matching latent {ok.B_r_Rb.median():.2f} '
                f'(range {ok.B_r_Rb.min():.2f}–{ok.B_r_Rb.max():.2f}); galaxies beating "age alone": '
                f'{", ".join(ok.galaxy[ok.B_r_Rb > ok.age_only]) or "none"}\n')
    HERE.with_name(f'summary_planB_{v}.txt').write_text(txt)
    print(txt)


if __name__ == '__main__':
    main()
