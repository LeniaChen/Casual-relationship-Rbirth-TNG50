# -*- coding: utf-8 -*-
"""
compute_moments.py — compress the inputs and truth of each galaxy into a covariance matrix (step 2 of the selection R_0_30_Z_30)

Why: with this selection each galaxy has 0.37–0.72 million stars, 6 GB in total. But the identifiability analysis only needs
  the covariances, means and number of stars of 8 quantities: the 6 inputs (logJZ, JPhi, ecc, age, FeH, MgH) and the 2 truths (R_bir, Rg):
  - the maximum-likelihood fit uses only the correlation matrix of the 6 inputs (paper Sec. 3.3; parameter_by_trek in the repo also uses only the covariance);
  - L1, L2 are linear combinations L = w·z of the 6 standardized inputs, so
        r(L, R) = cov(L, R) / [std(L)·std(R)] = (w·c) / sqrt(w·C·w · var(R))
    meaning: the correlation of L with the truth R follows directly from covariances, without going star by star.
      w: weights of L on the 6 standardized inputs;  c: covariances of the 6 standardized inputs with R (6 numbers);
      C: covariance matrix of the 6 standardized inputs (i.e. their correlation matrix);  var(R): variance of R
  This is exactly equal to the star-by-star result (checked on the earlier 7<R<10 data, see identifiability_moments.py).

Run (numpy, pandas):
    cd <repo root>
    $PY star_condition/R_0_30_Z_30/output_data/compute_moments.py
Outputs: output_data/moments_nonoise.json, moments_noise.json
      per galaxy: N (number of stars), names (names of the 8 quantities), mean (8 means), cov (8×8 covariance matrix, denominator N−1)
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
COND = HERE.parents[1]
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
TRUTH = ['R_bir', 'Rg']


def main():
    versions = sys.argv[1:] or ['nonoise', 'noise']
    for v in versions:
        out = {}
        for gid in IDS:
            e = pd.read_csv(COND / 'input_data' / v / f'{gid}_e.csv', usecols=COLS)[COLS]
            t = pd.read_csv(COND / 'true_value' / v / f'{gid}_truth.csv', usecols=TRUTH)[TRUTH]
            X = np.column_stack([e.to_numpy(float), t.to_numpy(float)])
            assert np.isfinite(X).all()
            out[gid] = dict(N=len(X), names=COLS + TRUTH, mean=X.mean(0).tolist(), cov=np.cov(X.T, ddof=1).tolist())
            print(v, gid, len(X), flush=True)
        json.dump(out, open(HERE.with_name(f'moments_{v}.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
