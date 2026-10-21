# -*- coding: utf-8 -*-
"""
linear_ceiling.py — linear ceiling: how well a least-squares linear regression on the 6 inputs predicts true Rbirth / Rg

Rationale: the latents L1, L2 computed in rlcd_A.py (hierarchical least squares) are linear combinations of the
6 standardized inputs, so for ANY set of edge coefficients |r(L1, Rbirth)| cannot exceed the correlation obtained by
regressing Rbirth directly on the 6 inputs (the multiple correlation coefficient). This number is therefore the best
the "frozen author graph, linear RLCD" approach could possibly do. Note that this regression uses the truth
(supervised); it is a reference only, not an RLCD result.

Run (numpy, pandas):
    $PY scripts/local/linear_ceiling.py
Output: results/rlcd_A/identifiability/linear_ceiling.csv, .txt
"""
import os
import numpy as np
import pandas as pd

IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
OUT = 'results/rlcd_A/identifiability'


def multiple_r(X, y):
    X = (X - X.mean(0)) / X.std(0)
    A = np.column_stack([X, np.ones(len(X))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return np.corrcoef(A @ coef, y)[0, 1]


def row(tag, x, rb, rg):
    X = x[COLS].to_numpy(float)
    return dict(tag=tag, n=len(x),
                r_ceiling_Rbir=multiple_r(X, rb), r_ceiling_Rg=multiple_r(X, rg),
                r_JPhi_Rg=abs(np.corrcoef(x.JPhi, rg)[0, 1]),
                r_FeH_Rbir=abs(np.corrcoef(x.FeH, rb)[0, 1]),
                r_age_Rbir=abs(np.corrcoef(x.age, rb)[0, 1]))


def main():
    rows = []
    g = 'Causal-Galactic-Archaeology/data/'
    if os.path.exists(g + 'g2.79e12_e.csv'):
        x = pd.read_csv(g + 'g2.79e12_e.csv')
        t = pd.read_csv(g + 'g2.79e12_multi_ele_action_Jin_SN.csv')
        rows.append(row('g279', x, t.R_bir.to_numpy(), t.Rg.to_numpy()))
    for v in ['nonoise', 'noise']:
        for gid in IDS:
            x = pd.read_csv(f'data/rlcd_input/{v}/{gid}_e.csv')
            t = pd.read_csv(f'data/rlcd_input/{v}/{gid}_truth.csv')
            rows.append(row(f'{v}_{gid}', x, t.R_bir.to_numpy(), t.Rg.to_numpy()))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, 'linear_ceiling.csv'), index=False)
    txt = ('Linear ceiling: correlation of a linear regression of the truth on the 6 inputs (supervised, reference only).\n'
           'L1 and L2 are linear combinations of the 6 inputs, so |r(L1,Rbirth)| <= r_ceiling_Rbir and |r(L2,Rg)| <= r_ceiling_Rg.\n\n'
           + tab.round(3).to_string(index=False) + '\n')
    open(os.path.join(OUT, 'linear_ceiling.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
