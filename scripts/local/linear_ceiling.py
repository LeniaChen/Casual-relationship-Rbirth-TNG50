# -*- coding: utf-8 -*-
"""
linear_ceiling.py —— 线性上限: 用 6 个输入量做最小二乘线性回归, 能把真 Rbirth / 真 Rg 预测到多好

原理: rlcd_A.py 求出的 L1、L2 (分层最小二乘) 都是 6 个标准化输入量的线性组合,
所以任何一组系数给出的 |r(L1, Rbirth)| 都不可能超过"直接用 6 个量对 Rbirth 做线性回归"的相关系数 (复相关系数)。
这个数给出了"冻结作者的图、线性 RLCD"最好能做到多好。注意这个回归用了真值 (有监督), 只作参照, 不是 RLCD 的结果。

运行 (numpy、pandas):
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/linear_ceiling.py
输出: results/rlcd_A/identifiability/linear_ceiling.csv, .txt
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
    txt = ('线性上限: 6 个输入量对真值做线性回归的相关系数 (有监督, 只作参照)。\n'
           'L1、L2 都是 6 个输入量的线性组合, 所以 |r(L1,Rbirth)| <= r_ceiling_Rbir, |r(L2,Rg)| <= r_ceiling_Rg。\n\n'
           + tab.round(3).to_string(index=False) + '\n')
    open(os.path.join(OUT, 'linear_ceiling.txt'), 'w').write(txt)
    print(txt)


if __name__ == '__main__':
    main()
