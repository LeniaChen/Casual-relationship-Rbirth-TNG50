# -*- coding: utf-8 -*-
"""
preview_selection.py —— 放宽选星条件之前的快速预览: 不同选星范围下, 出生半径的散布和它与各输入量的相关

运行 (h5py 或 pyfive 都可以; numpy):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/preview_selection.py            # 12 个星系
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/preview_selection.py 538905     # 指定星系
输出: 屏幕, 并追加到 results/selection_preview.txt

对每种选星条件 (基础条件与 make_rlcd_input.py 相同: 束缚、有出生半径且 < 20 kpc、J_z > 0、[Fe/H] > -0.5), 给出:
  N: 星数;  std(Rb): 真出生半径的标准差 (kpc)
  |r| age / FeH / MgH: 单个输入量与真出生半径的皮尔逊相关系数的绝对值
  linear ceiling: 6 个输入量 (log10 Jz, Jphi, e, age, [Fe/H], [Mg/H]) 对真出生半径做线性回归的相关系数,
                  即任何 L1 (6 个输入量的线性组合) 能达到的上限
f(L1) 的散布 = |r(L1, Rb)| × std(Rb), 所以 std(Rb) 决定"最宽能有多宽", r 决定"实际有多宽"。
"""
import sys
import numpy as np
try:
    import h5py as H
except ImportError:
    import pyfive as H

IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
SEL = {'A 7<R<10,|z|<0.3 (current)': (7, 10, 0.3), 'B 7<R<10,|z|<1': (7, 10, 1.0),
       'C 7<R<10,no z cut': (7, 10, 99), 'D 4<R<12,|z|<1': (4, 12, 1.0), 'E 4<R<12,no z cut': (4, 12, 99)}


def main():
    ids = sys.argv[1:] or IDS
    P = H.File('data/tng50_extract/data_obtain_all.hdf5', 'r')
    O = H.File('data/tng50_derived/orbits_all.hdf5', 'r')
    out = open('results/selection_preview.txt', 'a')
    for g in ids:
        s, o = P[g]['stars'], O[g]
        pos = s['pos'][:]
        R, z = np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2]
        age, feh, mgfe, rb = s['age'][:], s['FeH'][:], s['MgFe'][:], s['birth_R_2D'][:]
        hb = s['has_birth'][:].astype(bool)
        E, jz, jp, ecc = o['E'][:], o['Jz'][:], o['Jphi'][:], o['ecc'][:]
        base = (hb & np.isfinite(rb) & (rb < 20) & (E < 0) & np.isfinite(jz) & (jz > 0) & np.isfinite(ecc)
                & np.isfinite(feh) & (feh > -0.5) & np.isfinite(mgfe))
        X_all = np.column_stack([np.log10(np.where(jz > 0, jz, 1)), jp, ecc, age, feh, feh + mgfe])
        lines = [f'== {g}']
        for name, (r0, r1, zm) in SEL.items():
            m = base & (R > r0) & (R < r1) & (np.abs(z) < zm)
            X, y = X_all[m], rb[m]
            A = np.column_stack([(X - X.mean(0)) / X.std(0), np.ones(len(y))])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            cc = lambda a: abs(np.corrcoef(a, y)[0, 1])
            lines.append(f'  {name:28s} N={m.sum():7d}  std(Rb)={y.std():.2f}  |r| age={cc(X[:, 3]):.2f} '
                         f'FeH={cc(X[:, 4]):.2f} MgH={cc(X[:, 5]):.2f}  linear ceiling={np.corrcoef(A @ c, y)[0, 1]:.2f}')
        print('\n'.join(lines))
        out.write('\n'.join(lines) + '\n')
    out.close()


if __name__ == '__main__':
    main()
