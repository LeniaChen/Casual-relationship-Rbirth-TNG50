# -*- coding: utf-8 -*-
"""
preview_selection.py — quick preview before widening the star selection: spread of the birth radius and its
correlation with the inputs under different selections

Run (h5py or pyfive; numpy):
    cd <repo root>
    $PY scripts/local/preview_selection.py            # all 12 galaxies
    $PY scripts/local/preview_selection.py 538905     # selected galaxies
Output: screen, and appended to results/selection_preview.txt

For each selection (base cuts as in make_rlcd_input.py: bound, has a birth radius < 20 kpc, J_z > 0, [Fe/H] > -0.5):
  N: number of stars;  std(Rb): standard deviation of the true birth radius (kpc)
  |r| age / FeH / MgH: absolute Pearson correlation of a single input with the true birth radius
  linear ceiling: correlation of a linear regression of the true birth radius on the 6 inputs
                  (log10 Jz, Jphi, e, age, [Fe/H], [Mg/H]), i.e. the upper limit for any L1
                  (L1 is a linear combination of the 6 inputs)
The spread of f(L1) = |r(L1, Rb)| x std(Rb): std(Rb) sets how wide it can be, r sets how wide it actually is.
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
