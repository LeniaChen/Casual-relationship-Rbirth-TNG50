# -*- coding: utf-8 -*-
"""
rlcd_identifiability_summary.py — summary of the FIRST-ROUND identifiability analysis of the 12 galaxies (and the authors' g2.79e12); superseded by the strict second round (rlcd_identifiability_ci.py)

Run rlcd_identifiability.py for every galaxy first, then this script (numpy, pandas, matplotlib):
    $PY scripts/local/rlcd_identifiability_summary.py --version nonoise
    $PY scripts/local/rlcd_identifiability_summary.py --version noise

For each galaxy, among "equally good fits" (NLL within DNLL of the best NLL reached for that galaxy), report
  the minimum and maximum of |r(L1, true Rbirth)| and |r(L2, true Rg)|, using two sets of solutions:
    multistart: solutions reached by L-BFGS from 30 random starts
    profile:    L1->age fixed at 39 values in -0.95..0.95, the other 14 coefficients refitted
  This range is what the data alone allow: the wider it is, the less determined the meaning of L1.
  NOTE: DNLL = 0.01 is too loose for ~1e4 stars (it equals a likelihood-ratio statistic of 0.02 N); the strict
  criterion 2 N dNLL < 3.84 is used in rlcd_identifiability_ci.py.

⚠ With the authors' causal graph frozen, the 15 edge coefficients are constrained by only 14 independent data
  constraints (Jacobian rank = 14), so the correlation of L1 with the birth radius is not determined by the data. The paper's 0.64 is one point in this range.

Output: results/rlcd_A/identifiability/summary_<version>.csv, .txt, summary_<version>.png
"""
import os
import argparse
import numpy as np
import pandas as pd

D = 'results/rlcd_A/identifiability'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
DNLL = 0.01


def stats(tag):
    ms = pd.read_csv(os.path.join(D, f'{tag}_multistart.csv'))
    pf = pd.read_csv(os.path.join(D, f'{tag}_profile.csv'))
    ms, pf = ms[ms.valid], pf[pf.valid]
    best = min(ms.nll.min(), pf.nll.min())
    good = pd.concat([ms[ms.nll < best + DNLL][['abs_r_L1_Rbir', 'abs_r_L2_Rg']],
                      pf[pf.nll < best + DNLL][['abs_r_L1_Rbir', 'abs_r_L2_Rg']]])
    rank = None
    sfile = os.path.join(D, f'{tag}_summary.txt')
    nstar = None
    if os.path.exists(sfile):
        txt = open(sfile).read()
        import re
        m = re.findall(r'Jacobian rank \(random point \d\): (\d+) / (\d+)', txt)
        rank = '/'.join(sorted({f'{a}/{b}' for a, b in m}))
        m2 = re.search(r': (\d+) stars', txt)
        nstar = int(m2.group(1)) if m2 else None
    return dict(tag=tag, n_stars=nstar, jacobian_rank=rank, best_nll=best, n_equally_good=len(good),
                r1_min=good.abs_r_L1_Rbir.min(), r1_max=good.abs_r_L1_Rbir.max(),
                r2_min=good.abs_r_L2_Rg.min(), r2_max=good.abs_r_L2_Rg.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['nonoise', 'noise'])
    args = ap.parse_args()
    rows = []
    if os.path.exists(os.path.join(D, 'g279_multistart.csv')):
        rows.append(stats('g279'))
    for g in IDS:
        tag = f'{args.version}_{g}'
        if os.path.exists(os.path.join(D, f'{tag}_multistart.csv')):
            rows.append(stats(tag))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(D, f'summary_{args.version}.csv'), index=False)
    head = ("⚠ With the authors' graph frozen the edge coefficients are NOT identifiable (Jacobian rank 14/15). The table gives the range "
            f"of |r(L1, true Rbirth)| and |r(L2, true Rg)| among equally good fits (NLL within {DNLL} of the best).\n"
            "The paper (g2.79e12) reports r(L1, Rb) = 0.64, r(L2, Rg) ~ 0.98 -- those coefficients are not converged (NLL 1.59 vs best ~ -1.27).\n")
    txt = head + '\n' + tab.round(3).to_string(index=False) + '\n'
    with open(os.path.join(D, f'summary_{args.version}.txt'), 'w') as f:
        f.write(txt)
    print(txt)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(tab) + 1.6))
    y = np.arange(len(tab))[::-1]
    for yi, (_, r) in zip(y, tab.iterrows()):
        ax.plot([r.r1_min, r.r1_max], [yi + 0.12] * 2, color=BLUE, lw=5, solid_capstyle='round')
        ax.plot([r.r2_min, r.r2_max], [yi - 0.12] * 2, color=INK, lw=5, solid_capstyle='round')
    ax.axvline(0.64, color=MUTED, ls='--', lw=1)
    ax.text(0.645, y.max() + 0.55, 'paper r(L1,Rb)=0.64', color=MUTED, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace(f'{args.version}_', '') for t in tab.tag], fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel('|Pearson r| reachable by equally good fits')
    ax.plot([], [], color=BLUE, lw=5, label='|r(L1, true Rbirth)|')
    ax.plot([], [], color=INK, lw=5, label='|r(L2, true Rg)|')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f'NOT identifiable with the frozen author graph ({args.version}):\n'
                 'range of correlation with the truth among equally good fits', fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(D, f'summary_{args.version}.png'), dpi=130)


if __name__ == '__main__':
    main()
