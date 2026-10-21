# -*- coding: utf-8 -*-
"""
rlcd_identifiability_summary.py —— 汇总 12 个星系 (及作者的 g2.79e12) 的可识别性分析, 选项 (a) 的主结果

先对每个星系运行 rlcd_identifiability.py, 再运行本脚本 (numpy、pandas、matplotlib):
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/rlcd_identifiability_summary.py --version nonoise
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/rlcd_identifiability_summary.py --version noise

对每个星系, 在"同样好的拟合"里 (NLL 与该星系能达到的最优 NLL 相差 < DNLL) 统计:
  |r(L1, 真 Rbirth)| 和 |r(L2, 真 Rg)| 的最小值和最大值, 数据来自两部分:
    multistart: 30 个随机初值, L-BFGS 收敛到的解
    profile:    把 L1->age 固定在 -0.95..0.95 的 39 个值上, 其余 14 个系数重新拟合
  这个范围就是"数据本身允许的范围": 范围越宽, 说明 L1 的含义越不确定。

⚠ 问题标注: 作者的因果图冻结后, 15 条边系数只能由 14 个独立的数据约束确定 (雅可比矩阵秩 = 14),
  所以 L1 与出生半径的相关系数不是由数据唯一决定的。论文的 0.64 是这个范围里的一个点。

输出: results/rlcd_A/identifiability/summary_<version>.csv, .txt, summary_<version>.png
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
        m = re.findall(r'雅可比矩阵的秩 \(随机点 \d\): (\d+) / (\d+)', txt)
        rank = '/'.join(sorted({f'{a}/{b}' for a, b in m}))
        m2 = re.search(r': (\d+) 颗星', txt)
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
    head = ("⚠ 作者的因果图冻结后边系数不可识别 (雅可比矩阵秩 14/15)。下表给出在同样好的拟合 "
            f"(NLL 与最优相差 < {DNLL}) 里, |r(L1, 真 Rbirth)| 和 |r(L2, 真 Rg)| 能取到的范围。\n"
            "论文 (g2.79e12) 报告 r(L1, Rb) = 0.64, r(L2, Rg) ≈ 0.98 —— 那组系数并未收敛到最优 (NLL 1.59 vs 最优约 -1.27)。\n")
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
