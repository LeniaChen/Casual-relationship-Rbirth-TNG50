# -*- coding: utf-8 -*-
"""
rlcd_identifiability_ci.py —— 严格版的"可能范围": 剖面似然的 95% 置信区域 (rlcd_identifiability.py 的补充)

为什么需要这个脚本:
  rlcd_identifiability.py 用"NLL 与最优相差 < 0.01"界定"同样好的解"。但 NLL 是按每颗星平均的,
  对 7 千–12 万颗星来说, 相差 0.01 在统计上其实是分得开的。统计上"分不开"的标准是似然比检验:
      2·N·ΔNLL < 3.84
  意思是: 两组系数的拟合差别, 小到在 95% 的置信水平上不能说哪个更好。
      N: 恒星数;  ΔNLL: 与最优 NLL 的差 (NLL 按每颗星平均);  3.84: 自由度为 1 的 χ² 分布的 95% 分位数
  另外, 雅可比矩阵秩为 14 说明存在"拟合完全相同"的方向; 但最优解通常卡在"某个噪声方差 = 0"的边界上
  (因子分析里叫 Heywood 情形), 沿这个方向走很快就会让噪声方差变负, 所以实际能走多远, 要看拟合变差得多快。

做法:
  1. 把 <tag>_multistart.csv 里每个有效解 (30 个随机初值的结果) 用很严的收敛条件再优化一次 (polish),
     消除"优化器没收敛"造成的 NLL 差别; 可以用 --extra 再加更多随机初值 (g2.79e12 用了 --extra 200,
     因为它的似然面上有很多彼此很接近、但统计上分得开的局部最优);
  2. 保留与最优"统计上分不开" (2·N·ΔNLL < 3.84) 的解;
  3. 统计这些解的 |r(L1, 真 Rbirth)| 和 |r(L2, 真 Rg)| 的范围。
  这是"目前能找到的"范围: 初值有限, 如果还有没找到的更好的解, "最优"和这个范围都可能变 (g2.79e12 就是例子)。
  约束: 系数在 ±0.999 之内; 每个节点噪声方差 >= 0 (物理要求; rlcd_identifiability.py 用的是 >= 0.01 的软约束,
  对 TNG50 几乎没有影响, 但 g2.79e12 的最优解要求 [Fe/H] 的噪声方差 ≈ 0.0005, 所以这里用 >= 0)。

运行 (需要先跑完 rlcd_identifiability.py):
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/rlcd_identifiability_ci.py --version g279 --extra 200
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/rlcd_identifiability_ci.py --version nonoise --ids 538905
    ~/PYTHONPATH_01/venv/bin/python3 scripts/local/rlcd_identifiability_ci.py --version nonoise      # 12 个星系
输出: results/rlcd_A/identifiability/<tag>_ci.csv (polish 后的每个解), <tag>_ci.txt; 汇总 summary_ci_<version>.csv/.txt/.png
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rlcd_identifiability import COLS, DOT, OUT, load, setup, implied, nll, latents  # noqa: E402

IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
LR95 = 3.84      # 自由度 1 的 χ² 分布 95% 分位数
PEN = 1e7        # 噪声方差 < 0 的罚项系数
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['g279', 'nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--dot', default=DOT)
    ap.add_argument('--summary-only', action='store_true', help='只重新生成汇总表')
    ap.add_argument('--extra', type=int, default=0, help='另外再加多少个随机初值 (直接用严格条件优化)')
    args = ap.parse_args()
    ids = ['g279'] if args.version == 'g279' else (args.ids or IDS)
    if args.summary_only:
        write_summary(args.version)
        return
    names, idx, edges, order = setup(args.dot)
    n, ne = len(names), len(edges)
    xi = [idx['X_' + c] for c in COLS]
    lab = [f'{names[a]}->{names[b]}' for a, b in edges]

    for gid in ids:
        tag = gid if args.version == 'g279' else f'{args.version}_{gid}'
        X, T = load(args.version, gid)
        N = len(X)
        Semp = X[COLS].cov().to_numpy()
        ms = pd.read_csv(os.path.join(OUT, f'{tag}_multistart.csv'))
        ms = ms[ms.valid]

        # 1. polish (物理约束: 每个节点噪声方差 >= 0, 用很大的罚项近似硬约束)
        def nll_phys(th):
            _, _, om = implied(th, n, edges, order)
            return nll(th, n, edges, order, xi, Semp, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

        starts = [row[lab].to_numpy(float) for _, row in ms.iterrows()]
        starts += [np.tanh(np.random.default_rng(5000 + i).normal(0, 0.5, ne)) for i in range(args.extra)]
        pol = []
        for i_start, th_start in enumerate(starts):
            r = minimize(nll_phys, th_start, method='L-BFGS-B',
                         bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
            F, S, om = implied(r.x, n, edges, order)
            if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
                continue
            pol.append((nll(r.x, n, edges, order, xi, Semp, pen=0), r.x, i_start))
        pol.sort(key=lambda z: z[0])
        best = pol[0][0]

        rows = []
        for v, th, i_start in pol:
            F, S, om = implied(th, n, edges, order)
            lat = latents(X, names, idx, F)
            rows.append(dict(start=i_start, nll=v, lr=2 * N * (v - best), omega_min=om.min(), omega_min_node=names[int(np.argmin(om))],
                             abs_r_L1_Rbir=abs(np.corrcoef(lat['L1'], T['R_bir'])[0, 1]),
                             abs_r_L2_Rg=abs(np.corrcoef(lat['L2'], T['Rg'])[0, 1]),
                             **{l: val for l, val in zip(lab, th)}))

        ci = pd.DataFrame(rows).sort_values('nll')
        ci.to_csv(os.path.join(OUT, f'{tag}_ci.csv'), index=False)
        inn = ci[ci.lr < LR95]
        lines = [f'{tag}: {N} 颗星; polish 后有效解 {len(ci)} 个, 最优 NLL {best:.6f}',
                 f'与最优统计上分不开 (2·N·ΔNLL < {LR95}) 的解 {len(inn)} 个: '
                 f'L1->age 从 {inn["L1->X_age"].min():+.2f} 到 {inn["L1->X_age"].max():+.2f}',
                 f'它们的 |r(L1, Rbirth)| 范围 {inn.abs_r_L1_Rbir.min():.3f} – {inn.abs_r_L1_Rbir.max():.3f}; '
                 f'|r(L2, Rg)| 范围 {inn.abs_r_L2_Rg.min():.3f} – {inn.abs_r_L2_Rg.max():.3f}',
                 (f'只看第一轮 multistart 的 {len(ms)} 个有效解时的最优解: NLL {ci[ci.start < len(ms)].nll.iloc[0]:.6f}, '
                  f'|r(L1, Rbirth)| {ci[ci.start < len(ms)].abs_r_L1_Rbir.iloc[0]:.3f}, |r(L2, Rg)| {ci[ci.start < len(ms)].abs_r_L2_Rg.iloc[0]:.3f}'
                  if args.extra > 0 else '（没有另加初值）'),
                 f'最优解处噪声方差的最小值 {ci.omega_min.iloc[0]:.4f} ({ci.omega_min_node.iloc[0]}; ≈0 表示最优解卡在"噪声方差 ≥ 0"的边界上, 即某个观测量几乎完全由它的父节点决定)']
        txt = '\n'.join(lines)
        print(txt + '\n')
        open(os.path.join(OUT, f'{tag}_ci.txt'), 'w').write(txt + '\n')

    if args.version != 'g279':
        write_summary(args.version)


def write_summary(version):
    """读已有的 <tag>_ci.csv, 汇总成 summary_ci_<version>.csv/.txt (g2.79e12 放第一行)"""
    rows = []
    for tag in ['g279'] + [f'{version}_{g}' for g in IDS]:
        f = os.path.join(OUT, f'{tag}_ci.csv')
        if not os.path.exists(f):
            continue
        c = pd.read_csv(f).sort_values('nll')
        inn = c[c.lr < LR95]
        rows.append(dict(tag=tag, n_polished=len(c), n_indistinguishable=len(inn),
                         r1_min=inn.abs_r_L1_Rbir.min(), r1_max=inn.abs_r_L1_Rbir.max(),
                         r2_min=inn.abs_r_L2_Rg.min(), r2_max=inn.abs_r_L2_Rg.max(),
                         best_nll=c.nll.iloc[0], omega_min_at_best=c.omega_min.iloc[0],
                         omega_min_node=c.omega_min_node.iloc[0]))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, f'summary_ci_{version}.csv'), index=False)
    head = (f'严格标准: polish 后与最优统计上分不开 (2·N·ΔNLL < {LR95}, 95% 似然比检验; 噪声方差 >= 0) 的解, '
            f'|r(L1, 真 Rbirth)| 和 |r(L2, 真 Rg)| 的范围 (初值有限: 如果还有没找到的更好的解, 最优和范围都可能变)\n'
            f'omega_min_node: 最优解里噪声方差最小的节点 (≈0 表示它几乎完全由父节点决定)\n\n')
    txt = head + tab.round(3).to_string(index=False) + '\n'
    open(os.path.join(OUT, f'summary_ci_{version}.txt'), 'w').write(txt)
    print(txt)

    # 图: 每个星系一行; 粗蓝线 = 严格标准下 |r(L1,Rb)| 的范围, 细灰线 = 宽松标准 (ΔNLL<0.01) 的范围,
    #     黑三角 = 单用年龄的 |r(age,Rb)|, 叉 = 单用 [Fe/H], 空心圆 = 6 个量线性回归的上限
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
    loose = os.path.join(OUT, f'summary_{version}.csv')
    loose = pd.read_csv(loose).set_index('tag') if os.path.exists(loose) else None
    ceil = os.path.join(OUT, 'linear_ceiling.csv')
    ceil = pd.read_csv(ceil).set_index('tag') if os.path.exists(ceil) else None
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(tab) + 1.8))
    y = np.arange(len(tab))[::-1]
    for yi, (_, r) in zip(y, tab.iterrows()):
        if loose is not None and r.tag in loose.index:
            ax.plot([loose.loc[r.tag, 'r1_min'], loose.loc[r.tag, 'r1_max']], [yi] * 2, color=MUTED, lw=1.5)
        x0, x1 = r.r1_min, max(r.r1_max, r.r1_min + 0.006)
        ax.plot([x0, x1], [yi] * 2, color=BLUE, lw=6, solid_capstyle='round')
        if ceil is not None and r.tag in ceil.index:
            ax.plot(ceil.loc[r.tag, 'r_age_Rbir'], yi, marker='v', color=INK, ms=6, ls='none')
            ax.plot(ceil.loc[r.tag, 'r_FeH_Rbir'], yi, marker='x', color=INK, ms=6, ls='none')
            ax.plot(ceil.loc[r.tag, 'r_ceiling_Rbir'], yi, marker='o', mfc='none', mec=INK, ms=6, ls='none')
    ax.axvline(0.64, color=MUTED, ls='--', lw=1)
    ax.text(0.645, y.max() + 0.55, 'paper 0.64', color=MUTED, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace(f'{version}_', '') for t in tab.tag], fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xlabel('|Pearson r| with the true birth radius')
    ax.plot([], [], color=BLUE, lw=6, label='L1, statistically indistinguishable fits')
    ax.plot([], [], color=MUTED, lw=1.5, label='L1, loose criterion (dNLL<0.01)')
    ax.plot([], [], marker='v', color=INK, ls='none', label='age alone')
    ax.plot([], [], marker='x', color=INK, ls='none', label='[Fe/H] alone')
    ax.plot([], [], marker='o', mfc='none', mec=INK, ls='none', label='linear regression on all 6 inputs')
    ax.set_ylim(-2.6, len(tab) - 0.2)
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f'How well can L1 track the birth radius? ({version}, author graph frozen)', fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f'summary_ci_{version}.png'), dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
