# -*- coding: utf-8 -*-
"""
identifiability_moments.py —— 选星条件 R_0_30_Z_0.3 下的做法 A (冻结作者的因果图) + 严格的可识别性分析 (第 3 步)

方法与 scripts/local/rlcd_identifiability.py + rlcd_identifiability_ci.py 完全相同 (REPORT_RLCD_TNG50.md 第 3.9 节),
只是输入换成了 compute_moments.py 算好的协方差矩阵 (每个星系 数量见 input_data/selection.txt, 逐星计算没有必要):
  1. 雅可比矩阵的秩: 在 3 组随机系数处, 算"15 个边系数 → 15 个观测相关系数"的雅可比矩阵的秩 (作者的图应为 14)。
  2. 第一轮 30 个随机初值 (种子 1000+s), L-BFGS-B 最大似然, 噪声方差 >= 0.01 的软约束 (与原脚本相同);
  3. 严格收敛: 这 30 个解 + 200 个新随机初值 (种子 5000+i), 全部在"噪声方差 >= 0"(罚项 1e7) 下严格收敛 (polish);
  4. 与最优统计上分不开: 2·N·ΔNLL < 3.84 (95% 似然比检验);
       N: 星数;  ΔNLL: 与最优解的 NLL 之差 (NLL 按每颗星平均);  3.84: 自由度 1 的 χ² 分布的 95% 分位数
  5. 每个解的 |r(L1, 真 Rbirth)|、|r(L2, 真 Rg)|: L1、L2 用作者 notebook 的分层最小二乘求出, 它们是标准化输入 z 的
     线性组合 L = w·z, 所以
         r(L, R) = (w·c) / sqrt(w·C·w · var(R))
     意思是: 相关系数 = L 与 R 的协方差 ÷ (L 的标准差 × R 的标准差), 全部可以由协方差矩阵算出。
       w: L 对 6 个标准化输入的权重; c: 标准化输入与 R 的协方差; C: 标准化输入的相关矩阵; var(R): R 的方差
     (--selftest 会在旧的 7<R<10 数据上核对: 与逐星计算的结果一致)
  6. 参照线: 单用年龄、[Fe/H]、[Mg/H]、J_phi 与 Rbirth 的 |r|; 6 个输入量线性回归的上限 sqrt(c·C⁻¹·c / var(R))
     (意思是: 用 6 个输入量的最佳线性组合去预测出生半径, 能达到的相关系数; c、C、var(R) 同上)。
     ⚠ 这个条件包含整个盘 (0–30 kpc, |z|<0.3 kpc), 现在的角动量 J_phi (≈ 现在的半径) 与出生半径本身就强相关 (大多 0.8–0.9),
       所以 L1 必须和"单用 J_phi"比; 论文只用 7–10 kpc 的星, 正是为了去掉这部分"平凡"的信息。

运行 (numpy、scipy、pandas、matplotlib; 主 venv):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_0.3/output_data/identifiability_moments.py --version nonoise
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_0.3/output_data/identifiability_moments.py --version noise
    (每个星系约 5–10 分钟; 可以用 --ids 538905 只算一个; --extra 200 是新初值的个数)
输出 (本文件夹):
  <version>_<ID>_ci.csv      严格收敛后的全部有效解: start (初值编号, <30 为第一轮), nll, lr = 2·N·ΔNLL, omega_min,
                             omega_min_node (噪声方差最小的节点), abs_r_L1_Rbir, abs_r_L2_Rg, 15 个边系数
  <version>_<ID>_summary.txt 雅可比矩阵的秩、不可识别方向、统计上分不开的解的范围
  summary_ci_<version>.csv/.txt/.png   12 个星系的汇总 (严格范围 + 参照线)
"""
import os
import sys
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve()
COND = HERE.parents[1]
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'scripts' / 'local'))
from rlcd_identifiability import setup, implied, nll  # noqa: E402

DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LR95 = 3.84
PEN = 1e7
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def standardized_moments(m):
    """由 8×8 协方差得到: C (6 个输入的相关矩阵), c_R / c_G (标准化输入与 R_bir / Rg 的协方差), var_R, var_G"""
    cov = np.array(m['cov'])
    sd = np.sqrt(np.diag(cov)[:6])
    C = cov[:6, :6] / np.outer(sd, sd)
    c_R, c_G = cov[:6, 6] / sd, cov[:6, 7] / sd
    return C, c_R, c_G, cov[6, 6], cov[7, 7]


def latent_weights(coef):
    """作者 notebook 的分层最小二乘, 写成对 6 个标准化输入的权重。coef: {(父, 子): 系数}, 名字不带 X_"""
    w = {c: np.eye(6)[i] for i, c in enumerate(COLS)}

    def solve(L):
        num, den = np.zeros(6), 0.0
        for (p, c), a in coef.items():
            if p != L:
                continue
            r = w[c].copy()
            for (q, c2), b in coef.items():
                if c2 == c and q != L:
                    r = r - b * w[q]
            num, den = num + a * r, den + a * a
        return num / den
    w['L2'] = solve('L2')
    w['L1'] = solve('L1')
    return w['L1'], w['L2']


def corr_from_weights(wv, C, c, var):
    return abs(wv @ c) / np.sqrt((wv @ C @ wv) * var)


def analyse(tag, m, names, idx, edges, order, extra):
    n, ne = len(names), len(edges)
    xi = [idx['X_' + c] for c in COLS]
    iu = np.triu_indices(6, 1)
    lab = [f'{names[a]}->{names[b]}' for a, b in edges]
    N = m['N']
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    Semp = C
    lines = [f'{tag}: {N} 颗星; 冻结的图 {DOT.name}; 边 {ne} 条, 观测相关系数 {len(iu[0])} 个',
             f'理论最小 NLL = 0.5*logdet(S)+3 = {0.5*np.linalg.slogdet(Semp)[1]+3:.4f}']

    # 1. Jacobian rank
    cv = lambda th: implied(th, n, edges, order)[1][np.ix_(xi, xi)][iu]
    rng = np.random.default_rng(0)
    for trial in range(3):
        th = np.tanh(rng.normal(0, 0.5, ne))
        J = np.array([(cv(th + 1e-6 * np.eye(ne)[k]) - cv(th - 1e-6 * np.eye(ne)[k])) / 2e-6 for k in range(ne)]).T
        sv = np.linalg.svd(J, compute_uv=False)
        lines.append(f'雅可比矩阵的秩 (随机点 {trial+1}): {int(np.sum(sv > 1e-8 * sv[0]))} / {ne}')
    null = np.linalg.svd(J)[2][-1]
    lines.append('不可识别方向 (沿这个方向同时改系数, 观测相关矩阵不变), 分量从大到小:')
    for l, x in sorted(zip(lab, null), key=lambda z: -abs(z[1]))[:6]:
        lines.append(f'    {l:18s} {x:+.3f}')

    # 2. 第一轮 multistart (软约束, 与 rlcd_identifiability.py 相同)
    starts = []
    for s in range(30):
        th0 = np.tanh(np.random.default_rng(1000 + s).normal(0, np.sqrt(2 / n), ne))
        r = minimize(nll, th0, args=(n, edges, order, xi, Semp), method='L-BFGS-B',
                     bounds=[(-0.999, 0.999)] * ne, options={'maxiter': 5000})
        F, S, om = implied(r.x, n, edges, order)
        if om.min() > -1e-6 and np.linalg.slogdet(S[np.ix_(xi, xi)])[0] > 0:
            starts.append(r.x)
    n_first = len(starts)
    starts += [np.tanh(np.random.default_rng(5000 + i).normal(0, 0.5, ne)) for i in range(extra)]

    # 3. 严格收敛, 噪声方差 >= 0
    def nll_phys(th):
        _, _, om = implied(th, n, edges, order)
        return nll(th, n, edges, order, xi, Semp, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

    rows = []
    for i, th0 in enumerate(starts):
        r = minimize(nll_phys, th0, method='L-BFGS-B', bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
        F, S, om = implied(r.x, n, edges, order)
        if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
            continue
        coef = {(names[a].replace('X_', ''), names[b].replace('X_', '')): v for (a, b), v in zip(edges, r.x)}
        w1, w2 = latent_weights(coef)
        rows.append(dict(start=i, nll=nll(r.x, n, edges, order, xi, Semp, pen=0), omega_min=om.min(),
                         omega_min_node=names[int(np.argmin(om))],
                         abs_r_L1_Rbir=corr_from_weights(w1, C, c_R, var_R),
                         abs_r_L2_Rg=corr_from_weights(w2, C, c_G, var_G),
                         **{l: v for l, v in zip(lab, r.x)}))
    ci = pd.DataFrame(rows).sort_values('nll').reset_index(drop=True)
    ci.insert(2, 'lr', 2 * N * (ci.nll - ci.nll.iloc[0]))
    ci.to_csv(HERE.with_name(f'{tag}_ci.csv'), index=False)
    inn = ci[ci.lr < LR95]
    first = ci[ci.start < n_first]
    lines += [f'严格收敛后有效解 {len(ci)} 个 (第一轮 {n_first} 个 + 新初值 {extra} 个), 最优 NLL {ci.nll.iloc[0]:.6f}',
              f'与最优统计上分不开 (2·N·ΔNLL < {LR95}) 的解 {len(inn)} 个: '
              f'L1->age 从 {inn["L1->X_age"].min():+.2f} 到 {inn["L1->X_age"].max():+.2f}',
              f'它们的 |r(L1, Rbirth)| 范围 {inn.abs_r_L1_Rbir.min():.3f} – {inn.abs_r_L1_Rbir.max():.3f}; '
              f'|r(L2, Rg)| 范围 {inn.abs_r_L2_Rg.min():.3f} – {inn.abs_r_L2_Rg.max():.3f}',
              f'只看第一轮的最优解: NLL {first.nll.iloc[0]:.6f}, |r(L1, Rbirth)| {first.abs_r_L1_Rbir.iloc[0]:.3f}, '
              f'|r(L2, Rg)| {first.abs_r_L2_Rg.iloc[0]:.3f}',
              f'最优解处噪声方差的最小值 {ci.omega_min.iloc[0]:.4f} ({ci.omega_min_node.iloc[0]}; ≈0 表示该观测量被当成父节点的精确线性组合)']
    txt = '\n'.join(lines)
    print(txt + '\n', flush=True)
    open(HERE.with_name(f'{tag}_summary.txt'), 'w').write(txt + '\n')


def baselines(m):
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    r = lambda i: abs(c_R[i]) / np.sqrt(var_R)
    ceil = np.sqrt(c_R @ np.linalg.solve(C, c_R) / var_R)
    ceil_g = np.sqrt(c_G @ np.linalg.solve(C, c_G) / var_G)
    return dict(r_age_Rbir=r(3), r_FeH_Rbir=r(4), r_MgH_Rbir=r(5), r_JPhi_Rbir=r(1),
                r_best_single_Rbir=max(r(i) for i in range(6)), best_single=COLS[int(np.argmax([r(i) for i in range(6)]))],
                r_ceiling_Rbir=ceil, r_ceiling_Rg=ceil_g,
                r_JPhi_Rg=abs(c_G[1]) / np.sqrt(var_G), std_Rbir=np.sqrt(var_R))


def write_summary(version, moments):
    rows = []
    for gid in IDS:
        f = HERE.with_name(f'{version}_{gid}_ci.csv')
        if not f.exists():
            continue
        c = pd.read_csv(f)
        inn = c[c.lr < LR95]
        rows.append(dict(galaxy=gid, n_stars=moments[gid]['N'], n_valid=len(c), n_indistinguishable=len(inn),
                         r1_best=c.abs_r_L1_Rbir.iloc[0], r1_min=inn.abs_r_L1_Rbir.min(), r1_max=inn.abs_r_L1_Rbir.max(),
                         r2_best=c.abs_r_L2_Rg.iloc[0], r2_min=inn.abs_r_L2_Rg.min(), r2_max=inn.abs_r_L2_Rg.max(),
                         omega_min_at_best=c.omega_min.iloc[0], omega_min_node=c.omega_min_node.iloc[0],
                         **baselines(moments[gid])))
    tab = pd.DataFrame(rows)
    tab.to_csv(HERE.with_name(f'summary_ci_{version}.csv'), index=False)
    head = (f'选星条件 R_0_30_Z_0.3 (0<R<30 kpc, |z|<0.3 kpc), {version}; 作者的因果图冻结\n'
            f'r1/r2_best: 最大似然解的 |r(L1,Rbirth)| / |r(L2,Rg)|;  r1/r2_min–max: 与最优统计上分不开 (2·N·ΔNLL<{LR95}) 的解的范围\n'
            'r_age / r_FeH / r_MgH / r_JPhi _Rbir: 单个输入量与出生半径的 |r|;  r_best_single: 6 个输入量中单个最好的 (best_single 是哪个)\n'
            'r_ceiling: 6 个输入量线性回归的上限;  std_Rbir: 出生半径的标准差 [kpc]\n'
            '⚠ 这个选星范围包含整个盘 (0–30 kpc, |z|<0.3 kpc), 现在的角动量 J_phi (≈ 现在的半径) 本身就和出生半径强相关, '
            '所以要和"单用 J_phi"比, 而不只是和年龄、[Fe/H] 比\n\n')
    txt = head + tab.round(3).to_string(index=False) + '\n'
    open(HERE.with_name(f'summary_ci_{version}.txt'), 'w').write(txt)
    print(txt)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    INK, BLUE, MUTED = '#1f2937', '#2563eb', '#9ca3af'
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(tab) + 2.2))
    y = np.arange(len(tab))[::-1]
    for yi, (_, r) in zip(y, tab.iterrows()):
        ax.plot([r.r1_min, max(r.r1_max, r.r1_min + 0.006)], [yi] * 2, color=BLUE, lw=6, solid_capstyle='round')
        ax.plot(r.r1_best, yi, marker='|', color='white', ms=9, mew=2)
        ax.plot(r.r_age_Rbir, yi, marker='v', color=INK, ms=6, ls='none')
        ax.plot(r.r_FeH_Rbir, yi, marker='x', color=INK, ms=6, ls='none')
        ax.plot(r.r_JPhi_Rbir, yi, marker='s', mfc='#f59e0b', mec=INK, ms=6, ls='none')
        ax.plot(r.r_ceiling_Rbir, yi, marker='o', mfc='none', mec=INK, ms=6, ls='none')
    ax.axvline(0.64, color=MUTED, ls='--', lw=1)
    ax.text(0.645, y.max() + 0.55, 'paper 0.64', color=MUTED, fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(tab.galaxy, fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(-2.6, len(tab) - 0.2)
    ax.set_xlabel('|Pearson r| with the true birth radius')
    ax.plot([], [], color=BLUE, lw=6, label='L1, statistically indistinguishable fits (white tick: best fit)')
    ax.plot([], [], marker='v', color=INK, ls='none', label='age alone')
    ax.plot([], [], marker='x', color=INK, ls='none', label='[Fe/H] alone')
    ax.plot([], [], marker='s', mfc='#f59e0b', mec=INK, ls='none', label=r'$J_\phi$ alone (≈ present radius)')
    ax.plot([], [], marker='o', mfc='none', mec=INK, ls='none', label='linear regression on all 6 inputs')
    ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e5e7eb', lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f'R_0_30_Z_0.3 ({version}): how well can L1 track the birth radius?\n(author graph frozen)',
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(HERE.with_name(f'summary_ci_{version}.png'), dpi=130)
    plt.close(fig)


def selftest():
    """在旧的 7<R<10 数据上核对: 协方差法算出的 |r| 与逐星计算 (<tag>_ci.csv) 完全一致"""
    ci = pd.read_csv(ROOT / 'results' / 'rlcd_A' / 'identifiability' / 'nonoise_538905_ci.csv').sort_values('nll')
    e = pd.read_csv(ROOT / 'data' / 'rlcd_input' / 'nonoise' / '538905_e.csv')[COLS]
    t = pd.read_csv(ROOT / 'data' / 'rlcd_input' / 'nonoise' / '538905_truth.csv')
    X = np.column_stack([e.to_numpy(float), t[['R_bir', 'Rg']].to_numpy(float)])
    m = dict(N=len(X), cov=np.cov(X.T, ddof=1).tolist())
    C, c_R, c_G, var_R, var_G = standardized_moments(m)
    for k in range(3):
        row = ci.iloc[k]
        coef = {tuple(x.replace('X_', '') for x in col.split('->')): row[col] for col in ci.columns if '->' in col}
        w1, w2 = latent_weights(coef)
        print(f'解 {k}: 逐星 |r(L1,Rb)| {row.abs_r_L1_Rbir:.6f}  协方差法 {corr_from_weights(w1, C, c_R, var_R):.6f};  '
              f'逐星 |r(L2,Rg)| {row.abs_r_L2_Rg:.6f}  协方差法 {corr_from_weights(w2, C, c_G, var_G):.6f}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--extra', type=int, default=200)
    ap.add_argument('--summary-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    moments = json.load(open(HERE.with_name(f'moments_{args.version}.json')))
    if not args.summary_only:
        names, idx, edges, order = setup(str(DOT))
        for gid in args.ids or IDS:
            analyse(f'{args.version}_{gid}', moments[gid], names, idx, edges, order, args.extra)
    write_summary(args.version, moments)


if __name__ == '__main__':
    main()
