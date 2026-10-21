# -*- coding: utf-8 -*-
"""
planB_latent_vs_truth.py —— 做法 B: 每个 TNG50 星系"自己学出来"的因果图里, 隐变量与真出生半径、真导引半径的对比

运行 (numpy、pandas、matplotlib; 主 venv 即可), 先跑完 PlanB/output_data/rlcd_B.py 和 rlcd_B_eval.py:
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 "PlanB/真值对比/planB_latent_vs_truth.py"                  # 不加误差
    ~/PYTHONPATH_01/venv/bin/python3 "PlanB/真值对比/planB_latent_vs_truth.py" --version noise  # → planB_latent_vs_truth_noise.png
输出: 与本脚本同名、同目录的 planB_latent_vs_truth.png

做法 (每个星系用出现次数最多的结构; 该结构含无向边或没有隐变量时, 面板上注明并留空):
  1. 读 PlanB/input_data/<版本>/<ID>_e.csv 的 6 个输入, 每列减平均值、除以标准差, 得到 z_i;
  2. 每个隐变量 L = Σ_i w_i·z_i。w 来自 eval_<版本>_best.json, 是 rlcd_B_eval.py 用与作者 notebook 相同的分层最小二乘
     (从最下层的隐变量往上, 用它的子节点反推它) 在最大似然边系数下求出的权重;
  3. ⚠ 做法 B 学到的图里没有"哪个隐变量是出生半径"的标签 (编号是任意的)。这里直接挑与真出生半径 |r| 最大的隐变量
     当作"出生半径候选" L_Rb, 与真导引半径 |r| 最大的当作 L_Rg。这是"偷看了答案"的最好情况, 真实观测里做不到;
  4. 线性换算 (与作者 notebook 相同): f(L) = α·L + β, α、β 是让 Σ (真值 − α·L − β)² 最小的直线
     (α: 斜率, β: 截距); 横轴 f(L_Rb), 纵轴真出生半径; 黑线是 1:1;
  5. 面板里的参照: 单用年龄、[Fe/H]、J_φ 与真出生半径的 |r|, 6 个输入线性回归的上限 (任何"6 个输入的线性组合"
     与出生半径的 |r| 都不会超过它)。
  f(L) 的横向宽度 = |r| × 真值的宽度 (r 是 L 与真值的皮尔逊相关系数), 所以 |r| 小时点挤成窄竖条。
"""
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    v = args.version
    st = pd.read_csv(PLANB / 'output_data' / f'structures_{v}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    best = json.load(open(PLANB / 'output_data' / f'eval_{v}_best.json'))
    fig, axs = plt.subplots(6, 4, figsize=(17, 25))
    lines = []
    for k, gid in enumerate(IDS):
        axR, axG = axs.flat[2 * k], axs.flat[2 * k + 1]
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        sid, cnt = freq.index[0], freq.iloc[0]
        info = best.get(gid, {}).get(str(sid))
        e = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir', 'Rg'])
        z = ((e - e.mean()) / e.std()).to_numpy()
        rb = lambda c: abs(np.corrcoef(e[c], t.R_bir)[0, 1])
        base = f'age only {rb("age"):.2f}   [Fe/H] only {rb("FeH"):.2f}   $J_\\phi$ only {rb("JPhi"):.2f}'
        beta, *_ = np.linalg.lstsq(np.column_stack([z, np.ones(len(z))]), t.R_bir.to_numpy(), rcond=None)
        ceil = np.corrcoef(np.column_stack([z, np.ones(len(z))]) @ beta, t.R_bir)[0, 1]
        head = f'{gid} ({len(e):,} stars) — structure {sid} ({cnt}/{len(g)} seeds)'
        if info is None or not info['weights']:
            why = 'no latent variable' if info is not None else 'structure has undirected edges / not evaluated'
            for ax, nm in [(axR, 'birth radius'), (axG, 'guiding radius')]:
                ax.text(0.5, 0.5, why, ha='center', va='center', transform=ax.transAxes, fontsize=10, color='#6b7280')
                ax.set_title(f'{head}\n{nm}', fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            axR.text(0.04, 0.95, base + f'\nlinear ceiling {ceil:.2f}', transform=axR.transAxes, va='top', fontsize=7.5)
            lines.append(f'{gid}: 结构 {sid} — {why}')
            continue
        Lval = {L: z @ np.array(w) for L, w in info['weights'].items()}
        for ax, col, nm, lim in [(axR, 'R_bir', 'birth radius', 20), (axG, 'Rg', 'guiding radius', 20)]:
            y = t[col].to_numpy()
            rs = {L: abs(np.corrcoef(x, y)[0, 1]) for L, x in Lval.items()}
            L = max(rs, key=rs.get)
            a, b = np.polyfit(Lval[L], y, 1)
            f = a * Lval[L] + b
            r = np.corrcoef(f, y)[0, 1]
            ax.hexbin(f, y, bins='log', extent=(0, lim, 0, lim), gridsize=60, cmap='Blues', mincnt=1)
            ax.plot([0, lim], [0, lim], color='#1f2937', lw=1)
            ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal')
            ax.set_xlabel(f'$f({L})$ [kpc]', fontsize=9)
            ax.set_ylabel(('$R_b$' if col == 'R_bir' else '$R_g$') + ' (true) [kpc]', fontsize=9)
            txt = f'$r$ = {r:.2f}   ($f$ = {a:+.2f} {L} {b:+.2f})\nall latents |r|: ' + ', '.join(f'{k2} {v2:.2f}' for k2, v2 in sorted(rs.items()))
            if col == 'R_bir':
                txt += f'\n{base}\nlinear ceiling {ceil:.2f}'
            ax.text(0.04, 0.96, txt, transform=ax.transAxes, va='top', fontsize=7.3)
            ax.set_title(f'{head}\n{nm}: best-matching latent = {L}', fontsize=9)
            lines.append(f'{gid} {col}: 最匹配的隐变量 {L}, f = {a:+.3f}*{L} {b:+.3f}, r = {r:.3f}')
    fig.suptitle(f'Plan B — TNG50, 7<R<10 kpc, |z|<0.3 kpc ({v}): latents of the graph RLCD learned in each galaxy (most frequent of 30 seeds) vs truth\n'
                 'Latent names are arbitrary, so the latent with the highest |r| is chosen using the TRUE values (best case, not possible for real data).\n'
                 'f(L) = linear least-squares map of L onto the truth; spread of f = |r| × spread of the truth.', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = HERE.with_suffix('.png') if v == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')
    fig.savefig(out, dpi=100)
    print('\n'.join(lines))
    print(f'已保存 {out}')


if __name__ == '__main__':
    main()
