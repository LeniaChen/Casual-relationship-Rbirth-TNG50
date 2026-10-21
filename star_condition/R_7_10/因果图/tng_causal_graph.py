# -*- coding: utf-8 -*-
"""
tng_causal_graph.py —— TNG50 12 个星系的因果图 (做法 A: 冻结作者的图, 在每个星系上重新估计边系数; 仿照论文图 5)

运行 (numpy、pandas、matplotlib; 主 venv 即可):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 "results/TNG/因果图/tng_causal_graph.py"                  # 不加误差 (主版本)
    ~/PYTHONPATH_01/venv/bin/python3 "results/TNG/因果图/tng_causal_graph.py" --version noise  # 加误差 → tng_causal_graph_noise.png
输出: 与本脚本同名、同目录的 tng_causal_graph.png

用到的数据 (都已在 casuality 里):
  - 因果图的结构: 作者的最终图 results/rlcd_g279/author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot
    (与论文图 1 相同, 15 条边; 论文 4.2、4.3 节就是把这张图冻结后迁移到别的星系)
  - 边系数: results/rlcd_A/identifiability/<version>_<ID>_ci.csv 里 NLL 最小的那一行, 即
    scripts/local/rlcd_identifiability_ci.py 求出的最大似然解 (230 个初值、严格收敛、噪声方差 >= 0, 系数限制在 ±0.999)。
    这是"论文说的最大似然方法做到收敛"的结果。

⚠ 必须知道 (REPORT_RLCD_TNG50.md 第 4.3 节):
  - 这张图冻结后边系数不可识别: 统计上与最优分不开的解还有很多 (每个面板标题里的 n_eq), 它们的系数可以很不一样,
    这里画的只是其中 NLL 最小的一个;
  - 面板里标 "ω≈0" 的观测量: 最优解把它的噪声方差压到了 0 (模型认为它完全由父节点决定, 物理上不成立);
  - 系数等于 ±0.999 表示碰到了系数上限。
"""
import re
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # casuality/
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDENT = ROOT / 'results' / 'rlcd_A' / 'identifiability'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
LR95 = 3.84

# ---------------------------------------------------------------- 画图 (与 复现_Jin/因果图 脚本同一套样式)
NODES = {  # 名字: (x, y, 图上的文字, 是否隐变量)
    'L1': (0.55, 0.93, 'L1\nBirth Radius', True),
    'L2': (0.08, 0.62, 'L2\nGuiding Radius', True),
    'MgH': (0.36, 0.60, '[Mg/H]', False),
    'age': (0.74, 0.60, 'Stellar Age', False),
    'JPhi': (0.06, 0.12, r'$J_\phi$' + '\nAzimuthal Action', False),
    'ecc': (0.36, 0.12, r'$e$' + '\nEccentricity', False),
    'FeH': (0.70, 0.12, '[Fe/H]', False),
    'logJZ': (1.10, 0.60, r'$J_z$' + '\nVertical Action', False),
}
EW, EH = 0.105, 0.065          # 椭圆半宽、半高
BEND = {('MgH', 'logJZ'): -0.22, ('L1', 'logJZ'): 0.04}   # 弯曲量: 控制点沿边的左法向偏移 (负值向右/下弯)
LABEL_T = {('age', 'JPhi'): 0.22, ('MgH', 'FeH'): 0.55, ('age', 'ecc'): 0.30, ('MgH', 'ecc'): 0.55,
           ('L2', 'ecc'): 0.35, ('age', 'FeH'): 0.55, ('MgH', 'JPhi'): 0.62, ('MgH', 'logJZ'): 0.80,
           ('L1', 'FeH'): 0.30, ('L2', 'JPhi'): 0.5}
LATENT_FC, LATENT_EC = '#fde2e1', '#c2362f'
OBS_FC, OBS_EC = '#dbe8fb', '#2563eb'
INK, MUTED = '#1f2937', '#6b7280'


def _edge_point(center, toward):
    import numpy as np
    d = np.array(toward, float) - np.array(center, float)
    s = 1.0 / np.sqrt((d[0] / EW) ** 2 + (d[1] / EH) ** 2)
    return np.array(center, float) + d * s


def draw_graph(ax, coef, title=None, fontsize=8.5):
    """coef: {(父, 子): 系数 或 None}; 变量名不带 X_ 前缀"""
    import numpy as np
    from matplotlib.patches import Ellipse, FancyArrowPatch
    from matplotlib.path import Path as MPath
    ax.set_xlim(-0.07, 1.24)
    ax.set_ylim(0.03, 1.02)
    ax.set_aspect('equal')
    ax.axis('off')
    for name, (x, y, lab, lat) in NODES.items():
        ax.add_patch(Ellipse((x, y), 2 * EW, 2 * EH, fc=LATENT_FC if lat else OBS_FC,
                             ec=LATENT_EC if lat else OBS_EC, lw=1.4, zorder=3))
        ax.text(x, y, lab, ha='center', va='center', fontsize=fontsize - 0.5, color=INK, zorder=4)
    vals = [abs(v) for v in coef.values() if v is not None]
    vmax = max(vals) if vals else 1
    for (p, c), v in coef.items():
        P, Q = np.array(NODES[p][:2]), np.array(NODES[c][:2])
        d = Q - P
        C = (P + Q) / 2 + BEND.get((p, c), 0.0) * np.array([-d[1], d[0]]) / np.hypot(*d)
        P1, Q1 = _edge_point(P, C), _edge_point(Q, C)
        path = MPath([P1, C, Q1], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
        lw = 0.7 + 3.3 * (abs(v) / vmax if v is not None else 0.3)
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=10 + 2 * lw, lw=lw,
                                     color=INK, alpha=0.8, zorder=2))
        if v is not None:
            t = LABEL_T.get((p, c), 0.5)
            B = (1 - t) ** 2 * P1 + 2 * (1 - t) * t * C + t ** 2 * Q1
            ax.text(B[0], B[1], f'{v:.3f}' if abs(v) >= 0.995 else f'{v:.2g}', ha='center', va='center', fontsize=fontsize - 1, color=INK, zorder=5,
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.92))
    if title:
        ax.set_title(title, fontsize=fontsize + 1.5, color=INK)


def author_edges():
    txt = open(AUTHOR_DOT).read()
    return [(a.replace('X_', ''), b.replace('X_', '')) for a, b in re.findall(r'^\s*(\w+)\s*->\s*(\w+)', txt, re.M)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    out = HERE.with_suffix('.png') if args.version == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')
    edges = author_edges()
    col = {e: ('X_' if e[0] not in ('L1', 'L2') else '') + e[0] + '->' + ('X_' if e[1] not in ('L1', 'L2') else '') + e[1]
           for e in edges}

    fig, axs = plt.subplots(3, 4, figsize=(22, 15.5))
    for ax, gid in zip(axs.flat, IDS):
        ci = pd.read_csv(IDENT / f'{args.version}_{gid}_ci.csv').sort_values('nll')
        best = ci.iloc[0]
        coef = {e: float(best[col[e]]) for e in edges}
        n_eq = int((ci.lr < LR95).sum())
        draw_graph(ax, coef, fontsize=8.5)
        node = str(best.omega_min_node).replace('X_', '')
        if best.omega_min < 0.005 and node in NODES:
            x, y = NODES[node][:2]
            ax.text(x, y - EH - 0.035, 'ω≈0', ha='center', va='center', fontsize=8, color='#c2362f', fontweight='bold')
        ax.set_title(f'{gid}   (NLL {best.nll:.3f};  statistically indistinguishable fits n_eq = {n_eq})',
                     fontsize=10.5, color=INK)
    fig.suptitle(f'TNG50 ({args.version}): author\'s causal graph frozen, edge coefficients = maximum-likelihood fit per galaxy\n'
                 'ω≈0: noise variance of that observable driven to 0 at the optimum;  coefficients are NOT unique '
                 '(see REPORT_RLCD_TNG50.md §4.3)', fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=110)
    print(f'已保存 {out}')


if __name__ == '__main__':
    main()
