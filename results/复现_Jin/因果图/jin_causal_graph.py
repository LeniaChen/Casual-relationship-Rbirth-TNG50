# -*- coding: utf-8 -*-
"""
jin_causal_graph.py —— 复现 Jin et al. 论文图 1: g2.79e12 上 RLCD 学到的因果图 (带边系数)

运行 (必须用 RLCD 专用环境 venv_rlcd, 里面有 causal-learn、statsmodels 等):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 "results/复现_Jin/因果图/jin_causal_graph.py"
输出: 与本脚本同名、同目录的 jin_causal_graph.png (约 1 分钟)

做了什么:
  1. 结构学习 (真正重跑): 用仓库的 RLCD, 数据 data/g2.79e12_e.csv (6 列 logJZ, JPhi, ecc, age, FeH, MgH),
     参数与作者结果文件名相同: alpha0..3 = 1e-5, 1e-15, 1e-5, 1e-15, 第一阶段 GES (sparsity 2),
     与 rlcd/run_rlcd.py 的默认设置完全一样。
     ⚠ RLCD 的结果依赖 Python 的字符串哈希随机化 (问题 1): 只有 PYTHONHASHSEED=1 时才得到作者的图。
       本脚本启动时会自动以 PYTHONHASHSEED=1 重新启动自己; 想看别的种子: --seed 3。
       固定了种子, 结果就是确定的 (同一个种子每次都一样); "随机"指的是换种子会变。
       2026-09-27 扫了 PYTHONHASHSEED = 0–29 共 30 个种子: 作者的图出现 7 次 (1, 7, 8, 11, 14, 22, 25);
       最常见的是另一种 15 条边的图 (12 次, 如 seed 0、2、3); 其余 11 次是 18 条边的图 (4 种)。
     结构学习在内存里完成, 不写入仓库的 results_* 目录 (不会覆盖作者的文件)。
  2. 与作者提交的图 (results/rlcd_g279/author_result_..._ranktestscale1.dot) 逐条边比较
     (隐变量编号按连接方式自动对齐), 并检查论文 3.2 节的先验 (年龄没有观测父节点)。
  3. 边上的数字: 作者发布的参数学习结果 (…_param_by_trek_nll1.5861642734160428.dot, 即论文图 1 的数字)。
     ⚠ 这些数字不是本脚本算出来的, 是从作者发布的文件里直接读的, 所以当然和论文一模一样。
       参数学习用的是没设种子的随机初值 + SGD, 作者当时那次运行无法逐位重现 (作者仓库里同一张图的 13 次参数学习
       结果各不相同); 而且这组系数不是最大似然解 (问题 2)。自己重新求系数的结果见 results/rlcd_A/identifiability/g279_ci.*。
  4. 画图: 节点与论文图 1 相同 (布局略作调整以减少交叉); 红色 = 隐变量, 蓝色 = 观测量; 线宽 ∝ |系数|;
     图下方写明结构是否与作者的完全相同。
"""
import os
import re
import sys
import argparse
import itertools
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # casuality/
REPO = ROOT / 'Causal-Galactic-Archaeology'
DATA = REPO / 'data' / 'g2.79e12_e.csv'
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
COEF_DOT = REPO / 'rlcd' / 'results_causal_vertical_heating_g279_e' / \
    'result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1_param_by_trek_nll1.5861642734160428.dot'
OUT = HERE.with_suffix('.png')
ALPHA = (1e-5, 1e-15, 1e-5, 1e-15)

# ---------------------------------------------------------------- 画图 (与 TNG/因果图 脚本同一套样式)
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


def read_dot_edges(path, with_coef=False):
    """读 .dot 的有向边 → {(父, 子): 系数或 None}, 变量名去掉 X_"""
    txt = open(path).read()
    out = {}
    for m in re.finditer(r'^\s*(\w+)\s*->\s*(\w+)\s*(\[[^\]]*\])?', txt, re.M):
        a, b, attr = m.group(1).replace('X_', ''), m.group(2).replace('X_', ''), m.group(3) or ''
        lab = re.search(r'label="?([-+\deE.]+)', attr)
        out[(a, b)] = float(lab.group(1)) if (with_coef and lab) else None
    return out


def align_latents(ours, ref):
    """隐变量编号是任意的: 试遍所有换名方式, 取与参考图共同边最多的一种"""
    lat = sorted({v for e in ours for v in e if v.startswith('L')})
    best = None
    for perm in itertools.permutations(lat):
        mp = dict(zip(lat, perm))
        ren = {(mp.get(a, a), mp.get(b, b)) for a, b in ours}
        score = len(ren & set(ref))
        if best is None or score > best[0]:
            best = (score, ren, mp)
    return best[1], best[2]


# ---------------------------------------------------------------- 结构学习
def run_rlcd(seed):
    """在本进程里跑 RLCD (与 run_rlcd.py 相同的设置), 返回 (有向边集合, 无向边集合)"""
    os.chdir(REPO / 'rlcd')                    # 仓库代码按相对路径找文件
    sys.path.insert(0, str(REPO / 'rlcd'))
    import pandas as pd
    from DGM.DataModel import DataModel
    from utils.Chi2RankTest import Chi2RankTest
    from StructureLearning.RLCD.RLCD_alg import RLCD
    import logging
    logging.getLogger().setLevel(logging.WARNING)   # RLCD 的日志很多, 只保留警告

    df = pd.read_csv(DATA).dropna()                 # 与 scenarios.causal_vertical_heating_g279_e() 相同
    df = df.rename(columns={k: 'X_' + k.replace('-', '_') for k in df.columns})
    dgm = DataModel(df, df)
    df_x, _ = dgm.generate_data(N=-1, normalized=True)
    params = {'ranktest_method': Chi2RankTest(df_x.to_numpy(), 1), 'citest_method': None,
              'stage1_method': 'ges', 'alpha_dict': dict(enumerate(ALPHA)),
              'stage1_ges_sparsity': 2, 'stage1_partition_thres': 3}
    graph, _, _, _ = RLCD(1, dgm.xvars, df_x, params)
    strip = lambda v: v.replace('X_', '')
    directed = {(strip(a), strip(b)) for a, b in graph.dirEdges}
    undirected = {tuple(sorted(map(strip, e))) for e in graph.undirEdges}
    return directed, undirected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', default='1', help='PYTHONHASHSEED (作者的图对应 1)')
    args = ap.parse_args()
    if os.environ.get('PYTHONHASHSEED') != args.seed:          # 哈希种子只能在 Python 启动前设定
        env = dict(os.environ, PYTHONHASHSEED=args.seed)
        os.execve(sys.executable, [sys.executable, str(HERE)] + sys.argv[1:], env)

    ours, undirected = run_rlcd(args.seed)
    author = read_dot_edges(AUTHOR_DOT)
    ours_aligned, mapping = align_latents(ours, author)
    same = ours_aligned == set(author) and not undirected
    extra, missing = sorted(ours_aligned - set(author)), sorted(set(author) - ours_aligned)
    obs = {'JPhi', 'ecc', 'age', 'FeH', 'MgH', 'logJZ'}
    prior_ok = not any(c == 'age' and p in obs for p, c in ours_aligned)

    coef = read_dot_edges(COEF_DOT, with_coef=True)
    assert set(coef) == set(author), '系数文件的边与作者的图不一致'

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 7.6))
    if same:
        draw_graph(ax, coef, fontsize=10)
        note = (f'Structure re-learned with RLCD (PYTHONHASHSEED={args.seed}, alpha={ALPHA}): '
                f'IDENTICAL to the author\'s graph ({len(ours)} edges);\nage-has-no-observed-parent prior: '
                f'{"passed" if prior_ok else "FAILED"}.\n'
                'Edge numbers: READ from the author\'s released parameter-learning file (NLL 1.586), not re-computed '
                '(the paper\'s Fig. 1 numbers).')
    else:
        for L in ('L1', 'L2'):                     # 结构不同时, 隐变量不再有"出生半径/导引半径"的含义
            x, y, _, lat = NODES[L]
            NODES[L] = (x, y, L, lat)
        draw_graph(ax, {e: None for e in ours_aligned if e[0] in NODES and e[1] in NODES}, fontsize=10)
        note = (f'Structure re-learned with RLCD (PYTHONHASHSEED={args.seed}): DIFFERENT from the author\'s graph.\n'
                f'extra edges: {extra}   missing edges: {missing}   undirected: {sorted(undirected)}   '
                f'prior {"passed" if prior_ok else "FAILED"}')
    fig.text(0.5, 0.015, note, ha='center', va='bottom', fontsize=7.5, color=MUTED, wrap=True)
    fig.suptitle('Reproduction of Jin et al. Fig. 1 — causal graph of g2.79e12', fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    fig.savefig(OUT, dpi=150)

    print(f'PYTHONHASHSEED = {args.seed}; 隐变量对齐 {mapping}')
    print(f'我们的图 {len(ours)} 条有向边, 无向边 {len(undirected)} 条; 作者的图 {len(author)} 条')
    print('与作者的图完全相同' if same else f'与作者的图不同: 多出 {extra}, 缺少 {missing}')
    print(f'先验 (年龄没有观测父节点): {"通过" if prior_ok else "不通过"}')
    print(f'已保存 {OUT}')


if __name__ == '__main__':
    main()
