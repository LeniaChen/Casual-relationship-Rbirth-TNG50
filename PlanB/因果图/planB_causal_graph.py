# -*- coding: utf-8 -*-
"""
planB_causal_graph.py —— 做法 B: 每个 TNG50 星系"自己学出来"的因果图 (出现次数最多的那一种) 及其最大似然边系数

运行 (numpy、pandas、matplotlib; 主 venv 即可), 先跑完 PlanB/output_data/rlcd_B.py 和 rlcd_B_eval.py:
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 "PlanB/因果图/planB_causal_graph.py"                  # 不加误差
    ~/PYTHONPATH_01/venv/bin/python3 "PlanB/因果图/planB_causal_graph.py" --version noise  # 加误差 → planB_causal_graph_noise.png
    可选: --rank 2 画每个星系出现次数第 2 多的结构 (→ planB_causal_graph_rank2.png)
输出: 与本脚本同名、同目录的 planB_causal_graph.png

读的文件 (PlanB/output_data/):
  structures_<版本>.csv/.txt   每个星系 30 个随机种子 (PYTHONHASHSEED 0–29) 学到的图, 结构编号全体星系统一
  eval_<版本>_best.json        每个 (星系, 结构) 的最大似然解: 边系数、每个观测量的噪声方差
  eval_<版本>.csv              每个隐变量在最大似然解下与真出生半径、真导引半径的 |r|

图上的内容:
  - 红色椭圆 = 隐变量 (RLCD 自己发现的、没被观测到的量); 蓝色 = 6 个观测量;
  - 红色虚线 = 无向边 (RLCD 没能定下方向); 含无向边的图无法求边系数, 只画结构;
  - 隐变量的编号 L1, L2, ... 是任意的 (只是为了统一不同种子的结果, 按字典序最小的方式编号), 不代表"出生半径""导引半径";
    每个隐变量椭圆里写着它在最大似然解下与真值的相关系数: Rb = |r(L, 真出生半径)|, Rg = |r(L, 真导引半径)|;
  - 箭头上的数字 = 边系数 (标准化后, 即"父节点变化 1 个标准差, 子节点平均变化多少个标准差");
  - 面板标题: 星系、结构编号、在 30 个种子里出现几次、边数、雅可比秩 (= 边数时边系数才能由数据唯一确定)、
    "prior ✗" 表示违反论文 3.2 节的先验 (年龄有观测父节点);
  - 布局: 按因果顺序分层 (没有父节点的在最上面), 与作者图的位置无关。
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
OUT_DATA = PLANB / 'output_data'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LABEL = {'logJZ': r'$J_z$', 'JPhi': r'$J_\phi$', 'ecc': r'$e$', 'age': 'Age', 'FeH': '[Fe/H]', 'MgH': '[Mg/H]'}
EW, EH = 0.085, 0.052
LATENT_FC, LATENT_EC = '#fde2e1', '#c2362f'
OBS_FC, OBS_EC = '#dbe8fb', '#2563eb'
INK, MUTED = '#1f2937', '#6b7280'


def layout(nodes, edges):
    """分层布局: 层号 = 从没有父节点的节点出发的最长路径长度; 同层内按父节点的平均横坐标排序 (减少交叉)"""
    par = {v: [a for a, b in edges if b == v] for v in nodes}
    depth = {}

    def d(v):
        if v not in depth:
            depth[v] = 0 if not par[v] else 1 + max(d(p) for p in par[v])
        return depth[v]
    for v in nodes:
        d(v)
    nl = max(depth.values()) + 1
    layers = [[v for v in nodes if depth[v] == k] for k in range(nl)]
    pos = {}
    for k, lay in enumerate(layers):
        if k > 0:
            lay.sort(key=lambda v: np.mean([pos[p][0] for p in par[v]]))
        y = 0.9 - 0.8 * k / max(nl - 1, 1)
        for i, v in enumerate(lay):
            pos[v] = ((i + 1) / (len(lay) + 1), y)
    return pos


def boundary(p, q):
    """从椭圆 p 的中心朝 q 走, 与椭圆边界的交点"""
    dx, dy = q[0] - p[0], q[1] - p[1]
    t = 1 / np.sqrt((dx / EW) ** 2 + (dy / EH) ** 2)
    return p[0] + t * dx, p[1] + t * dy


def draw(ax, edges, coef, undirected, omega, ev, title):
    """edges: 有向边; coef: {边: 系数} (没有系数时为空); undirected: 无向边 (虚线, 方向没定下来); ev: 每个隐变量的 |r|"""
    lat = sorted({v for e in list(edges) + list(undirected) for v in e if v not in COLS})
    nodes = lat + COLS
    pos = layout(nodes, edges)
    vmax = max([abs(c) for c in coef.values()] + [1e-9]) if coef else 1
    nout = {}
    for (a, b) in edges:
        c = coef.get((a, b))
        j = nout.get(a, 0)
        nout[a] = j + 1
        pa, pb = pos[a], pos[b]
        same_layer = abs(pa[1] - pb[1]) < 1e-6
        rad = 0.35 if same_layer else 0.08
        s, e = boundary(pa, pb), boundary(pb, pa)
        lw = 0.6 + 2.4 * abs(c) / vmax if c is not None else 1.0
        ax.annotate('', xy=e, xytext=s, arrowprops=dict(arrowstyle='-|>', color=INK, lw=lw, shrinkA=0, shrinkB=0,
                    mutation_scale=9, connectionstyle=f'arc3,rad={rad}'))
        if c is None:
            continue
        # 标签位置: arc3 的控制点 = 中点 + rad·(dy, −dx), 取二次贝塞尔曲线上参数 t 处的点;
        # 同一个父节点的多条边, t 在 0.3 / 0.47 / 0.64 之间轮换, 避免标签挤在一起
        dx, dy = e[0] - s[0], e[1] - s[1]
        cpt = ((s[0] + e[0]) / 2 + rad * dy, (s[1] + e[1]) / 2 - rad * dx)
        tt = 0.3 + 0.17 * (j % 3)
        lx = (1 - tt) ** 2 * s[0] + 2 * tt * (1 - tt) * cpt[0] + tt ** 2 * e[0]
        ly = (1 - tt) ** 2 * s[1] + 2 * tt * (1 - tt) * cpt[1] + tt ** 2 * e[1]
        ax.text(lx, ly, f'{c:.2f}', fontsize=6.3, ha='center', va='center', color=INK, zorder=5,
                bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85))
    for (a, b) in undirected:
        s, e = boundary(pos[a], pos[b]), boundary(pos[b], pos[a])
        ax.annotate('', xy=e, xytext=s, arrowprops=dict(arrowstyle='-', color=LATENT_EC, lw=1.3, ls='--',
                    shrinkA=0, shrinkB=0, connectionstyle='arc3,rad=0.25'))
    for v in nodes:
        x, y = pos[v]
        is_lat = v not in COLS
        h = EH * (1.45 if is_lat and v in ev.index else 1)
        ax.add_patch(matplotlib.patches.Ellipse((x, y), 2 * EW * (1.12 if is_lat else 1), 2 * h,
                                                fc=LATENT_FC if is_lat else OBS_FC,
                                                ec=LATENT_EC if is_lat else OBS_EC, lw=1.1, zorder=6))
        if is_lat and v in ev.index:
            ax.text(x, y + 0.022, v, ha='center', va='center', fontsize=8.5, zorder=7, color=INK, fontweight='bold')
            ax.text(x, y - 0.02, f'Rb {ev.loc[v, "r_Rbir_best"]:.2f}\nRg {ev.loc[v, "r_Rg_best"]:.2f}',
                    ha='center', va='center', fontsize=6.3, zorder=7, color=LATENT_EC, linespacing=1.0)
        else:
            ax.text(x, y, v if is_lat else LABEL[v], ha='center', va='center', fontsize=8, zorder=7, color=INK)
        if not is_lat and omega.get('X_' + v, omega.get(v, 1)) < 1e-3:
            ax.text(x, y - EH - 0.012, 'ω≈0', ha='center', va='top', fontsize=6.5, color=LATENT_EC, fontweight='bold', zorder=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.03, 1.0)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=8.2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--rank', type=int, default=1, help='1 = 每个星系最常见的结构, 2 = 第二常见, ...')
    args = ap.parse_args()
    st = pd.read_csv(OUT_DATA / f'structures_{args.version}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    best = json.load(open(OUT_DATA / f'eval_{args.version}_best.json'))
    ev = pd.read_csv(OUT_DATA / f'eval_{args.version}.csv')
    ev['galaxy'] = ev.galaxy.astype(str)
    fig, axs = plt.subplots(3, 4, figsize=(20, 16))
    for ax, gid in zip(axs.flat, IDS):
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        if len(freq) < args.rank:
            ax.axis('off')
            ax.set_title(f'{gid}: only {len(freq)} structure(s) found', fontsize=8.2)
            continue
        sid, cnt = freq.index[args.rank - 1], freq.iloc[args.rank - 1]
        row = g[g.structure_id == sid].iloc[0]
        info = best.get(gid, {}).get(str(sid))
        head = f'{gid}   structure {sid}: {cnt}/{len(g)} seeds, {len(freq)} distinct structures'
        prior = '' if row.prior_ok else '   prior ✗'
        edges = [tuple(x) for x in json.loads(row.structure)]
        und = [tuple(x) for x in json.loads(row.undirected)]
        if info is None:          # 含无向边 (方向没定下来) 或没评估: 只画结构, 没有系数
            why = (f'{len(und)} undirected edge(s) (dashed red): direction not determined → no coefficients'
                   if row.n_undirected > 0 else 'not evaluated (freq < 20%)')
            draw(ax, edges, {}, und, {}, pd.DataFrame(), head + f'\n{why}{prior}')
            continue
        e = ev[(ev.galaxy == gid) & (ev.structure_id == sid)].set_index('latent')
        title = (head + f'\n{info["n_edges"]} edges, Jacobian rank {info["jacobian_rank"]}/{info["n_edges"]}'
                 f'{" (identifiable)" if info["jacobian_rank"] == info["n_edges"] else " (NOT identifiable)"}{prior}')
        coef = {tuple(k): c for k, c in zip(info['edges'], info['coef'])}
        draw(ax, edges, coef, [], info['omega'], e, title)
    rk = '' if args.rank == 1 else f' (rank-{args.rank} most frequent structure)'
    fig.suptitle(f'Plan B — TNG50, 7<R<10 kpc, |z|<0.3 kpc ({args.version}): causal graph learned from scratch by RLCD in each galaxy'
                 f'{rk}\nmost frequent structure over 30 hash seeds; edge numbers = maximum-likelihood coefficients; '
                 'inside each latent: |r| with the TRUE birth radius (Rb) and guiding radius (Rg).  Latent names are arbitrary.',
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    stem = HERE.stem + ('' if args.version == 'nonoise' else '_noise') + ('' if args.rank == 1 else f'_rank{args.rank}')
    out = HERE.with_name(stem + '.png')
    fig.savefig(out, dpi=110)
    print(f'已保存 {out}')


if __name__ == '__main__':
    main()
