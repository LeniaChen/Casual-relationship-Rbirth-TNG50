# -*- coding: utf-8 -*-
"""
rlcd_B_eval.py —— 做法 B 第 2 步: 对每个星系学到的因果图, 求边系数并检查每个隐变量与真出生半径、真导引半径的关系

先运行同一文件夹的 rlcd_B.py 得到 PlanB/output_data/structures_<版本>.csv, 再运行 (numpy、scipy、pandas; 主 venv 即可):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 PlanB/output_data/rlcd_B_eval.py --version nonoise
    ~/PYTHONPATH_01/venv/bin/python3 PlanB/output_data/rlcd_B_eval.py --version noise
    (--min-frac 0.2: 只评估在该星系出现频率 >= 20% 的结构, 且至少评估最常见的一种; --nstart 50: 随机初值个数)

对每个 (星系, 结构):
  1. 雅可比矩阵的秩: 边数 vs 独立的观测相关系数个数 (15), 判断边系数能否被数据唯一确定 (方法同 rlcd_identifiability.py);
  2. 最大似然: nstart 个随机初值, 噪声方差 >= 0, 严格收敛 (同 rlcd_identifiability_ci.py);
     统计上与最优分不开 (2·N·ΔNLL < 3.84) 的解都保留;
  3. 每个隐变量的值: 与作者 notebook 相同的分层最小二乘 (从最下层的隐变量往上), 它是 6 个标准化输入的线性组合,
     所以它与真值的相关系数可以直接由协方差矩阵算出 (见 star_condition/R_0_30_Z_30/output_data/identifiability_moments.py);
  4. 对每个隐变量记录 |r(L, 真出生半径)| 和 |r(L, 真导引半径)| 在"统计上分不开的解"里的范围;
     参照线: 单用年龄、[Fe/H]、J_phi 与出生半径的 |r|, 6 个输入量线性回归的上限。
输出 (PlanB/output_data/, 其他 --cond 在其子文件夹):
  eval_<版本>.csv   每个星系、结构、隐变量一行 (|r| 的范围、雅可比秩、参照线)
  eval_<版本>.txt   同上的文字版
  eval_<版本>_best.json  每个星系、结构的最大似然解: 边系数、每个隐变量对 6 个标准化输入的权重 w
                    (隐变量 L = Σ_i w_i·z_i, z_i 是第 i 个输入减平均值再除以标准差); 画图脚本用它算每颗星的隐变量
"""
import sys
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
sys.path.insert(0, str(ROOT / 'scripts' / 'local'))      # rlcd_identifiability.py (implied, nll)
sys.path.insert(0, str(HERE.parent))                     # 同一文件夹的 rlcd_B.py (放在最前, 以免用到 scripts/local 里的旧版)
from rlcd_identifiability import implied, nll  # noqa: E402
from rlcd_B import cond_paths, COLS, OBS, IDS  # noqa: E402

LR95 = 3.84
PEN = 1e7
TIGHT = {'maxiter': 20000, 'ftol': 1e-15, 'gtol': 1e-11, 'maxcor': 30}


def moments(cond, version, gid):
    """8×8 协方差 (6 个输入 + R_bir + Rg) 和星数"""
    src, _ = cond_paths(cond)
    if cond == 'main':
        e = pd.read_csv(PLANB / 'input_data' / version / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(PLANB / 'true_value' / version / f'{gid}_truth.csv', usecols=['R_bir', 'Rg'])
        X = np.column_stack([e.to_numpy(float), t[['R_bir', 'Rg']].to_numpy(float)])
        return np.cov(X.T, ddof=1), len(X)
    m = json.load(open(src / 'output_data' / f'moments_{version}.json'))[gid]
    return np.array(m['cov']), m['N']


def graph_setup(edges):
    names = COLS + sorted({v for e in edges for v in e if v not in OBS})
    idx = {v: i for i, v in enumerate(names)}
    E = [(idx[a], idx[b]) for a, b in edges]
    order, rem = [], list(range(len(names)))
    while rem:
        for k in rem:
            if all((p, k) not in E or p in order for p in range(len(names))):
                order.append(k)
                rem.remove(k)
                break
        else:
            raise ValueError('图里有环')
    return names, idx, E, order


def latent_weights(names, edges_named, coef):
    """所有隐变量对 6 个标准化输入的权重。
    与作者 notebook 相同的分层最小二乘: 从最下层开始, 每次取"子节点里的隐变量都已求出"的那一组隐变量 S,
    对 S 的所有子节点 c 先扣掉 S 以外的已知父节点的贡献 r_c = c − Σ a_qc·q, 再求使 Σ_c (r_c − Σ_{L∈S} a_Lc·L)² 最小的 S。
    S 只有一个隐变量时就是作者的公式 L = Σ a_Lc·r_c / Σ a_Lc²; 同一层有两个隐变量共享子节点时 (作者的图里没有这种情况)
    一起求解。更上层的 (还没求出的) 隐变量对子节点的贡献按 0 处理, 与作者逐层求解的做法一致。"""
    w = {c: np.eye(6)[i] for i, c in enumerate(COLS)}
    lat = [v for v in names if v not in OBS]
    todo = set(lat)
    while todo:
        S = [L for L in sorted(todo) if not any(c in todo for (p, c) in edges_named if p == L)]
        if not S:
            raise ValueError('隐变量之间有环')
        ch = sorted({c for (p, c) in edges_named if p in S})
        A = np.array([[coef.get((L, c), 0.0) for L in S] for c in ch])
        R = []
        for c in ch:
            r = w[c].copy()
            for (q, c2) in edges_named:
                if c2 == c and q not in S and q in w:
                    r = r - coef[(q, c2)] * w[q]
            R.append(r)
        Wm = np.linalg.pinv(A) @ np.array(R)            # |S| × 6
        for L, row in zip(S, Wm):
            w[L] = row
        todo -= set(S)
    return {L: w[L] for L in lat}


def evaluate(edges_named, cov8, N, nstart):
    names, idx, E, order = graph_setup(edges_named)
    n, ne = len(names), len(E)
    xi = [idx[c] for c in COLS]
    iu = np.triu_indices(6, 1)
    sd = np.sqrt(np.diag(cov8))
    C = cov8[:6, :6] / np.outer(sd[:6], sd[:6])
    cR, cG = cov8[:6, 6] / sd[:6], cov8[:6, 7] / sd[:6]
    vR, vG = cov8[6, 6], cov8[7, 7]
    cv = lambda th: implied(th, n, E, order)[1][np.ix_(xi, xi)][iu]
    th = np.tanh(np.random.default_rng(0).normal(0, 0.5, ne))
    J = np.array([(cv(th + 1e-6 * np.eye(ne)[k]) - cv(th - 1e-6 * np.eye(ne)[k])) / 2e-6 for k in range(ne)]).T
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(sv > 1e-8 * sv[0]))

    def f(t):
        _, _, om = implied(t, n, E, order)
        return nll(t, n, E, order, xi, C, pen=0) + PEN * np.sum(np.clip(-om, 0, None) ** 2)

    sols = []
    for s in range(nstart):
        t0 = np.tanh(np.random.default_rng(7000 + s).normal(0, 0.5, ne))
        r = minimize(f, t0, method='L-BFGS-B', bounds=[(-0.999, 0.999)] * ne, options=TIGHT)
        F, S, om = implied(r.x, n, E, order)
        if om.min() < -1e-6 or np.linalg.slogdet(S[np.ix_(xi, xi)])[0] <= 0:
            continue
        sols.append((nll(r.x, n, E, order, xi, C, pen=0), r.x, om.min(), names[int(np.argmin(om))]))
    if not sols:
        return rank, ne, 0, np.nan, pd.DataFrame(columns=['latent', 'nll', 'omega_min', 'omega_min_node', 'r_Rbir', 'r_Rg']), None
    sols.sort(key=lambda z: z[0])
    best = sols[0][0]
    out, best_info = [], None
    for v, t, ommin, omnode in sols:
        if 2 * N * (v - best) >= LR95:
            continue
        coef = {(names[a], names[b]): x for (a, b), x in zip(E, t)}
        W = latent_weights(names, edges_named, coef)
        if best_info is None:                       # sols 按 NLL 从小到大排, 第一个就是最大似然解
            _, _, om = implied(t, n, E, order)
            best_info = dict(nll=v, edges=[list(e) for e in edges_named], coef=[coef[e] for e in edges_named],
                             omega={names[k]: float(om[k]) for k in range(n)},
                             weights={L: [float(x) for x in wv] for L, wv in W.items()})
        for L, wv in W.items():
            q = wv @ C @ wv
            out.append(dict(latent=L, nll=v, omega_min=ommin, omega_min_node=omnode,
                            r_Rbir=abs(wv @ cR) / np.sqrt(q * vR) if q > 0 else np.nan,
                            r_Rg=abs(wv @ cG) / np.sqrt(q * vG) if q > 0 else np.nan))
    cols = ['latent', 'nll', 'omega_min', 'omega_min_node', 'r_Rbir', 'r_Rg']
    return rank, ne, len(sols), best, pd.DataFrame(out, columns=cols), best_info


def baselines(cov8):
    sd = np.sqrt(np.diag(cov8))
    R = cov8 / np.outer(sd, sd)
    C, cR = R[:6, :6], R[:6, 6]
    return dict(age=abs(R[3, 6]), FeH=abs(R[4, 6]), JPhi=abs(R[1, 6]), ceiling=np.sqrt(cR @ np.linalg.solve(C, cR)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cond', default='main', choices=['main', 'R_0_30_Z_30', 'R_0_30_Z_3', 'R_0_30_Z_0.3'])
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--min-frac', type=float, default=0.2)
    ap.add_argument('--nstart', type=int, default=50)
    ap.add_argument('--ids', nargs='*', default=None)
    args = ap.parse_args()
    _, out = cond_paths(args.cond)
    st = pd.read_csv(out / f'structures_{args.version}.csv')
    st['galaxy'] = st.galaxy.astype(str)
    rows, lines = [], [f'做法 B 评估 ({args.cond}, {args.version}): 每个星系评估出现频率 >= {args.min_frac:.0%} 的结构 (至少最常见的一种)',
                        '每个隐变量: 最大似然解的 |r(L, 出生半径)|、|r(L, 导引半径)|, 括号里是统计上与最优分不开 (2·N·ΔNLL<3.84) 的解的范围', '']
    bestjson = {}
    for gid in args.ids or IDS:
        g = st[st.galaxy == gid]
        if not len(g):
            continue
        cov8, N = moments(args.cond, args.version, gid)
        b = baselines(cov8)
        freq = g.structure_id.value_counts()
        pick = [sid for sid, c in freq.items() if c / len(g) >= args.min_frac] or [freq.index[0]]
        lines.append(f'== {gid} (N={N:,}); 参照: 单用年龄 {b["age"]:.2f}, [Fe/H] {b["FeH"]:.2f}, J_phi {b["JPhi"]:.2f}, 线性上限 {b["ceiling"]:.2f}')
        for sid in pick:
            gg = g[g.structure_id == sid]
            if gg.n_undirected.iloc[0] > 0:
                lines.append(f'  结构 {sid} ({len(gg)}/{len(g)} 次): 含无向边, 跳过')
                continue
            edges = [tuple(e) for e in json.loads(gg.structure.iloc[0])]
            rank, ne, nvalid, best, sol, info = evaluate(edges, cov8, N, args.nstart)
            lines.append(f'  结构 {sid} ({len(gg)}/{len(g)} 次, {gg.n_latent.iloc[0]} 个隐变量, {ne} 条边): 雅可比秩 {rank}/{ne}'
                         f'{" (可识别)" if rank == ne else " (不可识别)"}; 有效解 {nvalid}, 最优 NLL {best:.4f}')
            if info is not None:
                bestjson.setdefault(gid, {})[str(sid)] = dict(freq=len(gg) / len(g), n_seeds=len(g), jacobian_rank=rank,
                                                              n_edges=ne, N=N, **info)
            for L, s in sol.groupby('latent', sort=True):
                lines.append(f'    {L}: 最大似然解 |r(L,出生半径)| {s.r_Rbir.iloc[0]:.2f} (分不开的解 {s.r_Rbir.min():.2f}–{s.r_Rbir.max():.2f})'
                             f'   |r(L,导引半径)| {s.r_Rg.iloc[0]:.2f} ({s.r_Rg.min():.2f}–{s.r_Rg.max():.2f})')
                rows.append(dict(galaxy=gid, structure_id=sid, freq=len(gg) / len(g), n_latent=gg.n_latent.iloc[0], n_edges=ne,
                                 jacobian_rank=rank, latent=L, n_indistinguishable=len(s),
                                 r_Rbir_best=s.r_Rbir.iloc[0], r_Rbir_min=s.r_Rbir.min(), r_Rbir_max=s.r_Rbir.max(),
                                 r_Rg_best=s.r_Rg.iloc[0], r_Rg_min=s.r_Rg.min(), r_Rg_max=s.r_Rg.max(),
                                 omega_min_best=s.omega_min.iloc[0], omega_min_node=s.omega_min_node.iloc[0],
                                 base_age=b['age'], base_FeH=b['FeH'], base_JPhi=b['JPhi'], ceiling=b['ceiling']))
            print('\n'.join(lines[-(1 + sol.latent.nunique()):]), flush=True)
    pd.DataFrame(rows).to_csv(out / f'eval_{args.version}.csv', index=False)
    (out / f'eval_{args.version}.txt').write_text('\n'.join(lines) + '\n')
    json.dump(bestjson, open(out / f'eval_{args.version}_best.json', 'w'), indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
