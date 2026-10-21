# -*- coding: utf-8 -*-
"""
rlcd_A.py —— 做法 A: 冻结作者的因果图, 在每个星系上重新估计边系数, 算出每颗星的 L1、L2, 并与真值对比

必须用 RLCD 的环境运行 (需要 torch、pydot、pandas、matplotlib), 在 casuality 文件夹里:
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version g279            # 先用作者的数据验证整个流程
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version nonoise --ids 538905   # 试跑一个星系
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version nonoise          # 12 个星系, 无误差
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version noise            # 12 个星系, 加误差

流程 (与论文 3.3 节和作者的 notebook 相同):
  1. 读 6 个观测量, 标准化 (减平均值、除以标准差), 与仓库 DataModel.generate_data(normalized=True) 相同
  2. 因果图: 作者的 g2.79e12 最终图 (results/rlcd_g279/author_...dot), 结构固定不变
  3. 参数学习: 调用仓库 ParameterLearning/linear_parameter_identification.py 里的 parameter_by_trek
     (与 README 相同: lr=0.005, epochs=600, lr_decay=plateau)。它从随机初值开始优化, 结果与初值有关,
     所以重复 N_REPEAT 次 (torch 随机种子 0..N-1, 默认 50 次), 再从中挑一组:
       (a) 先去掉物理上不成立的解: 模型隐含的某个节点噪声方差 < 0, 或协方差矩阵不正定。
           仓库代码不检查这一点, 而且它算 NLL 时用的是行列式的绝对值, 不成立的解可能得到"很好"的 NLL
           (g2.79e12 上出现过 NLL = -7.03, 比理论最小值 -1.27 还小, 这是不可能的)
       (b) 在剩下的解里, 默认挑模型相关矩阵与实际相关矩阵之差 (RMSE) 最小的 (--select rmse);
           也可以挑 NLL 最小的 (--select nll)。
       理由: 作者 notebook 和论文用的系数是仓库里 NLL 第二小 (1.586) 的那组, 它的 RMSE 最小 (0.117);
       NLL 最小 (1.498) 的那组 RMSE 是 0.260, 而且 L1 几乎和化学量断开, L1 与出生半径的相关只有 0.01。
       NLL 被 [Fe/H] 与 [Mg/H] 的 0.997 强相关主导, 不是好的挑选标准。
     每次的 L1、L2 与真值的相关系数也记在 nll_repeats.csv 里, 只用于诊断, 不参与挑选
  4. 算每颗星的隐变量 (论文 3.3 节: 从最下层的隐变量开始, 逐层往上):
       对隐变量 L, 它的每个子节点 c 都满足  c = a_Lc * L + Σ(其他父节点 p) a_pc * p + 噪声
       令 r_c = c - Σ(其他父节点) a_pc * p  (去掉其他父节点贡献后, 剩下应由 L 解释的部分),
       使 Σ_c (r_c - a_Lc * L)^2 最小的 L 是   L = Σ_c a_Lc * r_c / Σ_c a_Lc^2
     这与作者 notebook 里逐颗星用 Nelder-Mead 最小化 chisq 的结果完全相同 (二次函数的最小值有解析解)
  5. 线性换算到物理单位: f(L1) = α1*L1 + β1 拟合真值 Rbirth, f(L2) = α2*L2 + β2 拟合真值 Rg (论文第 5-6 页)
     给出皮尔逊相关系数 r、残差的标准差、中位绝对偏差

输出: results/rlcd_A/<version>/<ID>/
        nll_repeats.csv      每次重复的 NLL 和所有边系数
        best_params.dot/.png NLL 最小那次的系数画在因果图上
        latents.csv          每颗星: iord (若有), L1, L2, f(L1), f(L2), R_bir, Rg
        fit.txt              换算公式和指标
      results/rlcd_A/<version>/summary.csv, summary.txt
      results/rlcd_A/<version>/L1_vs_Rbir.png, L2_vs_Rg.png  (每个星系一个小图)
"""
import os
import io
import sys
import time
import argparse
import contextlib
import numpy as np
import pandas as pd

REPO = 'Causal-Galactic-Archaeology'
sys.path.insert(0, os.path.join(REPO, 'rlcd', 'ParameterLearning'))
sys.path.insert(0, os.path.join(REPO, 'rlcd'))
import torch                                              # noqa: E402
import pydot                                              # noqa: E402
from linear_parameter_identification import (             # noqa: E402
    parameter_by_trek, get_F_bin_and_varnames_from_dotgraph, add_edge_coefficients_to_graph)

AUTHOR_DOT = 'results/rlcd_g279/author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def load(version, gid):
    """返回 (观测量 DataFrame, 真值 DataFrame)"""
    if version == 'g279':
        e = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_e.csv'))
        t = pd.read_csv(os.path.join(REPO, 'data', 'g2.79e12_multi_ele_action_Jin_SN.csv'))
    else:
        e = pd.read_csv(f'data/rlcd_input/{version}/{gid}_e.csv')
        t = pd.read_csv(f'data/rlcd_input/{version}/{gid}_truth.csv')
    assert len(e) == len(t)
    ok = np.isfinite(e[COLS].to_numpy()).all(axis=1)
    return e.loc[ok, COLS].reset_index(drop=True), t.loc[ok].reset_index(drop=True)


def estimate_latents(Xn, F, names):
    """按论文 3.3 节逐层求每颗星的隐变量 (解析最小二乘)"""
    val = {f'X_{c}': Xn[c].to_numpy() for c in COLS}
    lat = [v for v in names if not v.startswith('X_')]
    idx = {v: i for i, v in enumerate(names)}
    done = {}
    while len(done) < len(lat):
        progress = False
        for L in lat:
            if L in done:
                continue
            children = [c for c in names if F[idx[L], idx[c]] != 0]
            if any(c in lat and c not in done for c in children):
                continue                       # 先算下层的隐变量
            num = 0.0
            den = 0.0
            for c in children:
                xc = val[c] if c in val else done[c]
                r = xc.copy()
                for p in names:
                    if p != L and F[idx[p], idx[c]] != 0:
                        xp = val[p] if p in val else done.get(p)
                        if xp is None:
                            raise RuntimeError(f"{c} 的父节点 {p} 还没算出, 无法求 {L}")
                        r = r - F[idx[p], idx[c]] * xp
                a = F[idx[L], idx[c]]
                num = num + a * r
                den = den + a * a
            done[L] = num / den
            progress = True
        if not progress:
            raise RuntimeError('隐变量之间有环, 无法逐层求解')
    return done


def implied_cov(F, names):
    """由边系数算出模型隐含的协方差矩阵 (所有变量都标准化为方差 1, 与 parameter_by_trek 的 trek 规则等价),
    并给出每个节点的噪声方差 ω_k = 1 - Var(父节点的贡献)。ω_k < 0 说明这组系数在物理上不成立
    (某个变量被父节点"解释"得超过了它自己的全部方差), 这种解不能用。"""
    n = len(names)
    order, rem = [], list(range(n))
    while rem:
        for k in rem:
            if all(F[p, k] == 0 or p in order for p in range(n)):
                order.append(k)
                rem.remove(k)
                break
        else:
            raise RuntimeError('因果图里有环')
    S = np.zeros((n, n))
    omega = np.zeros(n)
    for pos, k in enumerate(order):
        par = [p for p in order[:pos] if F[p, k] != 0]
        for j in order[:pos]:
            S[j, k] = S[k, j] = sum(F[p, k] * S[j, p] for p in par)
        omega[k] = 1 - sum(F[p, k] * F[q, k] * S[p, q] for p in par for q in par)
        S[k, k] = 1.0
    return S, omega


def evaluate_solution(F, names, S_emp):
    """一组边系数的质量: 是否物理上成立、正确的 NLL、与实际相关矩阵的差 (RMSE)"""
    S, omega = implied_cov(F, names)
    xi = [i for i, v in enumerate(names) if v.startswith('X_')]
    Sx = S[np.ix_(xi, xi)]
    sign, logdet = np.linalg.slogdet(Sx)
    valid = bool(sign > 0 and omega.min() > -1e-6)
    nll = 0.5 * logdet + 0.5 * np.trace(np.linalg.solve(Sx, S_emp)) if sign > 0 else np.inf
    iu = np.triu_indices(len(xi), 1)
    rmse = float(np.sqrt(np.mean((Sx - S_emp)[iu] ** 2)))
    return dict(valid=valid, nll_true=float(nll), rmse=rmse, min_noise_var=float(omega.min()))


def fit_linear(L, truth):
    ok = np.isfinite(L) & np.isfinite(truth)
    a, b = np.polyfit(L[ok], truth[ok], 1)
    f = a * L + b
    res = truth[ok] - f[ok]
    r = np.corrcoef(f[ok], truth[ok])[0, 1]
    return a, b, f, r, float(np.std(res)), float(np.median(np.abs(res)))


def run_one(version, gid, dot_path, n_repeat, epochs, lr, lr_decay, outdir, select='rmse'):
    t0 = time.time()
    X, T = load(version, gid)
    Xn = (X - X.mean()) / X.std()
    df = Xn.rename(columns={c: f'X_{c}' for c in COLS})

    graph = pydot.graph_from_dot_file(dot_path)[0]
    F_bin, names = get_F_bin_and_varnames_from_dotgraph(graph)
    df = df[[v for v in names if v.startswith('X_')]]

    S_emp = df.cov().to_numpy()
    rows, sols = [], []
    for seed in range(n_repeat):
        torch.manual_seed(seed)
        np.random.seed(seed)
        with contextlib.redirect_stdout(io.StringIO()):
            F, nll = parameter_by_trek(df, F_bin, names, lr=lr, epochs=epochs, lr_decay=lr_decay)
        ev = evaluate_solution(F, names, S_emp)
        # 诊断用 (不参与挑选!): 这组系数给出的 L1、L2 与真值的相关系数
        lat_i = estimate_latents(Xn, F, names)
        diag_r1 = abs(np.corrcoef(lat_i['L1'], T['R_bir'])[0, 1])
        diag_r2 = abs(np.corrcoef(lat_i['L2'], T['Rg'])[0, 1])
        row = {'seed': seed, 'nll_code': float(nll), **ev,
               'diag_abs_r_L1_Rbir': diag_r1, 'diag_abs_r_L2_Rg': diag_r2}
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if F_bin[i, j]:
                    row[f'{a}->{b}'] = F[i, j]
        rows.append(row)
        sols.append((F, row))
        print(f"    [{gid}] 第 {seed+1:>2}/{n_repeat} 次: {'有效' if ev['valid'] else '无效'}  "
              f"NLL {ev['nll_true']:8.4f}  RMSE {ev['rmse']:.4f}   (诊断 |r(L1,Rb)| {diag_r1:.2f}, |r(L2,Rg)| {diag_r2:.2f})",
              flush=True)

    os.makedirs(outdir, exist_ok=True)
    tab = pd.DataFrame(rows)
    valid = [x for x in sols if x[1]['valid']]
    if not valid:
        raise RuntimeError(f"{gid}: {n_repeat} 次参数学习都没有物理上成立的解")
    key = 'rmse' if select == 'rmse' else 'nll_true'
    F, brow = min(valid, key=lambda x: x[1][key])
    nll, seed = brow['nll_true'], brow['seed']
    tab['selected'] = tab['seed'] == seed
    tab.sort_values(['valid', key], ascending=[False, True]).to_csv(os.path.join(outdir, 'nll_repeats.csv'), index=False)
    n_valid = len(valid)
    top = sorted(valid, key=lambda x: x[1][key])[:5]
    top_txt = ', '.join(f"{x[1]['diag_abs_r_L1_Rbir']:.2f}" for x in top)
    add_edge_coefficients_to_graph(graph, F_bin, F, names)
    graph.write_raw(os.path.join(outdir, 'best_params.dot'))
    try:
        graph.write_png(os.path.join(outdir, 'best_params.png'))
    except Exception as err:
        print(f"    (画 png 失败, 可能没装 graphviz: {err})")

    lat = estimate_latents(Xn, F, names)
    L1, L2 = lat['L1'], lat['L2']
    a1, b1, f1, r1, s1, m1 = fit_linear(L1, T['R_bir'].to_numpy())
    a2, b2, f2, r2, s2, m2 = fit_linear(L2, T['Rg'].to_numpy())

    out = pd.DataFrame({'L1': L1, 'L2': L2, 'f_L1': f1, 'f_L2': f2,
                        'R_bir': T['R_bir'], 'Rg': T['Rg'], 'age': X['age']})
    if 'iord' in T:
        out.insert(0, 'iord', T['iord'])
    out.to_csv(os.path.join(outdir, 'latents.csv'), index=False)

    txt = (f"星系 {gid} ({version}), 恒星 {len(X)}; {n_repeat} 次参数学习中物理上成立的 {n_valid} 次; "
           f"按 {key} 挑选: seed {seed}, NLL {nll:.4f}, RMSE {brow['rmse']:.4f}\n"
           f"  Rbirth = {a1:+.3f} * L1 {b1:+.3f}   皮尔逊 r = {r1:.3f}   残差标准差 {s1:.3f} kpc   中位绝对偏差 {m1:.3f} kpc\n"
           f"  Rg     = {a2:+.3f} * L2 {b2:+.3f}   皮尔逊 r = {r2:.3f}   残差标准差 {s2:.3f} kpc   中位绝对偏差 {m2:.3f} kpc\n"
           f"  (稳定性诊断: 按 {key} 排前 5 的有效解, |r(L1,Rb)| = {top_txt})\n"
           f"  用时 {time.time()-t0:.0f} 秒\n")
    with open(os.path.join(outdir, 'fit.txt'), 'w') as f:
        f.write(txt)
    print(txt, flush=True)
    return dict(ID=gid, n=len(X), n_valid=n_valid, best_nll=nll, best_rmse=brow['rmse'], best_seed=seed,
                a1=a1, b1=b1, r_L1_Rbir=r1, std_Rbir=s1, mad_Rbir=m1,
                a2=a2, b2=b2, r_L2_Rg=r2, std_Rg=s2, mad_Rg=m2), out


def plots(results, outroot, version):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(results)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    for key, truth, lab, rkey in [('f_L1', 'R_bir', 'Rbirth', 'r_L1_Rbir'), ('f_L2', 'Rg', 'Rg', 'r_L2_Rg')]:
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.2 * nrow), squeeze=False)
        axs = axs.ravel()
        for ax, (summ, out) in zip(axs, results):
            lim = (0, 16)
            ax.hexbin(out[key], out[truth], gridsize=45, extent=lim + lim, bins='log', cmap='Blues', mincnt=1)
            ax.plot(lim, lim, color='#1f2937', lw=0.8)
            ax.set_title(f"{summ['ID']}  r = {summ[rkey]:.2f}", fontsize=9)
            ax.set_xlabel(f'f({key[2:]}) [kpc]', fontsize=8)
            ax.set_ylabel(f'true {lab} [kpc]', fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in axs[n:]:
            ax.axis('off')
        fig.suptitle(f'{version}: f({key[2:]}) vs true {lab}', fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(outroot, f'{key[2:]}_vs_{truth}.png'), dpi=120)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['g279', 'nonoise', 'noise'])
    ap.add_argument('--ids', nargs='*', help='只算这些星系 (默认 12 个全算)')
    ap.add_argument('--n-repeat', type=int, default=50, help='参数学习重复次数')
    ap.add_argument('--select', default='rmse', choices=['rmse', 'nll'],
                    help='在物理上成立的解里, 按什么挑: rmse = 模型相关矩阵与实际相关矩阵之差最小 (默认); nll = 负对数似然最小')
    ap.add_argument('--epochs', type=int, default=600)
    ap.add_argument('--lr', type=float, default=0.005)
    ap.add_argument('--lr-decay', default='plateau')
    ap.add_argument('--dot', default=AUTHOR_DOT, help='冻结的因果图 (默认作者的 g2.79e12 图)')
    args = ap.parse_args()

    if not os.path.isdir(REPO):
        sys.exit('请在 casuality 文件夹里运行')
    ids = ['g279'] if args.version == 'g279' else (args.ids or IDS)
    outroot = os.path.join('results', 'rlcd_A', args.version)
    os.makedirs(outroot, exist_ok=True)
    print(f"做法 A: 冻结的因果图 {args.dot}\n版本 {args.version}, 星系 {ids}, 每个星系参数学习重复 {args.n_repeat} 次\n")

    results = []
    for gid in ids:
        summ, out = run_one(args.version, gid, args.dot, args.n_repeat, args.epochs, args.lr,
                            args.lr_decay, os.path.join(outroot, gid), select=args.select)
        results.append((summ, out))

    summ = pd.DataFrame([s for s, _ in results])
    tag = '' if (args.version == 'g279' or not args.ids) else '_' + '_'.join(ids)
    summ.to_csv(os.path.join(outroot, f'summary{tag}.csv'), index=False)
    with open(os.path.join(outroot, f'summary{tag}.txt'), 'w') as f:
        f.write(summ.to_string(index=False) + '\n')
    print(summ[['ID', 'n', 'n_valid', 'best_nll', 'best_rmse', 'r_L1_Rbir', 'std_Rbir', 'r_L2_Rg', 'std_Rg']].to_string(index=False))
    if args.version == 'g279':
        print("\n论文 (g2.79e12): Rb = -1.97 L1 + 9.31, r(L1, Rb) = 0.64;  Rg = -1.29 L2 + 9.52, r(L2, Rg) = 0.98")
        print("作者实际用的是 NLL 第二小的那组系数 (仓库 nll1.586 文件), 它的 RMSE 最小 (0.117)")
        print("注意: 隐变量的正负号是任意的, 所以斜率的符号可能和论文相反, 这不影响结果")
    plots(results, outroot, args.version)
    print(f"\n输出在 {outroot}/")


if __name__ == '__main__':
    main()
