# -*- coding: utf-8 -*-
"""
rlcd_B.py —— 做法 B 第 1 步: 在每个 TNG50 星系上, 用多个随机种子从头跑 RLCD 结构学习, 看数据自己会学出什么因果图

做法 B 第 1 步 (本脚本) → 第 2 步 rlcd_B_eval.py (同一文件夹) → 画图: PlanB/因果图/, PlanB/真值对比/

用 RLCD 专用环境运行 (需要 causal-learn 等):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 PlanB/output_data/rlcd_B.py --version nonoise
    ~/PYTHONPATH_01/venv_rlcd/bin/python3 PlanB/output_data/rlcd_B.py --version noise
    可选: --seeds 20 (PYTHONHASHSEED = 0..19), --ids 538905 (只算指定星系), --alpha 1e-5 1e-15 1e-5 1e-15,
          --cond (默认 main); 已算过的 (星系, 种子) 会跳过, 所以中断后重新运行会接着算

选星条件 (--cond):
  main           7 < R < 10 kpc, |z| < 0.3 kpc (与论文相同; 数据 PlanB/input_data/<版本>/)      → 输出 PlanB/output_data/
  R_0_30_Z_30    0 < R < 30 kpc, |z| < 30 kpc  (star_condition/R_0_30_Z_30/output_data/moments) → 输出 PlanB/output_data/R_0_30_Z_30/
  R_0_30_Z_3     0 < R < 30 kpc, |z| < 3 kpc   (同上)                                          → 输出 PlanB/output_data/R_0_30_Z_3/
  R_0_30_Z_0.3   0 < R < 30 kpc, |z| < 0.3 kpc (同上)                                          → 输出 PlanB/output_data/R_0_30_Z_0.3/

做法:
  1. RLCD 的设置与作者完全相同 (rlcd/run_rlcd.py 默认: 第一阶段 GES, sparsity 2; alpha 默认用作者结果文件名里的
     1e-5, 1e-15, 1e-5, 1e-15)。每个种子单独开一个 Python 进程并设置 PYTHONHASHSEED (哈希种子只能在启动时设),
     记录学到的图 (有向边、无向边、隐变量个数)。
  2. 为什么可以只用协方差矩阵: RLCD 用到数据的地方只有两处 —— 秩检验 (典型相关分析, 只用协方差矩阵和星数 N) 和
     第一阶段 GES 的 BIC 打分 (np.cov 和 N)。所以对 37–72 万颗星的大数据集, 本脚本用一个"样本协方差矩阵与真实数据
     完全相同、星数也相同"的高斯数据集代替原始数据, 结果与用原始数据完全一样 (--selftest 在 7<R<10 的 538905 上核对)。
     main 条件的数据本来就小, 直接用原始数据。
  3. 每张图与作者的图比较 (隐变量编号自动对齐), 并检查论文 3.2 节的先验 (年龄没有观测父节点)。
输出 (<输出目录>/):
  structures_<版本>.csv   每个星系、每个种子一行: 隐变量个数、有向边列表、无向边、是否与作者的图相同、先验是否通过、结构编号
  structures_<版本>.txt   每个星系出现过的各种结构及其频率
"""
import os
import sys
import json
import argparse
import itertools
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
REPO = ROOT / 'Causal-Galactic-Archaeology'
AUTHOR_DOT = ROOT / 'results' / 'rlcd_g279' / 'author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
OBS = set(COLS)


def cond_paths(cond):
    """返回 (输入所在目录, 输出目录)"""
    if cond == 'main':
        return PLANB / 'input_data', PLANB / 'output_data'
    return ROOT / 'star_condition' / cond, PLANB / 'output_data' / cond


def load_data(cond, version, gid):
    """返回 (6 列 DataFrame, 是否为'协方差相同的替身数据')"""
    src, _ = cond_paths(cond)
    if cond == 'main':
        return pd.read_csv(src / version / f'{gid}_e.csv')[COLS], False
    m = json.load(open(src / 'output_data' / f'moments_{version}.json'))[gid]
    return surrogate(np.array(m['cov'])[:6, :6], np.array(m['mean'])[:6], m['N']), True


def surrogate(cov, mean, N, seed=0):
    """造一个 N 行的高斯数据集, 它的样本平均值和样本协方差 (分母 N-1) 与给定的完全相同"""
    Z = np.random.default_rng(seed).standard_normal((N, len(mean)))
    Z -= Z.mean(0)
    Z = Z @ np.linalg.inv(np.linalg.cholesky(np.cov(Z.T, ddof=1))).T      # 样本协方差 = 单位矩阵
    X = Z @ np.linalg.cholesky(cov).T + mean
    return pd.DataFrame(X, columns=COLS)


def run_rlcd(df, alpha):
    """在本进程里跑 RLCD (与 run_rlcd.py 相同的设置), 返回 (有向边, 无向边), 变量名去掉 X_"""
    os.chdir(REPO / 'rlcd')
    sys.path.insert(0, str(REPO / 'rlcd'))
    import logging
    from DGM.DataModel import DataModel
    from utils.Chi2RankTest import Chi2RankTest
    from StructureLearning.RLCD.RLCD_alg import RLCD
    logging.getLogger().setLevel(logging.WARNING)
    df = df.rename(columns={k: 'X_' + k for k in df.columns})
    dgm = DataModel(df, df)
    df_x, _ = dgm.generate_data(N=-1, normalized=True)
    params = {'ranktest_method': Chi2RankTest(df_x.to_numpy(), 1), 'citest_method': None, 'stage1_method': 'ges',
              'alpha_dict': dict(enumerate(alpha)), 'stage1_ges_sparsity': 2, 'stage1_partition_thres': 3}
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):          # RLCD 打印很多中间信息
        graph, _, _, _ = RLCD(1, dgm.xvars, df_x, params)
    s = lambda v: v.replace('X_', '')
    return sorted((s(a), s(b)) for a, b in graph.dirEdges), sorted(tuple(sorted(map(s, e))) for e in graph.undirEdges)


def author_edges():
    import re
    txt = open(AUTHOR_DOT).read()
    return {(a.replace('X_', ''), b.replace('X_', '')) for a, b in re.findall(r'^\s*(\w+)\s*->\s*(\w+)', txt, re.M)}


def canonical(edges):
    """隐变量编号是任意的: 在所有换名方式里取字典序最小的边列表, 作为结构的"身份证" """
    lat = sorted({v for e in edges for v in e if v not in OBS})
    best = None
    for perm in itertools.permutations(range(1, len(lat) + 1)):
        mp = {L: f'L{p}' for L, p in zip(lat, perm)}
        key = tuple(sorted((mp.get(a, a), mp.get(b, b)) for a, b in edges))
        if best is None or key < best:
            best = key
    return best


def worker(args):
    df, _ = load_data(args.cond, args.version, args.ids[0])
    d, u = run_rlcd(df, args.alpha)
    print('RESULT ' + json.dumps({'directed': d, 'undirected': u}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cond', default='main', choices=['main', 'R_0_30_Z_30', 'R_0_30_Z_3', 'R_0_30_Z_0.3'])
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    ap.add_argument('--seeds', type=int, default=20)
    ap.add_argument('--ids', nargs='*', default=None)
    ap.add_argument('--alpha', nargs=4, type=float, default=[1e-5, 1e-15, 1e-5, 1e-15])
    ap.add_argument('--worker', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.worker:
        return worker(args)
    if args.selftest:
        return selftest(args)

    _, out = cond_paths(args.cond)
    out.mkdir(parents=True, exist_ok=True)
    csv = out / f'structures_{args.version}.csv'
    done = pd.read_csv(csv) if csv.exists() else pd.DataFrame()
    author = canonical(author_edges())
    rows = done.to_dict('records')
    for gid in args.ids or IDS:
        for s in range(args.seeds):
            if len(done) and ((done.galaxy.astype(str) == gid) & (done.seed == s)).any():
                continue
            env = dict(os.environ, PYTHONHASHSEED=str(s))
            cmd = [sys.executable, str(HERE), '--worker', '--cond', args.cond, '--version', args.version,
                   '--ids', gid, '--alpha', *map(str, args.alpha)]
            p = subprocess.run(cmd, env=env, capture_output=True, text=True)
            line = [l for l in p.stdout.splitlines() if l.startswith('RESULT ')]
            if not line:
                print(p.stderr[-2000:])
                raise RuntimeError(f'{gid} seed {s} 失败')
            res = json.loads(line[0][7:])
            d = [tuple(e) for e in res['directed']]
            key = canonical(d)
            lat = sorted({v for e in d for v in e if v not in OBS})
            rows.append(dict(galaxy=gid, seed=s, n_latent=len(lat), n_directed=len(d), n_undirected=len(res['undirected']),
                             same_as_author=(key == author and not res['undirected']),
                             prior_ok=not any(c == 'age' and p in OBS for p, c in d),
                             structure=json.dumps(key), undirected=json.dumps(res['undirected'])))
            print(f'{args.cond} {args.version} {gid} seed {s}: {len(lat)} 个隐变量, {len(d)} 条有向边'
                  f'{", 与作者的图相同" if rows[-1]["same_as_author"] else ""}', flush=True)
            pd.DataFrame(rows).to_csv(csv, index=False)
    summarize(out, args.version)


def summarize(out, version):
    df = pd.read_csv(out / f'structures_{version}.csv')
    codes = {k: i + 1 for i, k in enumerate(df.structure.value_counts().index)}   # 结构编号: 全部星系合计出现次数从多到少
    df['structure_id'] = df.structure.map(codes)
    df.to_csv(out / f'structures_{version}.csv', index=False)
    lines = [f'做法 B ({version}): 每个星系 {df.groupby("galaxy").size().max()} 个种子; 结构编号按所有星系合计的出现次数排序',
             f'与作者的图相同的次数: {int(df.same_as_author.sum())} / {len(df)}', '']
    for gid, g in df.groupby('galaxy'):
        lines.append(f'== {gid}')
        for sid, gg in g.groupby('structure_id'):
            r = gg.iloc[0]
            lines.append(f'  结构 {sid:>2}: {len(gg):>2} 次  隐变量 {r.n_latent}  有向边 {r.n_directed}  无向边 {r.n_undirected}  '
                         f'先验{"通过" if r.prior_ok else "不通过"}{"  ← 作者的图" if r.same_as_author else ""}  种子 {sorted(gg.seed.tolist())}')
    lines += ['', '各结构的边 (隐变量已统一编号):']
    for k, i in sorted(codes.items(), key=lambda x: x[1]):
        lines.append(f'  结构 {i}: ' + ', '.join(f'{a}->{b}' for a, b in json.loads(k)))
    (out / f'structures_{version}.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:3]))


def selftest(args):
    """核对: 协方差相同的替身数据与原始数据, RLCD 学到的图完全相同 (7<R<10 的 538905, 种子 1)"""
    gid = '538905'
    x = pd.read_csv(PLANB / 'input_data' / 'nonoise' / f'{gid}_e.csv')[COLS]
    sur = surrogate(np.cov(x.to_numpy().T, ddof=1), x.mean().to_numpy(), len(x))
    print('协方差最大差别', np.abs(np.cov(sur.to_numpy().T) - np.cov(x.to_numpy().T)).max())
    a = run_rlcd(x, args.alpha)
    b = run_rlcd(sur, args.alpha)
    print('原始数据:', a)
    print('替身数据:', b)
    print('完全相同' if a == b else '不同!')


if __name__ == '__main__':
    main()
