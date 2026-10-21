# -*- coding: utf-8 -*-
"""
summary_planB.py —— 做法 B 第 3 步: 12 个星系一行一个的汇总表, 并与做法 A (冻结作者的图) 对比

先运行同一文件夹的 rlcd_B.py 和 rlcd_B_eval.py, 再运行 (numpy、pandas):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 PlanB/output_data/summary_planB.py --version nonoise
    ~/PYTHONPATH_01/venv/bin/python3 PlanB/output_data/summary_planB.py --version noise
输出: PlanB/output_data/summary_planB_<版本>.csv 和 .txt

每个星系只看出现次数最多的结构 (main 条件: 7<R<10 kpc, |z|<0.3 kpc, 与论文和做法 A 主结果相同):
  n_structures   30 个种子学到几种不同的图;  top_freq 最常见的图出现几次
  n_latent / n_edges / n_undirected   隐变量个数 / 有向边 / 无向边 (有无向边就不能求系数)
  same_as_author 是否就是作者的图;  prior_ok 是否满足论文 3.2 节的先验 (年龄没有观测量父节点)
  jac_rank       雅可比秩 (= 边数时边系数可由数据唯一确定)
  B_best_L_Rb    与真出生半径 |r| 最大的隐变量 (用真值挑的, 最好情况) 和它在最大似然解下的 |r|
  B_best_L_Rg    同上, 导引半径
  A_r1_range     做法 A (作者的图冻结, results/rlcd_A/identifiability/summary_ci_<版本>.csv) 里 L1 的 |r(L1,出生半径)| 范围
  NLL_B / NLL_A  两种图的最大似然 NLL (平均每颗星的负对数似然, 越小越贴合数据; 两者用的是同一个观测相关矩阵, 可以直接比)
  dBIC_per_star  (NLL_B − NLL_A) + (边数_B − 15)·ln N / (2N): 按 BIC 把"边多占便宜"扣掉后的差 (N 为星数, 15 是作者图的边数);
                 负数 = 学到的图更好, 正数 = 作者的图更好
  age_only / JPhi_only / ceiling   单用年龄 / 单用 J_phi 与出生半径的 |r|; 6 个输入线性回归的上限
"""
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    v = ap.parse_args().version
    st = pd.read_csv(HERE.with_name(f'structures_{v}.csv'))
    st['galaxy'] = st.galaxy.astype(str)
    ev = pd.read_csv(HERE.with_name(f'eval_{v}.csv'))
    ev['galaxy'] = ev.galaxy.astype(str)
    best = json.load(open(HERE.with_name(f'eval_{v}_best.json')))
    A = pd.read_csv(ROOT / 'results' / 'rlcd_A' / 'identifiability' / f'summary_ci_{v}.csv').set_index('tag')
    rows = []
    for gid in IDS:
        g = st[st.galaxy == gid]
        freq = g.structure_id.value_counts()
        sid = freq.index[0]
        r = g[g.structure_id == sid].iloc[0]
        e = ev[(ev.galaxy == gid) & (ev.structure_id == sid)]
        info = best.get(gid, {}).get(str(sid))
        a = A.loc[f'{v}_{gid}']
        row = dict(galaxy=gid, n_structures=len(freq), top_structure=sid, top_freq=f'{freq.iloc[0]}/{len(g)}',
                   n_latent=r.n_latent, n_edges=r.n_directed, n_undirected=r.n_undirected,
                   same_as_author=bool(r.same_as_author), prior_ok=bool(r.prior_ok),
                   jac_rank=info['jacobian_rank'] if info else None)
        if len(e):
            iR, iG = e.r_Rbir_best.idxmax(), e.r_Rg_best.idxmax()
            row.update(B_best_L_Rb=f'{e.latent[iR]} {e.r_Rbir_best[iR]:.2f}', B_r_Rb=e.r_Rbir_best[iR],
                       B_best_L_Rg=f'{e.latent[iG]} {e.r_Rg_best[iG]:.2f}', B_r_Rg=e.r_Rg_best[iG])
        else:
            row.update(B_best_L_Rb='—', B_r_Rb=None, B_best_L_Rg='—', B_r_Rg=None)
        x = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')[COLS]
        y = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir']).R_bir.to_numpy()
        z = ((x - x.mean()) / x.std()).to_numpy()
        N = len(z)
        rb = lambda c: abs(np.corrcoef(x[c], y)[0, 1])
        X1 = np.column_stack([z, np.ones(N)])
        ceil = np.corrcoef(X1 @ np.linalg.lstsq(X1, y, rcond=None)[0], y)[0, 1]
        # BIC 差 (每颗星): (NLL_B − NLL_A) + (边数_B − 15)·ln N / (2N); 负数表示学到的图更好
        dbic = (info['nll'] - a.best_nll + (info['n_edges'] - 15) * np.log(N) / (2 * N)) if info else None
        row.update(A_r1_range=f'{a.r1_min:.2f}–{a.r1_max:.2f}', A_r2_range=f'{a.r2_min:.2f}–{a.r2_max:.2f}',
                   NLL_B=info['nll'] if info else None, NLL_A=a.best_nll, dBIC_per_star=dbic,
                   age_only=rb('age'), JPhi_only=rb('JPhi'), ceiling=ceil)
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab.to_csv(HERE.with_name(f'summary_planB_{v}.csv'), index=False)
    show = tab.drop(columns=['B_r_Rb', 'B_r_Rg'])
    txt = (f'做法 B 汇总 (main 条件 7<R<10 kpc, |z|<0.3 kpc; {v}); 每个星系只看 30 个种子里最常见的结构\n'
           'B_best_L_Rb / _Rg: 与真出生半径 / 真导引半径 |r| 最大的隐变量 (用真值挑的最好情况) 及其 |r| (最大似然解)\n'
           'A_r1_range / A_r2_range: 做法 A (冻结作者的图) 统计上分不开的解里 |r(L1,出生半径)| / |r(L2,导引半径)| 的范围\n'
           'NLL_B / NLL_A: 学到的图 / 作者的图的最大似然 NLL (越小越贴合数据); dBIC_per_star > 0 表示扣掉边数差别后作者的图仍更贴合\n\n'
           + show.round(3).to_string(index=False) + '\n')
    if tab.B_r_Rb.notna().any():
        ok = tab.dropna(subset=['B_r_Rb'])
        txt += (f'\n能求系数的星系 {len(ok)} 个: 最匹配隐变量的 |r(L,出生半径)| 中位数 {ok.B_r_Rb.median():.2f} '
                f'(范围 {ok.B_r_Rb.min():.2f}–{ok.B_r_Rb.max():.2f}); 超过"单用年龄"的星系: '
                f'{", ".join(ok.galaxy[ok.B_r_Rb > ok.age_only]) or "无"}\n')
    HERE.with_name(f'summary_planB_{v}.txt').write_text(txt)
    print(txt)


if __name__ == '__main__':
    main()
