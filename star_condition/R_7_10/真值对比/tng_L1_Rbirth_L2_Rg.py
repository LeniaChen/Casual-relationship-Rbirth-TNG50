# -*- coding: utf-8 -*-
"""
tng_L1_Rbirth_L2_Rg.py —— TNG50 12 个星系: f(L1) vs 真 Rbirth、f(L2) vs 真 Rg (对应论文图 3、图 6 的做法)

运行 (numpy、pandas、matplotlib; 主 venv 即可):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 "results/TNG/真值对比/tng_L1_Rbirth_L2_Rg.py"                  # 不加误差 (主版本)
    ~/PYTHONPATH_01/venv/bin/python3 "results/TNG/真值对比/tng_L1_Rbirth_L2_Rg.py" --version noise  # 加误差 → ..._noise.png
输出: 与本脚本同名、同目录的 tng_L1_Rbirth_L2_Rg.png

做法 (与 复现_Jin/真值对比/jin_L1_Rbirth_L2_Rg.py 完全相同, 只是数据和系数换成 TNG50 的):
  1. 数据: data/rlcd_input/<version>/<ID>_e.csv (6 列 logJZ, JPhi, ecc, age, FeH, MgH), 每列标准化;
     真值: 同目录 <ID>_truth.csv 的 R_bir (出生半径) 和 Rg (导引半径), 与 _e.csv 逐行对应。
  2. 边系数: 作者的因果图冻结 (论文 4.2、4.3 节的做法), 系数取
     results/rlcd_A/identifiability/<version>_<ID>_ci.csv 里 NLL 最小的解 (最大似然解, 与 TNG/因果图 里画的相同)。
  3. 求每颗星的 L2、L1 (先 L2 后 L1), 分层最小二乘:
         r_c = c − Σ_{p≠L} a_pc · p          (子节点 c 扣掉其他父节点 p 的贡献)
         L   = Σ_c a_Lc · r_c / Σ_c a_Lc²    (使 Σ_c (r_c − a_Lc·L)² 最小的 L)
       c: L 的子节点 (标准化后的值); p: c 的其他父节点; a_pc: p→c 的系数; a_Lc: L→c 的系数
  4. 线性换算: Rbirth ≈ α·L1 + β, Rg ≈ α'·L2 + β' (最小二乘, 用到真值, 只是换单位), 标注皮尔逊相关系数 r。
  5. 每个面板另外标注 (REPORT_RLCD_TNG50.md 第 5.6 节):
       range: 与最优统计上分不开的所有解 (2·N·ΔNLL < 3.84) 给出的 |r| 范围 —— 边系数不可识别, 画出的只是其中一个;
       age only / [Fe/H] only: 单用年龄或 [Fe/H] 与真 Rbirth 的 |r| (参照线);
       L1 resembles: L1 与 [Fe/H]、年龄的 |r|, 看 L1 实际上"像"哪个观测量。
  怎么看图: 很多 L1 面板是一条很窄的竖带, 这不是算错, 而是相关很弱的表现。线性换算 f = α·L + β 的斜率
  α = r·σ_R/σ_L (r: 相关系数, σ_R: 真值的标准差, σ_L: L 的标准差), 所以 f 的散布只有真值散布的 |r| 倍;
  r ≈ 0.1 时 f 全部挤在平均值 (约 7 kpc) 附近。
  已核对 L1 的计算 (2026-09-27): 与作者 notebook 的逐星 Nelder–Mead 最小化 χ² 结果相差 < 1e-7;
  作者自己的 g2.79e12 用同一段代码精确复现论文 (r = 0.640)。TNG50 的 L1 与出生半径相关弱, 是因为最大似然解让
  L1 主要等于 [Fe/H]、[Mg/H] (|r| ≈ 0.9), 而 TNG50 7–10 kpc 的星里 [Fe/H] 与出生半径几乎无关。
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # casuality/
IDENT = ROOT / 'results' / 'rlcd_A' / 'identifiability'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
LR95 = 3.84
INK = '#1f2937'


def coefficients(row):
    """<tag>_ci.csv 的一行 → {(父, 子): 系数}, 变量名去掉 X_"""
    out = {}
    for k, v in row.items():
        if '->' in str(k):
            a, b = k.split('->')
            out[(a.replace('X_', ''), b.replace('X_', ''))] = float(v)
    return out


def solve_latent(L, coef, val):
    num, den = 0.0, 0.0
    for (p, c), a in coef.items():
        if p != L:
            continue
        r = val[c].copy()
        for (q, c2), b in coef.items():
            if c2 == c and q != L:
                r = r - b * val[q]
        num, den = num + a * r, den + a * a
    return num / den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default='nonoise', choices=['nonoise', 'noise'])
    args = ap.parse_args()
    out = HERE.with_suffix('.png') if args.version == 'nonoise' else HERE.with_name(HERE.stem + '_noise.png')
    data = ROOT / 'data' / 'rlcd_input' / args.version

    fig, axs = plt.subplots(6, 4, figsize=(17, 25))
    lines = []
    for k, gid in enumerate(IDS):
        x = pd.read_csv(data / f'{gid}_e.csv')[COLS]
        t = pd.read_csv(data / f'{gid}_truth.csv')
        z = (x - x.mean()) / x.std()
        ci = pd.read_csv(IDENT / f'{args.version}_{gid}_ci.csv').sort_values('nll')
        coef = coefficients(ci.iloc[0])
        inn = ci[ci.lr < LR95]
        val = {c: z[c].to_numpy() for c in COLS}
        val['L2'] = solve_latent('L2', coef, val)
        val['L1'] = solve_latent('L1', coef, val)
        r_age = abs(np.corrcoef(x['age'], t['R_bir'])[0, 1])
        r_feh = abs(np.corrcoef(x['FeH'], t['R_bir'])[0, 1])
        row, c0 = divmod(k, 2)
        for j, (L, truth, lab, rng) in enumerate([('L1', 'R_bir', 'R_b', ('abs_r_L1_Rbir')),
                                                  ('L2', 'Rg', 'R_g', ('abs_r_L2_Rg'))]):
            ax = axs[row, 2 * c0 + j]
            y = t[truth].to_numpy()
            slope, icpt = np.polyfit(val[L], y, 1)
            f = slope * val[L] + icpt
            r = np.corrcoef(f, y)[0, 1]
            ax.hexbin(f, y, bins='log', extent=(0, 20, 0, 20), gridsize=60, cmap='Blues', mincnt=1)
            ax.plot([0, 20], [0, 20], color=INK, lw=1)
            ax.set_xlim(0, 20)
            ax.set_ylim(0, 20)
            ax.set_aspect('equal')
            ax.set_xlabel(f'$f(L_{L[1]})$ [kpc]', fontsize=9)
            ax.set_ylabel(f'${lab}$ (true) [kpc]', fontsize=9)
            ax.tick_params(labelsize=8)
            txt = (f'$r$ = {r:.2f}   ($f$ = {slope:+.2f} $L_{L[1]}$ {icpt:+.2f})\n'
                   f'range (indistinguishable fits): {inn[rng].min():.2f}–{inn[rng].max():.2f}')
            if L == 'L1':
                txt += (f'\nage only: {r_age:.2f}   [Fe/H] only: {r_feh:.2f}'
                        f'\nL1 resembles: |r(L1,[Fe/H])| {abs(np.corrcoef(val["L1"], x["FeH"])[0, 1]):.2f},'
                        f' |r(L1,age)| {abs(np.corrcoef(val["L1"], x["age"])[0, 1]):.2f}')
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=7.8, color=INK,
                    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.85))
            ax.set_title(f'{gid}  ({len(x):,} stars)  —  {"birth radius" if L == "L1" else "guiding radius"}',
                         fontsize=9.5, color=INK)
            lines.append(f'{gid} {truth} = {slope:+.3f}*{L} {icpt:+.3f}, r = {r:.3f}')
    fig.suptitle(f'TNG50 ({args.version}): author\'s causal graph frozen, maximum-likelihood edge coefficients per galaxy\n'
                 'f(L1) vs true birth radius, f(L2) vs true guiding radius.  Coefficients are not unique: '
                 '"range" = |r| over all statistically indistinguishable fits.\n'
                 'A narrow vertical band = weak correlation: the linear map f = a·L + b shrinks toward the mean '
                 '(spread of f = |r| × spread of the truth).', fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.968))
    fig.savefig(out, dpi=100)
    print('\n'.join(lines))
    print(f'已保存 {out}')


if __name__ == '__main__':
    main()
