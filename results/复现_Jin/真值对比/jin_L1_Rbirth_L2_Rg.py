# -*- coding: utf-8 -*-
"""
jin_L1_Rbirth_L2_Rg.py —— 复现 Jin et al. 论文图 3 的左右两栏: f(L1) vs 真 Rbirth, f(L2) vs 真 Rg (模拟星系 g2.79e12)

运行 (numpy、pandas、matplotlib; 主 venv 即可):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 "results/复现_Jin/真值对比/jin_L1_Rbirth_L2_Rg.py"
输出: 与本脚本同名、同目录的 jin_L1_Rbirth_L2_Rg.png

做法与作者 notebook (Causal-Galactic-Archaeology/L1L2calculation/GetL1_L2_279sim.ipynb) 完全相同:
  1. 数据: data/g2.79e12_multi_ele_action_Jin_SN.csv (119,236 颗星); MgH = MgFe + FeH; logjz = ln(JZ);
     每列减平均值、除以标准差 (标准化)。
  2. 边系数: 作者发布的参数学习结果
     rlcd/results_causal_vertical_heating_g279_e/result_..._param_by_trek_nll1.5861642734160428.dot
     (与 notebook 里写死的数字完全相同; 这就是论文图 1 上的数字)。
  3. 求每颗星的 L2、L1 (先 L2 后 L1): notebook 对每颗星用 Nelder-Mead 最小化 χ²,
     χ² 是二次函数, 所以最小值有精确解, 这里直接用精确解 (结果相同, 快几千倍):
         r_c = c − Σ_{p≠L} a_pc · p          (子节点 c 扣掉其他父节点 p 的贡献)
         L   = Σ_c a_Lc · r_c / Σ_c a_Lc²    (使 Σ_c (r_c − a_Lc·L)² 最小的 L)
       c: L 的子节点 (标准化后); p: c 的其他父节点; a_pc: p→c 的系数; a_Lc: L→c 的系数
  4. 线性换算 (notebook 的 fitpoints): Rbirth ≈ α·L1 + β, Rg ≈ α'·L2 + β' (最小二乘),
     f(L1) = α·L1 + β, f(L2) = α'·L2 + β'。
  5. 画图与 notebook 相同: hexbin (对数色标), 范围 0–20 kpc, 黑线为 1:1, 标注皮尔逊相关系数 r。
论文: Rb = −1.97 L1 + 9.31, r = 0.64;  Rg = −1.29 L2 + 9.52 (r ≈ 0.98)。

注意 (见 REPORT_RLCD_TNG50.md 第 4.2、4.3 节): 作者这组系数不是最大似然解 (NLL 1.593, 最优约 −1.272);
按论文自己说的最大似然方法重新拟合到收敛, 结果会不同。本脚本只负责"用作者的系数复现作者的图"。
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]                                   # casuality/
REPO = ROOT / 'Causal-Galactic-Archaeology'
DATA = REPO / 'data' / 'g2.79e12_multi_ele_action_Jin_SN.csv'
DOT = REPO / 'rlcd' / 'results_causal_vertical_heating_g279_e' / \
    'result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1_param_by_trek_nll1.5861642734160428.dot'
OUT = HERE.with_suffix('.png')

# 作者 notebook 里写死的系数 (用来核对 dot 文件)
NOTEBOOK = {('L2', 'JPhi'): -0.57, ('MgH', 'JPhi'): -0.082, ('age', 'JPhi'): -0.18,
            ('L2', 'ecc'): 0.085, ('MgH', 'ecc'): 0.098, ('age', 'ecc'): 0.4,
            ('L1', 'logJZ'): -0.099, ('age', 'logJZ'): 0.42, ('MgH', 'logJZ'): -0.13,
            ('L1', 'FeH'): 0.84, ('age', 'FeH'): 0.0032, ('MgH', 'FeH'): 0.13,
            ('L1', 'age'): -0.27, ('L1', 'MgH'): 0.83, ('L1', 'L2'): 0.43}


def read_coefficients(dot):
    """读 .dot 里每条边的系数, 返回 {(父, 子): 系数}, 变量名去掉 X_ 前缀"""
    coef = {}
    for a, b, lab in re.findall(r'(\w+)\s*->\s*(\w+)\s*\[[^\]]*label="?([-+\deE.]+)', open(dot).read()):
        coef[(a.replace('X_', ''), b.replace('X_', ''))] = float(lab)
    return coef


def solve_latent(L, coef, val):
    """分层最小二乘: 用 L 的所有子节点求 L (子节点里的其他父节点必须已知)"""
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
    coef = read_coefficients(DOT)
    assert coef == NOTEBOOK, f'dot 文件的系数与 notebook 不一致: {coef}'

    df = pd.read_csv(DATA)
    df['MgH'] = df['MgFe'] + df['FeH']
    df['logJZ'] = np.log(df['JZ'])
    z = (df - df.mean()) / df.std()
    val = {c: z[c].to_numpy() for c in ['JPhi', 'ecc', 'age', 'FeH', 'MgH', 'logJZ']}
    val['L2'] = solve_latent('L2', coef, val)          # 子节点: JPhi, ecc
    val['L1'] = solve_latent('L1', coef, val)          # 子节点: L2, age, MgH, FeH, logJZ

    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    lines = [f'g2.79e12: {len(df):,} stars']
    for ax, L, truth, lab in [(axs[0], 'L1', 'R_bir', 'R_b'), (axs[1], 'L2', 'Rg', 'R_g')]:
        y = df[truth].to_numpy()
        slope, icpt = np.polyfit(val[L], y, 1)
        f = slope * val[L] + icpt
        r = np.corrcoef(f, y)[0, 1]
        hb = ax.hexbin(f, y, bins='log', extent=(0, 20, 0, 20), gridsize=80, cmap='Blues', mincnt=1)
        ax.plot([0, 20], [0, 20], color='#1f2937', lw=1)
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 20)
        ax.set_aspect('equal')
        ax.set_xlabel(f'$f(L_{L[1]})$ [kpc]')
        ax.set_ylabel(f'${lab}$ (true) [kpc]')
        ax.text(0.04, 0.95, f'$r$ = {r:.3f}\n$f(L_{L[1]})$ = {slope:+.3f} $L_{L[1]}$ {icpt:+.3f}',
                transform=ax.transAxes, va='top', fontsize=11, color='#1f2937')
        fig.colorbar(hb, ax=ax, label='# of stars')
        lines.append(f'{truth} = {slope:+.3f}*{L} {icpt:+.3f}, Pearson r = {r:.3f}')
    axs[0].set_title('birth radius vs $L_1$ (paper Fig. 3 left: r = 0.64)', fontsize=11)
    axs[1].set_title('guiding radius vs $L_2$ (paper Fig. 3 right)', fontsize=11)
    fig.suptitle("Reproduction of Jin et al. Fig. 3 — g2.79e12, author's released edge coefficients (NLL 1.586)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT, dpi=130)
    print('\n'.join(lines))
    print(f'已保存 {OUT}')


if __name__ == '__main__':
    main()
