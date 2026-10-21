# -*- coding: utf-8 -*-
"""
make_input.py —— 选星条件 R_0_30_Z_3 (0 < R < 30 kpc, |z| < 3 kpc) 的 RLCD 输入和真值 (2026-09-27, 用户指定)

不重新从模拟数据选星, 而是直接从 star_condition/R_0_30_Z_30/ 的结果里筛出 |z| < 3 kpc 的星:
  - R_0_30_Z_30 的其余条件与这里完全相同 (0 < R < 30 kpc; 束缚且轨道量齐全; 有出生半径且 R_bir < 20 kpc;
    [Fe/H] > -0.5; J_z > 0; 不切 [Mg/Fe]), 唯一的区别就是 |z| 的上限 30 → 3 kpc;
  - true_value/<版本>/<ID>_truth.csv 里有每颗星现在的高度 z, 且与 input_data/<版本>/<ID>_e.csv 逐行对应;
  - 加误差 (noise) 版本直接沿用 R_0_30_Z_30 里已经加好的误差, 所以同一颗星在两个条件下的数值完全相同。
  用 make_rlcd_input.py 从头重跑也会得到同一批星, 只是 noise 版本抽到的随机误差不同 (随机数顺序不同)。

运行 (numpy、pandas):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_3/input_data/make_input.py            # 12 个星系
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_3/input_data/make_input.py 538905     # 只做指定星系
输出: input_data/<版本>/<ID>_e.csv, true_value/<版本>/<ID>_truth.csv (版本 = nonoise, noise), input_data/selection.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
COND = HERE.parents[1]                                   # star_condition/R_0_30_Z_3/
SRC = COND.parent / 'R_0_30_Z_30'
Z_MAX = 3.0
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    ids = sys.argv[1:] or IDS
    rows = {}
    for gid in ids:
        for v in ['nonoise', 'noise']:
            e = pd.read_csv(SRC / 'input_data' / v / f'{gid}_e.csv')
            t = pd.read_csv(SRC / 'true_value' / v / f'{gid}_truth.csv')
            assert len(e) == len(t)
            m = (np.abs(t['z'].to_numpy()) < Z_MAX)
            for sub, df in [('input_data', e), ('true_value', t)]:
                d = COND / sub / v
                d.mkdir(parents=True, exist_ok=True)
                df[m].to_csv(d / (f'{gid}_e.csv' if sub == 'input_data' else f'{gid}_truth.csv'), index=False)
            rows[gid] = (len(t), int(m.sum()))
        print(f'{gid}: R_0_30_Z_30 {rows[gid][0]:,} 颗 → |z|<{Z_MAX} 后 {rows[gid][1]:,} 颗', flush=True)
    rep = COND / 'input_data' / 'selection.txt'
    old = {}
    if rep.exists() and len(sys.argv) > 1:                       # 分批运行时合并
        for l in rep.read_text().splitlines():
            p = l.split()
            if p and p[0].isdigit():
                old[p[0]] = (int(p[1].replace(',', '')), int(p[2].replace(',', '')))
    old.update(rows)
    lines = [f'选星条件: 0 < R < 30 kpc, |z| < {Z_MAX} kpc; 其余同 R_0_30_Z_30 (束缚且轨道量齐全; has_birth 且 R_bir < 20; '
             '[Fe/H] > -0.5; J_z > 0; 不切 [Mg/Fe]); 从 R_0_30_Z_30 的文件中筛出',
             '', f"{'ID':>7} {'R_0_30_Z_30':>12} {'|z|<3':>10}"]
    for g in sorted(old):
        lines.append(f'{g:>7} {old[g][0]:>12,} {old[g][1]:>10,}')
    tot = [v[1] for v in old.values()]
    lines += ['', f'最终样本: 每个星系 {min(tot):,}–{max(tot):,} 颗, 合计 {sum(tot):,} 颗 ({len(old)} 个星系)']
    rep.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
