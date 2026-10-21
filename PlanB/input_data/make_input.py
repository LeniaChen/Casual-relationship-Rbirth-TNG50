# -*- coding: utf-8 -*-
"""
make_input.py —— 做法 B 的第 0 步: 准备输入数据和真值 (2026-09-27)

做法 B 先用与论文相同的主选星条件: 7 < R < 10 kpc, |z| < 0.3 kpc (其余: 束缚且轨道量齐全; 有出生半径且 R_bir < 20 kpc;
[Fe/H] > -0.5; J_z > 0), 即做法 A 主结果用的同一批星。这批星已经由 scripts/local/make_rlcd_input.py 选好,
放在 data/rlcd_input/<版本>/ 里; 本脚本只是把它们原样复制到 PlanB/ 下, 让做法 B 的所有文件都在 PlanB 里:
    data/rlcd_input/<版本>/<ID>_e.csv      → PlanB/input_data/<版本>/<ID>_e.csv     (RLCD 的 6 列输入)
    data/rlcd_input/<版本>/<ID>_truth.csv  → PlanB/true_value/<版本>/<ID>_truth.csv (真出生半径 R_bir、真导引半径 Rg 等)
两个文件逐行对应 (第 i 行是同一颗星)。版本: nonoise (不加观测误差), noise (加了观测误差)。

运行 (numpy、pandas):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 PlanB/input_data/make_input.py
输出: 上面两个文件夹里的 csv, 以及 PlanB/input_data/selection.txt (每个星系的星数)
"""
import shutil
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve()
PLANB = HERE.parents[1]
ROOT = PLANB.parent
SRC = ROOT / 'data' / 'rlcd_input'
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']


def main():
    lines = ['做法 B 选星条件 (与论文相同): 7 < R < 10 kpc, |z| < 0.3 kpc; 束缚且轨道量齐全; has_birth 且 R_bir < 20 kpc; '
             '[Fe/H] > -0.5; J_z > 0', f'来源: {SRC.relative_to(ROOT)}/<版本>/ (原样复制)', '',
             f"{'ID':>7} {'nonoise':>9} {'noise':>9}"]
    tot = 0
    for gid in IDS:
        n = {}
        for v in ['nonoise', 'noise']:
            for sub, name in [('input_data', f'{gid}_e.csv'), ('true_value', f'{gid}_truth.csv')]:
                d = PLANB / sub / v
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(SRC / v / name, d / name)
            e = pd.read_csv(PLANB / 'input_data' / v / f'{gid}_e.csv')
            t = pd.read_csv(PLANB / 'true_value' / v / f'{gid}_truth.csv', usecols=['R_bir'])
            assert len(e) == len(t), gid
            n[v] = len(e)
        tot += n['nonoise']
        lines.append(f"{gid:>7} {n['nonoise']:>9,} {n['noise']:>9,}")
        print(lines[-1], flush=True)
    lines += ['', f'合计 {tot:,} 颗 (12 个星系)']
    (PLANB / 'input_data' / 'selection.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
