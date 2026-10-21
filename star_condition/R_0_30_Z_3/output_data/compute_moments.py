# -*- coding: utf-8 -*-
"""
compute_moments.py —— 把每个星系的输入和真值压缩成"协方差矩阵" (选星条件 R_0_30_Z_3 的第 2 步)

为什么: 这个条件下每个星系有 37–70 万颗星, 数据共 6 GB。但可识别性分析需要的只是
  6 个输入量 (logJZ, JPhi, ecc, age, FeH, MgH) 和 2 个真值 (R_bir, Rg) 这 8 个量两两之间的协方差、平均值和星数:
  - 最大似然拟合只用 6 个输入量的相关矩阵 (论文 3.3 节; 仓库 parameter_by_trek 也只用协方差);
  - L1、L2 是 6 个标准化输入量的线性组合 L = w·z, 所以
        r(L, R) = cov(L, R) / [std(L)·std(R)] = (w·c) / sqrt(w·C·w · var(R))
    意思是: L 与真值 R 的相关系数, 可以直接由协方差算出来, 不需要逐颗星。
      w: L 对 6 个标准化输入量的权重;  c: 6 个标准化输入量与 R 的协方差 (6 个数);
      C: 6 个标准化输入量之间的协方差矩阵 (即相关矩阵);  var(R): R 的方差
  与逐星计算的结果完全相同 (已在 7<R<10 的旧数据上核对, 见 identifiability_moments.py)。

运行 (numpy、pandas):
    cd /Volumes/001HardDisk/MW/casuality
    ~/PYTHONPATH_01/venv/bin/python3 star_condition/R_0_30_Z_3/output_data/compute_moments.py
输出: output_data/moments_nonoise.json, moments_noise.json
      每个星系: N (星数), names (8 个量的名字), mean (8 个平均值), cov (8×8 协方差矩阵, 分母 N−1)
"""
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
COND = HERE.parents[1]
IDS = ['487742', '530330', '538905', '545437', '546870', '547293',
       '552581', '557721', '560421', '562029', '567607', '569251']
COLS = ['logJZ', 'JPhi', 'ecc', 'age', 'FeH', 'MgH']
TRUTH = ['R_bir', 'Rg']


def main():
    versions = sys.argv[1:] or ['nonoise', 'noise']
    for v in versions:
        out = {}
        for gid in IDS:
            e = pd.read_csv(COND / 'input_data' / v / f'{gid}_e.csv', usecols=COLS)[COLS]
            t = pd.read_csv(COND / 'true_value' / v / f'{gid}_truth.csv', usecols=TRUTH)[TRUTH]
            X = np.column_stack([e.to_numpy(float), t.to_numpy(float)])
            assert np.isfinite(X).all()
            out[gid] = dict(N=len(X), names=COLS + TRUTH, mean=X.mean(0).tolist(), cov=np.cov(X.T, ddof=1).tolist())
            print(v, gid, len(X), flush=True)
        json.dump(out, open(HERE.with_name(f'moments_{v}.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
