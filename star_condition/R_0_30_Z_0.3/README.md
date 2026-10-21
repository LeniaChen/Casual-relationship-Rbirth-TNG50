# 选星条件 R_0_30_Z_0.3（2026-09-27）

0 < R < 30 kpc，|z| < 0.3 kpc；其余条件与主流程相同（束缚且轨道量齐全、有出生半径且 < 20 kpc、[Fe/H] > −0.5、J_z > 0）。
每个星系 18.8–31.3 万颗星，共 279 万颗（从 R_0_30_Z_30 中筛出 |z| < 0.3 kpc 的星，见 input_data/make_input.py）。
|z| 上限与论文相同（论文用 7 < R < 10 kpc、|z| < 0.3 kpc），只是 R 放宽到整个盘。

| 文件夹 | 内容 | 怎么生成 |
|---|---|---|
| `input_data/` | `nonoise/`、`noise/<ID>_e.csv`（RLCD 的 6 列输入）；`selection.txt`（选星统计） | `input_data/make_input.py`（从 R_0_30_Z_30 筛选） |
| `true_value/` | `nonoise/`、`noise/<ID>_truth.csv`（出生半径 R_bir、导引半径 Rg 等真值） | 同上 |
| `output_data/` | `moments_<版本>.json`（每个星系 8 个量的协方差矩阵）；`<版本>_<ID>_ci.csv/_summary.txt`（最大似然解和统计上分不开的解）；`summary_ci_<版本>.csv/.txt/.png`（12 个星系汇总） | `compute_moments.py` → `identifiability_moments.py` |
| `因果图/` | 作者的因果图 + 每个星系的最大似然系数 | `causal_graph_R_0_30_Z_0.3.py` |
| `真值对比/` | f(L1)–出生半径、f(L2)–导引半径 | `L1_Rbirth_L2_Rg_R_0_30_Z_0.3.py` |

运行顺序（都在 casuality 目录下，用 `~/PYTHONPATH_01/venv/bin/python3`）：make_input.py → compute_moments.py → identifiability_moments.py --version nonoise（和 noise；建议设 OMP_NUM_THREADS=1，两个一起跑时快一倍多）→ 两个画图脚本（加 `--version noise` 得到 `_noise.png`）。

## 结果要点（nonoise）
- 最大似然解的 |r(L1, 出生半径)|：9 个星系 0.79–0.89，552581 为 0.45，569251 为 0.28，546870 为 0.17。
- 没有一个星系的最大似然解超过"单用 J_φ"（J_φ 单独与出生半径 0.71–0.92）；只有 538905、546870 在"统计上分不开的解"的上端略超。
- 统计上分不开的解仍然很多（2–79 个），L1 的 |r| 范围在 487742（0.22–0.79）、560421（0.19–0.80）、546870（0.06–0.63）等星系很宽 → 边系数仍不唯一（雅可比秩 14/15）。
- 与 |z|<3 相比：557721（0.58→0.89）、569251（0.05→0.28）、560421（0.69→0.79）的最优解变了，说明最大似然解对选星很敏感。
- noise 版本更不稳定：545437、547293、557721、560421 的最优解 |r| 掉到 0.10–0.41，但统计上分不开的解里又有 0.77–0.90 的。

⚠ 这个范围包含整个盘，现在的角动量 J_φ（≈ 现在的半径）本身就和出生半径强相关，所以 L1 要和"单用 J_φ"比。
