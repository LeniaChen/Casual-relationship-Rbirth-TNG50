# 选星条件 R_0_30_Z_3（2026-09-27）

0 < R < 30 kpc，|z| < 3 kpc；其余条件与主流程相同（束缚且轨道量齐全、有出生半径且 < 20 kpc、[Fe/H] > −0.5、J_z > 0）。
每个星系 37.1–69.9 万颗星，共 602 万颗（从 R_0_30_Z_30 中筛出 |z| < 3 kpc 的星，见 input_data/make_input.py）。

| 文件夹 | 内容 | 怎么生成 |
|---|---|---|
| `input_data/` | `nonoise/`、`noise/<ID>_e.csv`（RLCD 的 6 列输入）；`selection.txt`（选星统计） | `input_data/make_input.py`（从 R_0_30_Z_30 筛选） |
| `true_value/` | `nonoise/`、`noise/<ID>_truth.csv`（出生半径 R_bir、导引半径 Rg 等真值） | 同上 |
| `output_data/` | `moments_<版本>.json`（每个星系 8 个量的协方差矩阵）；`<版本>_<ID>_ci.csv/_summary.txt`（最大似然解和统计上分不开的解）；`summary_ci_<版本>.csv/.txt/.png`（12 个星系汇总） | `compute_moments.py` → `identifiability_moments.py` |
| `因果图/` | 作者的因果图 + 每个星系的最大似然系数 | `causal_graph_R_0_30_Z_3.py` |
| `真值对比/` | f(L1)–出生半径、f(L2)–导引半径 | `L1_Rbirth_L2_Rg_R_0_30_Z_3.py` |

运行顺序（都在 casuality 目录下，用 `~/PYTHONPATH_01/venv/bin/python3`）：make_input.py → compute_moments.py → identifiability_moments.py --version nonoise（和 noise）→ 两个画图脚本（加 `--version noise` 得到 `_noise.png`）。

⚠ 这个范围包含整个盘，现在的角动量 J_φ（≈ 现在的半径）本身就和出生半径强相关（多数星系 0.78–0.90），所以 L1 要和"单用 J_φ"比。论文只用 7–10 kpc 的星，就是为了去掉这部分信息。
