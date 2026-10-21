# 做法 B：在每个 TNG50 星系上从头学因果图（2026-09-27）

**问题**：做法 A 把作者在 g2.79e12 上学到的图"冻结"后搬到 TNG50。做法 B 反过来问：如果让 RLCD 在每个 TNG50 星系上自己从头学，会学出什么图？图里的隐变量是不是出生半径、导引半径？

**选星**：与论文、做法 A 主结果相同：7 < R < 10 kpc，|z| < 0.3 kpc（每个星系 0.73–3.5 万颗，共 24 万颗）。脚本也支持 `--cond R_0_30_Z_30 / R_0_30_Z_3 / R_0_30_Z_0.3`（输出到 `output_data/<条件>/`，尚未运行）。

| 文件夹 | 内容 | 怎么生成 |
|---|---|---|
| `input_data/` | `nonoise/`、`noise/<ID>_e.csv`（RLCD 的 6 列输入）；`selection.txt` | `input_data/make_input.py`（从 data/rlcd_input 原样复制） |
| `true_value/` | `nonoise/`、`noise/<ID>_truth.csv`（真出生半径 R_bir、真导引半径 Rg 等） | 同上 |
| `output_data/` | `structures_<版本>.csv/.txt`：每个星系 30 个种子学到的图；`eval_<版本>.csv/.txt/_best.json`：每个主要结构的雅可比秩、最大似然系数、每个隐变量与真值的 \|r\|；`summary_planB_<版本>.csv/.txt`：12 个星系汇总并与做法 A 对比 | `rlcd_B.py` → `rlcd_B_eval.py` → `summary_planB.py` |
| `因果图/` | 每个星系最常见的图 + 最大似然系数 | `planB_causal_graph.py`（`--version noise`；`--rank 2` 画第二常见的图） |
| `真值对比/` | 每个星系最匹配的隐变量 f(L) vs 真出生半径 / 真导引半径 | `planB_latent_vs_truth.py` |

运行顺序（在 casuality 目录下）：
1. `~/PYTHONPATH_01/venv/bin/python3 PlanB/input_data/make_input.py`
2. `~/PYTHONPATH_01/venv_rlcd/bin/python3 PlanB/output_data/rlcd_B.py --version nonoise --seeds 30`（和 `noise`；每次 RLCD 约 2 秒，12 个星系 × 30 个种子约 12 分钟；已算过的会跳过）
3. `~/PYTHONPATH_01/venv/bin/python3 PlanB/output_data/rlcd_B_eval.py --version nonoise`（和 `noise`；约 10 分钟；两个一起跑时设 `OMP_NUM_THREADS=1`）
4. `summary_planB.py`、两个画图脚本（都加 `--version noise` 得到加误差版本）

做法细节：
- 每个种子单独开一个 Python 进程并设 `PYTHONHASHSEED`（RLCD 的结果依赖集合遍历顺序）；RLCD 设置与作者相同（GES、sparsity 2、alpha 1e-5/1e-15/1e-5/1e-15）。隐变量编号按"字典序最小"统一，所以 L1、L2 的名字是任意的。
- 评估：每个星系出现频率 ≥ 20% 的结构，50 个随机初值求最大似然（噪声方差 ≥ 0），保留与最优统计上分不开的解（2·N·ΔNLL < 3.84；N 星数，ΔNLL 每颗星负对数似然比最优多多少）。隐变量用作者 notebook 的分层最小二乘从子节点反推。
- 核对：把作者的图放进 `rlcd_B_eval.py`，在 538905 上得到 NLL 0.4918、雅可比秩 14/15、L1/L2 的 \|r\| 范围落在做法 A 的范围内（与做法 A 一致）。
- Mac 上重跑 `summary_planB.py` 得到的 `summary_planB_nonoise.csv` 与云端只有末位浮点差别（.txt 完全相同）；其余文件 md5 一致。

## 结果（nonoise）
- **与作者的图完全相同：0/360 次**（noise 也是 0/360）。8 个星系 30 个种子全部学到同一张图，4 个有 2 种。
- 学到的图五花八门：隐变量 0–4 个，有向边 8–17 条。545437、569251 没有隐变量且有无向边，求不了系数。4 个星系违反论文先验（年龄有观测量父节点）。
- 10 个能求系数的星系里，只有 3 个的主要结构可识别（雅可比秩 = 边数）：546870、562029、567607，都只有 1 个隐变量；其余 7 个仍不可识别（边数多于能确定的量）。
- 用真值挑"最像出生半径"的隐变量（最好情况），\|r\| 中位数 0.51（0.24–0.77）；超过"单用年龄"的只有 530330（0.71 vs 0.68）、552581（0.53 vs 0.23）、562029（0.77 vs 0.69）。导引半径通常更好（多数 0.6–0.9）。
- 按 BIC 扣掉边数差别后（dBIC = (NLL_B − NLL_A) + (边数_B − 15)·ln N/(2N)，正数表示作者的图更好），**作者的图在每个能比较的星系都比 RLCD 在该星系自己学到的图更贴合数据**（dBIC 每颗星 +0.003 到 +0.42）——说明 RLCD 的结构学习（基于秩检验）找到的不是最大似然意义上的最好图。
- noise：结构更简单（多为 1–2 个隐变量），最匹配隐变量的 \|r(出生半径)\| 中位数只有 0.27；dBIC 除 487742（−0.007）外都是正的。

**结论**：在 TNG50 上，RLCD 不会自己学出"L1 = 出生半径、L2 = 导引半径"的图；学到的隐变量需要用真值才能认出哪个像出生半径，而即使这样也大多不比单用年龄好。
