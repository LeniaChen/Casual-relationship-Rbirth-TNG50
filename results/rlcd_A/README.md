# results/rlcd_A —— 做法 A（冻结作者的因果图）的结果

完整的解释见项目根目录的 `REPORT_RLCD_TNG50.md`（第 4、5 节）。

## ⚠ 读这些结果前必须知道的问题

**问题 3（最根本）：冻结作者的因果图后，边系数不能由数据可靠地确定。**
1. 数学上不可识别：作者的图有 15 条边，6 个观测量有 15 个相关系数，但"系数 → 相关系数"的雅可比矩阵秩只有 14，存在一个方向（主要是 L1→age 与 L1→MgH 此消彼长），沿着它改系数，观测相关矩阵完全不变。这是图的性质，与数据无关。
2. 真正求最大似然时，最优解几乎总是把某个观测量的噪声方差压到 0（24 个 TNG50 数据集中 22 个，多数是 J_φ；g2.79e12 是 [Fe/H] 或 e），似然面上有很多 NLL 只差万分之一左右的局部最优。
3. 后果：L1 的含义不稳定。TNG50 多数星系里，与最优统计上分不开（2·N·ΔNLL < 3.84）的解给出的 |r(L1, Rbirth)| 能从接近 0 取到 0.4–0.5；g2.79e12 上，初值从 30 个加到 230 个，最优解的 |r(L1, Rb)| 从 0.25 变成 0.62。**即使取最有利的解，L1 与 Rbirth 的相关也基本没有超过"单用年龄"或"单用 [Fe/H]"**（nonoise 12 个星系无一例外；noise 有 2 个在最有利的一端略高：0.58 对 0.57、0.71 对 0.67）。

→ 因此 `g279/`、`nonoise/` 里用 `rlcd_A.py` 挑出的"一组"系数给出的 r(L1, Rbirth)，只能作为示例，不能作为"方法有效/无效"的证据。
→ 主结果见 `identifiability/summary_ci_<version>.txt/.png`（严格标准）。

**更正（2026-09-26）**：第一版汇总（`identifiability/summary_<version>.*`）把"NLL 与最优相差 < 0.01"的解都算作"一样好"。NLL 是按每颗星平均的，对上万颗星来说这个标准太宽松（相当于似然比统计量 < 0.02·N）。这些文件保留作对比，结论以 `summary_ci_*` 为准。

**问题 2：参数学习的挑选。** 仓库的 `parameter_by_trek` 算 NLL 时用行列式的绝对值，而且不检查隐含的噪声方差是否 ≥ 0，所以物理上不成立的解可能得到"更好"的 NLL（g2.79e12 上出现过 −7.03，低于理论下限 −1.275）。作者论文用的系数是仓库里 NLL **第二小**（1.586）的那组，不是最小（1.498）的，论文没有说明；按正确公式重算它的 NLL 是 1.593，而最优约 −1.272，11.9 万颗星的似然比统计量约 68 万——那组系数离最大似然解很远。

**问题 1：结构学习有随机性。** 同一组 alpha，不同的 PYTHONHASHSEED 会得到不同的因果图；作者的图 5 次里出现 1 次（seed = 1）。见 `results/rlcd_g279/`。

以上三个问题论文和仓库都没有报告。

## 目录
- `g279/`：用作者的数据（g2.79e12）跑 `rlcd_A.py` v2 的结果（按 RMSE 挑中的解 r(L1,Rb) = 0.63，但前 5 个解从 0.01 到 0.63 不等）
- `nonoise/538905/`：TNG50 538905 的 `rlcd_A.py` v2 试跑（r(L1,Rb) = 0.06）
- `identifiability/`：**选项 (a) 的主结果**
  - 第一轮（宽松）：`scripts/local/rlcd_identifiability.py` → `<tag>_summary.txt`（雅可比矩阵的秩、不可识别方向）、`<tag>_multistart.csv`（30 个初值的解）、`<tag>_profile.csv/.png`（固定 L1→age 的剖面）；`rlcd_identifiability_summary.py` → `summary_<version>.csv/.txt/.png`
  - 第二轮（严格）：`scripts/local/rlcd_identifiability_ci.py` → `<tag>_ci.csv/.txt`（230 个初值严格收敛后的全部解）、`summary_ci_<version>.csv/.txt/.png`
  - 参照：`scripts/local/linear_ceiling.py` → `linear_ceiling.csv/.txt`（6 个量线性回归的上限、单用年龄 / [Fe/H] 的相关）
- `|r(...)|` 只用于诊断和给出范围，**没有被用来挑选任何解**
