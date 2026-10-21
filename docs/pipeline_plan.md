# TNG50 → Causal Galactic Archaeology (RLCD) 验证项目：完整流程与进度

> **这是本项目的主文档。** 新开一个对话时，先让 Claude 完整读这份文件（本地路径 `/Volumes/001HardDisk/MW/casuality/docs/pipeline_plan.md`，claude.ai 项目文档 `claude/pipeline_plan.md` 是同一内容的副本），就能接着当前进度继续。
> 最后更新：2026-09-27（REPORT 已补入做法 A 的全部结果并更新结论；做法 B 只在本文件记录，不进 REPORT）
> **给人看的完整汇报在项目根目录 `REPORT_RLCD_TNG50.md`**（项目是什么、每一步怎么做、问题、结果、结论、方法原理）；本文件偏重复现细节。
> **2026-10-21 起**：项目已打包成 GitHub 仓库（本机 `~/Project2026/Casual-relationship-Rbirth-TNG50`），以后代码和文档以仓库为准，新对话先读仓库根目录的 `STATUS.md` 和 `CLAUDE.md`。本文件里的相对路径都相对于“项目根目录”：原来是 `/Volumes/001HardDisk/MW/casuality`，现在是仓库根目录（命令里的 `cd /Volumes/001HardDisk/MW/casuality` 换成 `cd` 到仓库根目录即可）。大数据（约 17 GB）仍只在外置硬盘的 casuality 里，用 `bash tools/link_data.sh` 以软链接接回仓库，清单见 `DATA.md`。

---

## ⚠ 已发现的问题（论文和仓库都没有报告；用任何结果前先读这里）
1. **结构学习有随机性**：同一组 alpha，不同 PYTHONHASHSEED 会得到不同的因果图；作者的图 5 次里只出现 1 次（seed = 1）。来源：变量名集合的遍历顺序（受字符串哈希随机化影响），以及 `LatentGroups.py` 里没设种子的 `random.sample`。→ 第 6.3 节，`results/rlcd_g279/`
2. **参数学习的挑选有问题**：仓库 `parameter_by_trek` 用行列式的绝对值算 NLL，也不检查隐含的噪声方差是否 ≥ 0，所以物理上不成立的解可能得到"更好"的 NLL。作者论文用的系数是 NLL 第二小（1.586）的那组，不是最小（1.498）的，论文没有说明；那组系数也没有收敛到最优（最优约 −1.27）。→ 第 6.6 节
3. **（最根本）冻结作者的因果图后，边系数不能由数据可靠地确定**：
   - 雅可比矩阵秩 14 < 15（15 条边、15 个观测相关系数），存在拟合完全相同的一族系数（主要是 L1→age 与 L1→MgH 此消彼长）。图的性质，与数据无关。
   - 真正求最大似然（噪声方差 ≥ 0、严格收敛、230 个初值）时，最优解几乎总把某个观测量的噪声方差压到 0（24 个 TNG50 数据集中 22 个，多数是 J_φ；g2.79e12 是 [Fe/H] 或 e；Heywood 情形），似然面上有很多 NLL 只差 ~1e-4 的局部最优。
   - 后果：TNG50 多数星系里，与最优统计上分不开（2·N·ΔNLL < 3.84）的解给出的 |r(L1, Rbirth)| 能从接近 0 取到 0.4–0.5；g2.79e12 上初值从 30 个加到 230 个，最优解的 |r(L1,Rb)| 从 0.25 变成 0.62。**即使取最有利的解，L1 与 Rbirth 的相关也基本没超过单用年龄或单用 [Fe/H]**（nonoise 12/12；noise 有 2 个略高：0.58 对 0.57、0.71 对 0.67）。
   - 论文用的系数（NLL 按正确公式为 1.593）离最优（约 −1.272）很远：11.9 万颗星的似然比统计量约 68 万。
   - **更正**：第一版汇总（`summary_<version>.*`）用"ΔNLL < 0.01"算"一样好"，太宽松（相当于似然比 < 0.02·N）；以严格版 `summary_ci_<version>.*` 为准。
   → 第 6.6 节末尾，`results/rlcd_A/identifiability/`，`results/rlcd_A/README.md`，`REPORT_RLCD_TNG50.md` 第 4.3、5.6 节

**用户决定（2026-09-26）：选项 (a)**——把不可识别本身作为主要结果，对 12 个星系给出 |r(L1, Rbirth)| 和 |r(L2, Rg)| 在同样好的拟合里的可能范围，并在所有结果和文档里**明确标注这些问题**。

## 0. 给新对话的 Claude：开工前先做这几件事

1. **读完本文件**，特别是第 2 节"当前进度"和第 9 节"已知的坑"。
2. **申请文件夹访问**：`/Volumes/001HardDisk/MW/casuality`（本项目全部文件）。需要时再申请 `~/PYTHONPATH_01`（venv、AnastrisTNG、Agama 源码）、`/Volumes/001HardDisk_6/TNG_server/423_sample`（旧脚本和旧结果）、`~/Downloads`。
3. **确认最新状态**：
   - 看 `results/compute_orbits.log` 的末尾（先用 `tr '\r' '\n'` 把进度条的回车换成换行，再过滤掉"轨道积分"那些行），以及 `data/tng50_derived/` 里有哪些文件
   - 看 `results/` 里有没有新的日志或报告
4. **遵守工作约定**（第 1 节）。
5. **每次回复末尾都写"当前进度"和"下一步"**，并在推进后更新本文件的第 2 节。本地文件和项目文档要同时更新。

---

## 1. 工作约定（用户的要求）

- 用户用中文写就用中文回复。用户不熟悉数学符号：**每个公式都要配白话解释，并逐个解释每个符号的含义**。
- **每次回复末尾都要告诉用户：项目进行到哪里了、接下来做什么。**
- **本项目所有脚本、数据、结果都放在 `/Volumes/001HardDisk/MW/casuality`**，不要放到别处。
- **科研统一用 venv**：`~/PYTHONPATH_01/venv`（Python 3.11.7）。新包都装在这里，不要装进 anaconda base。用户终端默认处在 anaconda base，所以命令里要写 venv 的完整路径：`~/PYTHONPATH_01/venv/bin/python3`。
- 用户偏好把数据下载到本地计算，而不是在 TNG 服务器上算。
- 服务器端的提取在 **TNG JupyterLab** 上运行，那里的 AnastrisTNG 是**旧版**，不要按 2.0.0 的接口写服务器脚本。
- **Claude 不能直接在用户的 Mac 上运行 Python/AGAMA**（原因见第 9 节）。需要运行时，给用户完整的可复制命令，让他在终端里运行；输出存成日志放到 `results/`，Claude 通过读日志了解结果。

---

## 2. 当前进度（每次推进后更新）

### 第 1 阶段：准备 TNG50 数据（6 个输入量 + 2 个真值）
- [x] 1.1 读论文和仓库，确定输入输出格式（第 3 节）
- [x] 1.2 服务器提取 12 个星系：`scripts/server/data_obtain.py`，再用 `merge_hdf5.py` 合并成一个文件 → 下载为 `data/tng50_extract/data_obtain_all.hdf5`（0.96 GB）
- [x] 1.3 与旧的 `001_orbit` 结果逐颗恒星核对，一致（`scripts/local/check_vs_001_orbit.py`，报告 `results/check_vs_001_orbit.txt`）
- [x] 1.4a 在 venv 里装好 AGAMA 1.0（用 Homebrew libomp 做并行计算，安装方法见第 8 节）
- [x] 1.4b `compute_orbits.py --test` 测试通过（487742 的 5000 颗星，各项检查都通过）
- [x] 1.4c `compute_orbits.py` 跑完全部 12 个星系（2026-09-26，每个星系 3.8–8.0 分钟，都没有报错，日志 `results/compute_orbits.log`）→ `data/tng50_derived/orbits_all.hdf5`
  - 第一次运行时星系 530330 报错（ActionFinder 初始化失败），已修复：势的网格内边界改为 rmin = 0.2 kpc，版本号 `v2-rmin-grid`，详见第 6.2 节
  - 重跑时 487742 已按新设置重算完成，接着在算 530330
- [x] 1.5 检查 12 个星系的轨道结果：全部通过（`scripts/local/check_orbits.py` → `results/orbit_check/`，见第 6.4 节）
- [x] 1.6 选星并生成 RLCD 输入 csv（`scripts/local/make_rlcd_input.py` → `data/rlcd_input/{nonoise,noise}/<ID>_e.csv` 和 `<ID>_truth.csv`，见第 6.5 节）。**第 1 阶段完成（2026-09-26）**

### 第 2 阶段：运行 RLCD 并与真值对比
- [x] 2.1 安装 Causal-Galactic-Archaeology：venv `~/PYTHONPATH_01/venv_rlcd`（Python 3.12，来自 Homebrew），依赖清单 `scripts/local/rlcd_requirements_mac.txt`，仓库在 `casuality/Causal-Galactic-Archaeology/`，步骤见第 8.3 节
- [x] 2.2 用仓库自带的 g2.79e12 数据复现论文结果：跑了 4 组 alpha，对比见第 6.3 节和 `results/rlcd_g279/compare_4alpha.txt`。结论：安装正确；alpha0 = 1e-5 时能找到论文的隐变量结构（作者的 15 条边全部都在，另外多 3 条）；alpha0 更严时会变成 4 个隐变量。没有一组和作者的图完全相同
- [x] 2.2a 查明和作者结果不同的原因：**RLCD 每次运行的结果是随机的**。同一组 alpha、5 个不同的 PYTHONHASHSEED，只有 seed=1 和作者的图完全相同（见第 6.3 节）
- [x] 2.2b 用户决定（2026-09-26）：**做法 A 为主结果，做法 B 为补充**
  - A：冻结作者 g2.79 的因果图，只在每个 TNG50 星系上重新估计边系数（论文 4.2、4.3 节的做法）
  - B：在每个 TNG50 星系上用多个种子从头学因果图，统计各种结构出现的频率
- [x] **2.3A 做法 A，选项 (a)——12 个星系的可识别性分析（2026-09-26 完成）**，结果在 `results/rlcd_A/identifiability/`：
  - 第一轮（宽松）：`rlcd_identifiability.py`（每个星系：雅可比秩、30 个初值、L1→age 剖面）+ `rlcd_identifiability_summary.py` → `summary_{nonoise,noise}.*`
  - 线性上限：`linear_ceiling.py` → `linear_ceiling.*`
  - 第二轮（严格，主结果）：`rlcd_identifiability_ci.py --extra 200`（230 个初值、噪声方差 ≥ 0、严格收敛、似然比 95%）→ `<tag>_ci.*`、`summary_ci_{nonoise,noise}.*`
  - 这些都是在 Claude 云端环境里跑的（2 核）；在 Mac 上复现的命令见 `REPORT_RLCD_TNG50.md` 第 5.7 节
- [x] 2.4（做法 A 部分）L1 与 Rbirth、L2 与 Rg 对比 → 见上；**结论：按论文做法不可行**（`REPORT_RLCD_TNG50.md` 第 6 节）
- [x] 汇报文档 `REPORT_RLCD_TNG50.md`（项目根目录）。2026-09-27 更新：§0、§3.7（30 个种子：作者的图 7/30）、§3.11 出图脚本、§3.12 放宽选星（协方差方法）、§4.1、§4.3、§5.8（f(L1) 窄竖条的原因 σ_f = |r|·σ_R）、§5.9（三种宽选星的表和解读）、§6.1 第 4、6 条、§7、§8.4.4（只用协方差）。**用户决定：做法 B 结果不好，不写进 REPORT，也不在 REPORT 里专门提及，只在本文件正常记录**（已删去 REPORT 里原来提到做法 B 的几处）
- [x] 出图脚本（2026-09-27，用户要求：每个文件夹只放脚本，脚本输出同名 png 在同一目录）：
  - `results/复现_Jin/因果图/jin_causal_graph.py`（venv_rlcd）：在内存里重跑 RLCD 结构学习（自动以 PYTHONHASHSEED=1 重启自己，可用 `--seed` 换），与作者的图逐边比较，边上的数字用作者发布的系数（nll1.586），复现论文图 1。结果：与作者的图完全相同。不写入仓库的 results_* 目录
  - `results/复现_Jin/真值对比/jin_L1_Rbirth_L2_Rg.py`（venv）：作者的系数 + 分层最小二乘求 L1、L2 + 线性换算，复现论文图 3 左右两栏：Rb = −1.969 L1 + 9.309，r = 0.640；Rg = −1.289 L2 + 9.518，r = 0.976
  - `results/TNG/因果图/tng_causal_graph.py`（venv）：12 个星系，作者的图 + 每个星系的最大似然系数（`<version>_<ID>_ci.csv` 中 NLL 最小的一行），标出噪声方差≈0 的量；`--version noise` 输出 `_noise.png`
  - （2026-09-27 核对）TNG 的 L1 计算无误：与 notebook 式逐星 Nelder–Mead 相差 < 1e-7；用模型最优估计 E[L1|X] 得到的 |r| 也同样低。L1 面板是窄竖带是因为相关弱（f 的散布 = |r| × 真值散布）；最大似然解让 L1 ≈ [Fe/H]、[Mg/H]（|r|≈0.9），而 TNG50 里它们与出生半径几乎无关
  - （2026-09-27）用户想放宽选星再重跑。先用 `scripts/local/preview_selection.py` 预览了 5 种选星（输出 `results/selection_preview.txt`）：放宽到 4<R<12、不切 |z| 后，出生半径散布中位数只从 2.21 增到 2.46 kpc（作者 g2.79 是 3.67），[Fe/H] 与出生半径的 |r| 中位数仍约 0.15–0.17（作者 0.70），线性上限几乎不变（0.53–0.88）。结论：放宽选星能让 f(L1) 略宽，但改变不了"TNG 里化学量几乎不带出生半径信息"这一点。选星常数在 `make_rlcd_input.py` 第 47–50 行（R_MIN、R_MAX、Z_MAX、FEH_MIN），重跑时建议把 OUT_DIR 改成新目录以保留旧结果
  - `results/TNG/真值对比/tng_L1_Rbirth_L2_Rg.py`（venv）：同样的系数，12 个星系的 f(L1)–Rbirth、f(L2)–Rg，面板里标注统计上分不开的解的 |r| 范围和单用年龄/[Fe/H] 的 |r|
- [ ] 2.3B 做法 B（补充）：每个星系用多个 PYTHONHASHSEED 从头学结构，统计结构频率

- [x] **选星条件 R_0_30_Z_30（2026-09-27，用户要求）**：0<R<30 kpc、|z|<30 kpc，其余条件不变。全部文件在 `star_condition/R_0_30_Z_30/`（README.md 说明各子文件夹）。
  - 数据共 6 GB，太大无法整体传到云端，所以先在 Mac（device_bash 的虚拟机里跑，pyfive+numpy）上用 `output_data/compute_moments.py` 把每个星系压成 8×8 协方差矩阵，云端 `output_data/identifiability_moments.py` 只用协方差做同样的分析。L1、L2 是输入量的线性组合，r 可以直接由协方差算出，已在旧数据上核对到小数点后 6 位完全一致。`真值对比` 画图在虚拟机里逐星跑（matplotlib 可用）
  - 结果（nonoise）：|r(L1,Rbirth)| 多数星系 0.66–0.86（统计上分不开的范围）；但单用 J_φ 就有 0.78–0.90（整个盘里现在的半径 ≈ 出生半径）。L1 的最有利值仅在 4 个星系略超过最好的单个输入量（530330、538905、546870、557721）；546870、560421、569251 范围仍很宽（下限 0.03–0.13）。noise 版本的范围更宽。线性上限 0.74–0.93
  - 解读：L1 看起来"变好"，主要来自整个盘上"现在在哪 ≈ 出生在哪"这一平凡信息；论文只用 7–10 kpc 正是为了排除它
  - 注意：我之前的选星预览只放宽到 4<R<12，结论"帮助不大"不适用于全盘范围

- [x] **选星条件 R_0_30_Z_3（2026-09-27）**：从 R_0_30_Z_30 的文件里筛出 |z|<3 kpc（`input_data/make_input.py`，只少了约 3% 的星，noise 版本沿用同一批误差），其余流程与 R_0_30_Z_30 相同。结果与 R_0_30_Z_30 基本一致；但 546870（最优解 |r(L1,Rb)| 0.71 → 0.09）、569251（0.46 → 0.05）、557721（0.82 → 0.58）的最优解变了——同样说明最大似然解不稳定
- [x] **选星条件 R_0_30_Z_0.3（2026-09-27）**：从 R_0_30_Z_30 筛出 |z|<0.3 kpc（每个星系 18.8–31.3 万颗，共 279 万颗，约为 |z|<3 的一半），流程同上，全部文件在 `star_condition/R_0_30_Z_0.3/`（README.md 有结果要点）。nonoise 最大似然解的 |r(L1,Rb)|：9 个星系 0.79–0.89，552581 0.45、569251 0.28、546870 0.17；没有一个星系的最优解超过"单用 J_φ"（0.71–0.92）。与 |z|<3 相比 557721（0.58→0.89）、569251（0.05→0.28）的最优解又变了；noise 版本 545437/547293/557721/560421 的最优解掉到 0.10–0.41。统计上分不开的解仍多（雅可比秩 14/15）。经验：两个进程一起跑时要设 OMP_NUM_THREADS=1（否则 2 核机器上多线程互相抢，慢 2 倍多）
- [x] **做法 B（2026-09-27，主条件 7<R<10、|z|<0.3，与论文相同）**：全部文件在 `PlanB/`（input_data、true_value、output_data、因果图、真值对比；README.md 有运行顺序和结果）。脚本：`PlanB/input_data/make_input.py`（从 data/rlcd_input 复制）→ `PlanB/output_data/rlcd_B.py --seeds 30`（每个星系 30 个 PYTHONHASHSEED，venv_rlcd，每次约 2 秒）→ `rlcd_B_eval.py`（出现频率 ≥20% 的结构：雅可比秩、50 初值最大似然、分层最小二乘求隐变量权重，另存 eval_<版本>_best.json）→ `summary_planB.py`（与做法 A 对比，含 dBIC）→ `PlanB/因果图/planB_causal_graph.py`、`PlanB/真值对比/planB_latent_vs_truth.py`。已用作者的图核对 rlcd_B_eval（538905: NLL 0.4918、秩 14/15，与做法 A 一致）。注意 rlcd_B_eval.py 必须先 import 同文件夹的 rlcd_B（否则会用到 scripts/local 里输出到旧路径的版本）。
  结果：与作者的图相同 0/360（nonoise、noise 都是）；8 个星系 30 个种子学到同一张图；隐变量 0–4 个、有向边 8–17 条；545437、569251 无隐变量且有无向边（无法求系数）；4 个星系违反年龄先验；能求系数的 10 个里只有 3 个（546870、562029、567607，都是 1 个隐变量）可识别。用真值挑"最像出生半径"的隐变量（最好情况）|r| 中位数 0.51（0.24–0.77），只有 530330、552581、562029 超过单用年龄；noise 中位数 0.27。按 BIC 扣掉边数后，作者的图在所有能比较的星系都比 RLCD 自己学到的图更贴合（dBIC/星 +0.003–+0.42；noise 除 487742 外）。结论：RLCD 在 TNG50 上不会自己学出"L1=出生半径、L2=导引半径"的图。
- [ ] 可选后续：做法 B 用宽选星条件（`rlcd_B.py --cond R_0_30_Z_0.3` 等，用协方差替身数据；输出到 `PlanB/output_data/<条件>/`，画图脚本目前只读主条件）；等用户决定

### 第 3 阶段（以后）
- [ ] 如果 TNG50 验证可行，应用到银河系观测数据。需要决定是找作者要数据，还是自己从头做。

**下一步（具体，2026-09-26 晚更新）：** 做法 A 已完成，汇报文档已写好。等用户看完报告后决定：(1) 做法 B（每个 TNG50 星系多种子从头学图，检查最常见的图是否可识别）；(2) 让模型可识别（加 J_R 或更多元素、或非线性）；(3) 是否联系作者。建议先做 (1)。

**（更早的）下一步：** 第 1 阶段已完成。按用户要求，第 1 阶段做完才进入第 2 阶段后面的部分。下一步是第 2 阶段 2.2b：和用户确认 TNG50 用做法 A 还是 B。

**（旧）下一步：** 等 `compute_orbits.py` 跑完 → Claude 读 `results/compute_orbits.log`，检查每个星系的 v_c(8 kpc)、盘内恒星的中位数、有没有报错 → 写结果检查脚本（1.5）→ 和用户讨论选星条件 → 写生成 RLCD 输入的脚本（1.6）。

---

## 3. 研究目标与论文要点

### 3.1 目标
用 TNG50 模拟检验 Jin et al. 的 RLCD（Rank-based Latent Causal Discovery）方法：
- **输入 6 个可观测量**：Age、[Fe/H]、[Mg/H]、J_z、J_φ、e
- **程序输出**：两个隐变量 L1、L2
- **和 TNG50 真值对比**：L1 ↔ 出生半径 Rbirth，L2 ↔ 导引半径 Rg

### 3.2 论文（项目文件 `Causal_Galactic_Archaeology_ApJ.pdf`，Jin, Lu, Ting, Zheng, Buck）
- 模拟数据：NIHAO-UHD 的 g2.79e12（主要用）和 g8.26e11（检验迁移性）。选盘星：[Fe/H] > −1，现在位于 7–10 kpc。给年龄加 10% 误差，[Fe/H] 加 0.02 dex，[Mg/Fe] 加 0.06 dex。低 α 盘用 [Mg/Fe] < −0.01 划分，只研究低 α 盘。
- 轨道参数：AGAMA，轴对称多极展开势（包含暗物质、气体、恒星）。
- 观测数据：APOKASC-3。作用量用 galpy Stäckel fudge + gala MilkyWayPotential2022，v_c(R⊙) = 229 km/s。
- 发现的因果图：L1 是根节点，指向 [Fe/H]、[Mg/H]、Age、J_z 和 L2；L2 指向 J_φ 和 e。L1 对应 Rbirth（g2.79 上皮尔逊相关系数 0.64），L2 对应 Rg（非常紧）。
- 从 L 换算到物理量：Rb = −1.97 L1 + 9.31，Rg = −1.29 L2 + 9.52（g2.79）。
- 先验约束：Age 在观测变量里没有父节点。

### 3.3 仓库 https://github.com/ZehaoJin/Causal-Galactic-Archaeology 的实际数据格式
- RLCD 真正读入的是 `data/g2.79e12_e.csv`，6 列：`logJZ, JPhi, ecc, age, FeH, MgH`
  - `logJZ` = log₁₀(J_z)
  - `JPhi` 单位 kpc·km/s，模拟里是正值（观测文件里是负值）
  - `MgH` = MgFe + FeH
- 真值在 `data/g2.79e12_multi_ele_action_Jin_SN.csv`，列为：`age, iord, JR, JZ, JPhi, FeH, MgFe, ecc, Rg, R_bir, z_bir, R`
- **实际数据和论文描述不一致的地方**：
  - 实际选星范围是 R = 9.4–12.1 kpc（不是 7–10 kpc）
  - [Fe/H] 实际 > −0.5（不是 −1）
  - R_bir 最大约 20 kpc，看起来额外切了 R_bir < 20
  - csv 里看不出加过观测误差
  - Rg ≈ JPhi / 347（g2.79），也就是 Rg 满足 L_z = Rg · v_c(Rg)
- L1/L2 的计算在 `L1L2calculation/*.ipynb`，里面的**边系数是写死的**，换数据后要改成新学到的系数。那里用的是自然对数 `np.log(JZ)`，但因为先做了标准化，结果和 log₁₀ 等价。
- 运行命令（在 `rlcd/` 目录下）：
  - 结构学习：`python run_rlcd.py --s <scenario> --alpha0 1e-5 --alpha1 1e-15 --alpha2 1e-5 --alpha3 1e-15`
  - 系数学习：`cd ParameterLearning; python linear_parameter_identification.py --s <scenario> --dot_path <.dot> --lr 0.005 --epochs 600 --lr_decay plateau`
  - 新数据要在 `rlcd/scenarios.py` 里加一个 scenario 函数（读 csv，列名前加 `X_`）
- 依赖：Python 3.12+、Java（Tetrad）、Graphviz、`rlcd/py312_requirements.txt`

---

## 4. 文件夹结构（`/Volumes/001HardDisk/MW/casuality`）

```
casuality/
├── README.md
├── REPORT_RLCD_TNG50.md                   ← 汇报文档（给人看：项目、步骤、问题、结果、结论、方法）
├── docs/pipeline_plan.md                  ← 本文件（复现细节）
├── scripts/
│   ├── server/                            在 TNG JupyterLab 上运行
│   │   ├── data_obtain.py                 提取 12 个星系的粒子数据 + 出生半径 → ./data_obtain/<ID>.hdf5
│   │   └── merge_hdf5.py                  合并成 data_obtain_all.hdf5
│   ├── local/                             在用户 Mac 上用 venv 运行
│   │   ├── check_vs_001_orbit.py          核对新数据与旧 001_orbit 结果
│   │   ├── compute_orbits.py              AGAMA：引力势、作用量、Rg、偏心率
│   │   ├── check_orbits.py                检查 12 个星系的轨道结果 → results/orbit_check/
│   │   ├── make_rlcd_input.py             选星 + 生成 RLCD 输入 → data/rlcd_input/
│   │   ├── compare_rlcd_graphs.py         比较因果图（隐变量编号对齐 + 先验检查）
│   │   ├── rlcd_A.py                      做法 A：冻结作者的图 → 参数学习 → L1/L2 → 与真值对比（venv_rlcd）⚠ 结果不可识别
│   │   ├── rlcd_identifiability.py        可识别性检查：雅可比矩阵的秩、多初值最大似然、L1->age 剖面（numpy/scipy）
│   │   ├── rlcd_identifiability_summary.py 汇总第一轮（宽松标准 ΔNLL<0.01）的范围
│   │   ├── rlcd_identifiability_ci.py     第二轮严格标准：polish + 噪声方差≥0 + 似然比 95%（主结果，--extra 200）
│   │   ├── linear_ceiling.py              线性上限和单变量基线（有监督，只作参照）
│   │   └── rlcd_requirements_mac.txt      venv_rlcd 的依赖清单
│   └── reference/                         旧参考脚本（只读副本）：001_orbit.py、001_orbit_new.py、
│                                          TNG50_disk_trajectory.py、obtain_data_s98.py
├── Causal-Galactic-Archaeology/          RLCD 仓库（git clone），在其 rlcd/ 下运行
├── data/
│   ├── tng50_extract/data_obtain_all.hdf5 服务器提取的原始数据（不要改动）
│   ├── tng50_derived/
│   │   ├── orbits_all.hdf5                compute_orbits.py 的输出
│   │   ├── orbits_test.hdf5               --test 模式的输出（旧版势设置，可以删）
│   │   └── potentials/<ID>_halo.coef, <ID>_disk.coef   引力势系数
│   └── rlcd_input/{nonoise,noise}/<ID>_e.csv, <ID>_truth.csv   RLCD 输入和真值
└── results/
    ├── check_vs_001_orbit.txt             核对报告
    ├── compute_orbits.log                 轨道计算日志
    ├── rlcd_g279/                         RLCD 复现测试（作者的图、我们各次运行的图、比较报告）
    ├── orbit_check/                       轨道结果检查（summary + 图）
    ├── rlcd_input_selection.txt           选星每一步剩下的恒星数
    └── rlcd_A/<version>/                  做法 A 的结果（g279 验证、nonoise、noise）
```

**这个文件夹以外的相关位置：**
- `~/PYTHONPATH_01/venv`：科研用的 Python 环境
- `~/PYTHONPATH_01/venv_rlcd`：只用来跑 RLCD 的 Python 3.12 环境
- `~/PYTHONPATH_01/AnastrisTNG`：AnastrisTNG 源码，已更新到 2.0.0，在 venv 里是可编辑安装
- `~/PYTHONPATH_01/Agama`：AGAMA 源码，从这里编译安装
- `~/PYTHONPATH_01/galpos`：galpos 包（服务器上另有一份）
- `/Volumes/001HardDisk_6/TNG_server/423_sample/`：旧项目（423 个星系样本），包括：
  - `TNG50_disk_trajectory.hdf5`（375 MB，每个星系在每个快照的中心、速度、角动量）
  - `001_orbit/<ID>.txt`（旧的出生半径结果）
- `~/Downloads/001_orbit_s99_3D_z/<ID>.txt`：2026 年 4 月版的出生半径结果（多了 z 和三维距离两列）
- 服务器上：`~/MW/casuality/`（运行 data_obtain.py 的目录），`~/423_sample/`（轨道文件原件），模拟数据在 `/home/tnguser/sims.TNG/TNG50-1/output`

---

## 5. 计算环境

### 5.1 用户的 Mac（Apple Silicon，arm64）
- anaconda base：`/opt/anaconda3`，终端默认环境，`python3` 指向它。里面的 llvm-openmp 是 14.0.6，太旧，编译 AGAMA 时不能用。
- **venv `~/PYTHONPATH_01/venv`**（用 anaconda 的 Python 3.11.7 创建），已装：
  - AnastrisTNG 2.0.0（可编辑安装，指向 `~/PYTHONPATH_01/AnastrisTNG/src`）
  - pynbody 1.6.0、numpy 1.26.4、scipy 1.13.1、pandas 2.2.2、h5py 3.11、matplotlib 3.9、astropy 6.1、wheel 0.48
  - **agama 1.0**，从 `~/PYTHONPATH_01/Agama` 源码编译，OpenMP 用 Homebrew 的 libomp（`/opt/homebrew/opt/libomp`）
  - 没有 jupyter
- Homebrew 5.1.15、Xcode 命令行工具都已安装。

### 5.2 TNG JupyterLab 服务器
- 模拟数据：`/home/tnguser/sims.TNG/TNG50-1/output`
- 旧版 AnastrisTNG（约 1.x，接口与本地 2.0.0 不同）和 galpos
- 用户在 `~/MW/casuality/` 下运行 data_obtain.py。轨道文件需要软链接过来：`ln -s ~/423_sample/TNG50_disk_trajectory.hdf5 .`

### 5.3 RLCD（第 2 阶段）
- **已决定（用户同意）**：另建一个专用 venv `~/PYTHONPATH_01/venv_rlcd`，用 Homebrew 装的 Python 3.12 创建，只用来跑 RLCD。原因：仓库要求 Python 3.12，并固定了一批包的版本（scipy 1.11.4、pandas 2.3.1、matplotlib 3.10.5、torch 2.8.0 等），和主 venv 冲突；分开可以保证主 venv（AnastrisTNG、AGAMA）不受影响。
- 依赖清单 `scripts/local/rlcd_requirements_mac.txt`：由仓库的 `rlcd/py312_requirements.txt` 删去 nvidia-*、triton（只有 Linux + GPU 才有）以及 dataclasses、typing（旧版 Python 的补丁包）得到。
- Graphviz（`brew install graphviz`）用来把因果图画成 png。
- Java 只在 `--stage1_method fges` 时需要（Tetrad）；默认的 `ges` 用 causal-learn，不需要 Java。
- 仓库克隆到 `casuality/Causal-Galactic-Archaeology/`（约 260 MB，含预先算好的结果和 Tetrad jar）。RLCD 的命令都要在它的 `rlcd/` 子目录里运行，结果写在 `rlcd/results_<scenario>/`。

---

## 6. 流程详解

### 6.1 服务器端提取：`scripts/server/data_obtain.py`
**基础**：用户的 `001_orbit_new.py`，保留了其中所有修正。

**样本（12 个星系，snap 99 的 Subhalo ID）**：
487742, 530330, 538905, 545437, 546870, 547293, 552581, 557721, 560421, 562029, 567607, 569251
（原列表里的 607654 不在 423 样本中，轨道文件里没有它，已删除。）

**每个星系做的事**：
1. `snap99.load_particle(ID)` 加载 dm/star/gas/bh 四类粒子。恒星额外读取 `GFM_StellarFormationTime, BirthPos, BirthVel, GFM_InitialMass, GFM_Metals`。然后调用 `physical_units()`。
2. **必须在 face_on 之前**取出 BirthPos/BirthVel 等恒星数组，因为 face_on 会就地旋转所有三维数组。
3. 剔除风粒子，只保留满足 `0 < aform ≤ 1`、`tform > 0` 且两者都有限的恒星。
4. 丰度，按粒子数之比计算，太阳值用 Asplund 2009：
   - [Fe/H] = log₁₀((X_Fe / A_Fe) / (X_H / A_H)) + 4.50
     意思是：Fe 和 H 的原子个数之比，与太阳的同一比值相比，取对数。
     - X_Fe、X_H：Fe、H 的质量分数，来自 `GFM_Metals` 第 8、0 列
     - A_Fe = 55.845、A_H = 1.008：原子量，用来把质量比换成个数比
     - 4.50：太阳的 12 − logε_Fe = 12 − 7.50
   - [Mg/Fe] = log₁₀((X_Mg / A_Mg) / (X_Fe / A_Fe)) − 0.10
     - X_Mg：Mg 的质量分数，来自 `GFM_Metals` 第 6 列
     - A_Mg = 24.305：Mg 的原子量
     - 0.10：太阳的 logε_Mg − logε_Fe = 7.60 − 7.50
5. 出生时间夹取：出生时间比轨道终点晚、但晚不到 1e-3 Gyr 的恒星（z=0 附近出生，浮点舍入造成），把出生时间夹回轨道终点。
6. **出生半径**：用 galpos 的 `GalaxyPoseTrajectory`，读入轨道文件里的 `t, cpos, cvel, ang_mom`，`box_size = 35000/0.6774`。再用 `pynbody_decorate.make_star_birth(...)` 和 `align_with_galaxy(orientation_align=True)`，把每颗星的 BirthPos 换算到"它出生那一刻主前身星系的盘坐标系"。
   - `birth_R_2D` = √(x_b² + y_b²)，意思是出生位置在当时盘面上离中心的距离（柱坐标半径）；x_b、y_b 是出生时刻盘坐标系里的盘面坐标。**这就是 Rbirth 真值。**
   - `birth_R_3D` = √(x_b² + y_b² + z_b²)，意思是出生位置到当时星系中心的直线距离；z_b 是出生时离盘面的高度。
   - 出生时间早于轨道文件起点的恒星，出生半径为 NaN（`has_birth = False`）。这些星**保留**，因为算引力势需要它们的质量。
7. `face_on(rmax=8, alignwith='star')`：旋转整个子晕到 z=0 的盘坐标系，然后保存所有粒子。
8. 每个星系写一个 HDF5（float32，gzip 压缩）；先写临时文件再改名，中断后重跑会跳过已完成的星系。

**服务器上的运行命令**：
```bash
cd ~/MW/casuality
ln -s ~/423_sample/TNG50_disk_trajectory.hdf5 .        # 第一次运行前做一次
nohup python3 data_obtain.py > data_obtain.log 2>&1 &
python3 merge_hdf5.py                                   # 生成 data_obtain_all.hdf5
```
然后从 JupyterLab 下载，放到本地 `data/tng50_extract/`。
（注：下载文件的 header 属性 `script` 写的是 `001_orbit_causal.py`，这是改名前的文件名，内容相同。）

### `data_obtain_all.hdf5` 的结构
```
/<ID>/header   (只有属性) ID, snap=99, a≈1, h=0.6774, time_Gyr=13.80, face_on_rmax_kpc=8,
               face_on_alignwith='star', units, solar_scale, script
/<ID>/frame    center_pos(3), center_vel(3), rotation(3,3)     —— z=0 face_on 的返回值
/<ID>/stars    iord(u64), pos(N,3) kpc, vel(N,3) km/s, mass Msun, init_mass Msun, age Gyr, aform,
               FeH, MgFe, X_H, X_Mg, X_Fe, birth_pos(N,3) kpc [出生时刻盘坐标系],
               birth_R_2D kpc, birth_R_3D kpc, has_birth(int8, 0/1)
/<ID>/wind     pos, mass       (风粒子，只用于引力势)
/<ID>/dm       pos, mass
/<ID>/gas      pos, mass
/<ID>/bh       pos, mass       (每个星系 1 个黑洞)
```
pos、vel 都在 z=0 的盘坐标系里：原点是星系中心（ssc 方法定出），z 轴沿 8 kpc 内恒星的角动量方向。

### 各星系的粒子数（来自核对报告）
| ID | 恒星 | 风粒子 | 暗物质 | 气体 | 无出生半径 | 盘内 v_φ 中位数 (km/s) |
|---|---|---|---|---|---|---|
| 487742 | 608,237 | 15 | 3,011,552 | 1,869,429 | 31 | 158.8 |
| 530330 | 1,053,023 | 0 | 2,920,963 | 275,220 | 3 | 191.4 |
| 538905 | 893,618 | 0 | 2,294,758 | 507,938 | 1 | 181.7 |
| 545437 | 722,225 | 0 | 1,859,730 | 645,766 | 66 | 171.3 |
| 546870 | 739,723 | 0 | 2,323,462 | 171,851 | 24 | 158.9 |
| 547293 | 908,363 | 0 | 1,755,535 | 605,028 | 28 | 183.2 |
| 552581 | 779,508 | 1 | 1,771,153 | 419,905 | 55 | 182.4 |
| 557721 | 712,029 | 0 | 1,947,701 | 136,165 | 190 | 187.4 |
| 560421 | 526,387 | 2 | 1,296,156 | 455,951 | 7 | 166.0 |
| 562029 | 541,225 | 1 | 1,622,644 | 411,809 | 0 | 163.6 |
| 567607 | 515,043 | 27 | 1,494,362 | 388,236 | 156 | 173.2 |
| 569251 | 670,495 | 0 | 1,457,537 | 392,465 | 37 | 181.5 |

TNG50 的 [Mg/Fe] 中位数约 **+0.28 到 +0.33**，比 NIHAO 高得多，所以论文的低 α 切法不能照搬。[Fe/H] 中位数约 −0.21 到 −0.27。

### 核对结果（1.3）
- 新数据的年龄、质量、当前位置与旧的 `001_orbit` 和 `001_orbit_s99_3D_z` 逐颗恒星一致（差 < 1e-5）。
- **旧 `001_orbit/<ID>.txt` 的 `Current_R` 和 `Birth_R` 都是三维距离**（不是柱坐标半径），与新数据的 `birth_R_3D` 一致。
- 与 `001_orbit_s99_3D_z` 的出生半径差 ≤ 5 pc（约 3% 的恒星），可以忽略。
- 新数据比旧数据多出少数"z=0 附近出生"恒星的出生半径（旧数据里是 NaN），这是出生时间夹取的效果，属于预期。

### 6.2 本地轨道计算：`scripts/local/compute_orbits.py`
**运行**（在用户的 Mac 终端里）：
```bash
cd /Volumes/001HardDisk/MW/casuality
~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py --test      # 测试：第 1 个星系随机 5000 颗星 → orbits_test.hdf5
caffeinate -i ~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py 2>&1 | tee results/compute_orbits.log   # 全部
~/PYTHONPATH_01/venv/bin/python3 scripts/local/compute_orbits.py --ids 487742 530330   # 只算指定星系
```
- `caffeinate -i`：运行期间不让 Mac 休眠（合盖仍会停，要开盖、接电源）
- `2>&1 | tee ...`：把正常输出和报错一起显示在屏幕上，同时写进日志
- 已算完、且版本号一致的星系会跳过；版本号不同的会删掉重算。中断后重跑即可继续。

**单位**：`agama.setUnits(mass=1, length=1, velocity=1)`，即 Msun、kpc、km/s，时间单位 0.978 Gyr。

**引力势**（`build_potential`）：
- 盘部分用 **CylSpline**（轴对称柱坐标网格，30×30，R 和 z 的范围 rmin–60 kpc，第一个 z 格点 0.05 kpc）：r < 50 kpc 的恒星，加上 R < 30 kpc 且 |z| < 5 kpc 的气体
- 晕部分用 **Multipole**（轴对称，lmax = 8，径向 40 个格点，rmin–500 kpc）：全部暗物质、黑洞、风粒子，以及其余的恒星和气体
- 总势 = 两部分相加
- **网格内边界 rmin = 0.2 kpc**，约为 TNG50 在 z=0 的引力软化长度。最初让 AGAMA 自动定网格，结果铺到了 0.01 kpc，那里的粒子噪声让圆轨道角动量随半径不单调，星系 530330 的 ActionFinder 初始化因此报错。现在如果 0.2 kpc 仍然报错，会自动改用 0.5、1.0 kpc，实际用的值记在属性 `pot_rmin_kpc` 里。
- 两部分的系数导出到 `data/tng50_derived/potentials/<ID>_{halo,disk}.coef`，可以用 `agama.Potential(file=...)` 重新读入

**每颗星算的量**（只读 `/<ID>/stars/pos` 和 `vel`）：
- **E** = Φ(x) + ½|v|²
  意思是：总能量 = 所在位置的引力势能 + 动能（都按单位质量算）。
  - Φ(x)：引力势在这颗星所在位置的值
  - |v|：速度大小

  只处理 E < 0（被束缚）的恒星，其余所有轨道量记为 NaN（约占 0.6%，多是远处的晕星）。
- **J_φ** = x·v_y − y·v_x
  意思是：恒星绕 z 轴转动的角动量 L_z，正值表示和盘同向转。
  - x、y：盘面内的两个坐标
  - v_x、v_y：对应的两个速度分量

  对所有恒星都算。
- **J_R、J_z**：`agama.ActionFinder(pot, interp=False)`，Stäckel fudge 方法，返回三列 (J_R, J_z, J_φ)
- **Rg** = `pot.Rcirc(L=|J_φ|)`，即满足 |J_φ| = Rg · v_c(Rg) 的半径。
  意思是：带着这个角动量走正圆轨道时所在的半径。
  - v_c(R) = √(R · ∂Φ/∂R)：半径 R 处的圆周速度
  - ∂Φ/∂R：引力势随半径的变化率，也就是盘面上的径向引力
- **r_peri、r_apo、e、z_max**：用 `agama.orbit` 积分轨道，时长 = 15 × `pot.Tcirc(xv)`（和这颗星能量相同的正圆轨道的周期），均匀记录 1501 个点，积分精度 1e-7，每批 5000 条轨道
  - r_peri = 轨道上三维距离 r 的最小值，r_apo = 最大值
  - **e = (r_apo − r_peri) / (r_apo + r_peri)**
    意思是：轨道越扁，e 越接近 1；接近正圆，e 接近 0。r_apo 是远心距，r_peri 是近心距。
  - z_max = 轨道上 |z| 的最大值
  - **论文没写偏心率用三维 r 还是柱坐标 R**，这里用的是三维 r，尚待用户确认
- 旋转曲线：R = 0.5–30 kpc 共 60 个点，v_c = √(−R · F_R)，F_R 是径向引力

**`orbits_all.hdf5` 的结构**：
```
/<ID>/Jr, Jz, Jphi, E, Rg, rperi, rapo, ecc, zmax    float32，长度 = 该星系恒星数，
                                                     行顺序与 data_obtain_all.hdf5 的 /<ID>/stars 完全相同
/<ID>/rotation_curve/R, vc
/<ID> 属性: n_stars, test_mode, pot_version='v2-rmin-grid', pot_rmin_kpc, N_PERIOD=15, TRAJ_SIZE=1501, units
（测试模式另有 star_index，记录抽到了哪些恒星）
```

**测试结果（487742，5000 颗星，旧的自动网格；新设置下全量结果的数值几乎相同）**：
- 5000 颗中 4963 颗束缚
- 各项检查全部通过：J_φ 等于直接算的 L_z；现在的 r 都在 [r_peri, r_apo] 之间；|z| 都不超过 z_max；Rg · v_c(Rg) / |L_z| 的中位数为 1.000
- v_c(8 kpc) = 185 km/s，0.5–28 kpc 范围内 171–188 km/s，较平坦
- 盘内（4 < R < 12，|z| < 1）中位数：J_z = 18，J_φ = 1074，e = 0.25，Rg − R = −0.7 kpc
- 年轻星（< 3 Gyr）的 e 中位数约 0.20，比银河系薄盘（约 0.1）"热"。这是 TNG50 盘本身的特点，不是计算错误

**全量运行记录**：487742 用时 4.4 分钟（新设置下 v_c(8) = 185.4 km/s，结果与测试一致）；530330 约有 105 万颗星，约 10 分钟。

---

### 6.3 RLCD 复现测试（g2.79e12，作者数据）
运行命令见第 8.3 节。结果写到 `Causal-Galactic-Archaeology/rlcd/results_causal_vertical_heating_g279_e/result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot`（覆盖了作者提交的版本；作者的原版可以用 `git show HEAD:<路径>` 查看）。
- **作者的图（15 条边）**：L1 → {L2, age, FeH, MgH, logJZ}；L2 → {JPhi, ecc}；MgH → {FeH, JPhi, logJZ, ecc}；age → {JPhi, FeH, logJZ, ecc}。这和 `L1L2calculation/GetL1_L2_279sim.ipynb` 里写死的系数完全对应。
- **我们的图（18 条边）**：两个隐变量的名字对调了（我们的 L2 就是作者的 L1）。隐变量的结构相同：根隐变量 → {另一个隐变量, age, FeH, MgH, logJZ}，另一个隐变量 → {JPhi, ecc}。作者的 15 条边全部都在，另外多出 3 条：**MgH → age**、age → 导引半径隐变量、MgH → 导引半径隐变量。
- MgH → age 违反论文 3.2 节的先验（age 在观测变量里没有父节点）。论文的做法是网格搜索 alpha，得到多个候选图，再用这个先验筛选到唯一的图。
- **但要注意**：作者提交的图的文件名，正好就是我们用的这组 alpha（1e-5/1e-15/1e-5/1e-15）。也就是说，同一组参数，我们跑出的图和作者的不同。可能的原因：作者实际用的数据或代码版本和仓库里的不同；包的版本不同；或者作者按先验手动修改过这张图。
- 用户决定不跑全部 81 组，只再跑 3 组 alpha 对比。所有文件都在 `results/rlcd_g279/`（见下面的文件说明）。
- 用 `scripts/local/compare_rlcd_graphs.py` 比较（会自动按连接方式对齐隐变量编号，并检查先验；参考图可以写成 `git:<仓库目录>:<相对路径>`，读 git HEAD 里作者的原版）。报告在 `results/rlcd_g279/compare_4alpha.txt`。

| alpha (0/1/2/3) | 隐变量数 | 作者的 15 条边 | 多出的边 | 先验（age 没有观测父节点） |
|---|---|---|---|---|
| 1e-5 / 1e-15 / 1e-5 / 1e-15（作者文件名同款） | 2 | 15/15 | MgH→age、MgH→L2、age→L2 | 不通过 |
| 1e-5 / 1e-10 / 1e-5 / 1e-10 | 2 | 15/15 | age→MgH、MgH→L2、age→L2 | 通过 |
| 1e-10 / 1e-15 / 1e-10 / 1e-15 | 4 | 2/15 | 结构完全不同 | 通过 |
| 1e-15 / 1e-15 / 1e-15 / 1e-15 | 4 | 3/15 | 结构完全不同 | 通过 |

  - alpha0 = 1e-5 时：隐变量结构和论文一样（L1 → L2、age、FeH、MgH、logJZ；L2 → JPhi、ecc），作者的边全部都在，只多 3 条：age 和 MgH 之间一条（方向随 alpha 变化），以及 age、MgH 各一条指向 L2。
  - alpha0 = 1e-10 或 1e-15 时：出现 4 个隐变量，分成两对——一对指向 {age, FeH, MgH}，另一对指向 {JPhi, ecc, logJZ}——和论文的结构不同，但也满足先验。
  - 所以"只用先验筛选"并不能像论文说的那样收敛到唯一一张图；作者提交的最终图，比任何一组的直接输出都少几条边。可能是作者用了不同的代码或数据版本，或者手动删掉了一些边。安装本身没有问题。
- 第一阶段 GES 把 6 个变量都分进同一组（`result_stage1ges_sparsit2.dot`）。
- 注意：仓库自带的 `grid_search.py` 里调用的是 `python`，要先 `source ~/PYTHONPATH_01/venv_rlcd/bin/activate` 再运行，否则会用错 Python。

**`results/rlcd_g279/` 文件说明**（2026-09-26 整理）
- `author_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png`：作者的原图（从仓库的 git HEAD 取出）；`author_rerendered_same_style.png` 是用本机 Graphviz 重新画的同一张图
- `ours_seed{0..4}_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png` + `ours_seed{0..4}.log`：同一组 alpha、不同 PYTHONHASHSEED 的 5 次运行
- `ours_noseed_result_alpha<四个alpha>_ranktestscale1.dot/.png` + `ours_noseed_alpha<...>.log`：没有固定种子的 4 次运行（4 组不同 alpha）。其中 1e-5/1e-15/1e-5/1e-15 那张在仓库里被后来的运行覆盖了，dot 是按当时的输出重写的，png 是重新画的
- `ours_stage1ges_sparsit2.dot/.png`：第一阶段 GES 的结果（最后一次运行的）
- `compare_4alpha.txt`、`compare_5seeds.txt`：比较报告；`compare_author_vs_seed1.png`：作者的图和 seed 1 的图并排对比
- `compare_author_vs_seed3_diff.png`（及 `diff_author_vs_seed3_*.png`）：作者的图和 seed 3 的图，差异的边用颜色标出。两张图共有 11 条边，不同的 4 条是：作者 L1→FeH、L1→logJZ、L2→JPhi、L2→ecc；seed 3 L1→JPhi、L1→ecc、L2→FeH、L2→logJZ。也就是说，(JPhi, ecc) 和 (FeH, logJZ) 挂在哪个隐变量下面正好对调，第二个隐变量的物理含义因此从“导引半径”变成了别的东西。两张图乍看很像，要仔细看每条边
- 仓库 `Causal-Galactic-Archaeology/rlcd/results_causal_vertical_heating_g279_e/` 里被覆盖的两个文件已经用 `git checkout` 恢复成作者的原版；我们运行产生的新文件都移到了这里。**以后在仓库里运行 RLCD，结果会写回仓库的 results_* 目录，跑完要及时移出来，避免覆盖作者的文件。**

**随机性检验（同一组 alpha 1e-5/1e-15/1e-5/1e-15，重复 5 次）**
- 命令：在 `rlcd/` 下运行 `PYTHONHASHSEED=$s ~/PYTHONPATH_01/venv_rlcd/bin/python3 run_rlcd.py ...`，s = 0–4，每次运行后把 dot 复制出来（现存为 `results/rlcd_g279/ours_seed$s_result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot`）。比较报告：`results/rlcd_g279/compare_5seeds.txt`。
- 结果：
  - seed 1：**和作者的图完全相同**（15 条边）
  - seed 0、2、3：另一种 15 条边的结构——根隐变量 → {另一个隐变量, JPhi, ecc, MgH, age}，另一个隐变量 → {FeH, logJZ}。满足先验，但物理含义和论文不同
  - seed 4：seed 0 的结构再多 3 条边（X_age→X_MgH，以及 X_MgH、X_age → 非根隐变量）。满足先验
  - 最早那次没有固定种子的运行：作者的结构再多 3 条边，其中 MgH→age 不满足先验
- **原因**（读源码得到）：
  - (1) RLCD 大量用 Python 集合存变量名，按集合的顺序依次做检验、合并隐变量；而 Python 默认每次启动都随机设定字符串哈希，集合的遍历顺序因此每次不同。
  - (2) `StructureLearning/RLCD/LatentGroups.py` 第 815 行 `findRandomLatentChild` 里有 `random.sample(...)`，没有设种子（只在处理"非原子覆盖"时才会用到）。`misc.py` 第 211 行的 `bootStrapCovariance` 虽然也用随机数，但整个代码里没有调用它，不影响结果。
- **对论文 3.2 节先验的理解**：先验只是一个"筛子"——把 age 有观测父节点的候选图去掉，并不修改图本身。它能去掉最早那次运行的 18 边图，但 seed 0、2、3、4 的图都能通过先验，结构却不同于作者的。所以"用先验一筛就收敛到唯一的图"只是作者那几次运行碰巧如此；作者的图是 RLCD 可能输出的结果之一（5 次里出现 1 次）。
- **30 个种子的频率（2026-09-27，PYTHONHASHSEED = 0–29，同一组 alpha）**：作者的图 7 次（1、7、8、11、14、22、25，23%）；最常见的是 seed 0 那种 15 条边的图（12 次，40%）；其余 11 次是 4 种 18 条边的图。
- **以后运行 RLCD 必须**：固定 PYTHONHASHSEED 和随机数种子，并用多个种子重复，统计每种结构出现的频率。
- **作者有没有提到随机性**（2026-09-26 查过）：
  - 论文：没有提到随机性、种子或重复运行，只说用先验筛选后"收敛到唯一的图"。
  - Jin 等的仓库：README 和 notebook 里都没有提到。代码里有一个 `utils/utils.py: set_random_seed()`（会设置 random、numpy、torch 的种子和 PYTHONHASHSEED），但**从来没被调用**；`misc.py` 的 `reorderCovers` 注释写着 "Reduces some randomness"，说明 RLCD 原作者知道存在顺序带来的随机性，只做了部分处理。
  - RLCD 原始仓库（github.com/dongxinshuai/scm-identify，Dong et al.）：README 只提醒**参数学习**"may need multiple run to get ideal results due to random initialization"（优化器随机初始化），没有提结构学习的随机性。
  - 所以 Jin 等仓库里那 13 个 `*_param_by_trek_nll*.dot` 文件，就是对同一张图重复做参数学习的结果，文件名里的数字是负对数似然（NLL），越小拟合越好；最小的是 nll1.4979。**我们以后做参数学习也要重复多次，取 NLL 最小的那次。**

### 6.4 轨道结果检查（1.5）：`scripts/local/check_orbits.py`
运行：`~/PYTHONPATH_01/venv/bin/python3 scripts/local/check_orbits.py`（Claude 在虚拟机里用 `PYTHONPATH=$HOME/pylib python3 ...` 跑，用 pyfive 读 HDF5）。输出在 `results/orbit_check/`：`summary.txt`、`summary.csv`、`rotation_curves.png`、`ecc_disk.png`、`Rg_vs_R.png`、`logJz_disk.png`；另外 `MgFe_FeH.png` 是 4<R<12、|z|<3 的 [Mg/Fe]–[Fe/H] 分布图（选星讨论时画的）。
- 12 个星系全部通过：束缚比例 99.4–100%（487742 最低，为 99.42%）；束缚恒星中各量没有 NaN；|J_φ − L_z| < 0.004；r 都在 [r_peri, r_apo] 之间；|z| ≤ z_max；Rg·v_c(Rg)/|L_z| 中位数 = 1.000；e 在 [0, 1] 内。所有星系都用 rmin = 0.2 kpc（没有触发重试）。
- v_c(8 kpc) = 185–229 km/s，v_c 最大值 191–251 km/s。
- 盘内（4<R<12，|z|<1）中位数：J_z 9.9–18.5，J_φ 1074–1403，J_R 27–99，e 0.18–0.32（年轻星 < 3 Gyr 为 0.09–0.25），Rg − R 为 −0.36 到 −0.85 kpc，逆向转动 0.4–6.5%。
- Rg–R 图里 Rg ≈ 0 的水平带是几乎不转动的核球/晕星（L_z ≈ 0），属于正常现象。
- [Mg/Fe]：大多数星系**没有清楚的高/低 α 两条序列**，只有 545437、567607、560421、487742 在 [Mg/Fe] ≈ 0.25 附近有一团低 α 星。

### 6.5 选星与 RLCD 输入（1.6）：`scripts/local/make_rlcd_input.py`
运行：`~/PYTHONPATH_01/venv/bin/python3 scripts/local/make_rlcd_input.py`（约 2 分钟）。选星统计：`results/rlcd_input_selection.txt`。
**选星条件（2026-09-26 与用户确定）**：
1. 束缚（E < 0），且 J_R、J_z、J_φ、Rg、偏心率都有值
2. has_birth，且出生半径 R_bir < 20 kpc
3. **7 < R < 10 kpc**（论文正文的范围）
4. **|z| < 0.3 kpc**（用户指定）
5. **[Fe/H] > −0.5，不切 [Mg/Fe]**（TNG50 没有清楚的 α 双序列）
6. J_z > 0（要取对数）；逆向转动的星保留（和作者仓库数据一致）

**观测误差（noise 版本）**：年龄乘以 (1 + 0.10·N)，[Fe/H] 加 0.02·N，[Mg/Fe] 加 0.06·N（N 是标准正态随机数），[Mg/H] 用加误差后的 [Mg/Fe] + [Fe/H] 计算。随机种子 = 20260926 + 星系 ID。选星用的是不加误差的值，所以两个版本逐行对应同一批星。

**输出**（每个星系一对文件）：
- `data/rlcd_input/nonoise/<ID>_e.csv`、`noise/<ID>_e.csv`：6 列 `logJZ, JPhi, ecc, age, FeH, MgH`，格式和作者的 `data/g2.79e12_e.csv` 相同
- `<ID>_truth.csv`（和 `_e.csv` 逐行对应）：`iord, star_index`（在 data_obtain_all.hdf5 /stars 里的行号）`, R, z, Rg, R_bir, R_bir_3D, z_bir, age, FeH, MgFe, MgH, JR, JZ, JPhi, ecc, rperi, rapo, zmax, E`；noise 版本还多 `age_obs, FeH_obs, MgFe_obs, MgH_obs`

**每个星系的样本量**（最终）：487742: 13,190；530330: 24,154；538905: 34,970；545437: 16,222；546870: 11,484；547293: 32,525；552581: 13,573；557721: 25,326；560421: 19,696；562029: 7,299；567607: 12,496；569251: 29,035。合计 239,970。选星时样本主要被 R 范围和 |z| < 0.3 这两步大幅削减。
- 中位数：logJZ 0.48–0.93，J_φ 1356–1753，e 0.13–0.28，age 1.3–7.6 Gyr，R_bir 5.3–9.0 kpc（5–95% 约 2.5–13.5 kpc），Rg 7.1–8.3 kpc。
- 已核对：`10**logJZ` 等于 truth 里的 JZ；JPhi 一致；没有 NaN；noise 版本的年龄相对误差标准差 = 0.100，[Fe/H] 误差标准差 = 0.020。

### 6.6 做法 A：`scripts/local/rlcd_A.py`（用 venv_rlcd 运行）
```bash
cd /Volumes/001HardDisk/MW/casuality
~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version g279 2>&1 | tee results/rlcd_A_g279.log
~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version nonoise --ids 538905 2>&1 | tee results/rlcd_A_test538905.log
~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version nonoise 2>&1 | tee results/rlcd_A_nonoise.log
~/PYTHONPATH_01/venv_rlcd/bin/python3 scripts/local/rlcd_A.py --version noise 2>&1 | tee results/rlcd_A_noise.log
```
参数：`--n-repeat 10`（参数学习重复次数），`--epochs 600 --lr 0.005 --lr-decay plateau`（与仓库 README 相同），`--dot` 默认用作者的图 `results/rlcd_g279/author_...dot`。
1. **数据**：`data/rlcd_input/<version>/<ID>_e.csv`（g279 用仓库的 `data/g2.79e12_e.csv`，真值用 `g2.79e12_multi_ele_action_Jin_SN.csv`，两者逐行对应，已核对）。先标准化：每列减平均值、除以标准差（和仓库 `DataModel.generate_data(normalized=True)` 相同），列名前加 `X_`。
2. **参数学习**：直接调用仓库的 `parameter_by_trek`，它用观测量的协方差矩阵做线性高斯模型的最大似然估计（隐变量的噪声方差固定为 1），从随机初值开始。每次重复前设置 `torch.manual_seed(seed)`，seed = 0…N−1，取 NLL 最小的一次。
3. **求隐变量**：对隐变量 L 的每个子节点 c，先去掉 c 的其他父节点的贡献，r_c = c − Σ a_pc·p；再取 L = Σ a_Lc·r_c / Σ a_Lc²。意思是：找一个 L 值，让"由 L 预测出的子节点"和"实际子节点"差的平方和最小。
   - a_Lc：L 指向 c 的边系数；a_pc：c 的其他父节点 p 指向 c 的边系数
   - 从最下层的隐变量开始算（先 L2，再把 L2 当作 L1 的子节点算 L1）

   这和作者 notebook 里逐颗星用 Nelder-Mead 最小化 chisq 完全等价；已经用随机数据和 notebook 的系数核对过，结果一致到小数点后 6 位。
4. **换算到物理单位**：用 np.polyfit 做线性拟合，f(L1) = α·L1 + β 拟合 Rbirth，f(L2) 拟合 Rg；给出皮尔逊 r、残差标准差、中位绝对偏差。隐变量的正负号是任意的，所以斜率的符号不重要。

**输出**：`results/rlcd_A/<version>/<ID>/`，包括 `nll_repeats.csv`、`best_params.dot/.png`、`latents.csv`（iord, L1, L2, f_L1, f_L2, R_bir, Rg, age）、`fit.txt`；汇总在 `results/rlcd_A/<version>/summary.csv/.txt`；图是 `L1_vs_R_bir.png`、`L2_vs_Rg.png`。

**第一次试跑（v1 脚本：重复 10 次、取代码输出的 NLL 最小的一次）失败了，原因已查明（2026-09-26）：**
- g279 得到 r(L1, Rb) = 0.08（论文是 0.64）；538905 得到 r(L1, Rb) = 0.08、r(L2, Rg) = 0.25。
- **求隐变量和线性换算这两步没有问题**：用作者 notebook 里的系数和我们的公式，在 g279 上精确得到 Rb = −1.969·L1 + 9.309（论文 −1.97、9.31），r = 0.640；Rg = −1.289·L2 + 9.518（论文 −1.29、9.52），r = 0.976。
- **问题出在"挑哪次参数学习的结果"**：
  1. 仓库的 `parameter_by_trek` 算 NLL 时用 `torch.slogdet(...)[1]`，也就是行列式的**绝对值**的对数，而且不检查每个节点隐含的噪声方差是否 ≥ 0。所以物理上不成立的解也可能得到"很小"的 NLL。g279 的理论最小 NLL 是 −1.275（= ½·logdet(实际相关矩阵) + 3），但第 5 次得到 −7.03，显然不可能。v1 脚本挑中的正是这种无效解（g279 的 seed 4、538905 的 seed 5 都无效）。
  2. 即使只看有效解，NLL 最小也不是作者用的那组：作者 notebook 和论文用的是仓库里 **NLL 第二小（1.586）** 的系数，而不是最小（1.498）的。1.498 那组的 L1 几乎和化学量断开（L1→MgH = 0.012），L1 与 Rb 的相关只有 0.01。原因是 [Fe/H] 和 [Mg/H] 的相关系数高达 0.997，NLL 被这一项主导。**1.586 那组的 RMSE（模型相关矩阵与实际相关矩阵之差）最小，为 0.117**；1.498 那组是 0.260。
  3. 参数学习的结果对随机初值很敏感：同一个星系，不同的有效解给出的 |r(L1, Rb)| 从 0.01 到 0.6 都有。
- **v2 脚本的改法**（`rlcd_A.py` 现在的版本）：
  - 每次参数学习后，用 `implied_cov` 自己算模型隐含的相关矩阵和每个节点的噪声方差，去掉无效解（噪声方差 < 0 或矩阵不正定），并重新算出正确的 NLL；
  - 在有效解里，默认挑 RMSE 最小的（`--select rmse`，这个标准能挑中作者用的那组），也可以用 `--select nll`；
  - 默认重复 50 次（每次约 2 秒）；
  - 每次的 |r(L1,Rb)|、|r(L2,Rg)| 都写进 nll_repeats.csv，**只用于诊断，不参与挑选**；fit.txt 里给出排前 5 的有效解的 |r(L1,Rb)|，用来看结果稳不稳定。
- 这是论文和仓库都没有说明的第二个问题（第一个是结构学习的随机性）：参数学习有很多局部最优解，作者实际用的也不是 NLL 最小的解。

**v2 重跑结果与更深层的问题：作者的因果图冻结后，边系数不可识别（2026-09-26）**
- v2 重跑：g279 挑中 seed 31（RMSE 0.108），r(L1, Rb) = 0.626、r(L2, Rg) = 0.86；538905 挑中 seed 46（RMSE 0.078），r(L1, Rb) = 0.056、r(L2, Rg) = 0.78。但按 RMSE 排前 5 的有效解，|r(L1,Rb)| 分别是 0.63/0.01/0.31/0.01/0.11（g279）和 0.06/0.20/0.52/0.44/0.19（538905），极不稳定。
- 进一步用 `scripts/local/rlcd_identifiability.py` 检查（numpy/scipy，Claude 在云端跑的；结果在 `results/rlcd_A/identifiability/`）：
  1. **雅可比矩阵的秩 = 14 < 15**。作者的图有 15 条边，6 个观测量正好有 15 个相关系数；把"15 个边系数 → 15 个相关系数"这个映射对系数求导，得到的矩阵只有 14 个独立方向。也就是说，存在一个方向，沿着它同时改系数，观测到的相关矩阵**完全不变**。这个方向主要是 L1→age 和 L1→MgH 此消彼长，同时 age、MgH 指向 ecc、JPhi、FeH 的边做补偿。这是图的结构性质，和数据无关，g279 和 TNG50 都一样。
  2. **用真正收敛的优化器（L-BFGS，多个初值）做最大似然**：g279 有 14 个解都达到了理论最优 NLL（相差 < 0.01），它们的 |r(L1,Rb)| 从 0.06 到 0.61；538905 有 25 个这样的解，|r(L1,Rb)| 从 0.14 到 0.45。
  3. **剖面**：把 L1→age 固定在 −0.95…0.95 的一串值上，其余 14 个系数重新拟合。g279 上 NLL 几乎处处等于理论最小值 −1.275，而 |r(L1,Rb)| 从 0.00 连续变到 0.63（L1→age ≈ 0 附近 L1 最像 [Mg/H]，r 最大约 0.68，但那里 NLL 略差）。538905 上 |r(L1,Rb)| 从 0.07 变到 0.43。
- **结论**：冻结作者的图后，**数据本身无法确定 L1 是什么**。L1 与出生半径的相关系数取决于优化器碰巧停在这条"等价解"曲线的哪个位置，而不是由数据决定。论文的 r = 0.64 对应的是一个**没有完全收敛**的解（NLL 1.586，而理论最优是 −1.275），恰好停在 L1 ≈ [Mg/H] 的一端。做法 A 按论文的方式得到的 L1–Rbirth 结果，不能作为"方法有效"的证据。
- 这是本项目目前最重要的发现（第三个问题，也是最根本的一个）。下一步需要和用户讨论怎么处理（见第 7 节第 9 条）。

**12 个星系的结果与严格复查（2026-09-26 晚，用户选了 (a)）**
- 第一轮 `rlcd_identifiability.py`：13 个数据集（g279 + 12 个 nonoise；noise 同样）雅可比秩全部 14/15。按"ΔNLL < 0.01"算，|r(L1,Rb)| 范围多为 ~0 到 0.4–0.7（`summary_{nonoise,noise}.*`）。
- **这个标准太宽松**：NLL 按每颗星平均，统计上"分不开"应是 2·N·ΔNLL < 3.84。而且第一轮的解没有严格收敛（NLL 误差 ~1e-4）、噪声方差下限用的是 0.01。
- 第二轮 `rlcd_identifiability_ci.py --extra 200`：第一轮 30 个初值 + 200 个新初值，全部用"噪声方差 ≥ 0"（罚项 1e7）+ 严格收敛（ftol 1e-15，gtol 1e-11）重新优化，保留 2·N·ΔNLL < 3.84 的解。结果（`summary_ci_{nonoise,noise}.*`，数字见 `REPORT_RLCD_TNG50.md` 5.6 节）：
  - 最优解几乎都把一个观测量的噪声方差压到 0：nonoise 10/12 是 J_φ（546870、560421 的最优解在内部，最小噪声方差 0.04、0.03）；noise J_φ 9、[Fe/H] 2、e 1；g279 是 [Fe/H] 或 e（Heywood 情形 → 线性模型与数据对不上）。
  - TNG50：多数星系统计上分不开的解仍给出 |r(L1,Rb)| 从 ~0 到 0.4–0.5；546870 较窄。
  - g279：严格标准下只有 1 个解，但它随搜索变化：只用第一轮 30 个初值时最优解 |r(L1,Rb)| = 0.25、|r(L2,Rg)| = 0.76；230 个初值时最优解（NLL 再低 0.0009）是 0.62、0.55。
  - 作者论文的系数：按正确公式 NLL = 1.593（最优 −1.272），似然比统计量约 6.8e5。
  - `linear_ceiling.py`：6 个量线性回归的上限 0.53–0.85（TNG50），单用年龄 0.23–0.69、单用 [Fe/H] 0.00–0.46；g279 单用 [Fe/H] 0.70。**nonoise 所有 12 个星系里 L1 的最有利值都不超过单用年龄或 [Fe/H] 中较好的那个；noise 有 2 个例外（560421：0.58 对 0.57；567607：0.71 对 0.67）。**
- `results/rlcd_A/README.md` 已按上面更新（含更正说明）。


## 7. 待定的决策（需要和用户确认）

1. **偏心率**：用三维 r（当前做法）还是柱坐标 R？
2. **是否加独立检查**：沿轨道检查能量守恒；用轨道积分直接计算 J_z，与 Stäckel 近似的结果对比。
3. **（已定，见 6.5 节）选星条件**（1.6）：
   - 半径范围：论文写 7–10 kpc，但仓库数据实际是 9.4–12.1 kpc，可能按星系大小缩放过。TNG50 星系的 v_c 约 160–190 km/s，比银河系小。
   - 高度限制：比如 |z| < 1 kpc 或 z_max 的上限
   - 低 α 切法：TNG50 的 [Mg/Fe] 中位数约 +0.3，要看每个星系的 [Mg/Fe]–[Fe/H] 分布再定，或者先不分
   - [Fe/H] 下限（−0.5 或 −1）
   - 是否只要有出生半径的星（has_birth = True），以及 Rbirth 是否 < 20 kpc
   - 是否排除逆向转动（J_φ < 0）或未束缚的恒星
4. **（已定）观测误差**：两种都生成（nonoise 是主版本，noise 按论文加误差）。
5. **（已定，见 6.5 节）RLCD 输入格式**：和仓库相同的 6 列 `logJZ, JPhi, ecc, age, FeH, MgH`，每个星系一个 csv（比如 `data/rlcd_input/<ID>_e.csv`）；真值另存 `<ID>_truth.csv`，列为 `iord, R, z, Rg, R_bir, ...`。J_z = 0 时 log₁₀ 会变成 −inf，需要处理（加一个下限，或剔除）。
6. **RLCD 的 Python 版本**：3.12 新环境，还是用 3.11 的 venv。
7. **怎么检验**：每个星系单独做一次完整的 RLCD 结构学习；还是像论文那样，用一个星系学出结构，冻结后迁移到其他星系，只重新估计系数。

---

8. **（新，2.2b）TNG50 用什么因果图**，候选做法：
   - A. 按论文 4.2、4.3 节：冻结作者 g2.79 的最终图（15 条边），只在每个 TNG50 星系上重新估计边系数，再算 L1、L2（仓库里 `*_with_sim_dag` scenario 就是这么做的）。和论文的比较最直接，而且不依赖能不能精确复现结构学习。
   - B. 在每个 TNG50 星系上从头跑结构学习（比如 alpha0 = 1e-5），看能不能独立发现两个隐变量。更能检验方法本身，但结果会随 alpha 变化。
   - 也可以两个都做：A 作为主结果，B 作为补充。

9. **（已定：选 (a)）边系数不可识别，做法 A 怎么继续？** 候选：
   - (a) 把"不可识别"本身作为主要结果：对每个 TNG50 星系，用剖面给出 |r(L1, Rbirth)| 的可能范围（最好到最差），并指出论文的 0.64 只是这条曲线上的一点；
   - (b) 用一个额外的约束去掉这个自由度（比如固定 L1→age 或 L1→MgH 的某个值，或者要求 L1 与某个量的关系），但这个约束必须有物理理由，而且不能偷看真值；
   - (c) 增加观测量（比如再加一种元素丰度或 J_R），让相关系数多于边数，使模型可识别；这会改变因果图，更接近做法 B；
   - (d) 暂停做法 A，先做做法 B，看 TNG50 的数据自己会学出什么图，以及那张图是否可识别。

## 8. 安装记录（出问题时照着重做）

### 8.1 更新 AnastrisTNG（用户本地）
Claude 那边访问不了 GitHub（从用户电脑上的虚拟机），所以做法是：在云端 `git clone`，打包成 git bundle 传过去，再在本地 `git fetch <bundle>`，然后 `reset --hard origin/main`（这一步需要删除权限）。更新后在终端运行下面的命令重新安装（已完成：venv 里现在登记的是 anastristng 2.0.0）：
```bash
cd ~/PYTHONPATH_01/AnastrisTNG && ~/PYTHONPATH_01/venv/bin/pip install -e .
```

### 8.2 安装 AGAMA 到 venv（最终成功的做法）
```bash
xcode-select -p                                              # 确认有编译器
~/PYTHONPATH_01/venv/bin/python3 -m pip install wheel        # venv 原来缺 wheel
brew install libomp                                          # 新版 OpenMP
cd ~/PYTHONPATH_01 && git clone https://github.com/GalacticDynamics-Oxford/Agama.git && cd Agama
rm -f Makefile.local && make clean                           # 必须删 Makefile.local, 否则会沿用旧配置
env -u CONDA_EXE -u CONDA_PREFIX -u CONDA_DEFAULT_ENV \
  ~/PYTHONPATH_01/venv/bin/python3 -m pip install . --no-build-isolation --config-settings --build-option=--yes
grep -n omp Makefile.local                                   # 应看到 /opt/homebrew/opt/libomp
```
测试（先 `cd ~`，离开源码目录）：
```bash
~/PYTHONPATH_01/venv/bin/python3 -c "import agama; p=agama.Potential(type='Plummer'); print(p.potential([1,0,0]))"   # 约 -0.707
```
**踩过的坑**：
- `pip install agama`（PyPI 上的 1.3 kB 引导包）会报 `output.json` 找不到 → 改为从源码安装
- 安装脚本会优先链接 anaconda 的 `libiomp5.dylib`（llvm-openmp 14.0.6），新版 clang 需要的 `__kmpc_dispatch_deinit` 它没有，链接失败 → 用 `env -u CONDA_*` 让脚本找不到 anaconda，转而使用 Homebrew 的 libomp
- 重装时 `--yes` 会自动沿用已有的 `Makefile.local` → 必须先删掉它

### 8.3 安装 RLCD（进行中，用户在终端运行）
```bash
brew install python@3.12 graphviz
/opt/homebrew/bin/python3.12 -m venv ~/PYTHONPATH_01/venv_rlcd
~/PYTHONPATH_01/venv_rlcd/bin/python3 -m pip install --upgrade pip
cd /Volumes/001HardDisk/MW/casuality
git clone https://github.com/ZehaoJin/Causal-Galactic-Archaeology.git
~/PYTHONPATH_01/venv_rlcd/bin/python3 -m pip install -r scripts/local/rlcd_requirements_mac.txt
```
测试：用仓库自带的 g2.79e12 数据复现论文的因果图
```bash
cd /Volumes/001HardDisk/MW/casuality/Causal-Galactic-Archaeology/rlcd
~/PYTHONPATH_01/venv_rlcd/bin/python3 run_rlcd.py --s causal_vertical_heating_g279_e \
    --alpha0 1e-5 --alpha1 1e-15 --alpha2 1e-5 --alpha3 1e-15 2>&1 | tee ../../results/rlcd_g279/<日志名>.log
```
（这一步会覆盖仓库里预先算好的 `results_causal_vertical_heating_g279_e/result_alpha1e-05_1e-15_1e-05_1e-15_ranktestscale1.dot/.png`，可以用 `git diff` 或 `git status` 对比是否和作者的结果一致。）

---

## 9. 已知的坑（尤其是给 Claude 的）

- **Claude 的 `device_bash` 是跑在用户 Mac 上的一个 Linux 虚拟机**：能读写已授权的文件夹，但**不能运行用户 Mac 上的 Python/AGAMA**。虚拟机里只有 Python 3.10 + numpy + pandas，没有 h5py 和 scipy。
- **Claude 的环境（云端和虚拟机）都被组织网络策略挡住了 PyPI 和 conda**，pip 装不了包。云端能访问 GitHub，虚拟机不能。
- **在虚拟机里读 HDF5 的办法**：云端从 GitHub clone `NCAS-CMS/pyfive`（纯 Python 的 HDF5 读取库），打补丁后传到虚拟机的 `~/pylib`：
  - `h5d.py` 里 `import fsspec.utils` 改成 try/except
  - `__init__.py` 和 `h5d.py` 里的 `version("pyfive")` 改成固定字符串
  - 传输方法：先 commit 到 casuality，再 `mv` 到 `~/pylib` 并删除 casuality 里的副本
  - 用 `PYTHONPATH=$HOME/pylib python3 ...` 运行；pyfive 的 Dataset 没有 `len()`，要用 `.shape[0]`

  **每次新会话的虚拟机都是新的，需要重新传。** `check_vs_001_orbit.py` 已经做了"没有 h5py 就用 pyfive"的兼容。
- **控制电脑时，终端只允许点击、不允许打字**（系统限制，用户给权限也不行）。所以运行命令一律请用户在终端里操作。
- **`device_commit_files` 反复提交同一个 stagedPath 时，有时会写入旧内容**：每次都换一个新的暂存文件名，提交后用 `md5sum` 核对。
- `/Volumes/001HardDisk_6`（外置硬盘）有时在虚拟机里挂载不上，过一会儿再试，或者改用 `device_stage_files`。
- 在 casuality 里删除文件需要用户批准删除权限（每个会话一次）。
- `TNG50_disk_trajectory.hdf5` 早期快照的对齐会切换到 baryon/all，所以最老恒星的"盘面"定义不太可靠；早期研究可以参考 `birth_R_3D`。
- 服务器上运行脚本时，`TNG50_disk_trajectory.hdf5` 必须在当前目录（用软链接），否则每个星系都报 "No such file"。

---

- **RLCD 的结果带有随机性**（集合顺序 + 没设种子的随机数）：同样的参数，每次运行的因果图可能不同。运行时要设置 `PYTHONHASHSEED=<数字>`，并用多个种子重复。
- **`device_commit_files` 会重新编码 PNG**：提交后 png 的 md5 和本地不同，但像素完全一样（已核对），不是写错。其他文件 md5 应一致。
- **在云端用 `pgrep -f`/`pkill -f` 时小心匹配到自己**：bash 命令行里含有同样的字符串，会导致"等待"永远不结束或把自己杀掉。用 `pgrep -f "python3 -u scripts/..."` 这种只在 python 进程里出现的字符串。
- **Claude 云端环境只有 2 个 CPU 核**：12 个星系的可识别性分析每轮约 1–2 小时。
- **可识别性分析的"一样好"要用统计标准**：2·N·ΔNLL < 3.84，不要用固定的 ΔNLL 阈值；比较前先严格收敛（polish）。

## 10. 过程时间线（简要）
- 2026-09-26：读论文和仓库 → 更新本地 AnastrisTNG → 设计服务器提取 → 用户在服务器运行 data_obtain.py 和 merge_hdf5.py，下载数据 → 与旧结果核对，一致 → 安装 AGAMA（经过三次失败后成功）→ compute_orbits.py 测试通过 → 全量运行时 530330 报错 → 改网格内边界 rmin = 0.2 kpc → 重跑中 → 安装 RLCD（venv_rlcd）→ 4 组 alpha 对比 → 发现 RLCD 结果随机；5 个种子中 seed=1 完全复现作者的图 → 轨道计算完成 → 1.5 检查全部通过 → 1.6 选星（7–10 kpc，|z|<0.3，[Fe/H]>−0.5）生成 RLCD 输入 → **第 1 阶段完成** → 做法 A：rlcd_A v1 失败（挑到无效解）→ v2 → 发现不可识别（雅可比秩 14）→ 用户选 (a) → 12 个星系 nonoise/noise 第一轮 → 发现 ΔNLL<0.01 标准太宽松 → 严格复查（rlcd_identifiability_ci.py）→ 写汇报文档 `REPORT_RLCD_TNG50.md`
