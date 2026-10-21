# 给 Claude（或任何 AI 助手）的项目说明

本仓库：用 TNG50 模拟检验 Jin et al. 的 Causal Galactic Archaeology（RLCD）方法能否恢复恒星出生半径 Rbirth。用户是天文学博士生。

## 开工前先做这几件事
1. **读 `STATUS.md`**（进度、问题、下一步，一页）。
2. 需要细节时读 `docs/pipeline_plan.md`：第 0 节（开工清单）、第 2 节（进度）、第 6 节（每一步的做法和数据结构）、第 9 节（已知的坑）。要理解方法和结论时读 `REPORT_RLCD_TNG50.md`。
3. 确认数据是否就位：`bash tools/check_setup.sh`（缺数据就先 `bash tools/link_data.sh`）。
4. 每次回复末尾都告诉用户：**项目进行到哪里了、接下来做什么**。

## 与用户的工作约定
- 用户用中文写就用中文回复，用英文写就用英文。
- 用户不熟悉数学符号：**每个公式都要配白话解释，并逐个说明每个符号的含义**（文档里也是这样写的，保持这个风格）。
- 长时间的工作要主动汇报进展（在做什么、做到哪里），不要一直沉默到结束。
- **不要控制用户的电脑、不要读屏幕**，除非用户明确要求。需要在用户 Mac 上运行的东西（Python、AGAMA、RLCD），给出完整的可复制命令，让用户在终端里运行，输出存成日志放到 `results/`，再读日志了解结果。
- 用户终端默认在 anaconda base，所以命令里写 venv 的完整路径：
  - 主环境：`~/PYTHONPATH_01/venv/bin/python3`（Python 3.11.7；numpy、scipy、pandas、h5py、matplotlib、AGAMA）
  - RLCD 专用：`~/PYTHONPATH_01/venv_rlcd/bin/python3`（Python 3.12；作者仓库的依赖）
  - 新包装在对应的 venv 里，不要装进 anaconda base。
- 用户偏好把数据下载到本地算，而不是在 TNG 服务器上算。服务器（TNG JupyterLab）上的 AnastrisTNG 是**旧版**，服务器脚本不能按 2.0.0 的接口写。
- **推进后同步更新三处**：`STATUS.md`、`docs/pipeline_plan.md` 第 2 节、claude.ai 项目文档 `claude/pipeline_plan.md`（若能访问）。pipeline_plan 要详细到“一个没有记忆的新对话能照着复现并继续”。
- 用户的决定：做法 B 的结果**只记在 pipeline_plan 和 `PlanB/README.md`**，不写进 REPORT，也不在 REPORT 里专门提及。

## 仓库和数据的位置
- **本仓库是代码和文档的主副本**（GitHub 同步；本机路径 `~/Project2026/Casual-relationship-Rbirth-TNG50`）。所有脚本都用“相对仓库根目录”的路径（`ROOT = Path(__file__).parents[n]`），在仓库根目录下运行即可；脚本注释里写的 `cd /Volumes/001HardDisk/MW/casuality` 是打包前的旧位置，换成仓库根目录就行。
- **大数据（约 17 GB）不进 GitHub**，存在外置硬盘 `/Volumes/001HardDisk/MW/casuality`（打包前的原项目文件夹，也保留着打包时的代码副本，以后不再更新它的代码）。`tools/link_data.sh` 把它们以软链接接回仓库的同名位置；清单和校验和在 `DATA.md`、`data_manifest.tsv`。
- 作者的仓库 `Causal-Galactic-Archaeology/` 也不进 GitHub：本机由 link_data.sh 链接，新电脑用 `tools/setup_rlcd_repo.sh` 克隆（固定在 commit 8c0a5a4）。**在它的 `rlcd/` 下跑 RLCD 会覆盖作者提交的结果文件**，跑完把新文件移到 `results/`，作者原版用 `git checkout` 恢复。
- 新生成的大文件（> 50 MB 的 csv/hdf5）不要提交：放进 `.gitignore` 已覆盖的目录，或者在 `.gitignore`、`link_data.sh`、`data_manifest.tsv` 里补上。
- 不需要大数据就能重跑的部分：`star_condition/*/output_data/identifiability_moments.py`（只读已提交的 `moments_*.json`）、`results/复现_Jin/` 的两个脚本（只需要作者仓库）。

## 必须记住的技术要点
- RLCD 结构学习的结果是随机的：运行时一定设 `PYTHONHASHSEED=<数字>`，并用多个种子重复、统计结构频率。
- 估计边系数：要求噪声方差 ≥ 0，用能收敛的优化器（L-BFGS-B）和几百个初值，报告所有与最优“统计上分不开”（2·N·ΔNLL < 3.84）的解，不能只报一个；挑解时绝不看真值。
- 比较 |r| 时总是带上基线：单用年龄、单用 [Fe/H]、单用 J_φ、6 个量的线性回归上限。
- 两个分析进程同时跑时设 `OMP_NUM_THREADS=1`。
- Claude 自己的执行环境（云端或用户电脑上的 Linux 虚拟机）不能运行用户 Mac 上的 Python/AGAMA；虚拟机里没有 h5py/scipy，读 HDF5 可用纯 Python 的 pyfive（做法见 pipeline_plan 第 9 节）。

## 文档索引
| 文件 | 用途 |
|---|---|
| `STATUS.md` | 进度、问题、下一步（先读） |
| `README.md` | 给人看的使用说明：安装、数据、GitHub 同步、运行顺序 |
| `docs/pipeline_plan.md` | 复现细节：每个脚本的命令、数据结构、安装记录、已知的坑 |
| `REPORT_RLCD_TNG50.md` | 汇报：项目、步骤、问题、结果、结论、方法原理 |
| `DATA.md`、`data_manifest.tsv` | 不进 GitHub 的大数据：清单、大小、校验和、怎么重新生成 |
| `PlanB/README.md`、`star_condition/*/README.md`、`results/rlcd_A/README.md` | 各子项目的说明和结果要点 |
