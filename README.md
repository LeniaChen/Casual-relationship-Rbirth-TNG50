# TNG50 → Causal Galactic Archaeology（RLCD）验证项目

用 TNG50 模拟检验 Jin et al. 的因果发现方法 RLCD 能不能从 6 个可观测量（年龄、[Fe/H]、[Mg/H]、J_z、J_φ、偏心率）恢复恒星的**出生半径 Rbirth**（以及导引半径 Rg）。TNG50 里每颗星的真实出生半径都知道，可以直接对答案；如果可行，再用到银河系观测数据上。

- 论文：Jin, Lu, Ting, Zheng & Buck,《Causal Discovery of Latent Variables in Galactic Archaeology》，[arXiv:2507.00134](https://arxiv.org/abs/2507.00134)（PDF 本地放在 `docs/paper/`，不进 GitHub）
- 作者代码：<https://github.com/ZehaoJin/Causal-Galactic-Archaeology>（本项目用 commit 8c0a5a4）

> **结论摘要**：按论文的做法（冻结作者的因果图、只重估系数），在 TNG50 上**不能可靠地得到出生半径**——边系数不可识别，隐变量 L1 的含义随优化器初值变化，而且即使取最有利的解也基本不比“单用年龄”或“单用 [Fe/H]”更好。详见 `REPORT_RLCD_TNG50.md` 第 0、6 节。

---

## 0. 先看哪个文件

| 你想… | 看这个 |
|---|---|
| 知道项目进行到哪、有什么问题、下一步做什么 | **`STATUS.md`** |
| 让 Claude 接着做 | **`CLAUDE.md`**（给 AI 的说明；新对话里让它先读 `CLAUDE.md` 和 `STATUS.md`） |
| 理解项目、方法、结果和结论 | `REPORT_RLCD_TNG50.md`（汇报文档） |
| 复现某一步（命令、数据结构、安装记录、踩过的坑） | `docs/pipeline_plan.md` |
| 了解不进 GitHub 的大数据 | `DATA.md`、`data_manifest.tsv` |

## 1. 仓库里有什么、没有什么

**有**（约 60 MB）：全部脚本、文档、日志、汇总表、图、引力势系数、每个选星条件的协方差矩阵（`moments_*.json`）。

**没有**（约 17 GB，存在外置硬盘 `/Volumes/001HardDisk/MW/casuality`）：TNG50 提取数据和轨道数据（hdf5）、所有逐星的 RLCD 输入和真值 csv、作者的 RLCD 仓库。用 `tools/link_data.sh` 以“软链接”（快捷方式）接回仓库里的原位置，脚本不用改。详见 `DATA.md`。

## 2. 快速开始

### 2.1 本机（外置硬盘插着）
```bash
cd ~/Project2026/Casual-relationship-Rbirth-TNG50
git pull                       # 取回 GitHub 上的最新版本
bash tools/link_data.sh        # 把外置硬盘上的大数据接回来（只需做一次；换硬盘路径时重做）
bash tools/check_setup.sh      # 检查数据、作者仓库、两个 Python 环境、Graphviz 是否齐全
```
然后读 `STATUS.md`，或者打开 Claude 让它读 `CLAUDE.md` 和 `STATUS.md` 后继续。

### 2.2 新电脑
```bash
git clone https://github.com/<你的用户名>/Casual-relationship-Rbirth-TNG50.git
cd Casual-relationship-Rbirth-TNG50
bash tools/setup_rlcd_repo.sh                      # 克隆作者仓库（固定版本）
# 按第 5 节建两个 Python 环境
# 数据：把外置硬盘接上后  bash tools/link_data.sh /Volumes/<硬盘名>/MW/casuality
#      或把 casuality 拷到新电脑任意位置后  bash tools/link_data.sh <那个位置>
#      都没有时按 DATA.md 第 3 节重新生成
bash tools/check_setup.sh
```
不接数据也能做的事：读所有结果和图；重跑宽选星条件的可识别性分析（`star_condition/*/output_data/identifiability_moments.py`）；复现论文的图 1、图 3（`results/复现_Jin/`）。

## 3. 同步到 GitHub

### 3.1 第一次（只做一次）
1. 在 GitHub 网页上新建一个空仓库，名字 `Casual-relationship-Rbirth-TNG50`，**建议选 Private（私有）**：里面有还没发表的结果。不要勾选 “Add a README” 之类的选项。
2. 在终端里运行（把 `<你的用户名>` 换掉）：
```bash
cd ~/Project2026/Casual-relationship-Rbirth-TNG50
git init -b main
git add .
git status --short | head -50      # 看一眼：不应该出现 .hdf5、input_data/nonoise、true_value 之类的大数据
git commit -m "Initial import: TNG50 RLCD validation project"
git remote add origin https://github.com/<你的用户名>/Casual-relationship-Rbirth-TNG50.git
git push -u origin main
```
装了 GitHub 命令行工具 `gh` 的话，第 1 步和最后两行可以换成一行：`gh repo create Casual-relationship-Rbirth-TNG50 --private --source=. --push`。

> 注意：上一级文件夹 `~/Project2026` 本身也是一个 git 仓库。以后在 `~/Project2026` 里运行 `git add .` 时，git 会提示 “embedded git repository”（里面嵌着另一个仓库）。如果不想让上一级仓库管这个文件夹，在 `~/Project2026/.gitignore` 里加一行 `Casual-relationship-Rbirth-TNG50/`。

### 3.2 日常
```bash
cd ~/Project2026/Casual-relationship-Rbirth-TNG50
git pull                                  # 开始工作前：取回最新
git add -A && git commit -m "说明这次改了什么" && git push   # 工作结束后：上传
```
- 新生成的大文件（> 50 MB）不要提交。`git status` 里如果出现大的 csv/hdf5，把它的路径加进 `.gitignore`（以及 `tools/link_data.sh` 的列表、`data_manifest.tsv`）。
- 推进了研究进度后，记得更新 `STATUS.md` 和 `docs/pipeline_plan.md` 第 2 节（让 Claude 做也行）。

## 4. 文件夹结构
```
├── README.md  STATUS.md  CLAUDE.md  DATA.md  data_manifest.tsv
├── REPORT_RLCD_TNG50.md            汇报文档
├── docs/pipeline_plan.md           复现细节（主文档）；docs/paper/ 放论文 PDF（不进 GitHub）
├── tools/                          link_data.sh、check_setup.sh、setup_rlcd_repo.sh
├── setup/                          两个 Python 环境的依赖清单
├── scripts/
│   ├── server/                     在 TNG JupyterLab 上运行：data_obtain.py、merge_hdf5.py
│   ├── local/                      在 Mac 上运行：轨道、选星、做法 A、可识别性分析、线性上限、比较因果图……
│   └── reference/                  旧参考脚本（只读）
├── data/
│   ├── tng50_extract/              [不进 GitHub] data_obtain_all.hdf5
│   ├── tng50_derived/              potentials/*.coef（引力势系数）+ [不进 GitHub] orbits_all.hdf5
│   └── rlcd_input/{nonoise,noise}/ [不进 GitHub] 7–10 kpc 主选星的 RLCD 输入和真值
├── results/
│   ├── check_vs_001_orbit.txt、compute_orbits.log、rlcd_input_selection.txt、selection_preview.txt
│   ├── orbit_check/                轨道结果检查（汇总 + 图）
│   ├── rlcd_g279/                  复现作者的因果图（4 组 alpha、5 个种子）
│   ├── rlcd_A/                     做法 A；identifiability/ 是可识别性分析（主结果：summary_ci_*）
│   ├── 复现_Jin/{因果图,真值对比}/   复现论文图 1、图 3
│   └── TNG/{因果图,真值对比}/        TNG50 12 个星系的对应图
├── star_condition/                 放宽选星：R_0_30_Z_30、R_0_30_Z_3、R_0_30_Z_0.3（各有 README）；R_7_10（主选星的图）
├── PlanB/                          做法 B：每个星系 30 个种子从头学因果图（README 有运行顺序和结果）
└── Causal-Galactic-Archaeology/    [不进 GitHub] 作者的 RLCD 仓库
```

## 5. 计算环境（用户的 Mac，Apple Silicon）

| 环境 | 用途 | 安装 |
|---|---|---|
| `~/PYTHONPATH_01/venv`（Python 3.11.7） | 绝大多数脚本：轨道、选星、可识别性分析、画图 | `~/PYTHONPATH_01/venv/bin/python3 -m pip install -r setup/requirements_venv_main.txt`；AGAMA 从源码编译（`docs/pipeline_plan.md` 第 8.2 节） |
| `~/PYTHONPATH_01/venv_rlcd`（Python 3.12） | 运行作者的 RLCD：`rlcd_A.py`、`compare_rlcd_graphs.py`、PlanB 的 `rlcd_B.py`、`results/复现_Jin/因果图/` | 见下 |
| Graphviz | 把因果图画成图片 | `brew install graphviz` |

```bash
brew install python@3.12 graphviz
/opt/homebrew/bin/python3.12 -m venv ~/PYTHONPATH_01/venv_rlcd
~/PYTHONPATH_01/venv_rlcd/bin/python3 -m pip install --upgrade pip
~/PYTHONPATH_01/venv_rlcd/bin/python3 -m pip install -r setup/requirements_rlcd_mac.txt
```
终端默认在 anaconda base，所以命令里总是写 venv 的完整路径（例如 `~/PYTHONPATH_01/venv/bin/python3`），不要用 `python3`。

## 6. 主要流程（都在仓库根目录下运行）
```
TNG 服务器                   本地（venv）                                         本地（venv_rlcd）
data_obtain.py ─► compute_orbits.py ─► make_rlcd_input.py ─┬─► rlcd_identifiability_ci.py（做法 A，主结果）
(年龄/丰度/出生半径)  (势/作用量/Rg/e)    (选星 + 6 列输入 + 真值)  ├─► star_condition/*（放宽选星）
                                                          └─► PlanB/*（做法 B，需要 venv_rlcd）
```
每一步的完整命令、参数、输出说明见 `docs/pipeline_plan.md` 第 6 节和 `REPORT_RLCD_TNG50.md` 第 3、5.7 节；子项目的运行顺序见各自的 README。脚本注释里写的 `cd /Volumes/001HardDisk/MW/casuality` 是打包前的旧位置，现在换成仓库根目录即可（脚本都按相对路径找文件）。
