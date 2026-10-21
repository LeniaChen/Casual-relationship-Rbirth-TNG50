# 大数据说明（不进 GitHub 的部分）

GitHub 单个文件上限 100 MB，仓库也不适合放十几 GB 的数据，所以下面这些文件**只存在外置硬盘**的原项目文件夹 `/Volumes/001HardDisk/MW/casuality` 里（路径与本仓库完全相同）。仓库里只放代码、文档、日志、汇总结果和图（约 60 MB）。

- 完整清单：`data_manifest.tsv`（243 个文件；三列：相对路径、字节数、MD5 校验和）。MD5 是文件的“指纹”：文件内容哪怕只变一个字节，指纹就会不同，用来确认拷贝或下载的数据完好。
- 接回仓库：`bash tools/link_data.sh`（软链接，不复制）
- 检查齐全：`bash tools/check_setup.sh`，加 `--md5` 核对指纹

## 1. 清单（按重要性）

| 路径 | 文件数 | 大小 | 怎么来的 | 丢了能否重做 |
|---|---|---|---|---|
| `data/tng50_extract/data_obtain_all.hdf5` | 1 | 955 MB | TNG 服务器上 `scripts/server/data_obtain.py` + `merge_hdf5.py`，再下载 | **最难重做**：要 TNG JupyterLab 服务器、轨迹文件 `TNG50_disk_trajectory.hdf5`，跑几个小时。**建议另外备份一份** |
| `data/tng50_derived/orbits_all.hdf5` | 1 | 271 MB | 本地 `scripts/local/compute_orbits.py`（AGAMA，12 个星系约 1 小时） | 能，需要上一个文件和 AGAMA |
| `data/tng50_derived/orbits_test.hdf5` | 1 | 0.2 MB | `compute_orbits.py --test`（旧设置，可以删） | 不需要 |
| `data/rlcd_input/{nonoise,noise}/` | 48 | 245 MB | `scripts/local/make_rlcd_input.py`（约 2 分钟） | 能，需要上面两个 hdf5 |
| `PlanB/input_data/{nonoise,noise}/`、`PlanB/true_value/` | 48 | 245 MB | `PlanB/input_data/make_input.py`（从 `data/rlcd_input/` 原样复制） | 能 |
| `star_condition/R_0_30_Z_30/{input_data/*,true_value}/` | 48 | 6.3 GB | `star_condition/R_0_30_Z_30/input_data/make_input.py`（读两个 hdf5） | 能 |
| `star_condition/R_0_30_Z_3/...` | 48 | 6.0 GB | 同名 `make_input.py`（从 R_0_30_Z_30 筛出 \|z\| < 3 kpc） | 能，需要 R_0_30_Z_30 |
| `star_condition/R_0_30_Z_0.3/...` | 48 | 2.8 GB | 同名 `make_input.py`（从 R_0_30_Z_30 筛出 \|z\| < 0.3 kpc） | 能，需要 R_0_30_Z_30 |
| `Causal-Galactic-Archaeology/` | — | 257 MB | 作者仓库的 git clone（commit 8c0a5a4，未改动） | 能：`bash tools/setup_rlcd_repo.sh` |

合计约 16.9 GB（不含作者仓库）。

**说明**：
- `nonoise` = 不加观测误差（主版本）；`noise` = 按论文加误差：年龄 ×(1 + 0.10·N)、[Fe/H] + 0.02·N、[Mg/Fe] + 0.06·N。
  - 意思是：给每颗星的年龄加 10% 的相对误差，[Fe/H] 加 0.02 dex、[Mg/Fe] 加 0.06 dex 的绝对误差。
  - N：每颗星、每个量各抽一次的标准正态随机数（平均 0、标准差 1）；随机种子 = 20260926 + 星系 ID，所以重做会得到同一组误差。
- `<ID>_e.csv`：RLCD 的 6 列输入 `logJZ, JPhi, ecc, age, FeH, MgH`；`<ID>_truth.csv`：同一批星（逐行对应）的真值，出生半径是 `R_bir`、导引半径是 `Rg`。各列含义见 `docs/pipeline_plan.md` 第 6.5 节。
- `star_condition/R_0_30_Z_3`、`R_0_30_Z_0.3` 的 noise 版本沿用 R_0_30_Z_30 里已经加好的误差；如果改用 `make_rlcd_input.py` 从头选星，星是同一批，但抽到的误差会不同。

## 2. 哪些分析不需要大数据
- `star_condition/*/output_data/identifiability_moments.py`：只读已提交的 `moments_<版本>.json`（每个星系 8 个量的协方差矩阵），可直接重跑宽选星条件下的全部可识别性分析。
- `results/复现_Jin/` 两个脚本：只需要作者仓库。
- 所有汇总表、图、日志都已在仓库里，读结果不需要数据。

## 3. 完全没有数据时，按这个顺序重新生成
```bash
# 0) TNG JupyterLab 服务器上（旧版 AnastrisTNG），详见 docs/pipeline_plan.md 第 6.1 节
cd ~/MW/casuality
ln -s ~/423_sample/TNG50_disk_trajectory.hdf5 .            # 第一次运行前做一次
nohup python3 data_obtain.py > data_obtain.log 2>&1 &
python3 merge_hdf5.py                                        # → data_obtain_all.hdf5，下载到本仓库 data/tng50_extract/

# 1) 本地 Mac，在仓库根目录下
PY=~/PYTHONPATH_01/venv/bin/python3
caffeinate -i $PY scripts/local/compute_orbits.py 2>&1 | tee results/compute_orbits.log   # 需要 AGAMA，约 1 小时
$PY scripts/local/make_rlcd_input.py
$PY PlanB/input_data/make_input.py
$PY star_condition/R_0_30_Z_30/input_data/make_input.py
$PY star_condition/R_0_30_Z_3/input_data/make_input.py
$PY star_condition/R_0_30_Z_0.3/input_data/make_input.py

# 2) 核对：和原来的指纹比
bash tools/check_setup.sh --md5
```
重新生成的文件不一定和原来逐字节相同（例如 hdf5 里可能记有生成时间、浮点运算顺序不同），这时 `--md5` 会报“不符”；那就改看数值是否一致，例如重跑 `scripts/local/check_orbits.py`、比较 `results/rlcd_input_selection.txt` 里每一步的星数。
`caffeinate -i`：运行期间不让 Mac 休眠（合盖仍会停，要开盖、接电源）；`2>&1 | tee 文件`：屏幕上显示输出的同时写进日志。
