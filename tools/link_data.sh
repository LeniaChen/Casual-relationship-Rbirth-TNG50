#!/usr/bin/env bash
# link_data.sh —— 把不进 GitHub 的大数据（约 17 GB）以“软链接”的方式接回本仓库
#
# 为什么要这样做：所有脚本都按“仓库根目录下的相对路径”找数据（例如 data/rlcd_input/nonoise/），
# 而大数据只存在外置硬盘的 casuality 文件夹里。软链接就像一个“快捷方式”：仓库里的
# data/rlcd_input/nonoise 实际指向外置硬盘上的同名文件夹，脚本不用改就能运行，数据也不会复制一份。
#
# 用法（在终端里，任意目录都行）：
#   bash tools/link_data.sh                                   # 默认数据目录 /Volumes/001HardDisk/MW/casuality
#   bash tools/link_data.sh /path/to/另一个数据目录             # 换了硬盘或电脑时，指定数据所在的目录
#   CASUALITY_DATA_ROOT=/path/to/data bash tools/link_data.sh  # 也可以用环境变量指定
#
# 规则：
#   - 数据目录里没有的项目会跳过（打印“跳过”），不会报错；
#   - 仓库里已经是软链接的，会改成指向新的数据目录；
#   - 仓库里已经有真实文件/文件夹的，不覆盖（打印“已存在真实数据”）。
#   - 这个脚本只创建软链接，不删除、不移动任何数据。

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${1:-${CASUALITY_DATA_ROOT:-/Volumes/001HardDisk/MW/casuality}}"

if [ ! -d "$DATA_ROOT" ]; then
  echo "找不到数据目录：$DATA_ROOT"
  echo "外置硬盘插上了吗？或者用：bash tools/link_data.sh <数据所在目录>"
  exit 1
fi
DATA_ROOT="$(cd "$DATA_ROOT" && pwd)"
if [ "$(cd "$DATA_ROOT" && pwd -P)" = "$(cd "$REPO" && pwd -P)" ]; then
  echo "数据目录就是仓库本身，不需要链接。"
  exit 0
fi

ITEMS=(
  Causal-Galactic-Archaeology
  data/tng50_extract
  data/tng50_derived/orbits_all.hdf5
  data/tng50_derived/orbits_test.hdf5
  data/rlcd_input/nonoise
  data/rlcd_input/noise
  PlanB/input_data/nonoise
  PlanB/input_data/noise
  PlanB/true_value
)
for c in R_0_30_Z_30 R_0_30_Z_3 R_0_30_Z_0.3; do
  ITEMS+=("star_condition/$c/input_data/nonoise" "star_condition/$c/input_data/noise" "star_condition/$c/true_value")
done

n_link=0; n_skip=0; n_keep=0
for it in "${ITEMS[@]}"; do
  src="$DATA_ROOT/$it"
  dst="$REPO/$it"
  if [ ! -e "$src" ]; then
    echo "跳过（数据目录里没有）：$it"; n_skip=$((n_skip+1)); continue
  fi
  if [ -L "$dst" ]; then
    ln -sfn "$src" "$dst"; echo "更新链接：$it"; n_link=$((n_link+1)); continue
  fi
  if [ -e "$dst" ]; then
    echo "已存在真实数据，不覆盖：$it"; n_keep=$((n_keep+1)); continue
  fi
  mkdir -p "$(dirname "$dst")"
  ln -s "$src" "$dst"
  echo "链接：$it -> $src"; n_link=$((n_link+1))
done

echo
echo "完成：链接 $n_link 项，跳过 $n_skip 项，保留真实数据 $n_keep 项。"
echo "下一步可以运行：bash tools/check_setup.sh   （加 --md5 会逐个核对校验和，约 1 分钟）"
