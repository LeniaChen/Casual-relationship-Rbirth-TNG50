#!/usr/bin/env bash
# setup_rlcd_repo.sh —— 准备作者的 RLCD 仓库（Causal-Galactic-Archaeology，约 260 MB）
#
# 本机有外置硬盘时，tools/link_data.sh 已经把它链接进来了，不用运行这个脚本。
# 新电脑上没有外置硬盘时运行：从 GitHub 克隆，并切换到本项目用的版本（commit 8c0a5a4，2026-08-22）。
#
# 用法：bash tools/setup_rlcd_repo.sh

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DST="$REPO/Causal-Galactic-Archaeology"
COMMIT=8c0a5a46770024853a5d31decf35e15f3ab81ca2

if [ -e "$DST" ]; then
  echo "已存在：$DST（如果是软链接，说明已经链接到外置硬盘上的副本）"
  git -C "$DST" log -1 --format='当前版本：%h %ad %s' 2>/dev/null || true
  exit 0
fi
git clone https://github.com/ZehaoJin/Causal-Galactic-Archaeology.git "$DST"
git -C "$DST" checkout -q "$COMMIT" && echo "已切换到本项目用的版本 ${COMMIT:0:7}"
echo
echo "注意：在仓库的 rlcd/ 下运行 RLCD，结果会写回 rlcd/results_*/，可能覆盖作者提交的文件。"
echo "跑完要把新文件移到 results/ 下，作者的原版用 git -C Causal-Galactic-Archaeology checkout -- <路径> 恢复。"
