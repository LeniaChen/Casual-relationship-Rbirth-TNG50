#!/usr/bin/env bash
# check_setup.sh —— 检查这台电脑能不能直接接着做本项目：数据、作者仓库、两个 Python 环境、Graphviz
#
# 用法：
#   bash tools/check_setup.sh          # 快速检查（只看文件在不在、大小对不对）
#   bash tools/check_setup.sh --md5    # 另外逐个核对 data_manifest.tsv 里的校验和（约 17 GB，1 分钟左右）
#
# 只读检查，不会修改任何文件。

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
DO_MD5=0; [ "${1:-}" = "--md5" ] && DO_MD5=1

ok=0; bad=0
pass() { echo "  [OK]   $*"; ok=$((ok+1)); }
fail() { echo "  [缺少] $*"; bad=$((bad+1)); }

filesize() { wc -c < "$1" | tr -d ' '; }   # 文件字节数（macOS 和 Linux 通用）
md5of()    { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

echo "== 1. 大数据（清单：data_manifest.tsv）"
missing=0; wrong=0; checked=0
while IFS=$'\t' read -r path size md5sum_ref; do
  [ "$path" = "path" ] && continue
  if [ ! -f "$path" ]; then missing=$((missing+1)); continue; fi
  s=$(filesize "$path")
  if [ "$s" != "$size" ]; then echo "    大小不符：$path（应为 $size，实际 $s）"; wrong=$((wrong+1)); continue; fi
  if [ $DO_MD5 -eq 1 ]; then
    m=$(md5of "$path")
    [ "$m" != "$md5sum_ref" ] && { echo "    校验和不符：$path"; wrong=$((wrong+1)); continue; }
  fi
  checked=$((checked+1))
done < data_manifest.tsv
total=$(($(wc -l < data_manifest.tsv) - 1))
if [ $missing -eq 0 ] && [ $wrong -eq 0 ]; then
  pass "数据 $checked/$total 个文件齐全$([ $DO_MD5 -eq 1 ] && echo '，校验和一致')"
else
  fail "数据：缺 $missing 个、不符 $wrong 个（共 $total 个）。先运行 bash tools/link_data.sh；完全没有数据时见 DATA.md 第 3 节重新生成"
fi

echo "== 2. 作者的 RLCD 仓库"
if [ -d Causal-Galactic-Archaeology/rlcd ]; then
  c=$(git -C Causal-Galactic-Archaeology rev-parse --short HEAD 2>/dev/null || echo "?")
  pass "Causal-Galactic-Archaeology（commit $c；本项目用的是 8c0a5a4）"
else
  fail "Causal-Galactic-Archaeology：运行 bash tools/setup_rlcd_repo.sh"
fi

echo "== 3. Python 环境"
PY_MAIN="$HOME/PYTHONPATH_01/venv/bin/python3"
PY_RLCD="$HOME/PYTHONPATH_01/venv_rlcd/bin/python3"
if [ -x "$PY_MAIN" ]; then
  v=$("$PY_MAIN" -c 'import sys,numpy,scipy,pandas,h5py,matplotlib;print(sys.version.split()[0])' 2>/dev/null) \
    && pass "主环境 venv（Python $v，numpy/scipy/pandas/h5py/matplotlib 可用）" \
    || fail "主环境 venv 存在，但缺包：见 setup/requirements_venv_main.txt"
  "$PY_MAIN" -c 'import agama' 2>/dev/null && pass "AGAMA（只有重算轨道时才需要）" \
    || echo "  [提示] 主环境里没有 AGAMA：只有重算 compute_orbits.py 时才需要，安装见 docs/pipeline_plan.md 第 8.2 节"
else
  fail "主环境 $PY_MAIN：安装见 README 第 5 节"
fi
if [ -x "$PY_RLCD" ]; then
  "$PY_RLCD" -c 'import causallearn,torch' 2>/dev/null && pass "RLCD 环境 venv_rlcd（causal-learn、torch 可用）" \
    || fail "venv_rlcd 存在，但缺包：pip install -r setup/requirements_rlcd_mac.txt"
else
  fail "RLCD 环境 $PY_RLCD：安装见 README 第 5 节"
fi

echo "== 4. 其他工具"
command -v dot >/dev/null 2>&1 && pass "Graphviz（dot）" || fail "Graphviz：brew install graphviz"
command -v git >/dev/null 2>&1 && pass "git" || fail "git"

echo
echo "通过 $ok 项，缺少 $bad 项。"
[ $bad -eq 0 ] && echo "可以直接开始。先读 STATUS.md（进度、问题、下一步）。"
