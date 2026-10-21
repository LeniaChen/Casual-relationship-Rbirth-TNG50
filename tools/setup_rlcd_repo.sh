#!/usr/bin/env bash
# setup_rlcd_repo.sh — clone the authors' RLCD code (Causal-Galactic-Archaeology, ~260 MB)
# at the exact commit used in this project (8c0a5a4, 2026-08-22).
#
# Usage: bash tools/setup_rlcd_repo.sh

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DST="$REPO/Causal-Galactic-Archaeology"
COMMIT=8c0a5a46770024853a5d31decf35e15f3ab81ca2

if [ -e "$DST" ]; then
  echo "Already present: $DST"
  git -C "$DST" log -1 --format='current commit: %h %ad %s' 2>/dev/null || true
  exit 0
fi
git clone https://github.com/ZehaoJin/Causal-Galactic-Archaeology.git "$DST"
git -C "$DST" checkout -q "$COMMIT" && echo "Checked out ${COMMIT:0:7}"
echo
echo "Note: running RLCD inside its rlcd/ folder writes to rlcd/results_*/ and can overwrite the authors'"
echo "committed results. Move new outputs to results/ and restore originals with"
echo "  git -C Causal-Galactic-Archaeology checkout -- <path>"
