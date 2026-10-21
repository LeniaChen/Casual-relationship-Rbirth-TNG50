#!/usr/bin/env bash
# link_data.sh — link the full, non-public data (~17 GB) into this repo with symbolic links
#
# Only needed if you have the full data (TNG50 extraction, orbits, wide selections). Everyone else
# can use the public main-selection data: bash tools/unpack_data.sh
#
# Usage:
#   bash tools/link_data.sh /path/to/full/data      # folder with the same layout as this repo
#   (or set RBIRTH_DATA_ROOT in .env, see env.example)
#
# Rules: items missing in the data folder are skipped; existing symlinks are re-pointed;
# real files/folders in the repo are never overwritten. Nothing is deleted or moved.

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ -f "$REPO/.env" ] && set -a && . "$REPO/.env" && set +a
DATA_ROOT="${1:-${RBIRTH_DATA_ROOT:-}}"
DATA_ROOT="${DATA_ROOT/#\~/$HOME}"

if [ -z "$DATA_ROOT" ] || [ ! -d "$DATA_ROOT" ]; then
  echo "Data folder not found: '${DATA_ROOT}'"
  echo "Usage: bash tools/link_data.sh /path/to/full/data   (or set RBIRTH_DATA_ROOT in .env)"
  exit 1
fi
DATA_ROOT="$(cd "$DATA_ROOT" && pwd)"
if [ "$(cd "$DATA_ROOT" && pwd -P)" = "$(cd "$REPO" && pwd -P)" ]; then
  echo "The data folder is the repo itself; nothing to link."; exit 0
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
  src="$DATA_ROOT/$it"; dst="$REPO/$it"
  if [ ! -e "$src" ]; then echo "skip (not in data folder): $it"; n_skip=$((n_skip+1)); continue; fi
  if [ -L "$dst" ]; then ln -sfn "$src" "$dst"; echo "relink: $it"; n_link=$((n_link+1)); continue; fi
  if [ -e "$dst" ]; then echo "keep (real data already here): $it"; n_keep=$((n_keep+1)); continue; fi
  mkdir -p "$(dirname "$dst")"; ln -s "$src" "$dst"; echo "link: $it -> $src"; n_link=$((n_link+1))
done
echo
echo "Done: linked $n_link, skipped $n_skip, kept $n_keep."
echo "Next: bash tools/check_setup.sh   (add --md5 to verify checksums)"
