#!/usr/bin/env bash
# unpack_data.sh — unpack the public main-selection data shipped with this repo
#   data/release/rlcd_input_main_{nonoise,noise}.tar.gz  →  data/rlcd_input/{nonoise,noise}/
#
# Usage:  bash tools/unpack_data.sh
# Existing folders (or symlinks made by tools/link_data.sh) are left untouched.

set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

for v in nonoise noise; do
  tgz="data/release/rlcd_input_main_$v.tar.gz"
  out="data/rlcd_input/$v"
  if [ -e "$out" ]; then
    echo "skip   $out (already exists)"; continue
  fi
  [ -f "$tgz" ] || { echo "missing $tgz (did the clone finish?)"; exit 1; }
  tar -xzf "$tgz"
  echo "unpacked $out ($(ls "$out" | wc -l | tr -d ' ') files)"
done

# verify against the checksum list
md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
bad=0; n=0
while IFS=$'\t' read -r path size ref; do
  case "$path" in data/rlcd_input/*) ;; *) continue;; esac
  n=$((n+1))
  [ -f "$path" ] && [ "$(md5of "$path")" = "$ref" ] || { echo "checksum mismatch: $path"; bad=$((bad+1)); }
done < data_manifest.tsv
echo "checked $n files, $bad problems."
