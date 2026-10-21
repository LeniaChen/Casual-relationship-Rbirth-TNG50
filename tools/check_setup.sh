#!/usr/bin/env bash
# check_setup.sh — read-only check of data, the RLCD repo, Python environments and Graphviz
#
# Usage:
#   bash tools/check_setup.sh          # quick (existence + size)
#   bash tools/check_setup.sh --md5    # also verify checksums in data_manifest.tsv

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
[ -f .env ] && set -a && . ./.env && set +a
PY="${PY:-python3}"; PY="${PY/#\~/$HOME}"
PY_RLCD="${PY_RLCD:-$PY}"; PY_RLCD="${PY_RLCD/#\~/$HOME}"
DO_MD5=0; [ "${1:-}" = "--md5" ] && DO_MD5=1

ok=0; bad=0
pass() { echo "  [OK]      $*"; ok=$((ok+1)); }
fail() { echo "  [MISSING] $*"; bad=$((bad+1)); }
note() { echo "  [info]    $*"; }
filesize() { wc -c < "$1" | tr -d ' '; }
md5of()    { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

# $1 = path prefix regex, prints "present total wrong"
check_group() {
  local present=0 total=0 wrong=0
  while IFS=$'\t' read -r path size ref; do
    [ "$path" = "path" ] && continue
    echo "$path" | grep -Eq "$1" || continue
    total=$((total+1))
    [ -f "$path" ] || continue
    present=$((present+1))
    if [ "$(filesize "$path")" != "$size" ]; then wrong=$((wrong+1)); echo "    size mismatch: $path" >&2; continue; fi
    if [ $DO_MD5 -eq 1 ] && [ "$(md5of "$path")" != "$ref" ]; then wrong=$((wrong+1)); echo "    checksum mismatch: $path" >&2; fi
  done < data_manifest.tsv
  echo "$present $total $wrong"
}

echo "== 1. Public main-selection data (data/rlcd_input; needed for the main analysis)"
read -r p t w < <(check_group '^data/rlcd_input/')
if [ "$p" = "$t" ] && [ "$w" = 0 ]; then pass "data/rlcd_input: $p/$t files$([ $DO_MD5 -eq 1 ] && echo ', checksums OK')"
else fail "data/rlcd_input: $p/$t present, $w mismatched -> run: bash tools/unpack_data.sh"; fi

echo "== 2. Full data (not public; optional — see DATA.md)"
for g in '^data/tng50_' '^PlanB/' '^star_condition/R_0_30_Z_30/' '^star_condition/R_0_30_Z_3/' '^star_condition/R_0_30_Z_0\.3/'; do
  read -r p t w < <(check_group "$g")
  label=$(echo "$g" | sed 's/[\^\\]//g')
  if [ "$p" = "$t" ] && [ "$w" = 0 ]; then pass "$label: $p/$t files"
  elif [ "$p" = 0 ]; then note "$label: not present (optional)"
  else fail "$label: $p/$t present, $w mismatched"; fi
done

echo "== 3. Author's RLCD repository"
if [ -d Causal-Galactic-Archaeology/rlcd ]; then
  c=$(git -C Causal-Galactic-Archaeology rev-parse --short HEAD 2>/dev/null || echo "?")
  pass "Causal-Galactic-Archaeology (commit $c; this project used 8c0a5a4)"
else
  fail "Causal-Galactic-Archaeology -> run: bash tools/setup_rlcd_repo.sh"
fi

echo "== 4. Python environments (set PY / PY_RLCD in .env, see env.example)"
if v=$("$PY" -c 'import sys,numpy,scipy,pandas,h5py,matplotlib;print(sys.version.split()[0])' 2>/dev/null); then
  pass "analysis env PY=$PY (Python $v)"
else
  fail "analysis env PY=$PY: needs numpy scipy pandas h5py matplotlib (pip install -r setup/requirements.txt)"
fi
"$PY" -c 'import agama' 2>/dev/null && pass "AGAMA (only needed to recompute orbits)" \
  || note "AGAMA not found: only needed for scripts/local/compute_orbits.py (see docs/INSTALL_AGAMA.md)"
if "$PY_RLCD" -c 'import causallearn,torch' 2>/dev/null; then pass "RLCD env PY_RLCD=$PY_RLCD (causal-learn, torch)"
else fail "RLCD env PY_RLCD=$PY_RLCD: Python 3.12 + pip install -r setup/requirements-rlcd.txt"; fi

echo "== 5. Other tools"
command -v dot >/dev/null 2>&1 && pass "Graphviz (dot)" || fail "Graphviz: brew install graphviz / apt install graphviz"

echo
echo "$ok OK, $bad missing."
[ $bad -eq 0 ] && echo "Ready. Start with README.md, then STATUS.md."
