#!/bin/bash
# ============================================================================
# run_simulations.sh — local parallel smbh_binary sweep for macOS (M-series).
# Runs a joblist across N cores with xargs, collecting all CSV lines into one
# file under Results/. Uses per-worker temp files (no interleaving).
#
# Usage:
#   ./run_local.sh <binary> <joblist> <output.csv> [ncores]
# Example:
#   ./run_local.sh ./smbh_binary_PN joblist.txt results.csv 6
#
# Default ncores = 6 (the M5 Pro's performance cores), leaving the 12 efficiency
# cores + headroom free so the laptop stays responsive. Bump to 10 if idle.
# ============================================================================

set -euo pipefail

BIN="${1:?usage: run_local.sh <binary> <joblist> <output.csv> [ncores]}"
JOBLIST="${2:?need joblist}"
OUT="${3:?need output csv}"
NCORES="${4:-6}"

[ -x "$BIN" ]     || { echo "FATAL: binary $BIN not executable"; exit 1; }
[ -f "$JOBLIST" ] || { echo "FATAL: joblist $JOBLIST not found"; exit 1; }

# put outputs under Results/
mkdir -p Results
OUT="Results/$(basename "$OUT")"

TMPD=$(mktemp -d ./.tmp_local_XXXXXX)
trap 'rm -rf "$TMPD"' EXIT

# header once
"$BIN" --header | grep '^CSV,' > "$OUT"

NTOTAL=$(wc -l < "$JOBLIST" | tr -d ' ')
echo "running $NTOTAL sims on $NCORES cores -> $OUT"

# macOS xargs: no -d flag. Use -I{} with the whole line; -L1 = one line per call.
# Each worker writes its CSV line to a unique temp file (PID+RANDOM) -> no interleaving.
# -P $NCORES = parallelism.
cat "$JOBLIST" | xargs -P "$NCORES" -I{} bash -c '
    f="$2/r_$$_${RANDOM}.csv"
    "$0" $1 --report 0 | grep "^CSV," >> "$f"
' "$BIN" {} "$TMPD"

cat "$TMPD"/*.csv >> "$OUT" 2>/dev/null || true
NDONE=$(( $(wc -l < "$OUT" | tr -d ' ') - 1 ))
echo "done: $NDONE CSV rows in $OUT"
