#!/bin/bash
#SBATCH -J smbh_Newt_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 00:30:00
# Site-specific settings (log paths, partition, account) are NOT set here: the
# submit scripts pass them to sbatch from env.local.sh (see env.local.sh.example).
# Adjust --cpus-per-task / -t to your nodes if needed.

# TEST (Newt): first 1000 lines of joblist.txt on 40 cores (a few minutes).
# PREREQ: joblist.txt and logs/ exist in scratch; smbh_binary_Newt compiled.

echo "=== smbh Newt test $SLURM_JOB_ID started $(date) on $(hostname) ==="

: "${SMBH_ROOT:?SMBH_ROOT not set -- submit through the campaign submit_*.sh scripts}"
source "$SMBH_ROOT/env.sh"
RUNDIR="$SMBH_SCRATCH/smbh_starstar_collisions"
cd "$RUNDIR" || exit 1

BIN="$SMBH_BIN/smbh_binary_Newt"
JOBLIST="$RUNDIR/joblist.txt"
OUT="$RUNDIR/smbh_test_Newt_results.csv"
NLINES=1000

[ -f "$JOBLIST" ] || { echo "FATAL: joblist.txt missing"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: binary $BIN missing"; exit 1; }

TMPD=$(mktemp -d "$RUNDIR/.tmp_test_Newt_XXXXXX")
trap 'rm -rf "$TMPD"' EXIT

"$BIN" --header | grep '^CSV,' > "$OUT"

echo "running first $NLINES sims on ${SLURM_CPUS_PER_TASK} cores..."
head -n "$NLINES" "$JOBLIST" | \
xargs -P "${SLURM_CPUS_PER_TASK}" -d '\n' -I{} bash -c '
    f="$2/r_${BASHPID}_${RANDOM}.csv"
    "$0" $1 --report 0 | grep "^CSV," >> "$f"
' "$BIN" {} "$TMPD"

cat "$TMPD"/*.csv >> "$OUT" 2>/dev/null
echo "=== finished $(date): $(($(wc -l < "$OUT") - 1)) CSV rows in $OUT (expected $NLINES) ==="
