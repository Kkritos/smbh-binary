#!/bin/bash
#SBATCH -J smbh_PN_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH -t 1-00:00:00
# Site-specific settings (log paths, partition, account) are NOT set here: the
# submit scripts pass them to sbatch from env.local.sh (see env.local.sh.example).
# Adjust --cpus-per-task / -t to your nodes if needed.

# PN sweep: reads the shared joblist_shuffled.txt (made by submit_all.sh).
# Chunked (1 temp file/worker) + per-run timeout 180s (stuck runs logged, not categorized).

echo "=== smbh PN sweep $SLURM_JOB_ID started $(date) on $(hostname) ==="

: "${SMBH_ROOT:?SMBH_ROOT not set -- submit through the campaign submit_*.sh scripts}"
source "$SMBH_ROOT/env.sh"
RUNDIR="$SMBH_SCRATCH/smbh_starstar_collisions"
cd "$RUNDIR" || exit 1

BIN="$SMBH_BIN/smbh_binary_PN"
JOBLIST="$RUNDIR/joblist_shuffled.txt"
OUT="$RUNDIR/smbh_starstar_collisions_PN_results.csv"

[ -f "$JOBLIST" ] || { echo "FATAL: joblist_shuffled.txt missing (run submit_all.sh)"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: binary $BIN missing"; exit 1; }

TMPD=$(mktemp -d "$RUNDIR/.tmp_PN_XXXXXX")

"$BIN" --header | grep '^CSV,' > "$OUT"

NCORES=${SLURM_CPUS_PER_TASK}
TOTAL=$(wc -l < "$JOBLIST")
CHUNK=$(( (TOTAL + NCORES - 1) / NCORES ))

echo "running $TOTAL sims on $NCORES cores (chunked+shuffled, timeout 180s)..."

seq 0 $((NCORES-1)) | xargs -P "$NCORES" -I{} bash -c '
    k="$1"; bin="$2"; joblist="$3"; tmpd="$4"; chunk="$5"
    start=$(( k * chunk + 1 ))
    end=$(( start + chunk - 1 ))
    out="$tmpd/worker_${k}.csv"
    tout="$tmpd/timeouts_${k}.txt"
    sed -n "${start},${end}p" "$joblist" | while IFS= read -r line; do
        if ! timeout 180 "$bin" $line --report 0 | grep "^CSV,"; then
            echo "$line" >> "$tout"
        fi
    done >> "$out"
' _ {} "$BIN" "$JOBLIST" "$TMPD" "$CHUNK"

cat "$TMPD"/*.csv >> "$OUT" 2>/dev/null
NTO=$(cat "$TMPD"/timeouts_*.txt 2>/dev/null | wc -l)
echo "=== finished $(date): $(($(wc -l < "$OUT") - 1)) CSV rows, $NTO timeouts ==="
