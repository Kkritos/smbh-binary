#!/bin/bash
#SBATCH -J smbh_bhstar_arr
#SBATCH --partition=long-96core
#SBATCH -A <your-allocation>
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH -t 1-00:00:00
#SBATCH -e /path/to/scratch/smbh_bhstar_encounters/logs/err_%A_%a.e
#SBATCH -o /path/to/scratch/smbh_bhstar_encounters/logs/out_%A_%a.o

# BH-star + SMBH hyperbolic cross-section sweep (full PN) -- JOB ARRAY version.
# Splits the shared joblist_shuffled.txt into NARR contiguous shards (one per array
# task); the joblist is pre-shuffled, so contiguous shards are load-balanced.
# Each task runs its shard across 96 cores (per-run timeout 120s) and writes its own
# results_<task>.csv.  After all tasks finish, concatenate them (see footer).
#
# Set --array=0-(NARR-1) and NARR to the SAME number.  With --array=0-3 / NARR=4 the
# work is quartered -> ~1/4 the wall time, and a dead node loses only its shard.

NARR=8                                    # <-- must match --array=0-(NARR-1) above

: "${SMBH_ROOT:=/path/to/smbh-binary}"   # exported by submit_*.sh
source "$SMBH_ROOT/env.sh"
RUNDIR="$SMBH_SCRATCH/smbh_bhstar_encounters"
cd "$RUNDIR" || exit 1

BIN="$SMBH_BIN/smbh_binary_PN"
JOBLIST="$RUNDIR/joblist_shuffled.txt"
TASK=${SLURM_ARRAY_TASK_ID}
OUT="$RUNDIR/results_${TASK}.csv"

[ -f "$JOBLIST" ] || { echo "FATAL: joblist_shuffled.txt missing (run the submit script)"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: binary $BIN missing"; exit 1; }

echo "=== array task $TASK / $NARR  job $SLURM_ARRAY_JOB_ID  $(date) on $(hostname) ==="

# ---- this task's contiguous shard of the global joblist ----
TOTAL=$(wc -l < "$JOBLIST")
SHARD=$(( (TOTAL + NARR - 1) / NARR ))
START=$(( TASK * SHARD + 1 ))
END=$(( START + SHARD - 1 ))

TMPD=$(mktemp -d "$RUNDIR/.tmp_arr_${TASK}_XXXXXX")
sed -n "${START},${END}p" "$JOBLIST" > "$TMPD/shard.txt"
NLOCAL=$(wc -l < "$TMPD/shard.txt")

"$BIN" --header | grep '^CSV,' > "$OUT"     # each shard file carries its own header

NCORES=${SLURM_CPUS_PER_TASK}
CHUNK=$(( (NLOCAL + NCORES - 1) / NCORES ))
echo "task $TASK: lines ${START}-${END}  ($NLOCAL sims) on $NCORES cores, timeout 120s"

# fan the shard across the node's cores; 1 temp file per worker
seq 0 $((NCORES-1)) | xargs -P "$NCORES" -I{} bash -c '
    k="$1"; bin="$2"; shard="$3"; tmpd="$4"; chunk="$5"
    start=$(( k * chunk + 1 ))
    end=$(( start + chunk - 1 ))
    out="$tmpd/w_${k}.csv"
    tout="$tmpd/to_${k}.txt"
    sed -n "${start},${end}p" "$shard" | while IFS= read -r line; do
        if ! timeout 300 "$bin" $line | grep "^CSV,"; then
            echo "$line" >> "$tout"
        fi
    done >> "$out"
' _ {} "$BIN" "$TMPD/shard.txt" "$TMPD" "$CHUNK"

cat "$TMPD"/w_*.csv >> "$OUT" 2>/dev/null
NTO=$(cat "$TMPD"/to_*.txt 2>/dev/null | wc -l)
cat "$TMPD"/to_*.txt > "$RUNDIR/timeouts_${TASK}.txt" 2>/dev/null
echo "=== task $TASK done $(date): $(($(wc -l < "$OUT") - 1)) CSV rows, $NTO timeouts/crashes ==="

# ---------------------------------------------------------------------------
# AFTER all array tasks finish, merge the shards into one results file:
#   cd "$SMBH_SCRATCH/smbh_bhstar_encounters"
#   head -1 results_0.csv > smbh_bhstar_results.csv
#   for f in results_[0-9]*.csv; do tail -n +2 "$f" >> smbh_bhstar_results.csv; done
#   wc -l smbh_bhstar_results.csv
# ---------------------------------------------------------------------------
