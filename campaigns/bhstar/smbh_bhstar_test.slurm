#!/bin/bash
#SBATCH -J smbh_bhstar_test
#SBATCH --partition=short-96core-shared
#SBATCH -A <your-allocation>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 04:00:00
#SBATCH -e /path/to/scratch/smbh_bhstar_bmaxtest/logs/err_%j.e
#SBATCH -o /path/to/scratch/smbh_bhstar_bmaxtest/logs/out_%j.o

set -uo pipefail
: "${SMBH_ROOT:=/path/to/smbh-binary}"   # exported by submit_*.sh
source "$SMBH_ROOT/env.sh"
RUNDIR="${BHSTAR_RUNDIR:-$SMBH_SCRATCH/smbh_bhstar_bmaxtest}"   # override for quick tests (sbatch --export=ALL,BHSTAR_RUNDIR=...)
cd "$RUNDIR" || exit 1

BIN="$SMBH_BIN/smbh_binary_PN"
JOBLIST="$RUNDIR/joblist_shuffled.txt"
OUT="$RUNDIR/smbh_bhstar_results_test.csv"

[ -f "$JOBLIST" ] || { echo "FATAL: $JOBLIST missing (run the submit script)"; exit 1; }
[ -x "$BIN" ]     || { echo "FATAL: binary $BIN missing"; exit 1; }

echo "=== bmax test  job $SLURM_JOB_ID  $(date) on $(hostname) ==="

NLOCAL=$(wc -l < "$JOBLIST")
NCORES=${SLURM_CPUS_PER_TASK}
CHUNK=$(( (NLOCAL + NCORES - 1) / NCORES ))
TMPD=$(mktemp -d "$RUNDIR/.tmp_test_XXXXXX")

"$BIN" --header | grep '^CSV,' > "$OUT"           # header line
echo "running $NLOCAL sims on $NCORES cores, 120s/run timeout"

# fan the full joblist across the node's cores; 1 temp file per worker
seq 0 $((NCORES-1)) | xargs -P "$NCORES" -I{} bash -c '
    k="$1"; bin="$2"; shard="$3"; tmpd="$4"; chunk="$5"
    start=$(( k * chunk + 1 )); end=$(( start + chunk - 1 ))
    out="$tmpd/w_${k}.csv"; tout="$tmpd/to_${k}.txt"
    sed -n "${start},${end}p" "$shard" | while IFS= read -r line; do
        if ! timeout 120 "$bin" $line | grep "^CSV,"; then
            echo "$line" >> "$tout"
        fi
    done >> "$out"
' _ {} "$BIN" "$JOBLIST" "$TMPD" "$CHUNK"

cat "$TMPD"/w_*.csv >> "$OUT" 2>/dev/null
NTO=$(cat "$TMPD"/to_*.txt 2>/dev/null | wc -l)
cat "$TMPD"/to_*.txt > "$RUNDIR/timeouts_test.txt" 2>/dev/null
rm -rf "$TMPD"
echo "=== done $(date): $(($(wc -l < "$OUT") - 1)) CSV rows, $NTO timeouts/crashes ==="
