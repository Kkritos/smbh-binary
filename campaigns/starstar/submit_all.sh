#!/bin/bash
set -euo pipefail

# Master: generate joblist -> shuffle once (shared) -> submit both jobs.
# One shuffle so Newt row N and PN row N are the same initial condition.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/env.sh"   # sets SMBH_ROOT, SMBH_BIN, SMBH_SCRATCH, slurm_submit (reads env.local.sh)
RUNDIR="$SMBH_SCRATCH/smbh_starstar_collisions"
JOBDIR="$SMBH_ROOT/campaigns/starstar"
GEN="$JOBDIR/smbh_starstar_collisions_sweep.py"

mkdir -p "$RUNDIR/logs"
cd "$RUNDIR"

echo "generating joblist..."
python3 "$GEN" > joblist.txt
echo "  $(wc -l < joblist.txt) lines -> joblist.txt"

echo "shuffling (once, shared by both jobs)..."
shuf joblist.txt > joblist_shuffled.txt
echo "  $(wc -l < joblist_shuffled.txt) lines -> joblist_shuffled.txt"

echo "submitting jobs..."
slurm_submit "$JOBDIR/smbh_starstar_collisions_Newt_job.slurm" "$RUNDIR/logs" "$SMBH_SLURM_PARTITION"
slurm_submit "$JOBDIR/smbh_starstar_collisions_PN_job.slurm"   "$RUNDIR/logs" "$SMBH_SLURM_PARTITION"

echo "queue:"
squeue -u "$USER"
