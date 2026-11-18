#!/bin/bash
set -euo pipefail

# Test: generate the joblist and submit the two short test jobs (first 1000 runs, Newt + PN).

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/env.sh"   # sets SMBH_ROOT, SMBH_BIN, SMBH_SCRATCH, slurm_submit (reads env.local.sh)
RUNDIR="$SMBH_SCRATCH/smbh_starstar_collisions"
JOBDIR="$SMBH_ROOT/campaigns/starstar"
GEN="$JOBDIR/smbh_starstar_collisions_sweep.py"

mkdir -p "$RUNDIR/logs"
cd "$RUNDIR"

echo "generating joblist..."
python3 "$GEN" > joblist.txt
echo "  $(wc -l < joblist.txt) lines -> joblist.txt"

echo "submitting test jobs..."
slurm_submit "$JOBDIR/smbh_starstar_collisions_Newt_test.slurm" "$RUNDIR/logs" "$SMBH_SLURM_PARTITION_TEST"
slurm_submit "$JOBDIR/smbh_starstar_collisions_PN_test.slurm"   "$RUNDIR/logs" "$SMBH_SLURM_PARTITION_TEST"
squeue -u "$USER"
