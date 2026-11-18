#!/bin/bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/env.sh"   # sets SMBH_ROOT, SMBH_BIN, SMBH_SCRATCH, slurm_submit (reads env.local.sh)
JOBDIR="$SMBH_ROOT/campaigns/bhstar"
GEN="$JOBDIR/smbh_bhstar_sweep.py"
RUNDIR="$SMBH_SCRATCH/smbh_bhstar_bmaxtest"

mkdir -p "$RUNDIR/logs"
cd "$RUNDIR"

echo "generating b_max-test joblist (TEST grid, RP_TIDAL_FAC=10)..."
BHSTAR_TEST=1 BHSTAR_RP_TIDAL_FAC=10 python3 "$GEN" > joblist.txt
echo "  $(wc -l < joblist.txt) lines; grid_meta.csv written to $RUNDIR"

shuf joblist.txt > joblist_shuffled.txt

echo "submitting debug test job..."
slurm_submit "$JOBDIR/smbh_bhstar_test.slurm" "$RUNDIR/logs" "$SMBH_SLURM_PARTITION_TEST"
squeue -u "$USER"
