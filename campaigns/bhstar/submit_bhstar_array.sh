#!/bin/bash
set -euo pipefail

# Master (ARRAY version): generate the production joblist, shuffle, and submit the
# BH-star sweep as a 8-way job array (smbh_bhstar_array.slurm).  Writes to scratch.
#   - BHSTAR_PILOT=0 forces the full grid.
#   - shuffle makes the array's contiguous shards load-balanced.
#   - grid_meta.csv (b_max per grid point) lands in RUNDIR.

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/env.sh"   # sets SMBH_ROOT, SMBH_BIN, SMBH_SCRATCH, slurm_submit (reads env.local.sh)
JOBDIR="$SMBH_ROOT/campaigns/bhstar"
GEN="$JOBDIR/smbh_bhstar_sweep.py"
RUNDIR="$SMBH_SCRATCH/smbh_bhstar_encounters"

mkdir -p "$RUNDIR/logs"
cd "$RUNDIR"

echo "generating production joblist (PILOT off)..."
BHSTAR_PILOT=0 BHSTAR_RP_TIDAL_FAC=10 python3 "$GEN" > joblist.txt
echo "  $(wc -l < joblist.txt) lines -> joblist.txt"
echo "  grid_meta.csv written to $RUNDIR"

echo "shuffling (load balance across array shards + workers)..."
shuf joblist.txt > joblist_shuffled.txt
echo "  $(wc -l < joblist_shuffled.txt) lines -> joblist_shuffled.txt"

echo "submitting 8-way array..."
slurm_submit "$JOBDIR/smbh_bhstar_array.slurm" "$RUNDIR/logs" "$SMBH_SLURM_PARTITION"

squeue -u "$USER"
