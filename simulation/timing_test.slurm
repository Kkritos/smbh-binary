#!/bin/bash
#SBATCH --job-name=bss_timing
#SBATCH --partition=short-40core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=/path/to/scratch/smbh_binary_interactions/timing_%j.out

: "${SMBH_ROOT:=/path/to/smbh-binary}"   # exported by submit_*.sh
source "$SMBH_ROOT/env.sh"
cd "$SMBH_SCRATCH/smbh_binary_interactions"
BIN="$SMBH_BIN/smbh_binary_PN"

echo "=== gentle: beta=0.5, alpha=3 ==="
time "$BIN" --M_smbh 1e6 --m1 1 --type1 star --R1 1 --m2 1 --type2 star --R2 1 \
     --a_bin 0.0279028 --e_bin 0 --v_inf 25.2 --beta 0.5 \
     --mode classify --rtol 1e-9 --report 0

echo "=== deep: beta=10, alpha=3 ==="
time "$BIN" --M_smbh 1e6 --m1 1 --type1 star --R1 1 --m2 1 --type2 star --R2 1 \
     --a_bin 0.0279028 --e_bin 0 --v_inf 25.2 --beta 10 \
     --mode classify --rtol 1e-9 --report 0
