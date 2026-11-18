#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
SURV=scenes/survey; [ -d "$SURV" ] || SURV=scenes
sed -i 's/(wrt BH2)/(wrt SMBH)/' fallback_stage.py 2>/dev/null
mkdir -p refig_logs
n=0; ok=0
for gif in "$SURV"/*TDE-star2_by_bh1*.gif; do
  [ -e "$gif" ] || { echo "no micro-TDE scenes found under $SURV/"; break; }
  n=$((n+1))
  nm=$(basename "$gif"); nm=${nm%%_TDE-star2_by_bh1*}
  d="$SURV/$nm"
  raw=$(ls -t "$d"/s1/bss_*_full.txt 2>/dev/null | head -1)
  if [ -z "$raw" ]; then echo "[$nm] NO RAW traj -> needs a re-run"; continue; fi
  echo "[$nm] generating figures ..."
  if timeout 1800 python3 fallback_stage.py "$raw" "$d/fallback" \
        --Msmbh 1e6 --mbh 10 --mstar 1 --disruptor bh \
        --ndeb 4000 --ntime 1500 --nsnap 6 > "refig_logs/$nm.log" 2>&1; then
     echo "   OK -> $d/fallback.fallback.png , $d/fallback.energy.png"; ok=$((ok+1))
  else
     echo "   FAILED/timeout (see refig_logs/$nm.log)"
  fi
done
echo "done: $ok/$n micro-TDE scenes figured."
