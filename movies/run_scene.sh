#!/bin/bash
# ============================================================================
# run_scene.sh  --  one representative scene, end to end, with smbh_binary_PN.
#
#   ./run_scene.sh ejection|mutde|flyby  [--dry-run]
#
# 1) SEARCH  : random orientations in --mode classify until the target outcome.
# 2) STAGE 1 : re-run that exact IC with --mode full; the trajectory is written to a
#              file (bss_*_full.txt) in the run dir, which we parse to canonical CSV.
# 3) POST    : ejection -> RESTART the BH+SMBH to resolve the EMRI (continuous seam);
#              mutde    -> spawn+follow tidal debris (ballistic post-step);
#              flyby    -> nothing extra.
# 4) MOVIE   : convert to canonical and render the looping GIF.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
source ./config.sh
SCENE="${1:?usage: run_scene.sh ejection|mutde|flyby [--dry-run]}"
DRY=0; [ "${2:-}" = "--dry-run" ] && DRY=1
[ -x "$BIN" ] || { echo "FATAL: executable not found/exec: $BIN"; exit 1; }
D="$OUTDIR/$SCENE"; mkdir -p "$D"

# ---- per-scene physics, outer-orbit flags, and target outcome --------------
case "$SCENE" in
  ejection) A_BIN=0.1; E_BIN=0.0;  OUTER="--outer hyperbolic --beta 1.7 --v_inf 1.0 --n_dyn 2e4"; TARGET="EXCHANGE_bh1" ;;
  mutde)    A_BIN=0.1; E_BIN=0.85; OUTER="--outer hyperbolic --beta 1.5 --v_inf 1.0 --n_dyn 2e4"; TARGET="TDE:star2_by_bh1" ;;
  flyby)    A_BIN=0.3; E_BIN=0.0;  OUTER="--outer hyperbolic --beta 0.7 --v_inf 1.0 --n_dyn 2e4"; TARGET="FLYBY" ;;
  *) echo "unknown scene '$SCENE'"; exit 1 ;;
esac
echo "=== scene=$SCENE  a_bin=$A_BIN e_bin=$E_BIN  target=$TARGET ==="
BINSCALE="$($PY -c "print(4*$A_BIN*(1+$E_BIN))")"   # inset half-width (binary scale)

common() {   # $1 incl $2 Omega $3 omega $4 nu
  echo "--M_smbh $M_SMBH --m1 $M_BH --m2 $M_STAR --type1 $TYPE1 --type2 $TYPE2 \
--R1 $R1 --R2 $R2 --a_bin $A_BIN --e_bin $E_BIN $OUTER \
--incl $1 --Omega $2 --omega $3 --nu $4 --rtol $RTOL"
}
rnd() { $PY -c "import random,sys;print(f'{random.uniform(float(sys.argv[1]),float(sys.argv[2])):.6f}')" "$1" "$2"; }

# ---- 1. search for the target outcome --------------------------------------
WIN=""
for i in $(seq 1 "$N_SEARCH"); do
  A="$(common "$(rnd 0 3.14159)" "$(rnd 0 6.2832)" "$(rnd 0 6.2832)" "$(rnd 0 6.2832)")"
  if [ "$DRY" = 1 ]; then echo "DRY: $BIN $A --mode classify --report 0"; WIN="$A"; break; fi
  OC="$("$BIN" $A --mode classify --report 0 2>/dev/null | grep '^CSV,' | awk -F, -v f=$OUTCOME_FIELD '{print $f}')"
  if [ "$OC" = "$TARGET" ]; then echo "  found $TARGET at try $i"; WIN="$A"; break; fi
done
[ -z "$WIN" ] && { echo "no $TARGET in $N_SEARCH tries; adjust OUTER/E_BIN in run_scene.sh"; exit 1; }

# ---- 2. stage-1 trajectory --------------------------------------------------
# --mode full writes the trajectory to a FILE ($TRAJ_GLOB) in the run's cwd, so we
# run it inside its own dir ($D/s1) and pick that file up (stdout is just the summary).
echo "+ stage-1 trajectory"
RAW1=""
if [ "$DRY" = 1 ]; then
  echo "DRY: (cd $D/s1 && $BIN $WIN $TRAJ_FLAGS)   # trajectory -> $TRAJ_GLOB"
else
  mkdir -p "$D/s1"; rm -f "$D/s1"/$TRAJ_GLOB
  ( cd "$D/s1" && "$BIN" $WIN $TRAJ_FLAGS >stdout.log 2>&1 )
  RAW1="$(ls -t "$D/s1"/$TRAJ_GLOB 2>/dev/null | head -1)"
  [ -z "$RAW1" ] && { echo "FATAL: no trajectory file ($TRAJ_GLOB) in $D/s1"; tail -3 "$D/s1/stdout.log" 2>/dev/null; exit 1; }
  $PY parse_traj.py "$RAW1" "$D/stage1.csv"
fi

# ---- 3+4. scene-specific post + render -------------------------------------
case "$SCENE" in
  flyby)
    [ "$DRY" = 0 ] && cp "$D/stage1.csv" "$OUTDIR/flyby.csv" && cp "$D/stage1.events.json" "$OUTDIR/flyby.events.json" 2>/dev/null
    [ "$DRY" = 0 ] && $PY -c "import json;m={'center':'smbh','zoom_on':['bh'],'inset':{'center':'com','com':{'bh':10,'star':1},'show':['bh','star'],'hw':$BINSCALE,'label':'binary (zoom)'}};json.dump(m,open('$OUTDIR/flyby.events.json','w'))"
    [ "$DRY" = 0 ] && $PY animate.py "$OUTDIR/flyby.csv" -o "$OUTDIR/flyby.gif" --title "flyby — binary survives" --center smbh --zoom bh --fps $FPS ;;
  mutde)
    # debris seeded from the executable's EXACT disruption state ($RAW1): n=3 polytrope,
    # spin-up to breakup, ballistic in the BH+SMBH field. Sidecar drives center/zoom/inset.
    [ "$DRY" = 0 ] && $PY debris_follow.py "$RAW1" "$OUTDIR/mutde.csv" \
        --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR --ndeb $NDEB --norbits $MUTDE_N_ORBITS
    [ "$DRY" = 0 ] && $PY animate.py "$OUTDIR/mutde.csv" -o "$OUTDIR/mutde.gif" \
        --title "micro-TDE + debris fallback" --fps $FPS
    # extra: physical fallback dM/dt curve (same disruption state; long integration for the tail)
    [ "$DRY" = 0 ] && [ "${FALLBACK_ENABLE:-1}" = 1 ] && $PY fallback_stage.py "$RAW1" "$OUTDIR/mutde" \
        --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR --ndeb $NDEB --disruptor bh \
        --norbits $FALLBACK_N_ORBITS --norbits-bh $FALLBACK_N_ORBITS_BH --ntime $FALLBACK_NTIME \
        || echo "!! mutde fallback stage skipped/failed" ;;
  ejection)
    echo "+ EMRI restart (resolve the bound BH about the SMBH)"
    if [ "$DRY" = 0 ]; then
      $PY make_restart.py "$RAW1" "$D/restart.txt" --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR | tee "$D/restart.log"
      P=$(grep PERIOD_YR "$D/restart.log" | cut -d= -f2)
      TEND=$($PY -c "print(max(1e-3, $EMRI_N_ORBITS*float('$P')))")
      echo "  captured period ~ $P yr  ->  integrating $EMRI_N_ORBITS orbits (t_end=$TEND yr)"
      DABS="$(cd "$D" && pwd)"
      mkdir -p "$D/s2"; rm -f "$D/s2"/$TRAJ_GLOB
      ( cd "$D/s2" && "$BIN" --restart "$DABS/restart.txt" --t_end_yr "$TEND" --r_esc_AU "$EMRI_R_ESC_AU" $TRAJ_FLAGS >stdout.log 2>&1 )
      RAW2="$(ls -t "$D/s2"/$TRAJ_GLOB 2>/dev/null | head -1)"
      [ -z "$RAW2" ] && { echo "FATAL: no EMRI trajectory file in $D/s2"; tail -3 "$D/s2/stdout.log" 2>/dev/null; exit 1; }
      $PY parse_traj.py "$RAW2" "$D/stage2.csv"
      $PY concat_stages.py "$D/stage1.csv" "$D/stage2.csv" "$OUTDIR/ejection.csv"
      $PY -c "import json,os;f='$OUTDIR/ejection.events.json';m=json.load(open(f));m['inset']={'center':'com','com':{'bh':10,'star':1},'show':['bh','star'],'hw':$BINSCALE,'label':'binary (zoom)'};json.dump(m,open(f,'w'))"
      $PY animate.py "$OUTDIR/ejection.csv" -o "$OUTDIR/ejection.gif" --title "ejection + capture → EMRI" --center smbh --zoom bh --fps $FPS
    else
      echo "DRY: make_restart -> $BIN --restart restart.txt --t_end_yr <N*P> --r_esc_AU $EMRI_R_ESC_AU $TRAJ_FLAGS"
    fi ;;
esac
echo "=== done: $SCENE -> $OUTDIR/$SCENE.gif ==="
