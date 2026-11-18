#!/bin/bash
# ============================================================================
# run_sim.sh  --  ONE hyperbolic BH-star + SMBH encounter, dispatched by its
#                 ACTUAL outcome (not a pre-chosen scene).
#
#   ./run_sim.sh NAME A_BIN E_BIN "OUTER_FLAGS" [TARGET_OUTCOME]
#
# * TARGET given  -> random orientations are searched until the code reports it.
#   Special target FLYBY_ECC: search a (circular) binary for a FLYBY whose SURVIVING
#   binary eccentricity e_new >= FLYBY_ECC_MIN (keeps the highest-e flyby found).
# * TARGET empty  -> ONE random orientation is run and whatever happens is kept.
# Then the exact IC is re-run with the trajectory dump (and --report 1 so the
# human-readable SUMMARY is saved) and dispatched:
#   EXCHANGE_bh1       -> EMRI restart; HVS (star) v_esc + Peters t_merge recorded
#   EXCHANGE_star2     -> captured-star restart; runaway-BH v_esc recorded
#   TDE:star2_by_bh1   -> micro-TDE debris fallback (--disruptor bh)
#   TDE:star2_by_SMBH  -> SMBH-TDE debris stream    (--disruptor smbh)
#   IONIZATION         -> keep integrating (bodies fly apart) + time-lapse playback
#   FLYBY / BOUND_TRIPLE / ...  -> plain trajectory, labelled by outcome
# Output: scenes/survey/NAME_<OUTCOME>.gif (+ .mp4), per-sim summary.txt, and a row
# in scenes/survey/manifest.csv with ICs AND final conditions.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
source ./config.sh
NAME="${1:?usage: run_sim.sh NAME A_BIN E_BIN \"OUTER\" [TARGET] [\"incl Omega omega nu\"]}"
A_BIN="${2:?}"; E_BIN="${3:?}"; OUTER="${4:?}"; TARGET="${5:-}"; ORIENT="${6:-}"
[ -x "$BIN" ] || { echo "FATAL: executable not found: $BIN"; exit 1; }
OUT="$OUTDIR/survey"; D="$OUT/$NAME"; mkdir -p "$D"

rnd(){ $PY -c "import random,sys;print(f'{random.uniform(float(sys.argv[1]),float(sys.argv[2])):.6f}')" "$1" "$2"; }
pick_orient(){ INCL="$(rnd 0 3.14159)"; OMEGA="$(rnd 0 6.2832)"; SOMEGA="$(rnd 0 6.2832)"; NU="$(rnd 0 6.2832)"; }
common(){ echo "--M_smbh $M_SMBH --m1 $M_BH --m2 $M_STAR --type1 $TYPE1 --type2 $TYPE2 \
--R1 $R1 --R2 $R2 --a_bin $A_BIN --e_bin $E_BIN $OUTER \
--incl $INCL --Omega $OMEGA --omega $SOMEGA --nu $NU --rtol $RTOL"; }
classify_line(){ "$BIN" $1 --mode classify --report 0 2>/dev/null | grep '^CSV,' | tail -1; }
field(){ echo "$1" | awk -F, -v f="$2" '{print $f}'; }
outcome_of(){ field "$(classify_line "$1")" "$OUTCOME_FIELD"; }
ge(){ $PY -c "import sys;sys.exit(0 if float(sys.argv[1])>=float(sys.argv[2]) else 1)" "$1" "$2" 2>/dev/null; }
gt(){ $PY -c "import sys;sys.exit(0 if float(sys.argv[1])> float(sys.argv[2]) else 1)" "$1" "$2" 2>/dev/null; }

# ---- 1. pick an orientation / outcome (INCL/OMEGA/SOMEGA/NU are recorded) ----
# ORIENT given (replay): use that exact orientation, no search -> reproduces the sim.
WIN=""; OUTCOME=""; WANT_ECC=0; BASE_TARGET="$TARGET"
[ "$TARGET" = "FLYBY_ECC" ] && { WANT_ECC=1; BASE_TARGET="FLYBY"; }
if [ -n "$ORIENT" ]; then
  read -r INCL OMEGA SOMEGA NU <<< "$ORIENT"
  WIN="$(common)"; OUTCOME="$(outcome_of "$WIN")"
  echo "  [$NAME] locked orientation ($ORIENT) -> ${OUTCOME:-<none>}"
elif [ -n "$TARGET" ]; then
  BEST_E="-1"; BEST_A=""
  for i in $(seq 1 "$N_SEARCH"); do
    pick_orient; A="$(common)"; LINE="$(classify_line "$A")"
    [ -z "$LINE" ] && continue
    OC="$(field "$LINE" "$OUTCOME_FIELD")"
    if [ "$WANT_ECC" = "1" ]; then
      [ "$OC" != "FLYBY" ] && continue
      EE="$(field "$LINE" $((OUTCOME_FIELD+ENEW_OFF)))"
      if ge "$EE" "$FLYBY_ECC_MIN"; then
        WIN="$A"; OUTCOME="FLYBY"; echo "  [$NAME] FLYBY e_new=$EE >= $FLYBY_ECC_MIN at try $i"; break; fi
      if gt "$EE" "$BEST_E"; then BEST_E="$EE"; BEST_A="$A"; fi
    else
      [ "$OC" = "$BASE_TARGET" ] && { WIN="$A"; OUTCOME="$OC"; echo "  [$NAME] $BASE_TARGET at try $i"; break; }
    fi
  done
  if [ -z "$WIN" ] && [ "$WANT_ECC" = "1" ] && [ -n "$BEST_A" ]; then
    WIN="$BEST_A"; OUTCOME="FLYBY"; echo "  [$NAME] no FLYBY e>=$FLYBY_ECC_MIN in $N_SEARCH tries; using best e_new=$BEST_E"; fi
  [ -z "$WIN" ] && { echo "  [$NAME] no $TARGET in $N_SEARCH tries — skipped"; exit 2; }
else
  pick_orient; WIN="$(common)"; OUTCOME="$(outcome_of "$WIN")"; echo "  [$NAME] random -> ${OUTCOME:-<none>}"
fi
[ -z "$OUTCOME" ] && { echo "  [$NAME] classify failed — skipped"; exit 1; }

# ---- 2. stage-1 full trajectory (+ report 1 so the summary is saved) ---------
mkdir -p "$D/s1"; rm -f "$D/s1"/$TRAJ_GLOB
( cd "$D/s1" && "$BIN" $WIN $TRAJ_FLAGS --report 1 >stdout.log 2>&1 )
RAW1="$(ls -t "$D/s1"/$TRAJ_GLOB 2>/dev/null | head -1)"
[ -z "$RAW1" ] && { echo "  [$NAME] no trajectory file — skipped"; exit 1; }
$PY parse_traj.py "$RAW1" "$D/stage1.csv" >/dev/null
cp -f "$D/s1/stdout.log" "$D/summary.txt" 2>/dev/null || true

# ---- final conditions from the stage-1 CSV line (offsets from OUTCOME_FIELD) --
CSVLINE="$(grep '^CSV,' "$D/s1/stdout.log" 2>/dev/null | tail -1)"
A_NEW="nan"; E_NEW="nan"; V_ESC="nan"; TMERGE="nan"
if [ -n "$CSVLINE" ]; then
  A_NEW="$(field "$CSVLINE" $((OUTCOME_FIELD+ANEW_OFF)))"
  E_NEW="$(field "$CSVLINE" $((OUTCOME_FIELD+ENEW_OFF)))"
  V_ESC="$(field "$CSVLINE" $((OUTCOME_FIELD+VESC_OFF)))"
fi
if [ "$OUTCOME" = "EXCHANGE_bh1" ]; then
  TMERGE="$($PY peters.py "$A_NEW" "$E_NEW" --m1 $M_SMBH --m2 $M_BH 2>/dev/null || echo nan)"
  echo "  [$NAME] HVS (ejected star) v_esc=$V_ESC km/s | EMRI a=$A_NEW AU e=$E_NEW t_merge(Peters)=$TMERGE yr"
elif [ "$OUTCOME" = "EXCHANGE_star2" ]; then
  echo "  [$NAME] runaway BH v_esc=$V_ESC km/s | captured-star orbit a=$A_NEW AU e=$E_NEW"
fi

TITLE="$($PY - "$OUTCOME" <<'PY'
import sys; o=sys.argv[1]
m={"EXCHANGE_bh1":"ejection + capture → EMRI","TDE:star2_by_bh1":"micro-TDE + debris fallback",
   "FLYBY":"flyby — binary survives","IONIZATION":"ionization — binary unbound",
   "BOUND_TRIPLE":"bound triple (unresolved)","TDE:star2_by_SMBH":"SMBH TDE + debris stream",
   "EXCHANGE_star2":"star captured → orbits SMBH"}
print(m.get(o, ("plunge into SMBH" if o.startswith("PLUNGE") else
               "collision" if o.startswith("COLLISION") else
               "exchange / capture" if o.startswith("EXCHANGE") else o)))
PY
)"
GIF="$OUT/${NAME}_${OUTCOME//:/-}.gif"

# inset half-width: contains the binary at its initial OR final (post-encounter) size
BINSCALE="$($PY -c "
ab,eb=$A_BIN,$E_BIN
try: a=float('$A_NEW'); e=float('$E_NEW'); fin=a*(1+e) if a>0 and 0<=e<1 else 0.0
except Exception: fin=0.0
print(4.0*max(ab*(1+eb), fin))")"

render_traj(){   # $1 = inset half-width (AU); 0 = no inset
  cp "$D/stage1.csv" "$D/final.csv"
  $PY - "$D/final.events.json" "$1" <<'PY'
import json,sys
f,ihw=sys.argv[1],float(sys.argv[2])
m={"center":"smbh","zoom_on":["bh","star"]}
if ihw>0: m["inset"]={"center":"com","com":{"bh":10,"star":1},"show":["bh","star"],
                      "hw":ihw,"label":"binary (zoom)"}
json.dump(m,open(f,"w"))
PY
  $PY animate.py "$D/final.csv" -o "$GIF" --title "$TITLE" --center smbh --zoom bh,star --fps $FPS
}

# ---- record initial conditions AND final conditions -> manifest -------------
$PY - "$OUT/manifest.csv" "$NAME" "$A_BIN" "$E_BIN" "$OUTER" "$INCL" "$OMEGA" "$SOMEGA" "$NU" \
     "$OUTCOME" "$A_NEW" "$E_NEW" "$V_ESC" "$TMERGE" "$(basename "$GIF")" <<'PY'
import sys, os, csv
(man, name, ab, eb, outer, incl, Om, om, nu,
 outcome, a_new, e_new, v_esc, tmerge, gif) = sys.argv[1:16]
toks = outer.split()
gv = lambda fl: toks[toks.index(fl)+1] if fl in toks else "nan"
cols = ["name","a_bin","e_bin","v_inf_kms","beta","incl","Omega","omega","nu",
        "outcome","a_new_AU","e_new","v_esc_kms","t_merge_yr","gif"]
row  = [name, ab, eb, gv("--v_inf"), gv("--beta"), incl, Om, om, nu,
        outcome, a_new, e_new, v_esc, tmerge, gif]
rows = []
if os.path.exists(man):
    with open(man) as f: r = list(csv.reader(f))
    if r and r[0] == cols: rows = [x for x in r[1:] if x and x[0] != name]
rows.append(row)
with open(man, "w", newline="") as f:
    w = csv.writer(f); w.writerow(cols); w.writerows(rows)
PY
echo "  [$NAME] $OUTCOME  (ICs + final conditions logged to $OUT/manifest.csv; summary in $D/summary.txt)"

# ---- 3. dispatch by outcome -------------------------------------------------
case "$OUTCOME" in
  EXCHANGE_bh1|EXCHANGE_star2)
    # the captured body orbits the SMBH; restart and integrate it for several orbits
    if [ "$OUTCOME" = "EXCHANGE_star2" ]; then CAP=star; PKEY=PERIOD_STAR_YR; DROP=bh
    else CAP=bh; PKEY=PERIOD_BH_YR; DROP=star; fi
    $PY make_restart.py "$RAW1" "$D/restart.txt" --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR \
        | tee "$D/restart.log" >/dev/null
    P=$(grep "^$PKEY=" "$D/restart.log" | cut -d= -f2)
    TEND=$($PY -c "print(max(1e-3, $EMRI_N_ORBITS*float('$P')))")
    DABS="$(cd "$D" && pwd)"; mkdir -p "$D/s2"; rm -f "$D/s2"/$TRAJ_GLOB
    ( cd "$D/s2" && "$BIN" --restart "$DABS/restart.txt" --t_end_yr "$TEND" \
         --r_esc_AU "$EMRI_R_ESC_AU" $TRAJ_FLAGS >stdout.log 2>&1 )
    RAW2="$(ls -t "$D/s2"/$TRAJ_GLOB 2>/dev/null | head -1)"
    if [ -n "$RAW2" ]; then
      $PY parse_traj.py "$RAW2" "$D/stage2.csv" >/dev/null
      $PY concat_stages.py "$D/stage1.csv" "$D/stage2.csv" "$D/final.csv" --drop $DROP >/dev/null
      $PY -c "import json;f='$D/final.events.json';m=json.load(open(f));m['inset']={'center':'com','com':{'bh':10,'star':1},'show':['bh','star'],'hw':$BINSCALE,'label':'binary (zoom)'};json.dump(m,open(f,'w'))"
      $PY animate.py "$D/final.csv" -o "$GIF" --title "$TITLE" --center smbh --zoom $CAP --fps $FPS
    else
      echo "  [$NAME] restart wrote no file; rendering stage-1"; render_traj "$BINSCALE"
    fi ;;
  TDE:star2_by_bh1)
    $PY debris_follow.py "$RAW1" "$D/final.csv" --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR \
        --ndeb $NDEB --norbits $MUTDE_N_ORBITS --nframes $MUTDE_NFRAMES --disruptor bh
    $PY animate.py "$D/final.csv" -o "$GIF" --title "$TITLE" --fps $FPS
    # extra: physical fallback dM/dt curve (same ICs; long integration for the tail).
    # Wrapped in `timeout` so one slow encounter can never stall the whole survey.
    [ "${FALLBACK_ENABLE:-1}" = 1 ] && timeout "${FALLBACK_TIMEOUT:-900}" $PY fallback_stage.py "$RAW1" "$D/fallback" \
        --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR --ndeb $NDEB --disruptor bh \
        --norbits $FALLBACK_N_ORBITS --norbits-bh $FALLBACK_N_ORBITS_BH --ntime $FALLBACK_NTIME \
        || echo "  [$NAME] fallback stage skipped/failed/timed-out" ;;
  TDE:star2_by_SMBH)
    $PY debris_follow.py "$RAW1" "$D/final.csv" --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR \
        --ndeb $NDEB --norbits $MUTDE_N_ORBITS --nframes $MUTDE_NFRAMES --disruptor smbh
    $PY animate.py "$D/final.csv" -o "$GIF" --title "$TITLE" --fps $FPS
    # extra: physical fallback dM/dt curve (classic SMBH TDE -> expect t^-5/3)
    [ "${FALLBACK_ENABLE:-1}" = 1 ] && timeout "${FALLBACK_TIMEOUT:-900}" $PY fallback_stage.py "$RAW1" "$D/fallback" \
        --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR --ndeb $NDEB --disruptor smbh \
        --norbits $FALLBACK_N_ORBITS --ntime $FALLBACK_NTIME \
        || echo "  [$NAME] fallback stage skipped/failed/timed-out" ;;
  IONIZATION)
    # keep integrating the 3 now-unbound bodies so they visibly separate, then play
    # back with accelerated (time-lapse) time: slow through the disruption, fast apart.
    $PY make_restart.py "$RAW1" "$D/restart.txt" --Msmbh $M_SMBH --mbh $M_BH --mstar $M_STAR \
        | tee "$D/restart.log" >/dev/null
    T1="$(tail -1 "$D/stage1.csv" | cut -d, -f1)"
    TEND=$($PY -c "t='$T1'
try: t=float(t)
except Exception: t=0.0
print(max(1e-3, $IONIZ_EXTEND*t if t>0 else 1.0))")
    DABS="$(cd "$D" && pwd)"; mkdir -p "$D/s2"; rm -f "$D/s2"/$TRAJ_GLOB
    ( cd "$D/s2" && "$BIN" --restart "$DABS/restart.txt" --t_end_yr "$TEND" \
         --r_esc_AU "$IONIZ_ESC_AU" $TRAJ_FLAGS >stdout.log 2>&1 )
    RAW2="$(ls -t "$D/s2"/$TRAJ_GLOB 2>/dev/null | head -1)"
    if [ -n "$RAW2" ]; then
      $PY parse_traj.py "$RAW2" "$D/stage2.csv" >/dev/null
      $PY timewarp.py "$D/final.csv" "$D/stage1.csv" "$D/stage2.csv" \
          --nframes $IONIZ_NFRAMES --accel $IONIZ_ACCEL >/dev/null
    else
      echo "  [$NAME] ionization restart wrote no file; time-lapsing stage-1 only"
      $PY timewarp.py "$D/final.csv" "$D/stage1.csv" --nframes $IONIZ_NFRAMES --accel $IONIZ_ACCEL >/dev/null
    fi
    $PY -c "import json;json.dump({'center':'smbh','zoom_on':['bh','star'],'zoom_pct':100},open('$D/final.events.json','w'))"
    $PY animate.py "$D/final.csv" -o "$GIF" --title "$TITLE" --center smbh --zoom bh,star --trail 22 --fps $FPS ;;
  FLYBY|BOUND_TRIPLE)
    render_traj "$BINSCALE" ;;
  *)
    render_traj 0 ;;
esac
echo "  [$NAME] $OUTCOME -> $GIF"
