#!/bin/bash
# ============================================================================
# make_movies.sh  --  ONE BUTTON.  Builds all three seminar GIFs with the real
# smbh_binary_PN executable (search ICs -> integrate -> restart EMRI -> render).
#
#   ./make_movies.sh                # all three scenes, via SpaceHub
#   ./make_movies.sh ejection       # just one scene
#   ./make_movies.sh --probe        # run a tiny full-mode sim, print the dump format
#   ./make_movies.sh demo           # built-in integrator fallback (no executable)
#
# Everything runs inside this folder; the executable is found one folder up
# (or via config.sh).  Output GIFs land in scenes/.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
source ./config.sh
ARG="${1:-all}"

if [ "$ARG" = "--probe" ]; then
  echo "probe: a short full-mode run; the trajectory is written to a FILE ($TRAJ_GLOB),"
  echo "not stdout.  Running in a temp dir and showing that file's first lines."
  PD="$(mktemp -d)"
  ( cd "$PD" && "$BIN" --M_smbh $M_SMBH --m1 $M_BH --m2 $M_STAR --type1 $TYPE1 --type2 $TYPE2 \
         --R1 $R1 --R2 $R2 --a_bin 0.3 --e_bin 0.2 \
         --outer hyperbolic --beta 0.5 --v_inf 1.0 --n_dyn 2e2 \
         --incl 0.6 --Omega 1.0 --omega 2.0 --nu 3.0 --rtol $RTOL $TRAJ_FLAGS >stdout.log 2>&1 )
  RAW="$(ls -t "$PD"/$TRAJ_GLOB 2>/dev/null | head -1)"
  if [ -z "$RAW" ]; then
    echo "!! no trajectory file matched $TRAJ_GLOB in $PD"
    echo "   files written:"; ls -l "$PD"
    echo "   stdout/stderr:"; cat "$PD/stdout.log"
  else
    echo "--- trajectory file: $RAW"
    head -8 "$RAW"
    echo "..."
    echo "--- parse_traj.py check (-> canonical first rows):"
    $PY parse_traj.py "$RAW" "$PD/canon.csv" && head -4 "$PD/canon.csv"
  fi
  rm -rf "$PD"
  exit 0
fi

if [ "$ARG" = "survey" ]; then
  exec ./survey.sh "${2:-}"        # outcome gallery: many hyperbolic encounters
fi

if [ "$ARG" = "montage" ]; then    # tile the survey gifs into one grid video (gif + mp4)
  $PY montage.py "$OUTDIR/survey"/*.gif -o "$OUTDIR/survey/montage.gif" --fps 20
  echo "== montage: $OUTDIR/survey/montage.gif (+ .mp4) =="
  exit 0
fi

if [ "$ARG" = "demo" ]; then
  echo "== DEMO (built-in integrator, no executable) =="
  $PY demo_gen_eject.py && $PY animate.py scenes/ejection.csv --title "ejection + capture → GW inspiral" --center smbh --zoom bh --trail 60
  $PY mutde_gen.py       && $PY animate.py scenes/mutde.csv    --title "micro-TDE + debris fallback"      --center smbh --zoom bh,deb --trail 14
  $PY flyby_gen.py       && $PY animate.py scenes/flyby.csv    --title "flyby — binary survives"          --center smbh --zoom bh
  echo "== done =="; exit 0
fi

SCENES="ejection mutde flyby"
[ "$ARG" != "all" ] && SCENES="$ARG"
for s in $SCENES; do ./run_scene.sh "$s" || echo "!! $s failed"; done
echo "== GIFs: $(ls scenes/*.gif 2>/dev/null | tr '\n' ' ') =="
