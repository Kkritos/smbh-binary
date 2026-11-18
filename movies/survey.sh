#!/bin/bash
# ============================================================================
# survey.sh  --  a gallery of BH-star + SMBH encounters, always on a HYPERBOLIC
# approach, spanning outcomes as (a_bin, e_bin, v_inf, penetration) are varied.
# Each row is run with a random orientation; those with a TARGET search until that
# outcome occurs (guaranteeing the gallery contains it), the rest keep whatever the
# random orientation produces.  Dispatch + rendering is handled by run_sim.sh.
#
#   ./survey.sh                # run the whole gallery
#   ./survey.sh mutde          # run just one row by name
#
# Masses are fixed (m_BH=10, m_star=1, m_SMBH=1e6, in config.sh); edit the table
# freely — add rows, change parameters, set/clear targets.
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"

#     name        | a_bin | e_bin | outer flags (hyperbolic)                         | target (blank = natural outcome)
SIMS=(
 # ---- EXCHANGES (6): BH captured -> EMRI (HVS speed + Peters t_merge), or star captured ----
 "ejection       |0.1  |0.0  |--outer hyperbolic --beta 1.7 --v_inf 1.0 --n_dyn 2e4 |EXCHANGE_bh1"
 "ejection_ecc   |0.08 |0.3  |--outer hyperbolic --beta 2.0 --v_inf 1.0 --n_dyn 2e4 |EXCHANGE_bh1"
 "ejection_wide  |0.2  |0.1  |--outer hyperbolic --beta 1.5 --v_inf 1.5 --n_dyn 2e4 |EXCHANGE_bh1"
 "star_capture   |0.05 |0.6  |--outer hyperbolic --beta 2.0 --v_inf 2.0 --n_dyn 2e4 |EXCHANGE_star2"
 "star_capt_tight|0.04 |0.5  |--outer hyperbolic --beta 2.5 --v_inf 1.5 --n_dyn 2e4 |EXCHANGE_star2"
 "star_capt_ecc  |0.06 |0.7  |--outer hyperbolic --beta 1.8 --v_inf 2.0 --n_dyn 2e4 |EXCHANGE_star2"
 # ---- MICRO-TDEs (12): star disrupted by its BH companion, across a/e/v/beta ----
 "mutde_mid      |0.1  |0.85 |--outer hyperbolic --beta 1.5 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_tight    |0.05 |0.9  |--outer hyperbolic --beta 2.0 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_wide     |0.3  |0.8  |--outer hyperbolic --beta 1.2 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_slow     |0.1  |0.9  |--outer hyperbolic --beta 1.5 --v_inf 0.5 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_fast     |0.2  |0.85 |--outer hyperbolic --beta 1.5 --v_inf 5.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_deep     |0.1  |0.8  |--outer hyperbolic --beta 3.0 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_hecc     |0.08 |0.95 |--outer hyperbolic --beta 1.6 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_vwide    |0.5  |0.6  |--outer hyperbolic --beta 1.2 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_graze    |0.12 |0.8  |--outer hyperbolic --beta 1.0 --v_inf 2.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_vslow    |0.15 |0.9  |--outer hyperbolic --beta 2.2 --v_inf 0.3 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_close    |0.03 |0.7  |--outer hyperbolic --beta 2.5 --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_bh1"
 "mutde_mod      |0.1  |0.6  |--outer hyperbolic --beta 1.8 --v_inf 1.5 --n_dyn 2e4 |TDE:star2_by_bh1"
 # ---- a few other outcomes for variety (7): SMBH-TDE, ionization, flyby, triple ----
 "smbh_tde       |0.1  |0.0  |--outer hyperbolic --beta 12  --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_SMBH"
 "smbh_tde_ecc   |0.05 |0.3  |--outer hyperbolic --beta 15  --v_inf 1.0 --n_dyn 2e4 |TDE:star2_by_SMBH"
 "ioniz_a        |2.0  |0.0  |--outer hyperbolic --beta 1.5 --v_inf 100 --n_dyn 2e4 |IONIZATION"
 "ioniz_b        |4.0  |0.0  |--outer hyperbolic --beta 1.2 --v_inf 200 --n_dyn 2e4 |IONIZATION"
 "flyby_ecc      |0.3  |0.0  |--outer hyperbolic --beta 0.8 --v_inf 1.0 --n_dyn 2e4 |FLYBY_ECC"
 "flyby_dist     |0.3  |0.2  |--outer hyperbolic --beta 0.5 --v_inf 1.0 --n_dyn 2e4 |FLYBY"
 "triple         |0.2  |0.3  |--outer hyperbolic --beta 2.5 --v_inf 0.5 --n_dyn 2e4 |BOUND_TRIPLE"
)
ONLY="${1:-}"
MAN="scenes/survey/manifest.csv"
# Optional REPLAY table: name,incl,Omega,omega,nu .  If present (and REPLAY!=0), each row's
# orientation is LOCKED to these exact angles instead of being searched at random, so the
# survey reproduces a previous run bit-for-bit (angles to the precision stored in the file).
# Delete/rename orientations.csv, or run `REPLAY=0 ./survey.sh`, to go back to random search.
ORIENT_FILE="${ORIENT_FILE:-orientations.csv}"
REPLAY="${REPLAY:-1}"
[ "$REPLAY" = 1 ] && [ -f "$ORIENT_FILE" ] && echo "replay: locking orientations from $ORIENT_FILE"
[ -z "$ONLY" ] && { mkdir -p scenes/survey; rm -f "$MAN"; }   # fresh manifest for a full run
for row in "${SIMS[@]}"; do
  IFS='|' read -r name a e outer target <<< "$row"
  name="$(echo "$name" | xargs)"; a="$(echo "$a" | xargs)"; e="$(echo "$e" | xargs)"
  outer="$(echo "$outer" | sed 's/^ *//;s/ *$//')"; target="$(echo "$target" | xargs)"
  [ -n "$ONLY" ] && [ "$ONLY" != "$name" ] && continue
  ORI=""
  if [ "$REPLAY" = 1 ] && [ -f "$ORIENT_FILE" ]; then
    ORI="$(awk -F, -v n="$name" 'NR>1 && $1==n{printf "%s %s %s %s",$2,$3,$4,$5; exit}' "$ORIENT_FILE")"
  fi
  echo "=== $name  a_bin=$a e_bin=$e  [$outer]  target=${target:-<random>}${ORI:+  [locked: $ORI]} ==="
  ./run_sim.sh "$name" "$a" "$e" "$outer" "$target" "$ORI" || echo "  (skipped $name)"
done
echo "== survey GIFs: $(ls "scenes/survey"/*.gif 2>/dev/null | tr '\n' ' ') =="
if [ -f "$MAN" ]; then
  echo "== initial conditions -> outcome  (saved in $MAN) =="
  column -t -s, "$MAN" 2>/dev/null || cat "$MAN"
fi

# collect every TDE scene's direct-vs-frozen delay comparison into one table
CMP="scenes/survey/delay_comparison.csv"; first=1
for f in scenes/survey/*/fallback.compare.csv; do
  [ -f "$f" ] || continue
  nm="$(basename "$(dirname "$f")")"
  [ "$first" = 1 ] && { echo "name,$(head -1 "$f")" > "$CMP"; first=0; }
  echo "$nm,$(tail -1 "$f")" >> "$CMP"
done
if [ -s "$CMP" ]; then
  echo "== direct vs frozen-in Dt_delay, all TDE scenes  (saved in $CMP) =="
  column -t -s, "$CMP" 2>/dev/null || cat "$CMP"
fi
