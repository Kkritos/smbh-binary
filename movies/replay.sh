#!/bin/bash
# ============================================================================
# replay.sh  --  reproduce an EXACT gallery from a saved manifest.csv.
#
#   ./replay.sh survey_locked.csv            # replay every row
#   ./replay.sh survey_locked.csv mutde_mid  # replay one row by name
#
# survey.sh searches RANDOM, unseeded orientations and RESETS manifest.csv on every
# full run, so a re-run never reproduces a previous gallery.  This driver instead reads
# the LOCKED orientation (incl, Omega, omega, nu) + (a_bin, e_bin, v_inf, beta) recorded
# for each sim and re-runs that EXACT initial condition (masses are fixed in config.sh),
# so you get the same encounters back — now rendered with the current pipeline.
#
# To preserve a gallery you like before changing anything:
#     cp scenes/survey/manifest.csv survey_locked.csv
# then later:  ./replay.sh survey_locked.csv
#
# Expected columns (the manifest header):
#   name,a_bin,e_bin,v_inf_kms,beta,incl,Omega,omega,nu,outcome[,...]
# (extra trailing columns such as a_new_AU/e_new/v_esc_kms/t_merge_yr/gif are ignored).
# The outer orbit is reconstructed as hyperbolic with the recorded v_inf/beta (n_dyn only
# sets the output cadence, not the dynamics).
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
CSV="${1:?usage: replay.sh LOCKED_MANIFEST.csv [name]}"
ONLY="${2:-}"
[ -f "$CSV" ] || { echo "FATAL: manifest not found: $CSV"; exit 1; }
mkdir -p scenes/survey

tail -n +2 "$CSV" | while IFS=, read -r name a_bin e_bin vinf beta incl Om om nu outcome rest; do
  name="$(echo "$name" | xargs)"
  [ -z "$name" ] && continue
  [ -n "$ONLY" ] && [ "$ONLY" != "$name" ] && continue
  # guard against a manifest that lost the orientation columns
  case "$incl$Om$om$nu" in *[!0-9.eE+-]*|"") echo "  [$name] no orientation in manifest — skipped"; continue;; esac
  OUTER="--outer hyperbolic --beta $beta --v_inf $vinf --n_dyn 2e4"
  echo "=== replay $name  a_bin=$a_bin e_bin=$e_bin v_inf=$vinf beta=$beta  orient=($incl $Om $om $nu)  [was ${outcome:-?}] ==="
  ./run_sim.sh "$name" "$a_bin" "$e_bin" "$OUTER" "" "$incl $Om $om $nu" || echo "  (skipped $name)"
done

MAN="scenes/survey/manifest.csv"
echo "== replay GIFs: $(ls scenes/survey/*.gif 2>/dev/null | tr '\n' ' ') =="
if [ -f "$MAN" ]; then
  echo "== initial conditions -> outcome (reproduced; saved in $MAN) =="
  column -t -s, "$MAN" 2>/dev/null || cat "$MAN"
fi
