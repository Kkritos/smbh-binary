# ============================================================================
# config.sh  --  settings for the SpaceHub (smbh_binary_PN) movie pipeline.
# Sourced by run_scene.sh / make_movies.sh.  Defaults match Kostas' setup.
# ============================================================================

# ---- executable: built by `make -C simulation` in the repo root ---------------
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/env.sh"
BIN="${BIN:-$SMBH_BIN/smbh_binary_PN}"

# ---- how the full run dumps a trajectory -------------------------------------
# `--mode full` makes smbh_binary_PN write the full trajectory to a FILE (not
# stdout) via SpaceHub's DefaultWriter:  bss_<outer>_<type1><m1>_<type2><m2>_i<incl>_full.txt
# in the run's working directory.  stdout only carries the summary + classify CSV,
# so run_scene.sh runs each full-mode stage in its own dir and grabs that file.
# `--report 0` just silences the stdout summary (the file is written regardless).
TRAJ_FLAGS="--mode full --report 0"
TRAJ_GLOB="bss_*_full.txt"   # DefaultWriter output picked up after a full-mode run
OUTCOME_FIELD=22          # column # of 'outcome' in the classify CSV line (CSV,...)

# ---- fixed inner-binary bodies (BH + star) -----------------------------------
M_BH=10;  M_STAR=1;  TYPE1=bh;  TYPE2=star
R1=-1;    R2=-1                     # -1 = auto (point BH, 1 Rsun star)
M_SMBH=1e6
RTOL=1e-13

# ---- search: how many random orientations to try to hit each target outcome --
N_SEARCH=800              # orientations tried per target (higher helps rare outcomes)

# ---- EMRI restart stage (after ejection) -------------------------------------
EMRI_N_ORBITS=8           # captured-orbit periods to integrate (EMRI orbits, no merger)
EMRI_R_ESC_AU=1e8         # keep integrating after the star flies off (no early stop)

# ---- TDE debris stage (micro-TDE and SMBH-TDE) -------------------------------
NDEB=10000                # number of debris tracers drawn from an n=3 polytrope
MUTDE_N_ORBITS=18         # fallback periods to follow — long enough that bound debris RETURN
MUTDE_NFRAMES=420         # fallback frames (smoothness; the long baseline comes from N_ORBITS)

# ---- ionization: keep integrating past breakup, then time-lapse --------------
IONIZ_EXTEND=3.0          # keep integrating this * (stage-1 duration) so bodies fly apart
IONIZ_ESC_AU=1e12         # escape radius during the extension (disables the early stop)
IONIZ_NFRAMES=340         # frames in the time-lapsed movie
IONIZ_ACCEL=2.6           # time-warp exponent (>1: slow disruption, fast separation)

# ---- eccentric-flyby target (circular binary -> high-e surviving binary) -----
FLYBY_ECC_MIN=0.8         # search orientations for a FLYBY whose surviving e >= this

# ---- derived CSV field offsets (order AFTER 'outcome' is version-independent) -
# outcome(+0) exit_mode(+1) a_new_AU(+2) e_new(+3) v_esc_kms(+4) beta_actual(+5) ...
ANEW_OFF=2; ENEW_OFF=3; VESC_OFF=4   # add to OUTCOME_FIELD to read final conditions

# ---- output / tools ----------------------------------------------------------
OUTDIR="scenes"
PY=python3
FPS=30                    # smoother playback

# ---- self-heal: SeaWulf libexpat/pyexpat symbol mismatch ---------------------
# matplotlib imports plistlib -> xml.parsers.expat -> pyexpat.  If an old libexpat
# (< 2.4.0) is ahead on LD_LIBRARY_PATH (e.g. from a venv/module), pyexpat fails with
#   "undefined symbol: XML_SetBillionLaughsAttackProtectionMaximumAmplification"
# which breaks ALL matplotlib imports.  Fix: force the matching system libexpat via
# LD_PRELOAD.  We only adopt a candidate that actually makes `import expat` succeed.
if ! $PY -c "import xml.parsers.expat" >/dev/null 2>&1; then
  for _xpat in /usr/lib64/libexpat.so.1 /lib64/libexpat.so.1 \
               /usr/lib/x86_64-linux-gnu/libexpat.so.1 /usr/lib/libexpat.so.1; do
    if [ -e "$_xpat" ] && LD_PRELOAD="$_xpat${LD_PRELOAD:+:$LD_PRELOAD}" \
         $PY -c "import xml.parsers.expat" >/dev/null 2>&1; then
      export LD_PRELOAD="$_xpat${LD_PRELOAD:+:$LD_PRELOAD}"
      echo "note: worked around libexpat mismatch via LD_PRELOAD=$_xpat"
      break
    fi
  done
  if ! $PY -c "import xml.parsers.expat" >/dev/null 2>&1; then
    echo "WARNING: could not auto-fix the libexpat/pyexpat mismatch for $PY."
    echo "  matplotlib (rendering) will fail.  Diagnose with:"
    echo "    ldd \$($PY -c 'import pyexpat,os;print(pyexpat.__file__)') | grep expat"
    echo "    echo \$LD_LIBRARY_PATH"
    echo "  then LD_PRELOAD the libexpat.so.1 that has the missing symbol (>= 2.4.0)."
  fi
fi

# ---- debris FALLBACK-CURVE stage (physical dM/dt, long integration for the tail) ----
# Added as an extra, non-invasive step after each TDE scene's movie: seeds debris from
# the SAME disruption state + debris_ic (same rng seed) and integrates far longer than
# the movie to resolve the late-time fallback (full -> t^-5/3; moving potential -> ~t^-9/4,
# Coughlin & Nixon 2019).  The movies themselves are unchanged.
FALLBACK_ENABLE=1          # 1 = also compute the fallback dM/dt curve for each TDE scene
FALLBACK_N_ORBITS=600      # central (SMBH) fallback periods to follow  -> sets the tail length
FALLBACK_N_ORBITS_BH=50    # companion-BH mini-disk periods to follow (micro-TDE only)
FALLBACK_NTIME=1500        # log-time samples across the encounter integration (memory ~ N*NTIME)
FALLBACK_TIMEOUT=900       # max seconds for one scene's fallback stage; it is skipped past this
                           # (a safeguard so one slow encounter can't stall the whole survey)
