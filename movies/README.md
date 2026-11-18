# SMBH encounter movies — automated pipeline (smbh_binary_PN)

One button builds three looping GIFs **and mp4s** for talks, running the **real**
`smbh_binary_PN` executable end to end — set initial conditions, integrate, stop/restart
to resolve the EMRI, then render:

1. **`ejection.gif`** — Hills ejection + capture: star flung out as a hypervelocity
   star, BH left bound to the SMBH; a **restart run** then resolves the bound BH
   orbiting the SMBH for several precessing EMRI orbits (continuous across the seam).
2. **`mutde.gif`** — micro-TDE on a **near-parabolic** BH–SMBH encounter: the star (an
   **n=3 polytrope**) is driven to its tidal radius around its BH companion and disrupts;
   ~10⁴ debris tracers are seeded from the executable's exact disruption state (bulk
   velocity + spin-up to breakup, Evans & Kochanek 1989) and followed ballistically in
   the time-dependent BH+SMBH field. Because the encounter is ~parabolic the debris split
   roughly half **bound** (fall back to the SMBH over several orbits) / half **ejected**.
   Main panel = SMBH-centred fallback stream; inset = the receding BH + its mini-disk;
   debris coloured by specific energy wrt the SMBH (blue bound, red ejected) with a colorbar.
3. **`flyby.gif`** — a distant passage pumps the eccentricity of an initially circular
   binary, which survives.

Rendering (all scenes): a fixed encounter-type title, a white **scale bar** in both the
main window and the inset, no transient text, GIF **+ mp4** output. Orbits are resampled
by **angle swept around the SMBH** (not uniform time) so eccentric pericenter passages
are smooth and well resolved. In μTDE scenes the energy colorbar sits on the **left**
(the binary/debris come in from the right) and the inset is bottom-left.

Everything runs **inside this folder** (`movies/`); the executable is found via the
repo's `env.sh` automatically.

---

## One button

```bash
pip install numpy scipy matplotlib pillow imageio-ffmpeg   # once (imageio-ffmpeg = mp4s)
./make_movies.sh                 # all three scenes via smbh_binary_PN
```
`imageio-ffmpeg` bundles a static ffmpeg so the mp4s export with no system ffmpeg or
`module load`. Without it the GIFs still build and the mp4 is skipped with a note.
Output: `scenes/ejection.gif`, `scenes/mutde.gif`, `scenes/flyby.gif`.

Other entry points:
```bash
./make_movies.sh ejection        # just one scene
./make_movies.sh survey          # OUTCOME GALLERY: many hyperbolic encounters (see below)
./make_movies.sh montage         # tile scenes/survey/*.gif into one grid video (gif + mp4)
./make_movies.sh --probe         # short full-mode run; shows the bss_*_full.txt file + parse check
./make_movies.sh demo            # built-in integrator (no executable) — for styling
./run_scene.sh flyby --dry-run   # print every smbh_binary_PN command, run nothing
```

## Outcome gallery (`survey.sh`)

`./make_movies.sh survey` (or `./survey.sh`) fires a set of **hyperbolic** BH–star
encounters with the SMBH, spanning outcomes as `(a_bin, e_bin, v_inf, penetration β)`
are varied, each with a **random orientation**. Masses are fixed (m_BH=10, m_star=1,
m_SMBH=1e6). Each encounter is **dispatched by whatever actually happens**:

- `EXCHANGE_bh1` → ejection + capture, then an EMRI restart (the **BH** integrated for
  `EMRI_N_ORBITS` orbits about the SMBH); the ejected **star's** HVS speed `v_esc` and the
  EMRI **merger time** from Peters (1964), `t_merge`, are recorded
- `EXCHANGE_star2` → the **star** is captured; a restart integrates the **star's** orbit
  about the SMBH for several orbits (the BH is ejected); the runaway **BH** `v_esc` recorded
- `TDE:star2_by_bh1` → micro-TDE debris fallback (`debris_follow.py --disruptor bh`),
  integrated long enough (`MUTDE_N_ORBITS` fallback periods) that bound debris **return** —
  to the SMBH (main panel) and onto the **stellar BH** itself (the mini-disk in the inset;
  the run prints `bound(BH)` / `accreted(BH)`). A two-phase timestep (fine while the BH disk
  forms, coarse for the long fallback) keeps the long follow affordable.
- `TDE:star2_by_SMBH` → **SMBH-TDE** debris stream (`--disruptor smbh`): the star is
  disrupted by the SMBH itself (r_t ~ 0.46 AU), giving the classic ~half-bound/half-ejected
  t⁻⁵ᐟ³ fallback stream around the SMBH, with an inner-region inset
- `IONIZATION` → **keep integrating** the three now-unbound bodies so they visibly separate,
  then play back with **accelerated (time-lapse) time** — slow through the disruption, fast
  as they fly apart (`timewarp.py`)
- `FLYBY_ECC` → a FLYBY whose **surviving** binary eccentricity `e_new` ≥ `FLYBY_ECC_MIN`
  (search a **circular** binary until it is pumped eccentric but stays bound)
- `FLYBY` / `BOUND_TRIPLE` / `PLUNGE` / `COLLISION` / … → the trajectory, labelled by outcome

The default table has **25 encounters** (a tidy **5×5 montage**), weighted toward
**micro-TDEs and exchanges** while still spanning the outcome space: **12 micro-TDEs**
(across a/e/v/β), **6 exchanges** (3 BH-capture→EMRI + 3 star-capture), plus 2 SMBH-TDEs,
2 ionizations, 2 flybys (one circular→eccentric), and a bound triple. Every row is
**targeted**, so the composition is guaranteed (no dependence on a random draw). The whole
run is long (the μTDE/SMBH-TDE rows each integrate + render 10⁴ debris) — run one at a time
with `./survey.sh <name>`, or drop `NDEB` in `config.sh` for a quick pass.

Output: `scenes/survey/<name>_<OUTCOME>.gif` (+ `.mp4`), a per-sim human-readable
`scenes/survey/<name>/summary.txt` (the executable's `--report 1` block), **plus
`scenes/survey/manifest.csv`** recording initial conditions AND final conditions for every
encounter: `name, a_bin, e_bin, v_inf_kms, beta, incl, Omega, omega, nu, outcome,
a_new_AU, e_new, v_esc_kms, t_merge_yr, gif` (orientation angles + HVS speed, surviving/
captured orbit, and EMRI Peters time), printed as a table at the end of the run. Edit the `SIMS` table at the top
of `survey.sh` to add rows or change parameters; a row with a **target** outcome searches
orientations until that outcome occurs (so the gallery reliably contains it), a row with a
**blank** target keeps whatever the random orientation produces. Run one row with
`./survey.sh <name>`.

## Reproducing a gallery you like (`replay.sh`)

`survey.sh` searches **random, unseeded** orientations and **resets `manifest.csv`** on
every full run, so a re-run never reproduces a previous gallery. To keep one you like:

```bash
cp scenes/survey/manifest.csv survey_locked.csv   # BEFORE re-running anything
# ... later, after editing/upgrading the pipeline ...
./replay.sh survey_locked.csv                      # same ICs, re-rendered
./replay.sh survey_locked.csv mutde_mid            # just one row
```

`replay.sh` reads each row's recorded orientation (`incl, Omega, omega, nu`) and
`a_bin, e_bin, v_inf, beta` and re-runs that **exact** initial condition with **no
search** (masses are fixed in `config.sh`), so the encounters come back identical — now
rendered with the current pipeline. `run_sim.sh` also accepts a 6th argument
`"incl Omega omega nu"` to lock a single run's orientation directly.

## What each scene run does (`run_scene.sh`)

1. **Search** — random orientations in `--mode classify` until *your* classifier
   reports the target outcome (`EXCHANGE_bh1` / `TDE:star2_by_bh1` / `FLYBY`). The
   outcome is therefore guaranteed by the code, not assumed.
2. **Stage-1 trajectory** — re-run that exact IC with `--mode full`, which makes
   `smbh_binary_PN` write the full trajectory to a **file** (`bss_*_full.txt`, via
   SpaceHub's `DefaultWriter`) in the run's directory; `parse_traj.py` converts it to a
   canonical CSV. (The run is done inside its own subdir so the file is easy to pick up;
   stdout carries only the summary.)
3. **Post-step**
   - **ejection** → `make_restart.py` writes a restart file from stage-1's *exact
     final state* (3 body lines: SMBH, BH, star); `smbh_binary_PN --restart` then
     integrates the bound BH about the SMBH for a few orbital periods. The seam is
     continuous because stage 2 starts from stage 1's final state.
   - **mutde** → `debris_follow.py` reads the executable's exact 3-body disruption state
     (`bss_*_full.txt`), samples `NDEB` tracers from an n=3 polytrope scaled to R⋆, gives
     each the star's bulk velocity + solid-body spin-up to breakup Ω=√(Gm⋆/R⋆³) about the
     orbital angular momentum, then follows them as test particles in the time-dependent
     BH+SMBH field for `MUTDE_N_ORBITS` SMBH orbits (accretion-freezing tames the stiff
     bound orbits). The executable classifies the μTDE but can't emit 10⁴ particles, so
     this fallback is ballistic; the frozen-in energy spread arises self-consistently
     from each tracer's depth in the BH potential.
   - **flyby** → nothing extra.
4. **Render** — looping GIF + mp4 via `animate.py`.

## Setup — usually nothing to change

- **Executable**: `config.sh` sources the repo's `env.sh` and uses
  `$SMBH_BIN/smbh_binary_PN` (build it with `make -C simulation`). Override `BIN` if needed.
- **Trajectory output**: `--mode full` writes the trajectory to a **file**
  `bss_<outer>_<type1><m1>_<type2><m2>_i<incl>_full.txt` (SpaceHub `DefaultWriter`) in
  the run's cwd — **not** stdout (`--report` only toggles the stdout summary). The file
  is SpaceHub's dump: a header `time,id,mass,radius,px,py,pz,vx,vy,vz`, then one line per
  body per step (ids 0,1,2 = SMBH, BH, star). Units are SpaceHub-internal: positions are
  already **AU**, time is **internal** (`years = time/2π`), velocities are in units of
  `√(GM⊙/AU) = 29.78 km/s`. `parse_traj.py` handles all of this. Run
  `./make_movies.sh --probe` to see the file's first lines and a parse check.
- **Outcome column**: `OUTCOME_FIELD=22` in `config.sh` (position of `outcome` in the
  classify CSV line) — the same field your sweep uses.

Scene physics knobs live in the `case` block at the top of `run_scene.sh`
(a_bin, e_bin, outer-orbit flags, target outcome). The EMRI restart length is
`EMRI_N_ORBITS`; the debris count and follow time are `NDEB` / `MUTDE_N_ORBITS`
(all in `config.sh`).

## A note on the EMRI — orbits, not merger

The restart integrates the captured BH for `EMRI_N_ORBITS` orbital periods, showing the
bound, **precessing** EMRI orbit with the real PN forces, and stops short of coalescence.
Integrating to merger is **not** feasible directly: a realistic Hills capture
(a≈400–700 AU, e≈0.99) inspirals on a ~Myr–Gyr GW timescale — millions of orbits — so
reaching 4 r_g in tractable wall time isn't possible on the true orbit (physics, not a
code limit). The **merger time itself** is instead estimated analytically from the
captured orbit's (a_new, e_new) via Peters (1964) — `peters.py` — and written to
`manifest.csv` as `t_merge_yr` (typically ~10⁸–10⁹ yr for these captures).

## Files

| file | role |
|---|---|
| `make_movies.sh` | one-button driver (all scenes / one / `survey` / `--probe` / `demo`) |
| `survey.sh` | outcome gallery: table of hyperbolic encounters → `run_sim.sh` |
| `replay.sh` | reproduce an exact gallery from a saved `manifest.csv` (locked orientations) |
| `run_sim.sh` | one encounter, dispatched by its ACTUAL outcome (EMRI / debris / trajectory); optional locked orientation |
| `run_scene.sh` | the three preset scenes end-to-end against `smbh_binary_PN` |
| `config.sh` | executable path, trajectory flag, search/EMRI settings |
| `parse_traj.py` | `bss_*_full.txt` (SpaceHub dump) → canonical CSV (pivot by id, unit conversion) |
| `make_restart.py` | stage-1 final state → restart file (+ BH period) |
| `concat_stages.py` | glue stage-1 + stage-2 (EMRI) into one movie |
| `debris_follow.py` | TDE debris (`--disruptor bh`=μTDE / `smbh`=SMBH-TDE): n=3 polytrope + spin-up, ballistic follow → CSV + `.debris.npz` |
| `peters.py` | Peters (1964) GW merger time `t_merge(a,e,m1,m2)` for the EMRI |
| `timewarp.py` | stitch canonical CSVs + resample in accelerating time (ionization fly-apart) |
| `animate.py` / `scene_lib.py` | canonical CSV (+ debris npz) → looping GIF + mp4 (scale bars, auto-placed inset, energy colorbar) |
| `nbody.py`, `demo_gen_*.py`, `mutde_gen.py`, `flyby_gen.py` | built-in `demo` path |

## Canonical trajectory format

`t,smbh_x,smbh_y,bh_x,bh_y,star_x,star_y` — AU and yr, one row per frame, absent bodies
`nan`. Sidecar `<name>.events.json` carries `center`, `zoom_on`, the `inset` config, and
(μTDE) a `debris` block pointing at `<name>.debris.npz` (arrays `xy` [Nframes,Ndeb,2],
`energy` [Ndeb], `start_frame`). Anything that writes these can be animated.

## Requirements

Python 3.9+, `numpy scipy matplotlib pillow`; **ffmpeg** for the mp4s (if missing, the
GIF is still written and the mp4 is skipped with a note). Fast single-node runs — no
SLURM. Render on an interactive node if you build many at once (10⁴-tracer μTDE ≈ a
couple of minutes).

## Reproducing a previous run (orientation replay)

By default each survey row searches RANDOM orientations until it hits its target outcome, so
the exact realization (and the debris fractions) differ run-to-run even at fixed a/e/beta/v_inf.
To reproduce a specific past run, survey.sh reads `orientations.csv` (name,incl,Omega,omega,nu)
and LOCKS each row to those exact Euler angles via run_sim.sh's replay path (its 6th argument),
skipping the search. The shipped orientations.csv holds the angles of the 25-sim reference run.

  ./survey.sh                 # uses orientations.csv if present -> reproduces that run
  REPLAY=0 ./survey.sh        # ignore the file, go back to random search
  ORIENT_FILE=other.csv ./survey.sh      # use a different replay table

Angles are stored to the precision of the source (here ~6 sig figs); with the same executable
build and rtol this reproduces the setup and outcome. A 3-body scattering is formally chaotic,
so for bit-exact trajectories use each scene's raw s1/bss_*_full.txt (full-precision Cartesian
IC) instead. The angles were recovered from the per-run summary.txt CSV line (fields 17-20 =
incl,Omega,omega,nu; field 16 is the fixed value pi, not an angle).

## Debris fallback-rate stage (fallback_stage.py)

An extra, non-invasive step run after each TDE scene's movie (hooked into run_sim.sh and
run_scene.sh; toggle with FALLBACK_ENABLE in config.sh). It seeds debris from the SAME
disruption state and the SAME debris_ic / rng seed the movie uses. The movie output is untouched.

Why a separate stage: the movie's follow() batch-"settles" the BH1 mini-disk to stay cheap,
and the old Dt_delay used the Keplerian energy->period relation. In a time-varying potential
(a moving stellar BH near the SMBH) specific energy is not conserved, so that inversion breaks
down (Coughlin & Nixon 2019, ApJL 883 L17). Each parcel's fallback time is the time of its
FIRST return to pericenter about its host.

How it stays both correct and fast: the potential is time-varying only while the BH is close to
the SMBH (the "encounter", ~ r_peri / v_rel). We integrate SMBH+BH+ALL debris through that brief
encounter -- which captures the real energy redistribution and the fast mini-disk returns
directly -- and then, once the BH has departed and the field about each host is static, continue
any not-yet-returned parcel ANALYTICALLY along its now-fixed Kepler ellipse (exact once static;
no thousands of static orbits are integrated). So the Keplerian relation is used ONLY in the
regime where it holds, never during the time-varying phase. The analytic return-time helper is
validated against direct integration to ~1e-10 (test_kepler.py in the source scratchpad).

Two channels, two separate lines in the figure:
  * dM/dt onto BH1 (companion): the mini-disk -- fast returns, caught during the encounter.
  * dM/dt onto BH2 (central SMBH): the slow fallback -- returns continued analytically, so the
    late-time tail is resolved out to the longest-period bound debris.
The plot (<scene>.fallback.png) shows both curves (over a pale 'total'), vertical peak markers,
t^-5/3 and t^-9/4 reference slopes anchored at the central peak (full disruption vs
surviving-core / moving-potential), and a logarithmic-slope panel. Data -> <scene>.fallback.csv.

Accuracy of the frozen-in approximation (measured vs reduce_delay-style). The same figure
overlays, as thin DASHED curves in the channel colours, the OLD frozen-in Keplerian prediction:
each parcel assigned a fallback time t_K = 2 pi G M_host (2|eps_host,0|)^{-3/2} from its
DISRUPTION-INSTANT energy (a fixed Kepler orbit about a static host). So the solid (direct,
measured first-return) and dashed (frozen-in) dM/dt sit on the same axes. The title prints
Dt_delay direct vs frozen; a box prints the per-parcel ratio t_measured/t_Kepler (median,
16-84%) SEPARATELY for the BH and SMBH channels -- the SMBH (second-flare) channel is the
physically robust one; the BH mini-disk is softening-sensitive. The CSV now has 7 columns:
t_yr, dMdt_BH, dMdt_SMBH, dMdt_tot, dMdt_BH_Kfrozen, dMdt_SMBH_Kfrozen, dMdt_tot_Kfrozen,
so you can quantify the error offline and aggregate it across the survey's TDE scenes.

Scalar comparison outputs: each TDE scene also writes `<scene>.compare.csv` (one row:
t_peak and Dt_delay, direct vs frozen-in, their ratio, and the per-channel t_meas/t_Kepler
with 16-84%). After a full survey, survey.sh collects them all into
`scenes/survey/delay_comparison.csv` (one row per TDE scene) and prints the table -- that is
the micro-TDE "simulation vs frozen-in delay" comparison across the whole set.

Debris specific-energy distributions vs time (<scene>.energy.png / .energy.npz): eps_1 (wrt BH1,
blue) and eps_2 (wrt BH2, orange), same colour scheme as the fallback curves, as step-PDFs on a
symlog energy axis, at several snapshots following the disruption (time labelled in each panel).
Because the potential is time-varying, the distributions morph during the encounter and then
freeze once it is static -- the snapshots make that visible. Integration is bounded to the
encounter epoch for the same reason as above.

Robustness: the encounter window is set from the local dynamical time at the BH-SMBH
separation (never from sep/velocity, which diverges for a slow/near-marginal encounter and
used to balloon the window and hang the solver), and is hard-capped. The fallback subsamples
to <=4000 tracers and uses FALLBACK_NTIME time samples so memory stays bounded (the movie keeps
the full N). Each scene's fallback call is wrapped in `timeout FALLBACK_TIMEOUT` (default 900s),
so a pathological encounter is skipped and the survey continues instead of stalling.

Config knobs (config.sh): FALLBACK_ENABLE, FALLBACK_TIMEOUT (per-scene cap, s), FALLBACK_NTIME (time samples across the encounter
integration, default 5000). FALLBACK_N_ORBITS / FALLBACK_N_ORBITS_BH are still passed on the
command line but are now ignored (the window is set self-consistently by the encounter time).
The movie mini-disk settle time is --settle-orbits in debris_follow.py (default 40, so default
movies are unchanged); lower it to freeze/hold the inset disk after fewer orbits.

Standalone (any stage-1 trajectory):
  python3 fallback_stage.py <stage1_bss_full.txt> <out_prefix> --disruptor bh --ndeb 10000
writes <out_prefix>.fallback.{csv,png} and <out_prefix>.energy.{png,npz}.
