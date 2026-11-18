# SmbhBinary — binary + SMBH three-body encounters

Three-body scattering of a stellar binary (BH–star, star–star) with a supermassive (or
intermediate-mass) black hole, integrated with [SpaceHub](https://github.com/YihanWangAstro/SpaceHub)
including post-Newtonian terms (1PN + 2PN + 2.5PN). Science goals:

- **BH–star campaign** — outcome cross sections σ(outcome | M•, a_bin, e_bin, v∞/v_orb), with a
  focus on **micro-TDEs** (the star disrupted by its BH companion during the SMBH encounter):
  rates via the loss cone, the f_μ semi-analytic fit, debris fallback and the double-flare delay.
- **Star–star campaign** — collision fractions compared with **Yu & Lai (2024)**, Newtonian vs PN.
- **Movies** — talk-quality animations of each outcome class (ejection/EMRI, μTDE debris, flyby, …).

## Layout

```
env.sh                  source this: SMBH_ROOT, SPACEHUB, SMBH_BIN, SMBH_DATA, SMBH_SCRATCH
simulation/             the C++ code (smbh_binary.cpp), Makefile, local parallel runner, timing test
campaigns/bhstar/       BH-star sweep generator, SLURM array + submit scripts, bundle tools
campaigns/starstar/     star-star (Yu & Lai) sweep generator, SLURM jobs, submit script
movies/                 encounter-movie pipeline (see movies/README.md)
analysis/bhstar_mutde/  μTDE cross sections, branching, loss-cone rate, delay, debris (see its README)
analysis/debris_mutde/  mutde_lab.py: single-encounter μTDE + debris-fallback lab
analysis/yulai2024/     Newt vs PN collision fractions vs Yu & Lai (2024); digitised data
data/                   (gitignored) large data, same layout on Mac and SeaWulf
archive/                (gitignored, Mac only) superseded analyses and old outputs
```

`data/` layout:

```
data/bhstar/campaign2_rp10/   smbh_bhstar_bundle.h5 (~9.9 GB; rp_max = 10 r_t, ~1500 real./cell),
                              delay_per_event.csv, delay_cellsummary.csv
data/bhstar/campaign1_rp5/    smbh_bhstar_bundle_750.h5 (first campaign, rp_max = 5 r_t, 750/cell),
                              mutde_cellsummary.csv, grid_meta.csv, ecc_meta.csv
data/starstar/                smbh_starstar_collisions_{Newt,PN}_results.csv (312k rows each)
data/debris_mutde/            s1/ (16 GB SpaceHub dump of the mutde_lab encounter), still/, movie/
```

## Build

SpaceHub is an external header-only library (not in this repo). The Makefile finds it at
`~/codes/SpaceHub/src` (Mac) or `~/codes/SpaceHub/src` (SeaWulf); override with
`make SPACEHUB=/path/to/SpaceHub/src`.

```bash
make -C simulation          # -> simulation/bin/smbh_binary_PN and smbh_binary_Newt (-DNEWTONIAN_ONLY)
simulation/bin/smbh_binary_PN --help
simulation/bin/smbh_binary_PN --header     # CSV column names
```

One run = one encounter, one `CSV,...` line on stdout. `--mode full` also writes the
trajectory (`bss_*_full.txt`) and energy error (`bss_*_full.err`) to the working directory.
`--restart FILE` seeds a run from Cartesian phase space (e.g. the final-state columns of an
earlier run).

## Running campaigns

**SeaWulf (SLURM).** Each campaign has a `submit_*.sh` that generates the joblist in scratch,
shuffles it and submits the jobs. Scratch dirs: `$SMBH_SCRATCH/<name>` with
`SMBH_SCRATCH=/path/to/scratch/$USER`.

```bash
cd /path/to/smbh-binary && make -C simulation
campaigns/bhstar/submit_bhstar_test.sh      # small b_max test   -> $SMBH_SCRATCH/smbh_bhstar_bmaxtest
campaigns/bhstar/submit_bhstar_array.sh     # production, 8-way array -> $SMBH_SCRATCH/smbh_bhstar_encounters
campaigns/starstar/submit_all.sh            # Newt + PN jobs     -> $SMBH_SCRATCH/smbh_binary_interactions
```

SLURM scripts read `SMBH_ROOT` from the environment exported by the submit script (fallback:
`/path/to/smbh-binary`). The `#SBATCH -o/-e` log paths are literal
(`/path/to/scratch/...`) because SLURM does not expand variables there — edit them if
the scratch location changes.

After a BH-star run: merge the shards (footer of `smbh_bhstar_array.slurm`), then
`campaigns/bhstar/bundle_to_h5.py` (CSV → HDF5 bundle) and `reduce_bundle.py`
(→ `mutde_cellsummary.csv`). Put the bundle in `data/bhstar/<campaign>/`.

**Mac (local).** `simulation/run_local.sh <binary> <joblist> <out.csv> [ncores]` runs a joblist
in parallel with xargs and writes `Results/<out.csv>`.

## Analysis

- `analysis/bhstar_mutde/` — run order in its README. Notebooks run from `notebooks/` and read
  `data/bhstar/campaign2_rp10/` by relative path; override with `SMBH_BUNDLE` / `OLD_BUNDLE`.
  `pkg/*.py` scripts take the bundle path as an argument.
- `analysis/yulai2024/comparison_plots.ipynb` — run from its folder; reads `data/starstar/`.
- `analysis/debris_mutde/mutde_lab.py --help` — `--run` (needs `simulation/bin/smbh_binary_PN`),
  `--traj data/debris_mutde/s1/bss_..._full.txt` to reuse the stored dump, or `--demo`.

Python: numpy, scipy, matplotlib, h5py, pandas, scikit-learn, astropy, tqdm
(movies also: pillow, imageio-ffmpeg).

## Provenance and known gaps

- **Code versions.** Git tag `starstar-yulai-code` is the older `smbh_binary.cpp` (35-column CSV,
  r_start = max(100 r_t, 100 rp, 5 b)) that produced the star-star / Yu & Lai CSVs. The current
  code (from SeaWulf) adds `--restart`, 18 final phase-space columns, r_start = max(50 r_t, 10 rp)
  and an energy check in the escape criterion; it produced the BH-star bundles and movies.
- The raw BH-star scattering CSVs live only on SeaWulf scratch (if at all); the bundles are the
  archival product.
- `data/debris_mutde/s1` was made with the older binary, so re-running `mutde_lab.py --run`
  with the current one will not reproduce it bit for bit.
- `analysis/yulai2024/reproduce_fig5_Yu_Lai_2024.ipynb` is an early version whose input CSVs no
  longer exist; it is kept for its saved outputs.
- The SpaceHub commit on SeaWulf has not been checked against the Mac's (`190777d6`).
