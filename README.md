# smbh-binary — binary + SMBH three-body encounters

[![CI](https://github.com/Kkritos/smbh-binary/actions/workflows/ci.yml/badge.svg)](https://github.com/Kkritos/smbh-binary/actions/workflows/ci.yml)
[![version](https://img.shields.io/badge/version-1.0.0-blue)](CHANGELOG.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data DOI](https://img.shields.io/badge/data-10.5281%2Fzenodo.24372716-blue)](https://doi.org/10.5281/zenodo.24372716)

Code and analysis for **Kritos & Perna (2026), *Tidally-Induced Coda: Double Transients from Black Hole-Star Binary Encounters with Massive Black Holes*** (see [Citing](#citing)).

Three-body scattering of a stellar binary (BH–star, star–star) with a supermassive (or
intermediate-mass) black hole, integrated with [SpaceHub](https://github.com/YihanWangAstro/SpaceHub)
including post-Newtonian terms (1PN + 2PN + 2.5PN).

- **BH–star campaign** — outcome cross sections σ(outcome | M•, a_bin, e_bin, v∞/v_orb), with a
  focus on **micro-TDEs** (the star disrupted by its BH companion during the SMBH encounter):
  rates via the loss cone, the f_μ semi-analytic fit, debris fallback and the double-flare delay.
- **Star–star campaign** — collision fractions compared with **Yu & Lai (2024)**, Newtonian vs PN.
- **Movies** — animations of each outcome class (ejection/EMRI, μTDE debris, flyby, …).

The simulation outputs are archived on Zenodo: **[doi:10.5281/zenodo.24372716](https://doi.org/10.5281/zenodo.24372716)**.
With those files, every figure can be regenerated on a laptop without re-running the simulations.

---

## Reproducing the figures (quick start)

**1. Get the code and the Python packages** (Python ≥ 3.10):

```bash
git clone https://github.com/Kkritos/smbh-binary.git && cd smbh-binary
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

**2. Get the data** into `data/` (the downloader checks each file's checksum; files already
present are skipped):

```bash
python3 scripts/get_data.py --list                         # groups, sizes, what is present
python3 scripts/get_data.py starstar bhstar                # or: ... all
```

| group | size | needed for |
|---|---|---|
| `starstar` | ~130 MB | star–star / Yu & Lai figures |
| `bhstar` | ~10 GB | all BH–star / μTDE notebooks |
| `bhstar_rp5` | ~5 GB | only the rp5-vs-rp10 appendix of `analyze_branching` (skipped if absent) |
| `debris` | ~14 GB | the single-encounter debris movie and still (`mutde_lab.py`) |

Without the script: download the files from the Zenodo page and place them as in
[Data layout](#data-layout).

**3. Run the notebooks.** Each notebook reads `data/` by relative path, so start Jupyter from the
notebook's own folder (or run headless with `jupyter nbconvert --to notebook --execute --inplace <nb>`):

| folder | notebook (run in this order) | figures it writes (PDF + PNG, next to the notebook) | time\* |
|---|---|---|---|
| `analysis/yulai2024/` | `comparison_plots.ipynb` | `collision_fraction_{deep,shallow}` | 5 s |
| `analysis/bhstar_mutde/notebooks/` | 1. `analyze_mutde_cross_sections.ipynb` | `mutde_*`, `fmu_*`, `sigma_*`, `sigmav_rate`, `rate_gamma_fit_*`, `horizon6`, `convergence`; writes `fmu_fit.json`, `gamma6.npz` | 1.5 min |
| | 2. `mutde_losscone.ipynb` (reads `fmu_fit.json`) | `gammaM_fit`, `horizon_inset` | 15 s |
| | 3. `survey_delay.ipynb` | `delay_*` | 5 min |
| | 4. `analyze_branching.ipynb` | `branching_all`, `branching_strong`, `branching_convergence`, `branching_leakage_ea` | 1 min |
| | 5. `debris_fig1.ipynb` | `debris_fig1` | 2 min |

\* On an Apple M-series laptop. The BH–star notebooks stream the 10 GB bundle in chunks.

**4. The single-encounter debris movie and final still** (`analysis/debris_mutde/figures/`):

```bash
T=data/debris_mutde/s1/bss_hyperbolic_bh10_star1_i3.06878_full.txt
python3 analysis/debris_mutde/mutde_lab.py --traj $T --out data/debris_mutde/movie       # mutde.mp4 (+ gif), 1e4 debris
python3 analysis/debris_mutde/mutde_lab.py --traj $T --still --ndeb 100000 \
    --out data/debris_mutde/still                                                        # final_still.{png,pdf}, money_plot.png
```

Both reproduce the committed files bit for bit. Each run reads the whole 14 GB trajectory
into memory: **~27 GB RAM**, 3–8 min.

The committed figures in the repository are the outputs of exactly these steps, so you can
compare against them with `git diff` after a re-run (PDFs differ only in their creation date).

---

## Data layout

`data/` is not tracked by git; the Zenodo files go here:

```
data/starstar/                smbh_starstar_collisions_{Newt,PN}_results.csv   (312k encounters each)
data/bhstar/campaign2_rp10/   smbh_bhstar_bundle.h5    (main BH-star campaign: rp_max = 10 r_t,
                                                        ~1500 realizations per grid cell)
                              delay_per_event.csv, delay_cellsummary.csv  (from pkg/reduce_delay.py)
data/bhstar/campaign1_rp5/    smbh_bhstar_bundle_750.h5 (first campaign: rp_max = 5 r_t, 750/cell)
data/debris_mutde/s1/         bss_hyperbolic_bh10_star1_i3.06878_full.txt  (SpaceHub trajectory of
                                                        the mutde_lab encounter, SMBH + BH + star)
```

The notebooks also accept other locations through the environment variables `SMBH_BUNDLE`
(campaign 2 bundle) and `OLD_BUNDLE` (campaign 1 bundle); `SMBH_DATA` moves the whole `data/` root
for `get_data.py` and the scripts.

The delay tables can be rebuilt from the bundle (outputs land in the current directory):
`cd data/bhstar/campaign2_rp10 && python3 ../../../analysis/bhstar_mutde/pkg/reduce_delay.py smbh_bhstar_bundle.h5`
(see [`analysis/bhstar_mutde/README.md`](analysis/bhstar_mutde/README.md)).

---

## Repository layout

```
THEORY.md               tutorial on the physics and numerical methods behind the project
analysis/yulai2024/     Newt vs PN collision fractions vs Yu & Lai (2024); digitised Yu & Lai data
analysis/bhstar_mutde/  μTDE cross sections, branching, loss-cone rate, delay, debris figure
  notebooks/              the five notebooks and their figures
  pkg/                    shared Python: campaign conventions, plotting style, bundle reductions,
                          the ballistic debris tracer
analysis/debris_mutde/  mutde_lab.py: single-encounter μTDE + debris-fallback lab; figures/
simulation/             the C++ three-body code (smbh_binary.cpp) and its Makefile
campaigns/bhstar/       BH-star parameter sweep, SLURM scripts, CSV -> HDF5 bundle tools
campaigns/starstar/     star-star (Yu & Lai) sweep and SLURM scripts
movies/                 encounter-movie pipeline (see movies/README.md)
scripts/get_data.py     downloads the Zenodo data into data/
scripts/zenodo/         maintainers: publishing the data on Zenodo
tests/                  test suite (pytest), run by CI
env.sh                  sets SMBH_ROOT, SPACEHUB, SMBH_BIN, SMBH_DATA, SMBH_SCRATCH; reads env.local.sh
```

---

## Running the simulations

Only needed to produce new data; the figures need only the Zenodo files.

### Build

SpaceHub is an external header-only library. The results were produced with SpaceHub commit
`190777d6`:

```bash
git clone https://github.com/YihanWangAstro/SpaceHub.git external/SpaceHub
git -C external/SpaceHub checkout 190777d6
make -C simulation          # -> simulation/bin/smbh_binary_PN and smbh_binary_Newt (-DNEWTONIAN_ONLY)
```

The Makefile looks for SpaceHub in `external/SpaceHub/src` (gitignored); elsewhere, use
`make -C simulation SPACEHUB=/path/to/SpaceHub/src`. Needs a C++17 compiler.

```bash
simulation/bin/smbh_binary_PN --help
simulation/bin/smbh_binary_PN --header     # CSV column names
```

One run = one encounter, one `CSV,...` line on stdout. `--mode full` also writes the
trajectory (`bss_*_full.txt`) and energy error (`bss_*_full.err`) to the working directory.
`--restart FILE` seeds a run from Cartesian phase space (e.g. the final-state columns of an
earlier run).

### Campaigns

**Machine-specific settings** live in `env.local.sh` (gitignored), read by `env.sh`:

```bash
cp env.local.sh.example env.local.sh     # then replace the <...> placeholders
```

| variable | what | default |
|---|---|---|
| `SMBH_SCRATCH` | run directories, joblists, results and SLURM logs | `<repo>/scratch` |
| `SPACEHUB` | SpaceHub headers | `<repo>/external/SpaceHub/src` |
| `SMBH_SLURM_ACCOUNT` | allocation charged (`sbatch -A`) | none |
| `SMBH_SLURM_PARTITION` | partition for production jobs (~1 day, one full node) | cluster default |
| `SMBH_SLURM_PARTITION_TEST` | partition for short test jobs (≤ 4 h) | cluster default |

The `.slurm` files contain no paths, accounts or partitions; the submit scripts pass them to
`sbatch` (through the `slurm_submit` helper in `env.sh`). The core counts and wall times in
the `#SBATCH` lines (`--cpus-per-task`, `-t`) assume 96-core nodes; adjust them to your nodes.

The sweep generators (`campaigns/*/smbh_*_sweep.py`) write one command line per encounter.

- **Locally:** `simulation/run_local.sh <binary> <joblist> <out.csv> [ncores]` runs a joblist in
  parallel and writes `Results/<out.csv>`.
- **On a SLURM cluster:** each campaign has a `submit_*.sh` that builds the joblist in
  `$SMBH_SCRATCH/<name>`, shuffles it and submits the jobs:
  ```bash
  make -C simulation
  campaigns/bhstar/submit_bhstar_test.sh      # small b_max test        -> $SMBH_SCRATCH/smbh_bhstar_bmaxtest
  campaigns/bhstar/submit_bhstar_array.sh     # production, 8-way array -> $SMBH_SCRATCH/smbh_bhstar_encounters
  campaigns/starstar/submit_test.sh           # star-star tests (1000 runs, Newt + PN)
  campaigns/starstar/submit_all.sh            # star-star production    -> $SMBH_SCRATCH/smbh_starstar_collisions
  ```
  Logs go to `logs/` inside each run directory. A one-off job, e.g. the timing test:
  `source env.sh && slurm_submit simulation/timing_test.slurm "$SMBH_SCRATCH/logs" "$SMBH_SLURM_PARTITION_TEST"`.

After a BH-star run: merge the shards (footer of `smbh_bhstar_array.slurm`), then
`campaigns/bhstar/bundle_to_h5.py` (CSV → HDF5 bundle) and `reduce_bundle.py`. Put the bundle in
`data/bhstar/<campaign>/`.

---

## Provenance and known gaps

- **Code versions.** Git tag `starstar-yulai-code` is the older `smbh_binary.cpp` (35-column CSV,
  r_start = max(100 r_t, 100 rp, 5 b)) that produced the star–star / Yu & Lai CSVs. The current
  code adds `--restart`, 18 final phase-space columns, r_start = max(50 r_t, 10 rp) and an energy
  check in the escape criterion; it produced the BH–star bundles and the movies.
- The raw BH–star scattering CSVs are not archived; the HDF5 bundles are the archival product.
- The `debris` trajectory was made with the older binary, so re-running `mutde_lab.py --run` with
  the current one will not reproduce it bit for bit; use `--traj` as above.
- `analysis/yulai2024/reproduce_fig5_Yu_Lai_2024.ipynb` is an early version whose input CSVs no
  longer exist; it is kept for its saved outputs.

---

## Tests

```bash
pip install pytest && python3 -m pytest
```

The tests need no data: they check the campaign geometry conventions, the polytrope sampling of
the debris, the trajectory reader, the bundle tools on a small synthetic bundle, the consistency
of the notebooks and version strings, and — if `simulation/bin` is built — run short encounters
with both binaries. GitHub Actions builds the code against SpaceHub `190777d6` and runs them on
Linux and macOS for every push ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

## Versioning

Releases follow [semantic versioning](https://semver.org): the version is in [`VERSION`](VERSION),
changes are listed in [`CHANGELOG.md`](CHANGELOG.md), and each release is a git tag `vX.Y.Z`.
Version 1.0.0 is the code used for the paper. The tag `starstar-yulai-code` marks the older
`smbh_binary.cpp` that produced the star–star data.

## Citing

If you use this code or data, please cite the paper (GitHub's *Cite this repository* button uses
[`CITATION.cff`](CITATION.cff)):

> K. Kritos and R. Perna (2026), *Tidally-Induced Coda: Double Transients from Black Hole-Star Binary Encounters with Massive Black Holes*, in preparation.

and the data set: doi:[10.5281/zenodo.24372716](https://doi.org/10.5281/zenodo.24372716).
The integrator is [SpaceHub](https://github.com/YihanWangAstro/SpaceHub) (Wang et al. 2021, MNRAS).

## Acknowledgments

The packaging of this repository (documentation, tests, continuous integration and release
tooling) was prepared with the assistance of [Claude Code](https://claude.com/claude-code)
(Anthropic).

## License

Code: [MIT](LICENSE). Data on Zenodo: CC BY 4.0.
