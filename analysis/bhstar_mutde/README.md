# μTDE analysis & debris pipeline — v2 (re-run bundle)

Polished, bundle-adapted versions of the five notebooks and the debris pipeline, updated
for the **re-run scattering campaign** (`smbh_bhstar_bundle.h5`: ~1500 realizations/point,
`rp_max = 10·r_t,bin`, converged micro-TDEs). The physics and astrophysics are unchanged;
only the **sim-dependent inputs** were re-derived. Plotting conventions are unified across
every figure (and the debris figures/movies now match the `smbh_movies` layout).

## What changed from v1

| change | where | why |
|---|---|---|
| `rp_max = max(5·r_t, 8·r_g)` → `max(10·r_t, 8·r_g)` | `pkg/scales.py` (single source), and inline in `analyze_mutde_cross_sections`, `analyze_branching`, `mutde_losscone` | the re-run doubled the outer pericenter cut; `b_max` scales with it |
| f_μ fit coefficients `CF` **re-derived** and written to `fmu_fit.json` | `analyze_mutde_cross_sections` writes it; `mutde_losscone` reads it | the semi-analytic β_crit(ρ,e) law is refit on the new data |
| Δt_delay scaling fit **re-derived** | `survey_delay` on the new `delay_per_event.csv` | new events |
| branching "money plot": `b_max`-convergence panel → **appendix**, legend fills the 6th slot | `analyze_branching` | your instruction; the convergence test is a diagnostic, not a result panel |
| β panel floored at **β ≥ 0.1** (the rp10 tidal floor), without global-masking | `analyze_branching` | tidal encounters can't be shallower than `r_t/rp_max = 0.1`; β<0.1 events are only the r_g-dominated M•≳10⁹ cells and stay in the M/a/e/v panels |
| `mutde_cellsummary.csv` now **built from the bundle** | new `pkg/reduce_cellsummary.py` | the per-cell table the cross-section/loss-cone notebooks consume |
| science-region / loss-cone framing note | `analyze_mutde_cross_sections`, loss-cone | report μTDEs as a loss-cone **rate** (converged per-encounter P_μ × physical cut), not a bare flyby-normalized fraction |
| unified style + movies-tar layout on debris figs | `pkg/style.py`, `pkg/debris_tracer.py`, `debris_fig1` | CM serif, ticks-in all 4 sides, minor ticks incl. colorbar, Okabe-Ito/Wong; inset lower-left, title top, time lower-center, scale lower-right |

The debris **tracer physics** did not need adapting: it reads each μTDE event's actual
disruption-time phase space from the bundle, so it carries no rp5/rp10 assumption.

## Layout

```
pkg/
  scales.py             # THE campaign conventions: constants, cell_scales, b_max, rp10, grid_meta cross-check
  style.py              # shared plotting conventions (CM serif, ticks, palettes, inside titles)
  check_bundle.py       # run FIRST: schema + mu-TDE count + rp10 cross-check against grid_meta
  reduce_cellsummary.py # bundle -> mutde_cellsummary.csv   (per-(M,a,e,v) counts; chunked, laptop-safe)
  reduce_delay.py       # bundle -> delay_per_event.csv + delay_cellsummary.csv  (two-peak fallback delay)
  debris_tracer.py      # ballistic debris tracer (Lane-Emden seeds, vectorized integrator, classify, figs)
  run_bundle.py         # run the tracer on mu-TDE rows: trajectories + energy/fate + publication figure (+movie)
  make_test_bundle.py   # synthetic, schema-correct bundle on the REAL grid (for a dry run without the 9 GB file)
notebooks/
  analyze_mutde_cross_sections.ipynb   # sigma_mu(M,a,e,v), heatmaps, the f_mu fit -> fmu_fit.json
  survey_delay.ipynb                   # Delta t_delay maps + scaling fit
  analyze_branching.ipynb              # outcome branching "money plot" (+ convergence appendix)
  mutde_losscone.ipynb                 # loss-cone RATE: reads fmu_fit.json; R(<z), Gamma(M)
  debris_fig1.ipynb                    # Figure 1: one mu-TDE's debris, colored by eps_BH2 (+ 4 insets)
```

## Run order

```bash
cd pkg

# 0. sanity-check the real bundle (phase-space columns? rp10? mu-TDE count?)
python3 check_bundle.py ../../../data/bhstar/campaign2_rp10/smbh_bhstar_bundle.h5

# 1. reduce the bundle to the small tables the notebooks consume (run ONCE)
python3 reduce_cellsummary.py /path/to/smbh_bhstar_bundle.h5    # -> mutde_cellsummary.csv
python3 reduce_delay.py       /path/to/smbh_bhstar_bundle.h5    # -> delay_per_event.csv (+ cellsummary)
```

The notebooks read the bundle and `delay_per_event.csv` from
`data/bhstar/campaign2_rp10/` (relative to `notebooks/`; override the bundle with `SMBH_BUNDLE`).
Run the reductions from that data folder so their outputs land there, then run:

1. **`analyze_mutde_cross_sections.ipynb`** — refits `CF` and **writes `fmu_fit.json`**.
2. **`mutde_losscone.ipynb`** — reads `fmu_fit.json` (falls back to the v1 CF if absent) and the rp10 convention.
3. **`survey_delay.ipynb`** — re-derives the Δt_delay fit on the new events.
4. **`analyze_branching.ipynb`** — the money plot (reads the bundle directly).
5. **`debris_fig1.ipynb`** — set `SMBH_BUNDLE` or edit `BUNDLE`; pick a row, render Figure 1.
   (`run_bundle.py` is the batch/CLI version with trajectories + a movie.)

The notebooks find `pkg/` automatically (they add `../pkg`, `pkg`, `.` to `sys.path`). Run
them from the repo root or from `notebooks/` with `pkg/` alongside.

## One thing to verify on the real bundle

The debris notebooks (`debris_fig1`, `run_bundle`, `reduce_delay`) need the disruption-time
**phase-space columns** `x_smbh_AU … vz_b2_kms` (plus `m1, m2, R2, v_inf`). `check_bundle.py`
reports whether they are present. If the re-run scattering CSVs were written **without** the
final Cartesian state dump, those three products cannot run and the campaign would need to be
re-dumped with full final-state output — the cross-section / branching / loss-cone products do
**not** need the phase space and work regardless.

## Dry run (no 9 GB bundle)

```bash
cd pkg
python3 make_test_bundle.py smbh_bhstar_bundle.h5   # synthetic, real grid, schema-correct
python3 reduce_cellsummary.py smbh_bhstar_bundle.h5
python3 reduce_delay.py smbh_bhstar_bundle.h5
# then run the notebooks against these files
```

The synthetic bundle assigns outcomes by a known β_crit(ρ,e) law, so the cross-section
notebook's refit should ~recover it — a built-in self-consistency check. It is a shape/plumbing
stand-in only; science numbers come from the real bundle.

Requires: `numpy, scipy, h5py, matplotlib`, plus `scikit-learn` (survey_delay's random-forest
importance cell) and `pandas` if present. CM fonts fall back to DejaVu Serif if `cmr10` is absent.
