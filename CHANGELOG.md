# Changelog

All notable changes to this project are documented here. The project follows
[Semantic Versioning](https://semver.org): the version is in `VERSION`, and each release is
a git tag `vX.Y.Z`.

## [1.0.0] — unreleased

First public release, accompanying Kritos & Perna (2026), *Tidally Induced Coda: Double Transients from Black Hole–Star Binary Encounters with Massive Black Holes*.

- Three-body code `simulation/smbh_binary.cpp` (SpaceHub, 1PN + 2PN + 2.5PN), built against
  SpaceHub commit `190777d6`.
- BH–star and star–star campaign generators and SLURM scripts (site settings in `env.local.sh`).
- Analysis notebooks for the μTDE cross sections, loss-cone rate, fallback delay, branching
  ratios and debris figure; the Yu & Lai (2024) comparison; the single-encounter debris lab.
- `scripts/get_data.py` to fetch the simulation data from Zenodo; `scripts/zenodo/` to publish it.
- Test suite (`tests/`) and continuous integration (build + tests on every push).
- `THEORY.md`: a self-contained tutorial on the physics (orbits, black holes, PN gravity, stars, binaries,
  three-body encounters, tidal disruption, debris) and numerical methods (SpaceHub's AR-chain integrator,
  Monte Carlo cross sections, loss-cone rates) behind the project.

Earlier history: the tag `starstar-yulai-code` marks the older version of `smbh_binary.cpp`
that produced the star–star collision data.
