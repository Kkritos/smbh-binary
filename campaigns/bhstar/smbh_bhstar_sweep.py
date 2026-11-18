#!/usr/bin/env python3
"""
smbh_bhstar_sweep.py
--------------------
Joblist generator for the BH-star + SMBH hyperbolic (near-parabolic) scattering
campaign: outcome cross sections sigma(outcome | M_smbh, a_bin, e_bin, v_inf/v_orb).

Inner binary: BH (10 Msun) + star (1 Msun, 1 Rsun).  Full PN run (smbh_binary_PN).
Channels (all handled by the C++ pair_event; no C++ changes needed):
  - TDE:star2_by_bh1   -> mu-TDE: star disrupted by its BH companion
  - TDE:star2_by_SMBH  -> star disrupted by the SMBH (auto-off above the Hills mass)
  - PLUNGE / EXCHANGE / FLYBY / BOUND_TRIPLE / IONIZATION / COLLISION as usual.

ECCENTRICITY IS A GRID AXIS (not Monte-Carlo'd).  For each (M_smbh, a_bin, e_bin,
frac) cell we run N_b scattering experiments that differ only in the sampled impact
parameter and angles; e_bin is FIXED at a grid node.  This yields e-resolved cross
sections sigma(X | M,a,e,frac).  The thermal (or ANY) eccentricity average is done
at POST-PROCESSING by a weighted sum over the e nodes,
      sigma(X | M,a,frac) = sum_e w_e * sigma(X | M,a,e,frac),
with w_e proportional to the chosen p(e) (thermal p(e)=2e) integrated over the
SURVIVING range [0, e_max(a)].  Gridding + summing is stratified sampling: for the
same total runs it is at least as accurate as drawing e randomly, and it is
reweightable to any p(e) without re-running.

SURVIVAL CUT: a BH-star binary whose inner pericenter a(1-e) is inside the star's
tidal-disruption radius by the BH, r_t,inner = R_star (M_bin/m_star)^(1/3) ~ 0.0103 AU,
is "born disrupted" and not a real binary.  e_max(a) = 1 - r_t,inner/a; nodes with
e >= e_max(a) are SKIPPED (not simulated).  Set SURVIVAL_CUT=False to simulate them
anyway (they prompt-mu-TDE at t=0).

GRID:
  - M_smbh : 19 half-decade points, 10 ... 1e10   (comparable-mass -> IMBH -> SMBH)
  - a_bin  : 13 points, 0.03 ... 30 AU
  - frac   : 7 points, 0.001 ... 10               (hard -> soft/fast)
  - e_bin  : 12 nodes 0 ... 0.99 (survival-cut per a_bin); nodes denser near 1 so the
             thermal weight (largest at high e) is resolved
  - N_b    : 750 impact-parameter draws per surviving cell

Sampling per run: b uniform in b^2 in [0, b_max^2] -> sigma = f * pi * b_max^2;
cos(iota) ~ U(-1,1); Omega,omega ~ U(0,2pi); inner phase uniform in MEAN anomaly
(time) -> true anomaly nu (general for 0<=e<1, reduces to uniform nu at e=0).

Longer integration: n_dyn backstop 2e4 (slurm wall timeout 120 s).
Start radius fix (see C++): r_start = max(50 r_t_bin, 10 rp_out), which skips the
dynamically-empty infall that used to dominate the cost of slow encounters.

Outputs: stdout = one argument string per run; grid_meta.csv = b_max per (M,a,frac)
grid point (b_max is e-independent); ecc_meta.csv = e_max(a) survival limit per a_bin.

ENV CONTROLS:
  BHSTAR_PILOT=0          -> production grid (submit scripts set this).
  BHSTAR_TEST=1           -> b_max convergence test grid (few cells, overrides PILOT).
  BHSTAR_RP_TIDAL_FAC=20  -> enlarge rp_max = max(FAC*r_t_bin, 8 r_g) to push b_max out
                             until the mu-TDE (and bound-triple / plunge-star) cross
                             sections converge; default 5.0.  grid_meta.csv records the
                             resulting b_max, so post-processing stays consistent IF the
                             same factor is used there.
"""

import os
import sys
import math
import random

# ----------------------------------------------------------------------
# run mode: production / pilot / b_max-convergence test
# ----------------------------------------------------------------------
PILOT = os.environ.get("BHSTAR_PILOT", "1").strip().lower() not in ("0", "false", "no", "off")
TEST  = os.environ.get("BHSTAR_TEST",  "0").strip().lower() not in ("0", "false", "no", "off")

if TEST:
    # b_max convergence test: representative cells spanning M, a, frac, e.
    # Run with a large BHSTAR_RP_TIDAL_FAC; one pass contains the full b/b_max CDF.
    M_SMBH_GRID = [1e3, 1e6, 1e9]                    # Msun : low / mid / high
    A_BIN_GRID  = [0.1, 1.0, 10.0]                   # AU   : close / mid / wide
    FRAC_GRID   = [0.01, 1.0]                         # v_inf / v_orb : hard / ~v_orb
    E_GRID      = [0.0, 0.9]                          # circular and high-e
    N_PER_POINT = 2000                               # impact-parameter draws per surviving cell
elif PILOT:
    M_SMBH_GRID = [1e6]                              # Msun
    A_BIN_GRID  = [0.1]                              # AU
    FRAC_GRID   = [0.01]                             # v_inf / v_orb
    E_GRID      = [0.03, 0.5, 0.9]
    N_PER_POINT = 200                                # per (M,a,e,frac) cell
else:
    # 19 half-decade masses: 10, 10^1.5, 1e2, ... 1e10
    M_SMBH_GRID = [10.0 ** (1.0 + 0.5 * k) for k in range(19)]
    # 13 semimajor axes, 0.03 -> 30 AU (floor set by Roche lobe + mu-TDE radius)
    A_BIN_GRID  = [0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0]
    # 7 velocity ratios; 0.3/1/3/10 reach the SOFT / fast regime (v_inf >~ v_orb)
    FRAC_GRID   = [0.001, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0]
    # 12 eccentricity nodes, uniform to 0.9 then denser near 1 (thermal weight peaks high)
    E_GRID      = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    N_PER_POINT = 1500                                # impact-parameter draws per surviving cell

SEED = 42

# ----------------------------------------------------------------------
# integration backstop (passed as --n_dyn; C++ default is 1e4)
# ----------------------------------------------------------------------
N_DYN = 2.0e4

# ----------------------------------------------------------------------
# fixed inner-binary parameters: BH + star
# ----------------------------------------------------------------------
M_BH       = 10.0                    # Msun
M_STAR     = 1.0                     # Msun
R_STAR_RSUN = 1.0                    # Rsun (solar star)
M_BIN      = M_BH + M_STAR           # 11 Msun

V1    = 29.7847        # km/s  == sqrt(G Msun / AU)
C_KMS = 299792.458     # km/s
RSUN_KM = 695700.0
AU_KM   = 1.495978707e8
R_BH_RSUN = (2.0 * (V1 * V1 / (C_KMS * C_KMS)) * M_BH) * AU_KM / RSUN_KM  # 2 r_g(M_BH) in Rsun

RTOL = 1e-13

# ----------------------------------------------------------------------
# b_max control (rp_max encloses the tidal-interaction and plunge channels)
# RP_TIDAL_FAC is overridable: enlarge it until the tidally-driven channels
# (mu-TDE, bound triple, plunge-star) converge in the b/b_max CDF.
# ----------------------------------------------------------------------
RP_TIDAL_FAC  = float(os.environ.get("BHSTAR_RP_TIDAL_FAC", "5.0"))
RP_PLUNGE_FAC = 2.0

# ----------------------------------------------------------------------
# eccentricity survival cut: keep inner pericenter a(1-e) outside r_t,inner
# ----------------------------------------------------------------------
SURVIVAL_CUT = True          # False -> also simulate e >= e_max (born-disrupted) cells
E_MAX_CAP    = 0.995         # numerical ceiling on e_max (lets the e=0.99 node in for wide a)
RSUN_AU         = RSUN_KM / AU_KM
R_TIDE_INNER_AU = R_STAR_RSUN * RSUN_AU * ((M_BH + M_STAR) / M_STAR) ** (1.0 / 3.0)  # ~0.01034 AU


def e_max_for(a_bin_AU):
    """Largest inner eccentricity whose pericenter a(1-e) still clears r_t,inner."""
    if not SURVIVAL_CUT:
        return 1.0
    return max(0.0, min(1.0 - R_TIDE_INNER_AU / a_bin_AU, E_MAX_CAP))


def nu_from_uniform_time(rng, e):
    """Random inner phase uniform in mean anomaly (time) -> true anomaly."""
    Manom = rng.uniform(0.0, 2.0 * math.pi)
    E = Manom if e < 0.8 else math.pi          # starting guess
    for _ in range(100):
        dE = (E - e * math.sin(E) - Manom) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-13:
            break
    return 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(0.5 * E),
                            math.sqrt(1.0 - e) * math.cos(0.5 * E))


def v_orb_kms(a_bin_AU):
    return V1 * math.sqrt(M_BIN / a_bin_AU)


def r_tide_bin_AU(a_bin_AU, M_smbh):
    return a_bin_AU * (M_smbh / M_BIN) ** (1.0 / 3.0)


def r_g_AU(M_tot):
    return (V1 * V1 / (C_KMS * C_KMS)) * M_tot


def grid_point_scales(M_smbh, a_bin_AU, frac):
    m_tot   = M_smbh + M_BIN
    v_orb   = v_orb_kms(a_bin_AU)
    v_inf   = frac * v_orb
    r_t_bin = r_tide_bin_AU(a_bin_AU, M_smbh)
    r_g     = r_g_AU(m_tot)
    r_foc   = 2.0 * V1 * V1 * m_tot / (v_inf * v_inf)      # 2 G M_tot / v_inf^2, AU
    rp_max  = max(RP_TIDAL_FAC * r_t_bin, RP_PLUNGE_FAC * 4.0 * r_g)
    b_max2  = rp_max * rp_max + r_foc * rp_max
    b_max   = math.sqrt(b_max2)
    rp_floor = 2.0 * r_g
    return dict(m_tot=m_tot, v_orb=v_orb, v_inf=v_inf, r_t_bin=r_t_bin,
                r_g=r_g, r_foc=r_foc, rp_max=rp_max, b_max=b_max, b_max2=b_max2,
                rp_floor=rp_floor)


def beta_from_b2(b2, r_foc, r_t_bin, rp_floor):
    rp = (-r_foc + math.sqrt(r_foc * r_foc + 4.0 * b2)) / 2.0
    if rp < rp_floor:
        rp = rp_floor
    return r_t_bin / rp


def main():
    rng = random.Random(SEED)
    meta_rows = []
    n_total = 0
    n_cells = 0

    for M_smbh in M_SMBH_GRID:
        for a_bin in A_BIN_GRID:
            emax = e_max_for(a_bin)                       # survival limit (depends on a only)
            for frac in FRAC_GRID:
                s = grid_point_scales(M_smbh, a_bin, frac)
                # b_max is e-independent -> one meta row per (M,a,frac)
                meta_rows.append(
                    f"{M_smbh:g},{a_bin:g},{frac:g},{s['v_orb']:.6g},{s['v_inf']:.6g},"
                    f"{s['r_t_bin']:.6g},{s['r_g']:.6g},{s['rp_max']:.6g},{s['b_max']:.6g},{N_PER_POINT}"
                )
                for e_bin in E_GRID:
                    if e_bin >= emax:                     # born-disrupted -> skip
                        continue
                    n_cells += 1
                    for _ in range(N_PER_POINT):
                        b2    = rng.random() * s["b_max2"]
                        beta  = beta_from_b2(b2, s["r_foc"], s["r_t_bin"], s["rp_floor"])
                        cosi  = rng.uniform(-1.0, 1.0)
                        incl  = math.acos(cosi)
                        Omega = rng.uniform(0.0, 2.0 * math.pi)
                        omega = rng.uniform(0.0, 2.0 * math.pi)
                        nu    = nu_from_uniform_time(rng, e_bin)

                        print(
                            f"--M_smbh {M_smbh:g} --m1 {M_BH:g} --m2 {M_STAR:g} "
                            f"--type1 bh --type2 star "
                            f"--R1 {R_BH_RSUN:.6e} --R2 {R_STAR_RSUN:.6e} "
                            f"--a_bin {a_bin:g} --e_bin {e_bin:.10g} "
                            f"--outer hyperbolic --beta {beta:.10g} --v_inf {s['v_inf']:.10g} "
                            f"--incl {incl:.10g} --Omega {Omega:.10g} --omega {omega:.10g} --nu {nu:.10g} "
                            f"--n_dyn {N_DYN:g} "
                            f"--rtol {RTOL:g} --mode classify --report 0"
                        )
                        n_total += 1

    with open("grid_meta.csv", "w") as f:
        f.write("M_smbh,a_bin_AU,frac,v_orb_kms,v_inf_kms,r_t_bin_AU,r_g_smbh_AU,rp_max_AU,b_max_AU,N\n")
        f.write("\n".join(meta_rows) + "\n")

    with open("ecc_meta.csv", "w") as f:
        f.write("a_bin_AU,e_max_survival,n_e_nodes_used\n")
        for a_bin in A_BIN_GRID:
            emax = e_max_for(a_bin)
            nused = sum(1 for e in E_GRID if e < emax)
            f.write(f"{a_bin:g},{emax:.6g},{nused}\n")

    mode = "TEST" if TEST else ("PILOT" if PILOT else "PRODUCTION")
    sys.stderr.write(
        f"[smbh_bhstar_sweep] mode={mode}  RP_TIDAL_FAC={RP_TIDAL_FAC:g}  RP_PLUNGE_FAC={RP_PLUNGE_FAC:g}\n"
        f"[smbh_bhstar_sweep] (M,a,frac) points={len(meta_rows)}  "
        f"surviving (M,a,e,frac) cells={n_cells}  N_b/cell={N_PER_POINT}  total runs={n_total}\n"
        f"[smbh_bhstar_sweep] M_smbh: {len(M_SMBH_GRID)} pts {M_SMBH_GRID[0]:g}..{M_SMBH_GRID[-1]:g}  "
        f"a_bin: {len(A_BIN_GRID)} pts {A_BIN_GRID[0]:g}..{A_BIN_GRID[-1]:g} AU  "
        f"frac: {len(FRAC_GRID)} pts {FRAC_GRID[0]:g}..{FRAC_GRID[-1]:g}\n"
        f"[smbh_bhstar_sweep] e grid ({len(E_GRID)}): {E_GRID}  survival-cut e<e_max(a), cap={E_MAX_CAP}\n"
        f"[smbh_bhstar_sweep] r_t,inner={R_TIDE_INNER_AU:.5g} AU  "
        f"(e_max: a=0.03->{e_max_for(0.03):.3f}, a=0.1->{e_max_for(0.1):.3f}, a=1->{e_max_for(1.0):.3f})\n"
        f"[smbh_bhstar_sweep] n_dyn backstop={N_DYN:g}; thermal weight w_e ~ 2e applied at post-processing\n"
        f"[smbh_bhstar_sweep] wrote grid_meta.csv (b_max per M,a,frac) and ecc_meta.csv (e_max per a)\n"
    )


if __name__ == "__main__":
    main()
