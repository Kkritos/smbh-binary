#!/usr/bin/env python3
# =============================================================================
#  scales.py  --  the ONE place that knows the sim geometry conventions.
#  -----------------------------------------------------------------------------
#  Every analysis product (cross sections, branching, loss-cone rate) that
#  reconstructs the sampling geometry (b_max, r_foc, r_g, r_t, beta <-> b/b_max)
#  imports these helpers, so the campaign conventions live in exactly one file.
#
#  WHAT CHANGED FOR THE v2 (re-run) CAMPAIGN
#  -----------------------------------------------------------------------------
#  The first campaign used an outer cut  rp_max = max(5 r_t,bin, 8 r_g).
#  The re-run (higher statistics, converged mu-TDEs) uses  rp_max = max(10 r_t,bin, 8 r_g).
#  That single factor (5 -> 10) doubles b_max in the tidally-limited regime and
#  therefore rescales every flyby-normalised fraction and cross section.  It is
#  set here once, as RP_MAX_OVER_RT, and -- more robustly -- can be read straight
#  from the bundle's /grid_meta table (rp_max_AU, b_max_AU) with
#  bmax_from_grid_meta(), which is the source of truth when the bundle is present.
#
#  Units: G=1, [mass]=Msun, [length]=AU, [vel]=V1=sqrt(G Msun/AU)=29.7847 km/s.
# =============================================================================
import numpy as np

# ---- constants (shared by all notebooks) ------------------------------------
V1      = 29.7847                 # km/s  (= sqrt(G Msun / AU))
C_KMS   = 299792.458              # speed of light [km/s]
M_BIN   = 11.0                    # inner binary mass m1+m2 = 10 + 1  [Msun]
RSUN_AU = 695700.0 / 1.495978707e8
AU2_TO_PC2 = 1.0 / 4.2545e10      # 1 pc^2 = 4.2545e10 AU^2

# stellar tidal radius for the mu-TDE condition: r_TDE = R_star (m1/m_star)^(1/3)
#   = 1 Rsun * (10/1)^(1/3) = 0.01002 AU
R_TDE   = 1.0 * RSUN_AU * (10.0 / 1.0)**(1.0/3.0)

# the 7 canonical v_inf/v_orb nodes
FRAC_NODES = np.array([1e-3, 1e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0])

# ---- THE campaign convention (v2 re-run) ------------------------------------
RP_MAX_OVER_RT = 10.0             # v2: rp_max = max(10 r_t,bin, 8 r_g)   (was 5 in v1)
RG_FLOOR_FACTOR = 8.0             # pericenter interaction floor at 8 r_g
BETA_TIDAL_FLOOR = 1.0 / RP_MAX_OVER_RT   # = 0.1 : tidal encounters cannot be shallower than this


def cell_scales(M, a, frac, rp_max_over_rt=RP_MAX_OVER_RT):
    """Per-cell geometric scales, vectorised.  Returns (r_t, r_foc, rp_max, b_max^2),
    all in AU (b_max^2 in AU^2).  frac = v_inf/v_orb with v_orb = V1 sqrt(m_bin/a).

      r_t    = a (M/m_bin)^(1/3)                    binary tidal radius wrt the central BH
      v_inf  = frac * V1 * sqrt(m_bin/a)
      r_foc  = 2 V1^2 (M+m_bin) / v_inf^2           gravitational-focusing length
      r_g    = (V1/c)^2 (M+m_bin)                   gravitational radius of the central BH
      rp_max = max(rp_max_over_rt * r_t, 8 r_g)     outer pericenter cut  (v2: factor 10)
      b_max^2= rp_max^2 + r_foc rp_max              focused impact-parameter cut
    """
    M = np.asarray(M, float); a = np.asarray(a, float); frac = np.asarray(frac, float)
    m_tot = M + M_BIN
    r_t   = a * (M / M_BIN)**(1.0/3.0)
    v_inf = frac * V1 * np.sqrt(M_BIN / a)
    r_foc = 2.0 * V1 * V1 * m_tot / v_inf**2
    r_g   = (V1 * V1 / (C_KMS * C_KMS)) * m_tot
    rp_max = np.maximum(rp_max_over_rt * r_t, RG_FLOOR_FACTOR * r_g)
    b2    = rp_max**2 + r_foc * rp_max
    return r_t, r_foc, rp_max, b2


def b_max_AU(M, a, frac, rp_max_over_rt=RP_MAX_OVER_RT):
    """Convenience: b_max [AU] per cell (sqrt of the focused cut)."""
    _, _, _, b2 = cell_scales(M, a, frac, rp_max_over_rt)
    return np.sqrt(b2)


def frac_from_vinf(v_inf, a):
    """Snap v_inf [km/s] to the nearest canonical v_inf/v_orb node (log distance)."""
    v_inf = np.asarray(v_inf, float); a = np.asarray(a, float)
    v_orb = V1 * np.sqrt(M_BIN / a)
    r = v_inf / v_orb
    return FRAC_NODES[np.argmin(np.abs(np.log10(r)[..., None] - np.log10(FRAC_NODES)), axis=-1)]


def bb_from_beta(beta, M, a, frac, rp_max_over_rt=RP_MAX_OVER_RT):
    """Per-encounter b/b_max reconstructed from the stored penetration beta=r_t/r_p:
       r_p = r_t/beta, b^2 = r_p^2 + r_foc r_p, b/b_max = sqrt(b^2/b_max^2) in [0,1].
    This is the exact inverse of the b^2-uniform sampling, so (b/b_max)^2 ~ U[0,1]."""
    r_t, r_foc, _, b2max = cell_scales(M, a, frac, rp_max_over_rt)
    r_p = r_t / np.asarray(beta, float)
    b2  = r_p**2 + r_foc * r_p
    return np.sqrt(np.clip(b2 / b2max, 0.0, 1.0))


# ---- grid_meta: the source of truth when the bundle is present --------------
def read_grid_meta(bundle_path):
    """Return the /grid_meta table of the bundle as a dict of arrays.
    Columns (v2): M_smbh, a_bin_AU, frac, v_orb_kms, v_inf_kms, r_t_bin_AU,
                  r_g_smbh_AU, rp_max_AU, b_max_AU, N."""
    import h5py
    with h5py.File(bundle_path, "r") as f:
        g = f["grid_meta"]
        return {k: g[k][:] for k in g.keys()}


def check_rpmax_convention(bundle_path, tol=0.05):
    """Cross-check RP_MAX_OVER_RT against the bundle's stored grid_meta.
    Returns (ratio_median, ok): ratio = rp_max_AU / r_t_bin_AU over tidally-limited
    cells (where rp_max is set by r_t, not 8 r_g).  Warns if it disagrees with
    RP_MAX_OVER_RT by more than `tol`."""
    gm = read_grid_meta(bundle_path)
    rp, rt = gm["rp_max_AU"], gm["r_t_bin_AU"]
    rg = gm.get("r_g_smbh_AU", np.zeros_like(rp))
    tidal = rp > (RG_FLOOR_FACTOR * rg) * 1.001        # cells where r_t sets the cut
    ratio = np.median((rp / rt)[tidal]) if tidal.any() else np.median(rp / rt)
    ok = abs(ratio - RP_MAX_OVER_RT) <= tol * RP_MAX_OVER_RT
    if not ok:
        print(f"[scales] WARNING: grid_meta rp_max/r_t = {ratio:.3f} but "
              f"RP_MAX_OVER_RT = {RP_MAX_OVER_RT}; set the constant to match the bundle.")
    return ratio, ok
