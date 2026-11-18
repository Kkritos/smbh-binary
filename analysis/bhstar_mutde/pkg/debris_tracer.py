#!/usr/bin/env python3
# =============================================================================
#  debris_tracer.py
#  -----------------------------------------------------------------------------
#  Ballistic stellar-debris tracer for a micro-TDE (mu-TDE) in a BH-star + central-BH
#  encounter.  After the 1 Msun star is tidally disrupted by its 10 Msun companion
#  black hole (BH1), we replace the star with N massless "debris" test particles and
#  follow their purely ballistic motion in the time-dependent potential of the two
#  point masses (BH1 and the central massive BH, BH2).  We keep the FULL trajectory
#  of every particle and, at post-processing, work out which debris end up bound to
#  BH1, bound to BH2, or unbound.
#
#  Physical picture (Rees 1988, Nature 333, 523; Evans & Kochanek 1989, ApJ 346, L13):
#  ------------------------------------------------------------------------------
#  At the instant of disruption the star is torn apart and each fluid element's
#  orbit is "frozen in": it keeps the position and velocity it had at that moment.
#  Because different elements sit at slightly different depths in BH1's potential
#  (position offset dr across the stellar radius R*), they acquire a SPREAD in
#  specific orbital energy
#        Deps  ~  (G m_BH1 / r_disrupt^2) * R*          (the leading, "geometric" term)
#  This spread (NOT the stellar self-binding) is what sorts the debris into bound
#  and unbound material and, for a classic TDE, gives the ~flat dM/deps -> t^{-5/3}
#  fallback.  We realise it by (a) sampling each debris' POSITION inside the star
#  from the stellar density profile (a polytrope), and (b) giving each the bulk
#  stellar velocity PLUS a solid-body spin omega x dr, with omega the star's
#  instantaneous orbital angular velocity about BH1 at disruption (tidal locking).
#  BOTH ingredients matter: spin alone with every particle at the stellar centre
#  would give almost no energy spread.
#
#  Why the debris are MASSLESS test particles:
#  ------------------------------------------------------------------------------
#  The debris have no self-gravity in this model and do not back-react on the BHs,
#  so each one moves independently in the EXTERNAL field of BH1+BH2.  The two BHs
#  therefore follow their own 2-body orbit, identical regardless of the debris, and
#  we can advance ALL debris together in a single integration (no debris-debris
#  coupling).  (Caveat: the real debris carry ~1 Msun ~ 0.1 m_BH1, so neglecting
#  their pull on BH1 is a ~10% approximation; debris-stream self-gravity is also
#  neglected.  Both are standard for a first ballistic pass.)
#
#  Units (identical to the SpaceHub scattering code, so states can be swapped in):
#  ------------------------------------------------------------------------------
#        G        = 1
#        [mass]   = Msun
#        [length] = AU
#        [time]   = sqrt(AU^3 / (G Msun)) = 1 yr / (2*pi) ~ 0.159 yr
#        [vel]    = AU / [time] = V1 = 29.7847 km/s   (= sqrt(G Msun / AU))
#  A 1 Msun point mass at 1 AU has orbital period 2*pi in these units (= 1 yr).
# =============================================================================

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

# -----------------------------------------------------------------------------
#  Physical constants / unit conversions (for labelling and the stellar radius).
# -----------------------------------------------------------------------------
V1_KMS   = 29.7847            # velocity unit in km/s  (= sqrt(G Msun / AU))
AU_KM    = 1.495978707e8      # 1 AU in km
RSUN_AU  = 6.957e5 / AU_KM    # 1 Rsun in AU  (~0.004652 AU)
YR_PER_T = 1.0 / (2.0*np.pi)  # internal time unit in years  (T -> yr: multiply by this)
SEC_PER_YR = 3.155815e7
KMS_PER_V1 = V1_KMS           # velocity: internal -> km/s

# =============================================================================
#  1.  POLYTROPE (Lane-Emden) + debris initial conditions
# =============================================================================
def lane_emden(n, xi_max=20.0, n_pts=20000):
    """Integrate the Lane-Emden equation for polytropic index n and return the
    dimensionless radius xi, the solution theta(xi), and the first zero xi1
    (the dimensionless stellar surface, theta(xi1)=0).

    Lane-Emden:   (1/xi^2) d/dxi (xi^2 dtheta/dxi) = -theta^n ,
    with theta(0)=1, theta'(0)=0.  The density profile is rho ~ theta^n, so the
    enclosed-mass integrand is  dM/dxi ~ xi^2 theta^n.
    n=3 is the standard model for a Sun-like (radiative) star; n=1.5 for a fully
    convective / low-mass star.  Only n=0,1,5 have closed forms, so we integrate.
    """
    def rhs(xi, y):
        th, dth = y                       # y = [theta, dtheta/dxi]
        th_safe = max(th, 0.0)            # theta can dip <0 past the surface; clamp for theta^n
        # d(theta)/dxi = dth ;  d(dth)/dxi = -theta^n - 2/xi * dth
        return [dth, -(th_safe**n) - 2.0*dth/xi]

    # start just off xi=0 with the series expansion theta ~ 1 - xi^2/6 to avoid the 1/xi singularity
    xi0 = 1e-6
    y0  = [1.0 - xi0**2/6.0, -xi0/3.0]
    xi_grid = np.linspace(xi0, xi_max, n_pts)
    # stop the integration at the surface (first zero of theta):
    surf = lambda xi, y: y[0]; surf.terminal = True; surf.direction = -1
    sol = solve_ivp(rhs, (xi0, xi_max), y0, t_eval=xi_grid, events=surf,
                    rtol=1e-10, atol=1e-12, dense_output=True)
    xi    = sol.t
    theta = np.clip(sol.y[0], 0.0, None)          # clamp tiny negative overshoot to 0
    xi1   = sol.t_events[0][0] if sol.t_events[0].size else xi[-1]
    return xi, theta, xi1


def sample_polytrope_positions(N, n, R_star, rng):
    """Draw N debris positions isotropically inside a polytrope-n star of radius
    R_star (in AU), with radii distributed according to the enclosed-mass profile
    M(<r) so that EQUAL-MASS elements are sampled (each particle represents the
    same debris mass, which is what a mass-weighted dM/deps histogram wants).

    Returns dr : array [N,3] of offsets from the stellar centre (AU).
    """
    xi, theta, xi1 = lane_emden(n)
    # cumulative mass vs dimensionless radius:  M(<xi) ~ integral xi^2 theta^n dxi
    integrand = xi**2 * theta**n
    Mcum = cumulative_trapezoid(integrand, xi, initial=0.0)
    Mcum /= Mcum[-1]                               # normalise to 1 at the surface
    # inverse-CDF sampling: uniform u in [0,1] -> xi via M(<xi)=u, then r = R_star * xi/xi1
    u   = rng.random(N)
    xis = np.interp(u, Mcum, xi)
    r   = R_star * xis / xi1                       # physical radius of each particle (AU)
    # isotropic directions on the sphere:
    cost = rng.uniform(-1, 1, N); phi = rng.uniform(0, 2*np.pi, N)
    sint = np.sqrt(1 - cost**2)
    hat  = np.stack([sint*np.cos(phi), sint*np.sin(phi), cost], axis=1)
    return r[:, None] * hat                        # [N,3] offsets dr


def make_debris_ic(state, N, n=3.0, spin=True, spin_mode="tidal", seed=0):
    """Build debris initial conditions from a disruption 'state' (see
    build_representative_state / build_parabolic_state / read_bundle_row).

    Each debris particle k gets the whole star's bulk (solid-body) velocity plus the
    star's spin, imposed as a rigid rotation about the stellar centre:
        position  r_k = r_star + dr_k                 (dr_k sampled inside the polytrope)
        velocity  v_k = v_star + omega x dr_k         (bulk velocity + spin)

    Spin models (spin_mode), spin axis always along the orbital angular momentum
    L = r_rel x v_rel  (r_rel = r_star-r_BH1, v_rel = v_star-v_BH1):

      * "tidal"  (DEFAULT, penetration-dependent): the tidal torque that spins the star
        up scales as ~1/r^3 and is therefore impulsive at pericenter, so the acquired
        spin is approximately the ORBITAL ANGULAR VELOCITY AT PERICENTER of the star's
        orbit about BH1:  |omega| = h / r_peri^2 , h=|L|, with r_peri found from the
        orbit (energy eps = 1/2 v_rel^2 - G m_bin/r, m_bin=m1+m_star).  This GROWS with
        penetration beta = r_t,*/r_peri as ~beta^(3/2) and EXCEEDS break-up for
        beta >~ 0.8 (deeper plunges spin the star past break-up), exactly as expected.
        It is the leading (impulsive) result of integrating the tidal torque; a full
        affine/hydro torque integral would refine the O(1) prefactor only.
      * "breakup":  fixed at rotational break-up, |omega| = sqrt(G m_star/R_star^3)
        (surface speed = break-up speed), independent of penetration.
      * "orbital":  the INSTANTANEOUS orbital angular velocity |omega| = h/|r_rel|^2
        (equals "tidal" when the saved state is exactly at pericenter).
    spin=False zeros the spin (pure position-offset / frozen-in energy spread only).

    Returns (r_deb [N,3], v_deb [N,3], omega [3], info dict).
    """
    rng = np.random.default_rng(seed)
    R_star = state["R_star"]; m_star = state["m_star"]; m_bin = state["m_BH1"] + m_star
    dr = sample_polytrope_positions(N, n, R_star, rng)         # [N,3] offsets inside the star

    r_rel = state["r_star"] - state["r_BH1"]                   # star relative to BH1 at disruption
    v_rel = state["v_star"] - state["v_BH1"]
    Lvec  = np.cross(r_rel, v_rel); h = np.linalg.norm(Lvec)   # specific angular momentum
    Lhat  = Lvec/h if h > 0 else np.array([0.0, 0.0, 1.0])
    rr    = np.linalg.norm(r_rel)
    # solve the star's two-body orbit about BH1 for its pericenter r_peri:
    eps_orb = 0.5*np.dot(v_rel, v_rel) - m_bin/rr              # specific orbital energy (G=1)
    if abs(eps_orb) > 1e-9*(m_bin/rr):
        ecc_orb = np.sqrt(max(0.0, 1.0 + 2.0*eps_orb*h*h/m_bin**2))
        r_peri  = (-m_bin/(2.0*eps_orb))*(1.0 - ecc_orb)       # a(1-e); bound a>0, unbound a<0
    else:
        r_peri  = h*h/(2.0*m_bin)                              # parabolic limit
    r_peri = max(r_peri, 1e-12)
    Omega_break = np.sqrt(m_star/R_star**3)                    # break-up rate (G=1)
    if not spin:
        omega = np.zeros(3)
    elif spin_mode == "breakup":
        omega = Omega_break*Lhat
    elif spin_mode == "orbital":
        omega = Lvec/(rr*rr)
    else:                                                      # "tidal" (default): pericenter Omega
        omega = (h/r_peri**2)*Lhat

    r_deb = state["r_star"] + dr                               # [N,3]
    v_deb = state["v_star"] + np.cross(omega, dr)              # [N,3]  (omega x dr per particle)

    # diagnostics
    Deps_scale = state["m_BH1"] * R_star / r_peri**2           # ~ G m1 R*/r_peri^2 (frozen-in spread)
    v_break = np.sqrt(m_star/R_star)                           # break-up surface speed [V1]
    om_mag  = float(np.linalg.norm(omega))
    r_t_star = state.get("r_t_star", R_star*(state["m_BH1"]/m_star)**(1.0/3.0))
    info = dict(omega_mag=om_mag, sep_disrupt=rr, r_peri=r_peri,
                beta_pen=float(r_t_star/r_peri),              # penetration factor r_t,*/r_peri
                spin_over_breakup=om_mag/Omega_break,         # |omega| / Omega_breakup
                vspin_surface_kms=om_mag*R_star*V1_KMS,       # equatorial spin speed [km/s]
                v_break_kms=v_break*V1_KMS, Deps_scale=Deps_scale,
                spin_mode=(spin_mode if spin else "none"))
    return r_deb, v_deb, omega, info


# =============================================================================
#  2.  DISRUPTION STATE:  representative (constructed) OR read from the bundle
# =============================================================================
def build_representative_state(M2=1e6, m1=10.0, m_star=1.0, R_star_rsun=1.0,
                               a_bin=0.1, approach="radial"):
    """Construct a physically representative mu-TDE disruption state in code units
    (G=1, Msun, AU, V1).  This stands in for a real row of the scattering bundle
    while we develop/validate the tracer; swap in read_bundle_row() for science.

    Geometry (all positions in AU, velocities in V1):
      * BH2 (central, mass M2) at the origin, at rest (we work in BH2's frame).
      * BH1 (companion, mass m1) placed at the binary's tidal radius wrt BH2,
            r_t,bin = a_bin (M2/m_bin)^(1/3),  moving on a near-circular orbit
            about BH2 at the local circular speed (a representative pericentre).
      * the star at the stellar tidal radius from BH1,
            r_t,* = R* (m1/m_star)^(1/3),
            on a parabolic (marginally-bound) orbit about BH1 -> it is being
            disrupted right now.
    """
    m_bin   = m1 + m_star
    R_star  = R_star_rsun * RSUN_AU

    # --- BH1 placed near the binary's tidal radius wrt BH2, on a circular orbit ---
    r_t_bin = a_bin * (M2 / m_bin)**(1.0/3.0)                  # AU
    r_BH1   = np.array([r_t_bin, 0.0, 0.0])
    v_circ  = np.sqrt(M2 / r_t_bin)                            # G=1 -> v = sqrt(M2/r)
    v_BH1   = np.array([0.0, v_circ, 0.0])                     # orbiting BH2 in +y

    # --- star at the stellar tidal radius from BH1, on a parabolic orbit about BH1 ---
    r_t_star = R_star * (m1 / m_star)**(1.0/3.0)               # AU  (~0.01 AU for 1 Rsun, 10 Msun)
    v_par    = np.sqrt(2.0 * m1 / r_t_star)                    # parabolic relative speed about BH1
    # place the star offset from BH1 along +x (radially outward from BH2) and moving in +y
    # relative to BH1, so it is at pericentre of a prograde orbit about BH1:
    r_star = r_BH1 + np.array([r_t_star, 0.0, 0.0])
    v_star = v_BH1 + np.array([0.0, v_par, 0.0])

    return dict(M2=M2, m_BH1=m1, m_star=m_star, R_star=R_star, a_bin=a_bin,
                r_BH2=np.zeros(3), v_BH2=np.zeros(3),
                r_BH1=r_BH1, v_BH1=v_BH1,
                r_star=r_star, v_star=v_star,
                r_t_bin=r_t_bin, r_t_star=r_t_star,
                note="representative constructed state (not a real bundle row)")


def mutde_indices(bundle_path, mu_outcome="TDE:star2_by_bh1"):
    """Return the row indices in the bundle whose outcome is a mu-TDE."""
    import h5py
    with h5py.File(bundle_path, "r") as f:
        outcome = f["results/outcome"].asstr()[:]
    return np.where(outcome == mu_outcome)[0]


def read_bundle_row(bundle_path, row_index=None, match=None, mu_outcome="TDE:star2_by_bh1"):
    """Read ONE mu-TDE event's saved disruption-time phase space from the SeaWulf
    bundle and return a 'state' dict in the same units/format as
    build_representative_state (G=1, Msun, AU, V1).

    Column schema is exactly the smbh_binary.cpp --header line:
      central BH (M_smbh):  x_smbh_AU,y_smbh_AU,z_smbh_AU, vx_smbh_kms,vy_smbh_kms,vz_smbh_kms
      body1 = the 10 Msun BH (bh1):  x_b1_AU...vz_b1_kms     (mass column 'm1')
      body2 = the star (star2):      x_b2_AU...vz_b2_kms     (mass 'm2', radius 'R2' [Rsun])
    Positions are in AU, velocities in km/s (we convert km/s -> V1 units), in the
    system COM frame.  For the 'TDE:star2_by_bh1' outcome the saved state is taken
    at the terminal TDE event, i.e. the moment of disruption — exactly our seed point.

    Select the event with EITHER row_index (index into the mu-TDE subset, 0-based)
    or match=dict(M_smbh=..., a_bin=..., e_bin=...) to pick the nearest cell in log.
    """
    import h5py
    POS_SM = ("x_smbh_AU","y_smbh_AU","z_smbh_AU"); VEL_SM = ("vx_smbh_kms","vy_smbh_kms","vz_smbh_kms")
    POS_B1 = ("x_b1_AU","y_b1_AU","z_b1_AU");       VEL_B1 = ("vx_b1_kms","vy_b1_kms","vz_b1_kms")
    POS_B2 = ("x_b2_AU","y_b2_AU","z_b2_AU");       VEL_B2 = ("vx_b2_kms","vy_b2_kms","vz_b2_kms")
    with h5py.File(bundle_path, "r") as f:
        res = f["results"]
        need = POS_SM+VEL_SM+POS_B1+VEL_B1+POS_B2+VEL_B2+("M_smbh","m1","m2","R2")
        missing = [c for c in need if c not in res]
        if missing:
            raise KeyError(f"bundle is missing columns {missing}; available: {sorted(res.keys())}")
        outcome = res["outcome"].asstr()[:]
        mu = np.where(outcome == mu_outcome)[0]
        if mu.size == 0:
            raise ValueError(f"no {mu_outcome} rows in this bundle")
        if match is not None:
            M = res["M_smbh"][:][mu]; a = res["a_bin"][:][mu]; e = res["e_bin"][:][mu]
            score = (np.log10(M)-np.log10(match["M_smbh"]))**2 \
                  + (np.log10(a)-np.log10(match["a_bin"]))**2 \
                  + (e-match.get("e_bin", e))**2
            sel = int(mu[np.argmin(score)])
        else:
            sel = int(mu[row_index if row_index is not None else 0])
        g = lambda keys: np.array([res[k][sel] for k in keys], float)
        r_bh2, v_bh2 = g(POS_SM), g(VEL_SM)/V1_KMS     # central BH  (km/s -> V1)
        r_bh1, v_bh1 = g(POS_B1), g(VEL_B1)/V1_KMS     # 10 Msun BH (companion / disruptor)
        r_star, v_star = g(POS_B2), g(VEL_B2)/V1_KMS   # the star (debris progenitor)
        M2 = float(res["M_smbh"][sel]); m1 = float(res["m1"][sel])
        m_star = float(res["m2"][sel]); R2_rsun = float(res["R2"][sel])
    return dict(M2=M2, m_BH1=m1, m_star=m_star, R_star=R2_rsun*RSUN_AU,
                r_BH2=r_bh2, v_BH2=v_bh2, r_BH1=r_bh1, v_BH1=v_bh1,
                r_star=r_star, v_star=v_star,
                r_t_star=np.linalg.norm(r_star-r_bh1),
                note=f"bundle row {sel} ({mu_outcome})")


# =============================================================================
#  3.  VECTORIZED NEWTONIAN INTEGRATOR  (massless debris in the BH1+BH2 field)
# =============================================================================
def _accelerations(P, m1, M2, soft2):
    """Newtonian accelerations (G=1) for the packed position array P of shape
    [(2+N),3]: index 0 = BH1 (mass m1), index 1 = BH2 (mass M2), 2.. = massless
    debris.  Debris feel BH1 and BH2 but NOT each other; the BHs feel only each
    other (the massless debris exert no force).  'soft2' is a squared softening
    length that keeps the acceleration finite for debris that pass very close to a
    BH (physically, such material would be swallowed; the softening is a numerical
    regulator, set ~ R_star)."""
    A = np.zeros_like(P)
    # --- the two BHs: a clean 2-body problem (debris do not enter) ---
    d12 = P[1] - P[0]; r12 = (d12 @ d12)**1.5
    A[0] =  M2 * d12 / r12            # BH1 pulled toward BH2
    A[1] = -m1 * d12 / r12            # BH2 pulled toward BH1
    # --- debris (vectorized over all of them) ---
    rp = P[2:]                                   # [N,3]
    d1 = P[0] - rp; d2 = P[1] - rp               # vectors from each debris to BH1, BH2
    inv1 = (np.einsum("ij,ij->i", d1, d1) + soft2)**-1.5
    inv2 = (np.einsum("ij,ij->i", d2, d2) + soft2)**-1.5
    A[2:] = m1 * d1 * inv1[:, None] + M2 * d2 * inv2[:, None]
    return A


def integrate_debris(state, r_deb, v_deb, t_end, n_frames=300, soft=None,
                     rtol=1e-9, atol=1e-12, method="DOP853"):
    """Advance the 2 BHs + N massless debris from the disruption state to t_end
    (code-time units), returning dense fixed-cadence output:

        t        : [n_frames]                 times
        traj_bh  : [2, n_frames, 6]           BH1,BH2 phase space  (x,y,z,vx,vy,vz)
        traj_deb : [N, n_frames, 6]           every debris particle's phase space
        dE_bh    : [n_frames]                 fractional 2-body energy drift (accuracy check)

    The integration is a single ODE solve over the packed state y = [positions,
    velocities]; dy/dt = [velocities, accelerations].  Because the debris are
    massless and mutually non-interacting, one solve advances them all.
    """
    N = r_deb.shape[0]
    m1, M2 = state["m_BH1"], state["M2"]
    if soft is None:
        soft = state["R_star"]                   # default softening ~ one stellar radius
    soft2 = soft*soft

    # pack [BH1, BH2, debris] positions then the same ordering for velocities
    P0 = np.vstack([state["r_BH1"], state["r_BH2"], r_deb])     # [(2+N),3]
    Vv0 = np.vstack([state["v_BH1"], state["v_BH2"], v_deb])
    nbod = N + 2
    y0 = np.concatenate([P0.ravel(), Vv0.ravel()])

    def rhs(t, y):
        P = y[:3*nbod].reshape(nbod, 3)
        Vv = y[3*nbod:].reshape(nbod, 3)
        A = _accelerations(P, m1, M2, soft2)
        return np.concatenate([Vv.ravel(), A.ravel()])

    t_eval = np.linspace(0.0, t_end, n_frames)
    sol = solve_ivp(rhs, (0.0, t_end), y0, t_eval=t_eval,
                    method=method, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError("integration failed: " + sol.message)

    Y = sol.y.T                                  # [n_frames, 6*nbod]
    pos = Y[:, :3*nbod].reshape(-1, nbod, 3)     # [n_frames, nbod, 3]
    vel = Y[:, 3*nbod:].reshape(-1, nbod, 3)
    phase = np.concatenate([pos, vel], axis=2)   # [n_frames, nbod, 6]
    phase = np.transpose(phase, (1, 0, 2))       # [nbod, n_frames, 6]
    traj_bh, traj_deb = phase[:2], phase[2:]

    # accuracy monitor: the 2-BH orbit must conserve its own energy (debris are massless)
    rb = traj_bh[:, :, :3]; vb = traj_bh[:, :, 3:]
    dsep = np.linalg.norm(rb[1]-rb[0], axis=1)
    Ekin = 0.5*m1*np.sum(vb[0]**2,1) + 0.5*M2*np.sum(vb[1]**2,1)
    Epot = -m1*M2/dsep
    Etot = Ekin + Epot
    dE_bh = np.abs(Etot - Etot[0]) / np.abs(Etot[0])
    return dict(t=sol.t, traj_bh=traj_bh, traj_deb=traj_deb, dE_bh=dE_bh,
                soft=soft, N=N, state=state)


# =============================================================================
#  4.  ENERGIES & CLASSIFICATION  (done at post-processing on the saved trajectory)
# =============================================================================
def classify(result, frame=-1):
    """At output 'frame' (default last), label each debris particle by its fate,
    using a HIERARCHICAL (not "most bound") criterion.  This matters because BH2's
    potential well is enormously deeper than BH1's: near BH1 essentially every
    particle has eps_BH2<0 (BH1 itself orbits BH2), so a naive "most bound" test
    would call everything BH2-bound.  The physically meaningful question is nested:

        1. is it still bound to its LOCAL disruptor BH1?  eps_BH1 < 0
              -> label 'BH1'  (stays as a mini-disk / accretes onto the 10 Msun BH)
        2. else, has it escaped BH1 but stays bound to the central BH2?  eps_BH2 < 0
              -> label 'BH2'  (spread into the nuclear region around the central BH)
        3. else positive energy wrt both
              -> label 'unbound'  (ejected from the system)

    Specific energy wrt body b:  eps_b = 1/2 |v - v_b|^2 - G m_b / |r - r_b|  (G=1).
    These are INSTANTANEOUS two-body energies at the chosen frame; evaluate at late
    time (after the debris have separated from BH1 and the BH1-BH2 tidal exchange
    has played out) for the asymptotic fate.  Returns eps_BH1, eps_BH2, integer
    label per particle, and the population fractions.
    """
    st = result["state"]; m1, M2 = st["m_BH1"], st["M2"]
    bh = result["traj_bh"][:, frame, :]          # [2,6]  BH1,BH2 phase space
    deb = result["traj_deb"][:, frame, :]        # [N,6]
    r, v = deb[:, :3], deb[:, 3:]
    def eps(mb, rb, vb):
        dr = r - rb; dv = v - vb
        return 0.5*np.einsum("ij,ij->i", dv, dv) - mb/np.linalg.norm(dr, axis=1)
    eps1 = eps(m1, bh[0, :3], bh[0, 3:])         # wrt companion BH1
    eps2 = eps(M2, bh[1, :3], bh[1, 3:])         # wrt central BH2
    label = np.full(len(r), 2, int)              # start all 'unbound' (2)
    label[eps2 < 0] = 1                          # bound to BH2 (overwritten by BH1 next)
    label[eps1 < 0] = 0                          # bound to BH1 takes precedence (local disruptor)
    frac = dict(BH1=np.mean(label == 0), BH2=np.mean(label == 1),
                unbound=np.mean(label == 2))
    return dict(eps_BH1=eps1, eps_BH2=eps2, label=label, frac=frac)


# =============================================================================
#  5.  SAVE FULL TRAJECTORIES TO HDF5  (no classification baked in)
# =============================================================================
def save_h5(result, path):
    """Write the complete run to one HDF5 file: BH and debris trajectories, the
    time array, and the disruption state as attributes.  Classification is left to
    the reader (classify()), per the 'save everything' requirement."""
    import h5py
    with h5py.File(path, "w") as f:
        f.create_dataset("t", data=result["t"])                              # [n_frames]
        f.create_dataset("traj_bh",  data=result["traj_bh"],  compression="gzip")   # [2,n_frames,6]
        f.create_dataset("traj_deb", data=result["traj_deb"], compression="gzip")   # [N,n_frames,6]
        f.create_dataset("dE_bh", data=result["dE_bh"])
        g = f.create_group("state")
        for k, val in result["state"].items():
            g.attrs[k] = val if np.isscalar(val) or isinstance(val, str) else np.asarray(val)
        f.attrs["units"] = "G=1, mass=Msun, length=AU, vel=V1(29.7847km/s), time=yr/2pi"
        f.attrs["columns"] = "x,y,z,vx,vy,vz"
        f.attrs["soft_AU"] = result["soft"]; f.attrs["N_debris"] = result["N"]


# =============================================================================
#  6.  EXPORT INITIAL CONDITIONS FOR THE SpaceHub C++ DRIVER (consistency check)
# =============================================================================
def write_spacehub_ic(state, r_deb, v_deb, t_end, n_frames, path,
                      soft=None, m_tiny=1e-12):
    """Write a plain-text IC file for debris_spacehub.cpp.  Format (all in the same
    units the restart mode of smbh_binary.cpp uses: Msun, Rsun, AU, km/s):

        line 1 : N_bodies  t_end_yr  n_frames  soft_AU
        then one line per body:  mass_Msun  radius_Rsun  x_AU y_AU z_AU  vx_kms vy_kms vz_kms
        body order: BH1 (10 Msun), BH2 (central), then the N debris.

    Debris are given a TINY mass m_tiny (not exactly 0): SpaceHub's default
    algorithmic-regularization method weights its time transformation by masses, so
    a strictly-massless particle can upset the chain; 1e-12 Msun is dynamically
    negligible (<<1e-10 of BH1) yet keeps the regularization well-defined.  Velocities
    are converted internal V1 -> km/s on write so the C++ reads km/s like --restart.
    """
    if soft is None: soft = state["R_star"]
    N = r_deb.shape[0]
    def kms(v): return np.asarray(v) * V1_KMS
    with open(path, "w") as f:
        f.write(f"{N+2} {t_end*YR_PER_T:.10e} {n_frames} {soft:.10e}\n")
        # BH1 (companion, 10 Msun), radius ~ 2 r_g (point mass; use 0 radius here)
        f.write(f"{state['m_BH1']:.10e} 0 "
                f"{state['r_BH1'][0]:.10e} {state['r_BH1'][1]:.10e} {state['r_BH1'][2]:.10e} "
                f"{kms(state['v_BH1'])[0]:.10e} {kms(state['v_BH1'])[1]:.10e} {kms(state['v_BH1'])[2]:.10e}\n")
        # BH2 (central)
        f.write(f"{state['M2']:.10e} 0 "
                f"{state['r_BH2'][0]:.10e} {state['r_BH2'][1]:.10e} {state['r_BH2'][2]:.10e} "
                f"{kms(state['v_BH2'])[0]:.10e} {kms(state['v_BH2'])[1]:.10e} {kms(state['v_BH2'])[2]:.10e}\n")
        # debris (tiny mass, zero radius)
        vk = kms(v_deb)
        for i in range(N):
            f.write(f"{m_tiny:.3e} 0 "
                    f"{r_deb[i,0]:.10e} {r_deb[i,1]:.10e} {r_deb[i,2]:.10e} "
                    f"{vk[i,0]:.10e} {vk[i,1]:.10e} {vk[i,2]:.10e}\n")
    return path


# =============================================================================
#  7.  MOVIE:  two panels — (left) system view with BOTH BHs + debris moving,
#              (right) co-moving zoom on BH1 showing the mini-disk.
#      Debris are coloured by their FINAL fate (bound BH1 / bound BH2 / unbound),
#      so you can watch which material ends up where.
# =============================================================================
def make_movie(result, path, labels=None, zoom_rt=60.0, fps=25, stride_lim=30):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.lines import Line2D
    try:
        import style as _sty; _sty.use_paper_style(base=12)
    except Exception:
        pass
    tb, td, t, st = result["traj_bh"], result["traj_deb"], result["t"], result["state"]
    nfr, N = len(t), td.shape[0]
    # Colour each debris by its FROZEN-IN specific energy wrt BH1 (the fallback
    # energy set at disruption) -- a continuous quantity, fixed for the whole movie,
    # shown with a colorbar.  eps<0 = bound to BH1 (will fall back / mini-disk),
    # eps>0 = unbound from BH1.  (The categorical bound-BH1/BH2/unbound split is in
    # the static debris_energy figure instead.)
    cval = classify(result, frame=0)["eps_BH1"]                # [N]  frozen-in eps wrt BH1 [V1^2]
    clim = float(np.nanpercentile(np.abs(cval), 98)) + 1e-30   # symmetric colour range about 0
    CMAP = "coolwarm"

    # --- WIDE view limits: enclose both BH orbits + the debris over the whole run ---
    bh_xy  = tb[:, :, :2].reshape(-1, 2)
    sub    = td[:, ::max(1, nfr//stride_lim), :2].reshape(-1, 2)   # subsample debris for limits
    allxy  = np.concatenate([bh_xy, sub], axis=0)
    xmn, ymn = np.nanpercentile(allxy, 0.5, axis=0)
    xmx, ymx = np.nanpercentile(allxy, 99.5, axis=0)
    pad = 0.08*max(xmx-xmn, ymx-ymn) + 1e-9
    # --- ZOOM half-width: a few tens of stellar tidal radii about BH1 ---
    rt  = st.get("r_t_star", float(np.linalg.norm(td[:,0,:3]-tb[0,0,:3], axis=1).mean()))
    span = zoom_rt*rt

    fig, (axW, axZ) = plt.subplots(1, 2, figsize=(13.4, 6.0))
    # left: system frame (BHs move)
    scW = axW.scatter(td[:,0,0], td[:,0,1], s=3, c=cval, cmap=CMAP, vmin=-clim, vmax=clim)
    # black dots for the BHs: BH1 (10 Msun) small, BH2 (central) larger
    b1W, = axW.plot([tb[0,0,0]], [tb[0,0,1]], "o", ms=7,  color="k")
    b2W, = axW.plot([tb[1,0,0]], [tb[1,0,1]], "o", ms=14, color="k")
    axW.set_xlim(xmn-pad, xmx+pad); axW.set_ylim(ymn-pad, ymx+pad); axW.set_aspect("equal")
    axW.set_xlabel("x [AU]"); axW.set_ylabel("y [AU]"); axW.set_title("system view — both BHs + debris")
    axW.legend(handles=[Line2D([0],[0],marker='o',ls='',color='k',ms=5,label='BH1 (10 M$_\\odot$)'),
                        Line2D([0],[0],marker='o',ls='',color='k',ms=9,label='BH2 (central)')],
               fontsize=8, loc='upper right', framealpha=0.9)
    # right: co-moving zoom on BH1
    rel0 = td[:,0,:2] - tb[0,0,:2]
    scZ = axZ.scatter(rel0[:,0], rel0[:,1], s=4, c=cval, cmap=CMAP, vmin=-clim, vmax=clim)
    axZ.plot([0],[0], "o", ms=9, color="k")
    axZ.set_xlim(-span, span); axZ.set_ylim(-span, span); axZ.set_aspect("equal")
    axZ.set_xlabel(r"x $-$ x$_{\rm BH1}$ [AU]"); axZ.set_ylabel(r"y $-$ y$_{\rm BH1}$ [AU]")
    axZ.set_title("zoom on BH1 — mini-disk")
    # shared colorbar: the frozen-in fallback energy wrt BH1
    cb = fig.colorbar(scW, ax=[axW, axZ], fraction=0.046, pad=0.02)
    cb.set_label(r"frozen-in $\varepsilon_{\rm BH1}$  [$V_1^2$]   (blue: bound to BH1  $\cdot$  red: unbound)")
    cb.ax.tick_params(which="both", direction="in")
    # movies-tar on-canvas layout: title TOP-centre, time LOWER-centre, scale LOWER-right.
    fig.suptitle("micro-TDE debris", y=0.98)
    tlabel = fig.text(0.5, 0.045, "", ha="center", va="bottom", fontsize=11)
    # scale bar, lower-right of the zoom panel (a round number of AU)
    nice = 10.0**np.floor(np.log10(span/2)); Lb = nice*max(1, round((span/2)/nice))
    axZ.plot([span-0.08*span-Lb, span-0.08*span], [-span+0.10*span, -span+0.10*span],
             "k-", lw=2.4, solid_capstyle="butt")
    axZ.text(span-0.08*span-Lb/2, -span+0.13*span, f"{Lb:g} AU",
             ha="center", va="bottom", fontsize=9)

    def upd(i):
        scW.set_offsets(td[:,i,:2]); b1W.set_data([tb[0,i,0]],[tb[0,i,1]]); b2W.set_data([tb[1,i,0]],[tb[1,i,1]])
        scZ.set_offsets(td[:,i,:2] - tb[0,i,:2])
        tlabel.set_text(f"t = {t[i]*YR_PER_T*SEC_PER_YR/86400:.2f} days   (frame {i+1}/{nfr})")
        return scW, b1W, b2W, scZ, tlabel
    an = animation.FuncAnimation(fig, upd, frames=nfr, interval=1000/fps, blit=True)
    try:
        an.save(path, writer=animation.FFMpegWriter(fps=fps, bitrate=3000))
    except Exception:
        path = path.replace(".mp4", ".gif"); an.save(path, writer=animation.PillowWriter(fps=min(fps,20)))
    plt.close(fig)
    return path


# =============================================================================
#  8.  PUBLICATION FIGURE:  BH1 orbit + a representative subsample of debris
#      trajectories (in BH2's rest frame, so the central BH is fixed at the
#      origin), coloured by frozen-in fallback energy, with an inset showing the
#      full (all-N) debris energy distribution.  Saves both PNG and PDF.
# =============================================================================
def paper_figure(result, path, n_curves=48, seed=0, frame0=0, frame1=None, figsize=(7.4, 7.2)):
    # Publication single panel, styled to match debris Fig 1 and the smbh_movies
    # layout: energy inset at LOWER-LEFT, panel title INSIDE top-left in a white
    # bbox, a scale bar LOWER-RIGHT, CM-serif fonts and ticks-in on all four sides.
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    try:
        import style as _sty; _sty.use_paper_style(base=12)
        _inside = _sty.inside_title
    except Exception:
        _inside = lambda ax, s, **k: ax.set_title(s)
    tb, td, t, st = result["traj_bh"], result["traj_deb"], result["t"], result["state"]
    f1 = frame1 if frame1 is not None else len(t)

    # work in BH2's frame: subtract the central BH's position at every time so BH2
    # sits still at the origin (BH1 and the debris are drawn relative to it).
    rbh2 = tb[1, frame0:f1, :2]                       # [nf,2]
    bh1  = tb[0, frame0:f1, :2] - rbh2                # BH1 trajectory rel. to BH2

    # colour / select by the FROZEN-IN specific energy wrt BH1 (fallback energy):
    cval = classify(result, frame=0)["eps_BH1"]       # [N]
    clim = float(np.nanpercentile(np.abs(cval), 98)) + 1e-30
    norm = mpl.colors.Normalize(-clim, clim); cmap = mpl.cm.coolwarm
    # a REPRESENTATIVE subsample: evenly spaced across the sorted energy, so the
    # drawn curves span the full bound->unbound range rather than a random clump.
    order = np.argsort(cval)
    idx = order[np.linspace(0, len(order)-1, min(n_curves, len(order))).astype(int)]

    fig, ax = plt.subplots(figsize=figsize)
    for k in idx:                                     # thin coloured debris trajectories
        xy = td[k, frame0:f1, :2] - rbh2
        ax.plot(xy[:, 0], xy[:, 1], lw=0.7, alpha=0.85, color=cmap(norm(cval[k])))
    ax.plot(bh1[:, 0], bh1[:, 1], "k-", lw=1.8, label=r"BH1 orbit (10 M$_\odot$)")
    ax.plot([0], [0], "o", ms=14, color="k", label="BH2 (central)")
    ax.plot([bh1[0, 0]], [bh1[0, 1]], "o", ms=6, color="k")        # BH1 at disruption
    ax.set_xlabel("x [AU]"); ax.set_ylabel("y [AU]"); ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    _inside(ax, f"{len(idx)} representative debris trajectories", loc="upper left", fontsize=10)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"frozen-in $\varepsilon_{\rm BH1}$ [$V_1^2$]  (blue: bound $\cdot$ red: unbound)")
    cb.ax.tick_params(which="both", direction="in")

    # scale bar, LOWER-RIGHT (movies convention): a round number of AU
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    span = x1 - x0
    nice = 10.0**np.floor(np.log10(span/4)); L = nice*max(1, round((span/4)/nice))
    xb = x1 - 0.06*span - L; yb = y0 + 0.08*(y1 - y0)
    ax.plot([xb, xb+L], [yb, yb], "k-", lw=2.2, solid_capstyle="butt")
    ax.text(xb+L/2, yb+0.02*(y1-y0), f"{L:g} AU", ha="center", va="bottom", fontsize=9)

    # energy inset, LOWER-LEFT (movies convention)
    axin = ax.inset_axes([0.065, 0.085, 0.34, 0.26])
    axin.hist(cval, bins=60, range=(-clim, clim), color="0.45")
    axin.axvline(0, color="k", ls="--", lw=0.8)
    axin.set_xlabel(r"$\varepsilon_{\rm BH1}$ [$V_1^2$]", fontsize=7)
    axin.set_ylabel(r"$dN/d\varepsilon$", fontsize=7)
    axin.tick_params(labelsize=6, which="both", direction="in")
    axin.text(0.5, 1.02, f"all {len(cval)} debris", transform=axin.transAxes,
              ha="center", va="bottom", fontsize=7)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.rsplit(".", 1)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return path


# =============================================================================
#  9.  REPRESENTATIVE STATE WITH A PARABOLIC COM ORBIT ABOUT BH2
#      The BH1+star pair approaches BH2 on a MARGINALLY-bound (parabolic, E_COM=0)
#      orbit and the mu-TDE happens at PERICENTER.  This is the physically honest
#      stand-in for a real encounter (unlike the bound-circular placeholder): the
#      debris energy kick now straddles the BH2 escape boundary, so material comes
#      out bound-to-BH1, bound-to-BH2, AND unbound (ejected) -- all three fates.
# =============================================================================
def build_parabolic_state(M2=1e6, m1=10.0, m_star=1.0, R_star_rsun=1.0,
                          a_bin=0.1, beta_pen=1.0, r_peri=None):
    """beta_pen = penetration factor r_t,*/r_peri,inner (>=1): how deep inside the
    stellar tidal radius the inner-binary pericenter plunges.  beta_pen=1 is a grazing
    disruption (pericenter = r_t,*); larger beta_pen means a deeper plunge, a stronger
    disruption, and a faster (super-break-up) stellar spin."""
    m_bin  = m1 + m_star
    R_star = R_star_rsun * RSUN_AU
    # pericenter of the parabola = binary's tidal radius wrt BH2 (strong encounter):
    if r_peri is None:
        r_peri = a_bin * (M2 / m_bin)**(1.0/3.0)
    # --- COM on the parabola at pericenter: position along +x, velocity tangential (+y),
    #     speed = local escape speed sqrt(2 G M2 / r_peri)  (E_COM = 0 exactly) ---
    r_com = np.array([r_peri, 0.0, 0.0])
    v_com = np.array([0.0, np.sqrt(2.0*M2/r_peri), 0.0])
    # --- the BOUND BH1-star BINARY, caught at the pericenter of its (eccentric) inner orbit ---
    # Binary semimajor axis a_bin; the star is disrupted at an inner pericenter
    #   r_peri_in = r_t,* / beta_pen   (beta_pen>=1: deeper plunge inside the tidal radius),
    # a bound, eccentric orbit (e_in = 1 - r_peri_in/a_bin).  Kepler pericenter speed
    #   v_peri = sqrt(G m_bin (2/r_peri_in - 1/a_bin)),  G=1,  perpendicular to the separation.
    # Separation oriented tangential (|| v_com), relative velocity radial (_|_ v_com) so the
    # debris energy spread wrt BH2 straddles the escape boundary -> all three fates (a
    # representative orientation; the real bundle rows carry the true one).
    r_t_star  = R_star * (m1 / m_star)**(1.0/3.0)          # stellar tidal radius
    r_peri_in = r_t_star / beta_pen                        # actual inner pericenter (deeper if beta>1)
    e_in      = 1.0 - r_peri_in/a_bin                      # inner-binary eccentricity at disruption
    v_peri    = np.sqrt(m_bin*(2.0/r_peri_in - 1.0/a_bin)) # bound-orbit pericenter speed (G=1)
    sep   = np.array([0.0, r_peri_in, 0.0])                # separation (= inner pericenter), tangential
    v_rel = np.array([v_peri, 0.0, 0.0])                   # relative velocity at pericenter, radial
    # split the binary COM into BH1 and the star (reduced-mass partition):
    r_BH1  = r_com - (m_star/m_bin)*sep;  v_BH1  = v_com - (m_star/m_bin)*v_rel
    r_star = r_com + (m1/m_bin)*sep;      v_star = v_com + (m1/m_bin)*v_rel
    return dict(M2=M2, m_BH1=m1, m_star=m_star, R_star=R_star, a_bin=a_bin, e_in=e_in,
                r_BH2=np.zeros(3), v_BH2=np.zeros(3),
                r_BH1=r_BH1, v_BH1=v_BH1, r_star=r_star, v_star=v_star,
                r_t_bin=r_peri, r_t_star=r_t_star, r_peri_in=r_peri_in, beta_pen=beta_pen,
                note="bound BH1-star binary (COM on a parabola about BH2), disrupted at inner pericenter")
