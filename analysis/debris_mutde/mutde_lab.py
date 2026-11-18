#!/usr/bin/env python3
# =============================================================================
#  mutde_lab.py
#  -----------------------------------------------------------------------------
#  Single-case micro-TDE ("mu-TDE") + debris-fallback laboratory, for a laptop.
#
#  WHAT IT DOES  (one representative encounter, start to finish)
#  ------------------------------------------------------------
#   STAGE 1 (SpaceHub) : the ACCURATE incoming 3-body orbit -- SMBH + stellar BH
#                        + star -- is integrated by your compiled SpaceHub
#                        executable `smbh_binary_PN` until the star is disrupted
#                        by its BH companion (the mu-TDE).  We need SpaceHub here
#                        for the fidelity of the binary evolution and the deep
#                        approach to the SMBH.
#
#   STAGE 2 (ballistic): at the disruption instant we replace the star by N
#                        massless test particles (debris) sampled from an n=3
#                        polytrope, following Evans & Kochanek (1989): each
#                        element keeps the position/velocity it had at disruption,
#                        i.e. bulk stellar velocity + solid-body spin  omega x dr,
#                        with omega the star's ORBITAL angular velocity about the
#                        BH at disruption.  The debris then move purely
#                        ballistically in the time-dependent field of the two
#                        moving BHs (no self-gravity, no back-reaction, no
#                        debris-debris force).  The two BHs move on their own
#                        mutual two-body orbit.
#
#   MULTI-TIMESCALE    : the problem is stiff -- the BH mini-disk orbital time
#                        (~hours) and the SMBH return time (~months) differ by
#                        ~1000x.  We therefore integrate in TWO phases:
#                          phase A: small dt, resolve disruption + mini-disk
#                                   formation (shown in slow motion);
#                          FREEZE : once the BH-bound debris have phase-mixed
#                                   into the mini-disk, we FREEZE them (they then
#                                   ride rigidly with the BH) and
#                          phase B: large dt, resolve the SMBH-bound debris
#                                   returning to the SMBH (shown sped up).
#
#   OUTPUTS            : (1) an animation (gif + mp4) with the slow->fast time
#                        warp, debris coloured by specific energy wrt the SMBH
#                        (erg/g, symmetric about 0: blue = bound, red = unbound),
#                        and a zoom inset on the stellar-BH mini-disk;
#                        (2) the NUMERICAL mass-fallback curve dM/dt (measured
#                        from first-return times, NOT the frozen-in estimate),
#                        with the frozen-in curve overplotted for comparison;
#                        (3) the "money plot" (Fig. 1) -- the final snapshot with
#                        the energy distribution and dM/dt insets.
#
#  USAGE
#  -----
#    python3 mutde_lab.py --run                 # run SpaceHub, then everything
#    python3 mutde_lab.py --traj bss_..._full.txt   # reuse an existing bss dump
#    python3 mutde_lab.py --demo                 # synthesize stage-1 (no SpaceHub)
#    options:  --ndeb 100000   --out outdir   --exe path/to/smbh_binary_PN   --no-movie   --no-money
#
#  Everything you would normally change lives in the CONFIG dataclass below.
# =============================================================================

import os, sys, glob, time, shutil, argparse, subprocess, tempfile

# repo root = two levels up from analysis/debris_mutde/
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dataclasses import dataclass, field
import numpy as np

# -----------------------------------------------------------------------------
#  UNITS.  Internally we use  AU, year, Msun  (so G = 4*pi^2).  These match the
#  SpaceHub dump (positions in AU; its internal speed unit = 2*pi AU/yr).
# -----------------------------------------------------------------------------
G      = 4.0 * np.pi**2                 # gravitational constant [AU^3 yr^-2 Msun^-1]
C_AUYR = 63241.077                      # speed of light          [AU/yr]
RSUN_AU = 695700.0 / 1.495978707e8      # 1 Rsun in AU            = 0.0046491
TWO_PI  = 2.0 * np.pi

# cgs conversions (for the energy colourbar in erg/g = cm^2/s^2)
_AU_CM   = 1.495978707e13               # 1 AU in cm
_YR_S    = 365.25636042 * 86400.0       # 1 yr in s
_AUYR_CMS = _AU_CM / _YR_S              # 1 AU/yr in cm/s  (~4.74e5)
SPEC_E_AUYR2_TO_ERGG = _AUYR_CMS**2     # (AU/yr)^2 -> erg/g  (~2.25e11)
KMS_PER_AUYR = 4.740572                  # 1 AU/yr in km/s   (for v_inf -> AU/yr)
DAY_PER_YR   = 365.25636042             # days per year     (time axes shown in days)


# =============================================================================
#  CONFIG -- edit here, or override the initial conditions on the command line.
#  Defaults = the "mutde_slow" representative hard-binary mu-TDE.
# =============================================================================
@dataclass
class Config:
    # ---- SpaceHub executable (stage 1) --------------------------------------
    #      default: $SMBH_BIN/smbh_binary_PN, else <repo>/simulation/bin/smbh_binary_PN
    bss_exe: str = os.path.join(os.environ.get("SMBH_BIN", os.path.join(REPO_ROOT, "simulation", "bin")),
                                "smbh_binary_PN")

    # ---- initial conditions (the encounter) ---------------------------------
    M_smbh: float = 1.0e6      # SMBH mass                       [Msun]
    m_bh:   float = 10.0       # stellar-mass BH (the disruptor) [Msun]
    m_star: float = 1.0        # star mass                       [Msun]
    R_star: float = 1.0        # star radius                     [Rsun]
    a_bin:  float = 0.10       # inner-binary semimajor axis     [AU]
    e_bin:  float = 0.85       # inner-binary eccentricity
    beta:   float = 1.5        # penetration  r_t,bin / r_p
    v_inf:  float = 1.0        # asymptotic approach speed       [km/s]  (SpaceHub --v_inf)
    n_dyn:  float = 2.0e4      # SpaceHub outer-orbit length     [t_dyn units]
    # locked orientation (radians) -- "mutde_mid": a STABLE binary (pericentre
    # a_bin*(1-e)=0.015 AU > r_t=0.010 AU) that is disrupted by the SMBH tide near
    # its pericentre r_p = r_t,bin/beta ~ 3 AU, NOT by the companion in isolation.
    # (The condition a_bin*(1-e_bin) > R_star*(m_bh/m_star)^(1/3) MUST hold, else the
    #  binary self-disrupts far out and no debris ever returns to the SMBH.)
    incl:  float = 3.06878
    Omega: float = 0.445981
    omega: float = 0.347833
    nu:    float = 3.60568
    rtol:  float = 1e-13

    # ---- debris sampling ----------------------------------------------------
    n_deb:   int  = 10_000     # number of debris test particles (1e4 default; 1e5 = sharper streams)
    poly_n:  float = 3.0       # polytropic index of the star (n=3 ~ solar)
    seed:    int  = 7

    # ---- two-phase ballistic integration ------------------------------------
    # phase A: resolve disruption + mini-disk;  phase B: resolve SMBH return.
    steps_per_Pdisr: int   = 240      # phase-A timestep = P_disr / this  (mini-disk resolution)
    relax_orbits:    float = 8.0      # integrate the mini-disk this many BH-orbits before freezing
    disk_capture_rt: float = 5.0      # a parcel is "mini-disk" if bound to BH and within this * r_t
    return_margin:   float = 1.6      # integrate phase B to this * (first SMBH return time)
    steps_phaseB:    int   = 1200     # number of phase-B steps across the SMBH return

    # ---- movie --------------------------------------------------------------
    # stage 1 (read from SpaceHub): the incoming SMBH+BH+star motion up to disruption.
    frames_S: int = 110       # animation frames of the incoming binary (approach from the right)
    zoom_bin_fac: float = 1.8  # stage-1 zoom-inset half-width = this * a_bin (resolves the binary)
    # frames_A: HALF are concentrated in the first ~2 disruption orbits (true slow-mo
    # of the tear-apart + first fallback); the rest cover the mini-disk relaxation.
    frames_A: int = 320       # animation frames during disruption/mini-disk (slow motion)
    slowmo_orbits: float = 2.0  # fraction of frames_A packed into the first this-many P_disr
    frames_B: int = 170       # animation frames during the SMBH return (sped up)
    fps:      int = 30

    # ---- output -------------------------------------------------------------
    out: str = os.path.join(os.environ.get("SMBH_DATA", os.path.join(REPO_ROOT, "data")), "debris_mutde", "out")
    use_tex: bool = True       # render text/labels with LaTeX (falls back to CM mathtext
                               # if a latex+dvipng install is not found); set False to force CM


# =============================================================================
#  STAGE 1 -- SpaceHub incoming orbit (or a synthetic stand-in for --demo)
# =============================================================================
def run_spacehub(cfg: Config, workdir: str) -> str:
    """Run `smbh_binary_PN --mode full` for the configured encounter and return
    the path to the `bss_*_full.txt` trajectory it writes.  SpaceHub drops the
    file in its working directory, so we run inside `workdir`."""
    if not (os.path.isfile(cfg.bss_exe) and os.access(cfg.bss_exe, os.X_OK)):
        sys.exit(f"[stage1] SpaceHub executable not found/executable:\n  {cfg.bss_exe}\n"
                 f"        build it with `make -C simulation`, pass --exe, or use --traj <file> / --demo.")
    cmd = [cfg.bss_exe,
           "--M_smbh", repr(cfg.M_smbh), "--m1", repr(cfg.m_bh), "--m2", repr(cfg.m_star),
           "--type1", "bh", "--type2", "star", "--R1", "-1", "--R2", repr(cfg.R_star),
           "--a_bin", repr(cfg.a_bin), "--e_bin", repr(cfg.e_bin),
           "--outer", "hyperbolic", "--beta", repr(cfg.beta),
           "--v_inf", repr(cfg.v_inf), "--n_dyn", repr(cfg.n_dyn),
           "--incl", repr(cfg.incl), "--Omega", repr(cfg.Omega),
           "--omega", repr(cfg.omega), "--nu", repr(cfg.nu),
           "--rtol", repr(cfg.rtol), "--mode", "full", "--report", "1"]
    print("[stage1] SpaceHub:", " ".join(cmd))
    for f in glob.glob(os.path.join(workdir, "bss_*_full.txt")):
        os.remove(f)
    with open(os.path.join(workdir, "stdout.log"), "w") as log:
        subprocess.run(cmd, cwd=workdir, stdout=log, stderr=subprocess.STDOUT, check=True)
    hits = sorted(glob.glob(os.path.join(workdir, "bss_*_full.txt")), key=os.path.getmtime)
    if not hits:
        sys.exit("[stage1] SpaceHub produced no bss_*_full.txt (see stdout.log)")
    return hits[-1]


def parse_bss_full(path: str):
    """Parse a SpaceHub DefaultWriter dump.  Format: optional header
    `time,id,mass,radius,px,py,pz,vx,vy,vz`, then one line PER BODY
    (ids 0,1,2 = SMBH, BH, star), three rows per output time.
    Returns t[yr], pos[nt,3,3] in AU, vel[nt,3,3] in AU/yr."""
    rows = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            try:
                vals = [float(x) for x in s.replace(",", " ").split()]
            except ValueError:
                continue                       # header line
            if len(vals) >= 10:
                rows.append(vals[:10])
    a = np.asarray(rows)                        # (3*nt, 10): time,id,mass,radius,px,py,pz,vx,vy,vz
    nt = a.shape[0] // 3
    a = a[:3 * nt].reshape(nt, 3, 10)
    order = np.argsort(a[0, :, 1])              # sort bodies by id (0,1,2) using first step
    a = a[:, order, :]
    t_internal = a[:, 0, 0]                     # SpaceHub internal time (year = 2*pi)
    pos = a[:, :, 4:7]                          # AU (already)
    vel = a[:, :, 7:10] * TWO_PI               # internal speed -> AU/yr  (1 unit = 2*pi AU/yr)
    t_yr = t_internal / TWO_PI
    return t_yr, pos, vel


def synth_stage1(cfg: Config, nt: int = 400):
    """DEMO ONLY: synthesize a short INCOMING trajectory (SMBH + BH + star) that ENDS
    at the disruption state, so the movie's stage-1 animation (incoming binary +
    intrinsic orbit) can be exercised WITHOUT SpaceHub.  We integrate the outer
    (SMBH-BH) and inner (BH-star) two-body orbits backward from disruption over half
    an outer period.  Not physically faithful -- use --run (which reads the real
    SpaceHub dump) for science.  Returns t[yr], pos[nt,3,3] AU, vel[nt,3,3] AU/yr."""
    rt_bin = cfg.a_bin * (cfg.M_smbh / (cfg.m_bh + cfg.m_star))**(1.0 / 3.0)
    Rp = rt_bin / cfg.beta
    Mout = cfg.M_smbh + cfg.m_bh
    e_out = 0.9                                          # marginally-bound outer orbit
    a_out = Rp / (1.0 - e_out)
    vperi = np.sqrt(G * Mout * (1 + e_out) / Rp)         # BH speed rel. SMBH at pericentre
    # --- disruption state at t = 0 (same geometry as before) ---
    r_sm0 = np.zeros(3);           v_sm0 = np.zeros(3)
    r_bh0 = np.array([Rp, 0., 0.]); v_bh0 = np.array([0., vperi, 0.])
    r_t = cfg.R_star * RSUN_AU * (cfg.m_bh / cfg.m_star)**(1.0 / 3.0)
    vrel = np.sqrt(G * cfg.m_bh * (1 + cfg.e_bin) / r_t)
    r_star0 = r_bh0 + np.array([0., 0., r_t])
    v_star0 = v_bh0 + np.array([vrel, 0., 0.])
    # --- integrate BOTH orbits backward from t=0 over half an outer period ---
    P_out = TWO_PI*np.sqrt(a_out**3/(G*Mout))
    tg_back = np.linspace(0.0, -0.5*P_out, nt)           # decreasing: t=0 is disruption
    rsm_t, rbh_t, _ = twobody_track(r_bh0, v_bh0, r_sm0, v_sm0, cfg.M_smbh, cfg.m_bh, tg_back)
    # inner binary: relative star-BH orbit (reuse twobody_track with BH as the "primary")
    ibh_t, istar_t, _ = twobody_track(r_star0, v_star0, r_bh0, v_bh0, cfg.m_bh, cfg.m_star, tg_back)
    rel_in = istar_t - ibh_t                             # star-BH separation(t)
    rstar_t = rbh_t + rel_in                             # ride the inner orbit on the BH track
    order = np.argsort(tg_back)                          # back to chronological order
    pos = np.stack([rsm_t, rbh_t, rstar_t], axis=1)[order]
    vel = np.tile(np.array([v_sm0, v_bh0, v_star0]), (nt, 1, 1))   # only [-1] is used downstream
    return tg_back[order], pos, vel


# =============================================================================
#  STAGE 2a -- debris initial conditions (Evans & Kochanek 1989 frozen-in)
# =============================================================================
def sample_polytrope_n3(N, R_star_AU, poly_n, rng):
    """Draw N positions inside a polytrope of index `poly_n` (n=3 ~ solar), scaled
    to radius R_star.  We sample the enclosed-mass CDF numerically from the
    Lane-Emden density, then assign isotropic directions."""
    # --- Lane-Emden theta(xi) for index n, to get rho(r) ~ theta^n ---
    xi = np.linspace(1e-6, 10.0, 4000)
    th = np.ones_like(xi); dth = np.zeros_like(xi)
    for i in range(1, len(xi)):
        h = xi[i] - xi[i-1]
        thc = max(th[i-1], 0.0)
        ddth = -thc**poly_n - 2.0/xi[i-1]*dth[i-1]       # Lane-Emden: th'' = -th^n - (2/xi)th'
        dth[i] = dth[i-1] + ddth*h
        th[i]  = th[i-1] + dth[i]*h
        if th[i] <= 0:                                   # first zero = stellar surface xi_1
            xi = xi[:i+1]; th = th[:i+1]; break
    th = np.clip(th, 0.0, None)
    rho = th**poly_n
    r   = xi / xi[-1] * R_star_AU                        # scale xi in [0,xi_1] -> [0,R_star]
    Menc = np.concatenate([[0.0], np.cumsum(0.5*(rho[1:]*r[1:]**2 + rho[:-1]*r[:-1]**2)
                                            * np.diff(r)) * 4*np.pi])
    Menc /= Menc[-1]
    u = rng.random(N)                                    # inverse-CDF sample of radius
    rr = np.interp(u, Menc, r)
    cth = rng.uniform(-1, 1, N); sth = np.sqrt(1-cth**2); ph = rng.uniform(0, TWO_PI, N)
    return np.column_stack([rr*sth*np.cos(ph), rr*sth*np.sin(ph), rr*cth])


def make_debris(cfg: Config, r_star, v_star, r_bh, v_bh, rng):
    """Replace the star by N debris: positions dr ~ polytrope about the star's
    centre; velocity = v_star + omega x dr, with omega the star's ORBITAL angular
    velocity about the BH at disruption (Evans & Kochanek 1989)."""
    dr = sample_polytrope_n3(cfg.n_deb, cfg.R_star * RSUN_AU, cfg.poly_n, rng)
    r_rel = r_star - r_bh
    v_rel = v_star - v_bh
    omega = np.cross(r_rel, v_rel) / np.dot(r_rel, r_rel)     # instantaneous orbital ang. velocity
    dpos = r_star[None, :] + dr
    dvel = v_star[None, :] + np.cross(omega[None, :], dr)     # bulk velocity + solid-body spin
    return dpos, dvel


# =============================================================================
#  STAGE 2b -- the ballistic machinery
# =============================================================================
def twobody_track(r_bh, v_bh, r_sm, v_sm, M_smbh, m_bh, t_grid):
    """Propagate the SMBH+BH two-body problem across t_grid (seconds of code: yr),
    in the inertial COM frame.  Returns r_sm(t), r_bh(t), v_bh(t) arrays.  The BHs
    move on their own orbit regardless of the (massless) debris."""
    from scipy.integrate import solve_ivp
    Mtot = M_smbh + m_bh
    # COM frame
    r_com = (M_smbh*r_sm + m_bh*r_bh)/Mtot
    v_com = (M_smbh*v_sm + m_bh*v_bh)/Mtot
    rel0 = np.concatenate([(r_bh-r_sm), (v_bh-v_sm)])         # relative BH-SMBH state
    def rhs(t, y):
        d = y[:3]; r3 = (d@d)**1.5
        return np.concatenate([y[3:], -G*Mtot*d/r3])
    sol = solve_ivp(rhs, (t_grid[0], t_grid[-1]), rel0, t_eval=t_grid,
                    method="DOP853", rtol=1e-11, atol=1e-13)
    rel = sol.y[:3].T                                        # (nt,3) BH-SMBH separation vector
    vrel = sol.y[3:].T
    r_sm_t = -(m_bh/Mtot)*rel + (r_com[None,:] + np.outer(t_grid - t_grid[0], v_com))
    r_bh_t =  (M_smbh/Mtot)*rel + (r_com[None,:] + np.outer(t_grid - t_grid[0], v_com))
    v_bh_t =  (M_smbh/Mtot)*vrel + v_com[None,:]
    return r_sm_t, r_bh_t, v_bh_t


def accel(dpos, r_sm, r_bh, M_smbh, m_bh, SS2, SB2):
    """Ballistic acceleration on each debris parcel from the two (softened) point
    masses.  Softening keeps the mini-disk finite and the integration stable."""
    a = np.zeros_like(dpos)
    ds = dpos - r_sm; rs3 = (np.einsum("ij,ij->i", ds, ds) + SS2)**1.5
    db = dpos - r_bh; rb3 = (np.einsum("ij,ij->i", db, db) + SB2)**1.5
    a -= G*M_smbh*ds/rs3[:,None]
    a -= G*m_bh  *db/rb3[:,None]
    return a


def integrate_debris(cfg: Config, r_sm0, v_sm0, r_bh0, v_bh0, dpos, dvel):
    """Two-phase velocity-Verlet integration of the debris in the time-dependent
    field of the moving BHs, with the stiff mini-disk FROZEN between phases.

    Returns a dict with: movie frames (times, debris positions in the SMBH frame,
    specific energy wrt the SMBH), the host classification, and the measured
    first-return times for the dM/dt curve."""
    r_t = cfg.R_star*RSUN_AU*(cfg.m_bh/cfg.m_star)**(1.0/3.0)   # mu-TDE tidal radius
    Rg_sm = G*cfg.M_smbh/C_AUYR**2
    SS2 = (3.0*Rg_sm)**2                 # SMBH softening (gravitational radius scale)
    SB2 = (0.3*r_t)**2                   # BH softening (a fraction of the tidal radius)
    N   = dpos.shape[0]
    m_parcel = cfg.m_star/N

    # characteristic times: mini-disk orbital period, and an estimate of the delay
    P_disr = TWO_PI*np.sqrt(r_t**3/(G*cfg.m_bh))
    t_relax = cfg.relax_orbits*P_disr
    dtA = P_disr/cfg.steps_per_Pdisr
    nA  = int(np.ceil(t_relax/dtA))

    # frozen-in energies at disruption (host assignment + the frozen-in dM/dt later)
    eps_bh0 = 0.5*np.einsum("ij,ij->i", dvel-v_bh0, dvel-v_bh0) - G*cfg.m_bh /np.linalg.norm(dpos-r_bh0, axis=1)
    eps_sm0 = 0.5*np.einsum("ij,ij->i", dvel-v_sm0, dvel-v_sm0) - G*cfg.M_smbh/np.linalg.norm(dpos-r_sm0, axis=1)

    # --- time grids for the two phases (phase B length set by the predicted delay) ---
    beta = cfg.beta
    Dt_pred = TWO_PI*np.sqrt(cfg.a_bin**3/(G*(cfg.m_bh+cfg.m_star)))*(cfg.M_smbh/(cfg.m_bh+cfg.m_star))**0.5/(8*beta)**0.75
    t_endB = cfg.return_margin*max(Dt_pred, 5*P_disr)
    dtB = (t_endB - t_relax)/max(cfg.steps_phaseB, 1)

    # one two-body track covering BOTH phases (dense; interpolated at every step)
    t_all = np.concatenate([np.arange(nA+1)*dtA, t_relax + (np.arange(1, cfg.steps_phaseB+1))*dtB])
    r_sm_t, r_bh_t, v_bh_t = twobody_track(r_bh0, v_bh0, r_sm0, v_sm0, cfg.M_smbh, cfg.m_bh, t_all)

    # frame bookkeeping: SLOW MOTION of the disruption -- pack half of the phase-A
    # frames into the first `slowmo_orbits` mini-disk periods (where the star tears
    # apart and the first debris fall back), then spread the rest over the relaxation.
    n_fast = min(int(round(cfg.slowmo_orbits*cfg.steps_per_Pdisr)), nA)
    fA = np.unique(np.concatenate([
        np.linspace(0, n_fast, cfg.frames_A//2).astype(int),
        np.linspace(n_fast, nA, cfg.frames_A - cfg.frames_A//2).astype(int)]))
    fB = np.unique(np.linspace(nA+1, len(t_all)-1, cfg.frames_B).astype(int))
    frame_idx = np.concatenate([fA, fB])
    frames = {"t": [], "pos": [], "eps": [], "eps_bh": [], "r_sm": [], "r_bh": [], "phase": []}

    # state
    x = dpos.copy(); v = dvel.copy()
    frozen = np.zeros(N, bool)          # mini-disk parcels, rigidly riding the BH
    off = np.zeros((N, 3))              # their offset from the BH at freeze time
    tret_bh = np.full(N, np.nan)        # first return-to-pericentre about the BH
    tret_sm = np.full(N, np.nan)        # first return-to-pericentre about the SMBH
    prev_rb = np.linalg.norm(x - r_bh0, axis=1); rising_b = np.zeros(N, bool)
    prev_rs = np.linalg.norm(x - r_sm0, axis=1); rising_s = np.zeros(N, bool)

    def record(i, phase):
        rs = r_sm_t[i]; rb = r_bh_t[i]; vb = v_bh_t[i]
        xx = x.copy(); vv = v.copy()
        if frozen.any():
            # the mini-disk is RIGIDLY locked to the BH: it comoves with the BH, so for
            # energy purposes its velocity IS the BH velocity (not its frozen orbital v --
            # that internal motion is around the BH and must not inflate 0.5|v|^2 wrt the SMBH).
            xx[frozen] = rb + off[frozen]; vv[frozen] = vb
        rel = xx - rs                              # debris position in the SMBH frame
        # INSTANTANEOUS specific energies (they evolve while the potential is time-varying)
        eps_sm = 0.5*np.einsum("ij,ij->i", vv, vv) - G*cfg.M_smbh/np.linalg.norm(rel, axis=1)
        eps_bh = 0.5*np.einsum("ij,ij->i", vv-vb, vv-vb) - G*cfg.m_bh/np.linalg.norm(xx-rb, axis=1)
        frames["t"].append(t_all[i]); frames["pos"].append(rel.copy())
        frames["eps"].append(eps_sm); frames["eps_bh"].append(eps_bh)
        frames["r_sm"].append(rs.copy()); frames["r_bh"].append(rb - rs); frames["phase"].append(phase)

    acc = accel(x, r_sm_t[0], r_bh_t[0], cfg.M_smbh, cfg.m_bh, SS2, SB2)
    fset = set(frame_idx.tolist())
    for i in range(1, len(t_all)):
        dt = t_all[i] - t_all[i-1]
        live = ~frozen
        # velocity-Verlet (KDK) for the live (non-frozen) parcels
        v[live] += 0.5*dt*acc[live]
        x[live] += dt*v[live]
        acc_new = accel(x, r_sm_t[i], r_bh_t[i], cfg.M_smbh, cfg.m_bh, SS2, SB2)
        v[live] += 0.5*dt*acc_new[live]
        acc = acc_new

        # ---- first-return detection (pericentre about each host = local min of r) ----
        rb = np.linalg.norm(x - r_bh_t[i], axis=1)
        rs = np.linalg.norm(x - r_sm_t[i], axis=1)
        rising_b |= (rb > prev_rb)      # must go out before coming back
        hit_b = live & rising_b & (rb < prev_rb) & np.isnan(tret_bh) & (eps_bh0 < 0)
        tret_bh[hit_b] = t_all[i]
        rising_s |= (rs > prev_rs)
        hit_s = live & rising_s & (rs < prev_rs) & np.isnan(tret_sm) & (eps_sm0 < 0)
        tret_sm[hit_s] = t_all[i]
        prev_rb, prev_rs = rb, rs

        # ---- FREEZE the mini-disk at the end of phase A ----
        if i == nA:
            disk = live & (eps_bh0 < 0) & (rb < cfg.disk_capture_rt*r_t)
            off[disk] = x[disk] - r_bh_t[i]
            frozen |= disk
            print(f"[stage2] froze {disk.sum()} mini-disk parcels at t={t_all[i]:.3e} yr "
                  f"({cfg.relax_orbits:g} BH orbits); {(~frozen).sum()} parcels continue to the SMBH.")

        if i in fset:
            record(i, 0 if i <= nA else 1)

    # ---- classification at DISRUPTION (frozen-in energies, Evans-Kochanek snapshot) ----
    cls = np.where(eps_bh0 < 0, 0, np.where(eps_sm0 < 0, 1, 2))   # 0=BH-bound,1=SMBH-bound,2=ejected

    # ---- SETTLED classification from the FINAL-frame energies (after the BH has receded
    #      and the orbits have relaxed -- captures the post-disruption energy exchange the
    #      frozen-in estimate misses).  Mini-disk parcels comove with the BH (vv=vb). ----
    rs_f, rb_f, vb_f = r_sm_t[-1], r_bh_t[-1], v_bh_t[-1]
    xx = x.copy(); vv = v.copy()
    if frozen.any():
        xx[frozen] = rb_f + off[frozen]; vv[frozen] = vb_f
    eps_sm_f = 0.5*np.einsum("ij,ij->i", vv, vv) - G*cfg.M_smbh/np.linalg.norm(xx-rs_f, axis=1)
    eps_bh_f = 0.5*np.einsum("ij,ij->i", vv-vb_f, vv-vb_f) - G*cfg.m_bh/np.linalg.norm(xx-rb_f, axis=1)
    cls_final = np.where(eps_bh_f < 0, 0, np.where(eps_sm_f < 0, 1, 2))

    return dict(frames=frames, cls=cls, cls_final=cls_final,
                eps_bh0=eps_bh0, eps_sm0=eps_sm0, eps_bh_f=eps_bh_f, eps_sm_f=eps_sm_f,
                tret_bh=tret_bh, tret_sm=tret_sm, r_t=r_t, m_parcel=m_parcel,
                P_disr=P_disr, t_relax=t_relax, Dt_pred=Dt_pred)


# =============================================================================
#  dM/dt -- numerical (first-return times) vs frozen-in (Keplerian) comparison
# =============================================================================
def fallback_curves(cfg: Config, res, nbins=40):
    """Measured dM/dt from first-return times (BH + SMBH channels), and the
    frozen-in Keplerian estimate for the same parcels (for comparison)."""
    m = res["m_parcel"]
    tb, ts = res["tret_bh"], res["tret_sm"]
    tmeas = np.concatenate([tb[np.isfinite(tb)], ts[np.isfinite(ts)]])
    # frozen-in: t_K = 2*pi*G*M*(2|eps|)^-3/2 for the bound parcels about each host
    eb, es = res["eps_bh0"], res["eps_sm0"]
    tK = np.full(len(eb), np.nan)
    mb = eb < 0; ms = (~mb) & (es < 0)
    tK[mb] = TWO_PI*(G*cfg.m_bh)  *(2*np.abs(eb[mb]))**-1.5
    tK[ms] = TWO_PI*(G*cfg.M_smbh)*(2*np.abs(es[ms]))**-1.5
    allt = np.concatenate([tmeas, tK[np.isfinite(tK)]])
    lo, hi = np.nanpercentile(allt, [0.3, 99.7])
    bins = np.logspace(np.log10(lo), np.log10(hi), nbins+1)
    cen  = np.sqrt(bins[:-1]*bins[1:]); dt = np.diff(bins)
    h = lambda arr: m*np.histogram(arr[np.isfinite(arr)], bins=bins)[0]/dt
    return dict(t=cen,
                bh=h(tb), sm=h(ts), tot=h(np.concatenate([tb, ts])),
                K_bh=h(np.where(mb, tK, np.nan)), K_sm=h(np.where(ms, tK, np.nan)),
                K_tot=h(tK))


# =============================================================================
#  (movie + money plot live in mutde_lab_plots.py-style functions below)
# =============================================================================
def _apply_latex(plt, want=True):
    """Use real LaTeX for all text when a latex+dvipng toolchain is present; otherwise
    fall back to Computer-Modern mathtext with a serif family (same look, no TeX needed)."""
    import shutil
    ok = bool(want) and bool(shutil.which("latex")) and bool(shutil.which("dvipng"))
    plt.rcParams.update({
        "text.usetex": ok,
        "font.family": "serif",
        # NB: avoid cmr10 as the text font -- it lacks the Unicode minus (renders as a box).
        # DejaVu Serif has full coverage; math still uses the CM fontset for the LaTeX look.
        "font.serif": ["CMU Serif", "DejaVu Serif", "Times New Roman"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
    })
    return ok


def _sym_norm(eps_ergg):
    import matplotlib.colors as mcolors
    v = np.nanpercentile(np.abs(eps_ergg), 98) or 1.0
    return mcolors.TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)


C_BH, C_SM = "#E69F00", "#0072B2"       # orange = BH channel, blue = SMBH channel (used everywhere)


def _symlog_bins(vmax, lt, nside=26):
    """Symmetric log-spaced bin edges on [-vmax, vmax] with a linear core |x|<lt."""
    pos = np.logspace(np.log10(lt), np.log10(vmax), nside)
    return np.concatenate([-pos[::-1], [0.0], pos])


def _smooth_logy(xc, y, size=3):
    """Light smoothing of a binned dM/dt curve in log-space, over the populated bins only,
    so the Poisson bin-to-bin scatter is damped without shifting the overall shape."""
    g = y > 0
    xs = xc[g]; ys = y[g]
    if ys.size < 3:
        return xs, ys
    from scipy.ndimage import uniform_filter1d
    return xs, 10.0**uniform_filter1d(np.log10(ys), size=size, mode="nearest")


def render_movie(cfg: Config, res, path_gif, path_mp4=None, stage1=None, render_max=20000, still_base=None):
    """Animation in a FIXED, SMBH-centred frame.  Two acts:
      STAGE 1 (from SpaceHub): the incoming SMBH+BH+star motion -- the binary enters
               from the RIGHT and spirals in to disruption; the inset resolves the
               intrinsic binary orbit (star about the BH);
      STAGE 2 (ballistic): the disruption -> mini-disk -> SMBH-return of the debris.
    Panels: (left) debris/bodies coloured by specific energy wrt the SMBH + zoom inset;
            (top-R) energy distributions dN/deps (eps_BH orange, eps_SMBH blue);
            (bot-R) DIRECT dM/dt (solid) + frozen-in (dashed), building up in time.
    `stage1` = (t[yr], pos[nt,3,3], vel[nt,3,3]) for bodies (SMBH, BH, star), as read
    from the SpaceHub dump (or synthesised in --demo)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.ticker import MaxNLocator, SymmetricalLogLocator
    from matplotlib.animation import FuncAnimation, PillowWriter
    _apply_latex(plt, cfg.use_tex)
    fr = res["frames"]; nf = len(fr["t"])
    pos     = fr["pos"]
    epsSraw = fr["eps"]        # specific energy wrt SMBH  [(AU/yr)^2]
    epsBraw = fr["eps_bh"]     # specific energy wrt BH    [(AU/yr)^2]
    r_bh = fr["r_bh"]; phase = fr["phase"]; tarr = fr["t"]
    N    = pos[0].shape[0]
    # the movie subsamples the on-screen scatter for speed/size; a still draws ALL points.
    pick = np.arange(N) if (still_base is not None or N <= render_max) \
        else np.random.default_rng(0).choice(N, render_max, False)

    # ---- frozen-in energy spread (Evans-Kochanek): Delta_eps = G m_bh R_star / r_t^2 with
    #      r_t = R_star (m_bh/m_star)^(1/3).  Normalises BOTH the colour scale and the
    #      energy-distribution x-axis. ----
    r_t = res["r_t"]
    d_eps = G*cfg.m_bh*(cfg.R_star*RSUN_AU)/r_t**2 or 1.0   # [(AU/yr)^2]
    cvmax = np.nanpercentile(np.abs(np.concatenate(epsSraw))/d_eps, 99.5) or 1.0
    cnorm = mcolors.TwoSlopeNorm(vmin=-cvmax, vcenter=0.0, vmax=cvmax)
    colS = [e/d_eps for e in epsSraw]                      # eps_bullet / Delta_eps per frame

    # ---- fixed, SMBH-centred frame half-width (covers the whole returning stream) ----
    Rmax = max(np.nanpercentile(np.linalg.norm(pp[:, :2], axis=1), 97) for pp in pos)
    Hv = 1.2*max(Rmax, 10.0)

    # ---- stage-1 incoming trajectory: select the part within the frame, rotate so the
    #      binary approaches from the RIGHT (+x), and resample to frames_S frames ----
    def _rotm(th): c, s = np.cos(th), np.sin(th); return np.array([[c, -s], [s, c]])
    S = None; Rot = np.eye(2)
    if stage1 is not None:
        t1, pos1, vel1 = stage1
        rsm1 = pos1[:, 0, :2]; rbh1 = pos1[:, 1, :2]; rst1 = pos1[:, 2, :2]
        dist = np.linalg.norm(rbh1 - rsm1, axis=1)
        within = np.where(dist <= Hv)[0]
        i0 = int(within[0]) if within.size else 0
        idxS = np.unique(np.linspace(i0, len(t1)-1, cfg.frames_S).astype(int))
        b0 = rbh1[idxS[0]] - rsm1[idxS[0]]
        Rot = _rotm(-np.arctan2(b0[1], b0[0]))          # map the first shown BH dir -> +x
        S = dict(idx=idxS, rsm=rsm1, rbh=rbh1, rst=rst1, t=t1, tdisr=t1[-1])
    nS = len(S["idx"]) if S else 0
    rot = lambda a2: np.asarray(a2) @ Rot.T              # apply the same rotation everywhere

    # ---- energy-distribution x-axis: eps/Delta_eps (same frozen-in spread as the colour scale) ----
    epsS_n = [e/d_eps for e in epsSraw]; epsB_n = [e/d_eps for e in epsBraw]
    xmax = np.nanpercentile(np.abs(np.concatenate([np.concatenate(epsS_n), np.concatenate(epsB_n)])), 99.5)
    elt  = max(xmax/1e3, 1e-6)
    ebins = _symlog_bins(xmax, elt); ecen = 0.5*(ebins[1:]+ebins[:-1])

    # ---- dM/dt axis (log-t in DAYS); bins span measured + frozen-in return times ----
    tb, ts = res["tret_bh"], res["tret_sm"]; m = res["m_parcel"]
    eb, es = res["eps_bh0"], res["eps_sm0"]
    mb = eb < 0; ms = (~mb) & (es < 0)                         # BH-bound / SMBH-bound parcels
    tK = np.full(len(eb), np.nan)                              # frozen-in t_K = 2piGM (2|eps|)^-3/2
    tK[mb] = TWO_PI*(G*cfg.m_bh)  *(2*np.abs(eb[mb]))**-1.5
    tK[ms] = TWO_PI*(G*cfg.M_smbh)*(2*np.abs(es[ms]))**-1.5
    tKb = np.where(mb, tK, np.nan); tKs = np.where(ms, tK, np.nan)
    fin = lambda arr: arr[np.isfinite(arr)]
    allret = np.concatenate([fin(tb), fin(ts), fin(tK)])
    tlo, thi = np.nanpercentile(allret, [0.5, 99.5]); tlo = max(tlo, 1e-6)
    tbins = np.logspace(np.log10(tlo), np.log10(thi), 28); tcc = np.sqrt(tbins[1:]*tbins[:-1]); dtw = np.diff(tbins)
    tcc_d = tcc*DAY_PER_YR                                    # bin centres in DAYS (for the x-axis)
    def mdot_upto(ret, tnow):
        sub = ret[np.isfinite(ret) & (ret <= tnow)]
        return m*np.histogram(sub, bins=tbins)[0]/dtw        # [Msun/yr]
    ymax_mdot = max(1e-30,
                    (m*np.histogram(np.concatenate([fin(tb), fin(ts)]), bins=tbins)[0]/dtw).max(),
                    (m*np.histogram(fin(tK), bins=tbins)[0]/dtw).max())

    # ---------- figure / static styling ----------
    # NB: axis labels use axes.labelcolor (NOT text.color) -- set it white or they render
    # black and vanish on the black background.
    plt.rcParams.update({"text.color": "w", "axes.edgecolor": "w", "axes.labelcolor": "w"})
    fig = plt.figure(figsize=(13.4, 7.4)); fig.patch.set_facecolor("black")
    # explicit margins so the axis labels are never clipped (animation has no tight-bbox)
    gs  = fig.add_gridspec(2, 2, width_ratios=[1.62, 1.0], height_ratios=[1, 1],
                           left=0.07, right=0.965, bottom=0.095, top=0.95, wspace=0.34, hspace=0.36)
    ax  = fig.add_subplot(gs[:, 0]); axE = fig.add_subplot(gs[0, 1]); axF = fig.add_subplot(gs[1, 1])
    def _style(a, lblsize=9):
        """white major+minor ticks (pointing in) on a black panel."""
        a.set_facecolor("black")
        a.tick_params(which="major", colors="w", labelcolor="w", direction="in",
                      length=5, width=0.9, top=True, right=True, labelsize=lblsize)
        a.tick_params(which="minor", colors="w", direction="in",
                      length=2.8, width=0.6, top=True, right=True)
        for s in a.spines.values(): s.set_color("w")
    for a in (ax, axE, axF): _style(a)
    ax.minorticks_on()      # linear scatter: show minor ticks between labels
    sc = ax.scatter([], [], s=2.5, c=[], cmap="coolwarm", norm=cnorm, lw=0)
    trailln, = ax.plot([], [], "-", color="0.55", lw=0.9, alpha=0.8)      # stage-1 BH approach track
    bhdot, = ax.plot([], [], "o", mfc="black", mec="w", ms=7,  mew=1.3)   # stellar BH
    stardot, = ax.plot([], [], "*", mfc="gold", mec="darkorange", ms=12, mew=0.6)  # the star (stage 1)
    smdot, = ax.plot([], [], "o", mfc="black", mec="w", ms=11, mew=1.6)   # SMBH (both are black holes)
    txt = ax.text(0.03, 0.97, "", transform=ax.transAxes, color="w", fontsize=12, va="top")
    ax.set_xlabel(r"$x-x_{\bullet}$ [AU]"); ax.set_ylabel(r"$y-y_{\bullet}$ [AU]")
    ax.set_xlim(-Hv, Hv); ax.set_ylim(-Hv, Hv); ax.set_aspect("equal")   # fixed frame (set ONCE)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(r"$\varepsilon_{\bullet}/\Delta\varepsilon_{\rm frozen}$")
    cb.ax.yaxis.set_tick_params(color="w"); plt.setp(plt.getp(cb.ax, "yticklabels"), color="w")
    # zoom inset (BH frame): binary orbit in stage 1, mini-disk in stage 2 -- ticks + labels,
    # placed well inside so its tick labels do not run into the main panel.
    axi = ax.inset_axes([0.075, 0.10, 0.29, 0.29]); _style(axi, 5); axi.minorticks_on()
    axi.xaxis.set_major_locator(MaxNLocator(3)); axi.yaxis.set_major_locator(MaxNLocator(3))
    axi.tick_params(pad=1.5)
    axi.set_title(r"mini-disk (zoom)", color="w", fontsize=8)
    axi.set_xlabel(r"$x-x_{\rm BH}$ [AU]", fontsize=6, labelpad=1)
    axi.set_ylabel(r"$y-y_{\rm BH}$ [AU]", fontsize=6, labelpad=1)
    sci = axi.scatter([], [], s=2, c=[], cmap="coolwarm", norm=cnorm, lw=0)
    starnow, = axi.plot([], [], "*", mfc="gold", mec="darkorange", ms=9, mew=0.5)  # star (no trail)
    bhzoom, = axi.plot([], [], "o", mfc="black", mec="w", ms=6, mew=1.0)
    # text labels next to the black holes (SMBH fixed at origin; BH position updated per frame)
    smbh_lbl = ax.text(0.02*Hv, 0.02*Hv, "SMBH", color="w", fontsize=9, ha="left", va="bottom")
    bh_lbl   = ax.text(0, 0, "BH", color="w", fontsize=9, ha="right", va="bottom")
    bh_lbl_i = axi.text(0, 0, "BH", color="w", fontsize=6, ha="left", va="bottom")

    def _placeholder(a, msg):                               # right panels during stage 1
        a.cla(); _style(a, 7); a.set_xticks([]); a.set_yticks([])
        a.text(0.5, 0.5, msg, transform=a.transAxes, color="0.55", ha="center", va="center", fontsize=9)

    def upd(k):
        smdot.set_data([0], [0]); bhzoom.set_data([0], [0])

        # ===================== STAGE 1: incoming binary (from SpaceHub) =====================
        if S and k < nS:
            j = S["idx"][k]
            sc.set_offsets(np.empty((0, 2)))                # no debris yet
            bh2 = rot(S["rbh"][j] - S["rsm"][j]); st2 = rot(S["rst"][j] - S["rsm"][j])
            bhdot.set_data([bh2[0]], [bh2[1]]); stardot.set_data([st2[0]], [st2[1]])
            bh_lbl.set_position((bh2[0]-0.015*Hv, bh2[1]+0.015*Hv))
            jj = S["idx"][:k+1]
            tr = rot(S["rbh"][jj] - S["rsm"][jj]); trailln.set_data(tr[:, 0], tr[:, 1])   # approach track
            rel = rot(S["rst"][j] - S["rbh"][j])            # star position about the BH (no orbit trail)
            starnow.set_data([rel[0]], [rel[1]])
            sci.set_offsets(np.empty((0, 2)))
            hw = cfg.zoom_bin_fac*cfg.a_bin; axi.set_xlim(-hw, hw); axi.set_ylim(-hw, hw)
            bh_lbl_i.set_position((0.1*hw, 0.1*hw))
            txt.set_text(r"$t = %+.1f$ days" % ((S['t'][j]-S['tdisr'])*DAY_PER_YR) + "\nincoming binary")
            _placeholder(axE, "(pre-disruption)"); _placeholder(axF, "(pre-disruption)")
            return sc,

        # ===================== STAGE 2: ballistic debris =====================
        k2 = k - nS
        p = pos[k2][pick]; cs = colS[k2][pick]
        sc.set_offsets(rot(p[:, :2])); sc.set_array(cs)
        rb2 = rot(r_bh[k2][:2]); bhdot.set_data([rb2[0]], [rb2[1]]); stardot.set_data([], []); trailln.set_data([], [])
        bh_lbl.set_position((rb2[0]-0.015*Hv, rb2[1]+0.015*Hv))
        # inset: ALL debris inside the BH-centred box (square cut, not a circle), so parcels
        # stay visible right up to the inset edges.
        hw = 6*res["r_t"]; dxy = rot(p[:, :2]-r_bh[k2][:2])
        inbox = (np.abs(dxy[:, 0]) < hw) & (np.abs(dxy[:, 1]) < hw)
        sci.set_offsets(dxy[inbox] if inbox.any() else np.empty((0, 2)))
        sci.set_array(cs[inbox] if inbox.any() else np.array([]))
        starnow.set_data([], []); bh_lbl_i.set_position((0.1*hw, 0.1*hw))
        axi.set_xlim(-hw, hw); axi.set_ylim(-hw, hw)
        ph = "disruption / mini-disk" if phase[k2] == 0 else "return to SMBH"
        txt.set_text(r"$t = %+.2f$ days" % (tarr[k2]*DAY_PER_YR) + "\n" + ph)

        # --- energy distributions (morphing in time), x = 2*eps/v_inf^2 ---
        axE.cla()
        hB, _ = np.histogram(epsB_n[k2], bins=ebins); hS, _ = np.histogram(epsS_n[k2], bins=ebins)
        axE.step(ecen, hB+0.1, where="mid", color=C_BH, lw=1.5, label=r"$\varepsilon_{\rm BH}$")
        axE.step(ecen, hS+0.1, where="mid", color=C_SM, lw=1.5, label=r"$\varepsilon_{\bullet}$")
        axE.axvline(0, color="0.5", lw=0.8, ls="--")
        axE.set_xscale("symlog", linthresh=elt); axE.set_yscale("log")
        axE.set_xlim(-xmax, xmax); axE.set_ylim(0.5, N)
        # label the |x|>=1 decades (+0): keep +/-10^0, +/-10^1, ... but drop the
        # crowded +/-10^-1; minor ticks fill in across the whole symlog x-axis.
        ks = list(range(0, int(np.floor(np.log10(xmax)))+1)) or [0]
        axE.set_xticks([0.0] + [s*10.0**kk for kk in ks for s in (-1, 1)])
        axE.xaxis.set_minor_locator(SymmetricalLogLocator(base=10.0, linthresh=elt, subs=np.arange(2, 10)))
        axE.set_xlabel(r"$\varepsilon/\Delta\varepsilon_{\rm frozen}$"); axE.set_ylabel(r"$dN/d\varepsilon$")
        axE.legend(fontsize=8, frameon=False, labelcolor="w", loc="upper left")
        _style(axE, 7)

        # --- dM/dt (DIRECT solid + frozen-in dashed), x in DAYS ---
        axF.cla()
        yB = mdot_upto(tb, tarr[k2]);  yS = mdot_upto(ts, tarr[k2])      # direct (measured returns)
        yKb = mdot_upto(tKb, tarr[k2]); yKs = mdot_upto(tKs, tarr[k2])   # frozen-in (Keplerian)
        # light log-space smoothing damps the finite-N Poisson wiggles in the building curves
        drew = False
        for yy, col, lsty, lww, lab in [(yB, C_BH, "-", 1.7, r"$dM/dt_{\rm BH}$"),
                                        (yS, C_SM, "-", 1.7, r"$dM/dt_{\bullet}$"),
                                        (yKb, C_BH, "--", 1.2, r"BH frozen-in"),
                                        (yKs, C_SM, "--", 1.2, r"$\bullet$ frozen-in")]:
            xs, ys = _smooth_logy(tcc_d, yy)
            if xs.size:
                axF.loglog(xs, ys, color=col, lw=lww, ls=lsty, alpha=(0.85 if lsty == "--" else 1.0), label=lab)
                drew = True
        axF.axvline(tarr[k2]*DAY_PER_YR, color="w", lw=1.0, ls=":")      # current time marker (days)
        axF.set_xlim(tlo*DAY_PER_YR, thi*DAY_PER_YR); axF.set_ylim(ymax_mdot*1e-5, ymax_mdot*3)
        axF.set_xlabel(r"$t$ [days]"); axF.set_ylabel(r"$dM/dt$ [$M_\odot\,{\rm yr}^{-1}$]")
        if drew:
            axF.legend(fontsize=6.5, frameon=False, labelcolor="w", loc="upper right", ncol=2)
        _style(axF, 7)
        return sc,

    # --- STILL mode: render only the final frame (full return) at paper resolution ---
    if still_base is not None:
        upd(nS+nf-1)
        fig.savefig(still_base+".png", dpi=300, facecolor="black", bbox_inches="tight")
        fig.savefig(still_base+".pdf", facecolor="black", bbox_inches="tight")
        print("[still]", still_base+".png", "+", still_base+".pdf"); plt.close(fig); return

    anim = FuncAnimation(fig, upd, frames=nS+nf, blit=False)
    anim.save(path_gif, writer=PillowWriter(fps=cfg.fps)); print("[movie]", path_gif)
    if path_mp4 and shutil.which("ffmpeg"):
        try:
            anim.save(path_mp4, fps=cfg.fps, dpi=130); print("[movie]", path_mp4)
        except Exception as e:
            print("[movie] mp4 skipped:", e)
    plt.close(fig)


def money_plot(cfg: Config, res, fb, path):
    """The Fig.-1 'money plot': final snapshot in the SMBH frame coloured by specific
    energy (left), and the direct vs frozen-in two-peak dM/dt (right).  Two SEPARATE
    panels (gridspec) so nothing overlaps regardless of the debris extent."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    _apply_latex(plt, cfg.use_tex)
    fr = res["frames"]; p = fr["pos"][-1]; e = fr["eps"][-1]           # (AU/yr)^2
    r_t = res["r_t"]; d_eps = G*cfg.m_bh*(cfg.R_star*RSUN_AU)/r_t**2 or 1.0   # frozen-in spread
    cvmax = np.nanpercentile(np.abs(e)/d_eps, 99.5) or 1.0
    cnorm = mcolors.TwoSlopeNorm(vmin=-cvmax, vcenter=0.0, vmax=cvmax)
    fig = plt.figure(figsize=(12.6, 6.2))
    gs  = fig.add_gridspec(1, 2, width_ratios=[1.22, 1.0], wspace=0.30)
    ax  = fig.add_subplot(gs[0, 0]); axd = fig.add_subplot(gs[0, 1])

    # ---- left: debris snapshot, coloured by eps_bullet / Delta_eps_frozen ----
    sc = ax.scatter(p[:, 0], p[:, 1], s=3, c=e/d_eps, cmap="coolwarm", norm=cnorm, lw=0)
    ax.plot(0, 0, "o", mfc="k", mec="0.4", ms=11, mew=1.3)                          # SMBH
    bx, by = fr["r_bh"][-1][0], fr["r_bh"][-1][1]
    ax.plot(bx, by, "o", mfc="k", mec="0.4", ms=6, mew=1.1)                         # stellar BH
    ax.set_xlabel(r"$x-x_{\bullet}$ [AU]"); ax.set_ylabel(r"$y-y_{\bullet}$ [AU]"); ax.set_aspect("equal")
    # clip to the bulk so a few far-flung ejecta do not squash the interesting region
    half = 1.15*np.nanpercentile(np.linalg.norm(p[:, :2], axis=1), 97) + 1e-6
    ax.set_xlim(-half, half); ax.set_ylim(-half, half)
    ax.text(0.02*half, 0.02*half, "SMBH", fontsize=9, ha="left", va="bottom")
    ax.text(bx-0.02*half, by+0.02*half, "BH", fontsize=9, ha="right", va="bottom")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(r"$\varepsilon_{\bullet}/\Delta\varepsilon_{\rm frozen}$")
    # fractions from the SETTLED (final-energy) classification, so they match the colours
    clf = res.get("cls_final", res["cls"]); fBH=(clf==0).mean(); fSM=(clf==1).mean(); fEJ=(clf==2).mean()
    ax.set_title(r"$\mu$TDE debris (settled):  $f(\varepsilon_{\rm BH}<0)=%.2f$,  "
                 r"$f(\varepsilon_{\bullet}<0)=%.2f$,  $f(\varepsilon_{\bullet}>0)=%.2f$"
                 % (fBH, fSM, fEJ), fontsize=10)

    # ---- right: two-peak mass fallback (DAYS), direct (solid) vs frozen-in (dashed) ----
    td = fb["t"]*DAY_PER_YR
    for y, c, l in [(fb["tot"], "k", r"total"), (fb["bh"], C_BH, r"to BH ($dM/dt_{\rm BH}$)"),
                    (fb["sm"], C_SM, r"to SMBH ($dM/dt_{\bullet}$)")]:
        xs, ys = _smooth_logy(td, y)
        if xs.size: axd.loglog(xs, ys, c=c, lw=1.6, label=l)
    xs, ys = _smooth_logy(td, fb["K_tot"])
    if xs.size: axd.loglog(xs, ys, "k--", lw=1.0, alpha=0.7, label=r"frozen-in (total)")
    axd.set_xlabel(r"$t$ [days]"); axd.set_ylabel(r"$dM/dt$ [$M_\odot\,{\rm yr}^{-1}$]")
    axd.set_title(r"mass fallback rate (direct)", fontsize=10)
    axd.legend(fontsize=8, frameon=False); axd.grid(True, which="both", alpha=0.15)

    fig.savefig(path, dpi=160, bbox_inches="tight"); print("[money]", path); plt.close(fig)


# =============================================================================
#  main
# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="single-case micro-TDE + debris fallback lab")
    ap.add_argument("--run", action="store_true", help="run SpaceHub for stage 1")
    ap.add_argument("--traj", help="reuse an existing bss_*_full.txt instead of running SpaceHub")
    ap.add_argument("--demo", action="store_true", help="synthesize stage 1 (no SpaceHub)")
    ap.add_argument("--ndeb", type=int, help="number of debris particles (override config)")
    ap.add_argument("--out", help="output directory (default: data/debris_mutde/out)")
    ap.add_argument("--exe", help="path to smbh_binary_PN (default: simulation/bin/smbh_binary_PN)")
    ap.add_argument("--no-movie", action="store_true"); ap.add_argument("--no-money", action="store_true")
    ap.add_argument("--still", action="store_true",
                    help="save the final-frame composite as PNG+PDF (all points, 300 dpi) instead of the movie")
    a = ap.parse_args()

    cfg = Config()
    if a.ndeb: cfg.n_deb = a.ndeb
    if a.out:  cfg.out = a.out
    if a.exe:  cfg.bss_exe = a.exe
    os.makedirs(cfg.out, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)

    # ---- STAGE 1 ----
    if a.traj:
        t, pos, vel = parse_bss_full(a.traj)
    elif a.demo:
        t, pos, vel = synth_stage1(cfg)
    else:
        work = os.path.join(cfg.out, "s1"); os.makedirs(work, exist_ok=True)
        t, pos, vel = parse_bss_full(run_spacehub(cfg, work))
    r_sm, r_bh, r_star = pos[-1]            # disruption = final snapshot
    v_sm, v_bh, v_star = vel[-1]
    d_smbh  = np.linalg.norm(r_bh - r_sm)              # BH-SMBH separation at disruption
    rt_bin  = cfg.a_bin*(cfg.M_smbh/(cfg.m_bh+cfg.m_star))**(1.0/3.0)   # binary tidal radius wrt SMBH
    print(f"[stage1] disruption state: |r_bh-r_sm|={d_smbh:.4g} AU  (= {d_smbh/rt_bin:.2f} r_t,bin), "
          f"|r_star-r_bh|={np.linalg.norm(r_star-r_bh):.4g} AU")
    # SANITY GUARD: a genuine micro-TDE is driven by the SMBH tide, so the star must be
    # disrupted DEEP in the SMBH potential (|r_bh-r_sm| ~ r_p = r_t,bin/beta, i.e. a few
    # r_t,bin at most).  If instead it happened many tidal radii out, the binary
    # self-disrupted on its own pericentre and NO debris will return to the SMBH.
    if d_smbh > 10.0*rt_bin:
        print("\n*** WARNING: the star was disrupted {:.0f} tidal radii from the SMBH. ***".format(d_smbh/rt_bin))
        peri_bin = cfg.a_bin*(1.0 - cfg.e_bin)
        r_t_star = cfg.R_star*RSUN_AU*(cfg.m_bh/cfg.m_star)**(1.0/3.0)
        print("    The binary self-disrupted far from the SMBH, so no debris will return")
        print("    and you will NOT get a double flare.  Cause: the inner binary is at/inside")
        print("    its Roche limit in isolation -- pericentre a_bin*(1-e_bin) = {:.4g} AU".format(peri_bin))
        print("    vs stellar tidal radius r_t = {:.4g} AU.".format(r_t_star))
        print("    FIX: lower e_bin so a_bin*(1-e_bin) > r_t  (here e_bin < {:.3f}).\n".format(
            1.0 - r_t_star/cfg.a_bin))

    # ---- STAGE 2 ----
    dpos, dvel = make_debris(cfg, r_star, v_star, r_bh, v_bh, rng)
    res = integrate_debris(cfg, r_sm, v_sm, r_bh, v_bh, dpos, dvel)
    fb = fallback_curves(cfg, res)
    cls = res["cls"]; clf = res["cls_final"]
    print(f"[result] frozen-in (at disruption):  mini-disk = {(cls==0).mean():.3f}   "
          f"SMBH-bound = {(cls==1).mean():.3f}   ejected = {(cls==2).mean():.3f}")
    print(f"[result] settled   (final energies): mini-disk = {(clf==0).mean():.3f}   "
          f"SMBH-bound = {(clf==1).mean():.3f}   ejected = {(clf==2).mean():.3f}  "
          f"<- matches the figure colours")
    print(f"[result] predicted Dt_delay = {res['Dt_pred']*12:.2f} months")

    np.savez_compressed(os.path.join(cfg.out, "fallback.npz"),
                        t=fb["t"], bh=fb["bh"], sm=fb["sm"], tot=fb["tot"], K_tot=fb["K_tot"])
    if not a.no_movie:
        if a.still:                                   # paper still of the final frame (PNG + PDF)
            render_movie(cfg, res, None, None, stage1=(t, pos, vel),
                         still_base=os.path.join(cfg.out, "final_still"))
        else:                                         # the full animation (gif + mp4)
            render_movie(cfg, res, os.path.join(cfg.out, "mutde.gif"), os.path.join(cfg.out, "mutde.mp4"),
                         stage1=(t, pos, vel))     # the full SMBH+BH+star trajectory read from SpaceHub
    if not a.no_money:
        money_plot(cfg, res, fb, os.path.join(cfg.out, "money_plot.png"))
    print("[done] outputs in", cfg.out)


if __name__ == "__main__":
    main()
