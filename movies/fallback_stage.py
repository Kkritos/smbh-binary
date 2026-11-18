#!/usr/bin/env python3
"""fallback_stage.py -- physical debris FALLBACK RATE dM/dt onto the stellar BH and
the SMBH, measured directly from a long ballistic integration of the SAME debris the
movie uses (same disruption state from smbh_binary_PN, same debris_ic, same rng seed).

Why this exists: the movie's debris_follow.follow() batch-"settles" the inner disk to
stay cheap, and the old Dt_delay used the Keplerian energy->period relation.  In a
time-varying potential (a moving stellar BH near the SMBH) specific energy is NOT
conserved, so that inversion breaks down (Coughlin & Nixon 2019, ApJL 883 L17).  Here
each debris parcel's fallback time is the time of its FIRST return to pericenter about
its host.  We integrate SMBH+BH+all debris THROUGH the (brief) time-varying encounter --
capturing the real energy redistribution and the fast mini-disk returns directly -- and
then, once the BH has departed and the potential about each host is static, continue any
not-yet-returned parcel ANALYTICALLY along its now-fixed Kepler ellipse (exact, and cheap:
no thousands of static orbits are integrated).  Kepler is used only where it is valid, so
the late-time tail is resolved (full disruption -> t^-5/3; surviving-core / moving-potential
-> ~t^-9/4).  dM/dt is the mass histogram of the first-return times.

  python3 fallback_stage.py <stage1_bss_full.txt> <out_prefix> \
      [--Msmbh 1e6 --mbh 10 --mstar 1 --Rstar 1 --ndeb 10000 --seed 7] \
      [--disruptor bh|smbh] [--ntime 4000]   (--norbits/--norbits-bh accepted but unused)

Writes <out_prefix>.fallback.csv  (t[yr], dMdt_BH, dMdt_SMBH, dMdt_tot [Msun/yr])
   and <out_prefix>.fallback.png  (dM/dt curves + logarithmic-slope panel).
Reuses debris_follow.debris_ic / sample_polytrope_n3 so the debris are byte-identical
to the movie's; the MOVIE output is not touched.
"""
import argparse, numpy as np
from scipy.integrate import solve_ivp
import parse_traj
import debris_follow as df

G = 4.0*np.pi**2                      # AU, yr, Msun  (matches debris_follow)
C_AUYR = 63241.077
RSUN_AU = df.RSUN_AU
KMS_TO_AUYR = df.KMS_TO_AUYR


def _rhs(t, y, nbod, m_sm, m_bh, SS2, SB2):
    """[SMBH, BH, debris...]; SMBH<->BH exact 2-body, debris feel both (softened)."""
    P = y[:3*nbod].reshape(nbod, 3); V = y[3*nbod:].reshape(nbod, 3)
    A = np.zeros_like(P)
    d = P[1] - P[0]; r3 = (d @ d)**1.5
    A[0] = G*m_bh*d/r3; A[1] = -G*m_sm*d/r3
    q = P[2:]
    rs = P[0] - q; rb = P[1] - q
    rs3 = (np.einsum("ij,ij->i", rs, rs) + SS2)**1.5
    rb3 = (np.einsum("ij,ij->i", rb, rb) + SB2)**1.5
    A[2:] = G*m_sm*rs/rs3[:, None] + G*m_bh*rb/rb3[:, None]
    return np.concatenate([V.ravel(), A.ravel()])


def _first_return_time(logt, r):
    """First return-to-pericenter time (code units via exp(logt)) for one parcel:
    earliest local minimum of r reached after the parcel was farther out than it."""
    if r.size < 3:
        return np.nan
    dr = np.diff(r)
    loc_min = np.where((dr[:-1] < 0) & (dr[1:] >= 0))[0] + 1
    run_max = np.maximum.accumulate(r)
    for i in loc_min:
        if run_max[i-1] > r[i]*1.001:
            y0, y1, y2 = r[i-1], r[i], r[i+1]
            denom = (y0 - 2*y1 + y2)
            xmin = logt[i] + (0.5*(y0 - y2)/denom)*(logt[i+1]-logt[i]) if denom != 0 else logt[i]
            return np.exp(xmin)
    return np.nan


def _background(r_sm, v_sm, r_bh, v_bh, m_sm, m_bh, t_eval, rtol=1e-10, atol=1e-13):
    """Integrate the SMBH+BH two-body problem alone (cheap) and return their phase space
    at t_eval: P[2,3,nt], V[2,3,nt] with index 0=SMBH, 1=BH."""
    y0 = np.concatenate([r_sm, r_bh, v_sm, v_bh])
    def rhs(t, y):
        d = y[3:6] - y[0:3]; r3 = (d @ d)**1.5
        return np.concatenate([y[6:9], y[9:12], G*m_bh*d/r3, -G*m_sm*d/r3])
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), y0, t_eval=t_eval, method="DOP853",
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError("background integration failed: " + sol.message)
    return sol.y[:6].reshape(2, 3, -1), sol.y[6:].reshape(2, 3, -1)


def _kepler_t_to_peri(dr, dv, GM):
    """Time to the NEXT pericentre passage for a bound two-body state (dr, dv) about a
    fixed point mass GM, from analytic Kepler.  Valid once the potential is static (the
    perturbing BH has departed), so energy is conserved and the orbit no longer evolves.
    Returns code time (yr); nan if the parcel is unbound (no return)."""
    r = np.sqrt(dr @ dr); v2 = dv @ dv
    eps = 0.5*v2 - GM/r
    if eps >= -1e-30:
        return np.nan                                   # unbound/parabolic -> never returns
    a = -GM/(2.0*eps); n = np.sqrt(GM/a**3)             # semimajor axis, mean motion
    h = np.cross(dr, dv)
    e_vec = np.cross(dv, h)/GM - dr/r                    # Laplace-Runge-Lenz
    e = np.sqrt(e_vec @ e_vec)
    if e < 1e-8:
        return 2*np.pi/n                                # ~circular: pericentre ill-defined
    cosE = np.clip((1.0 - r/a)/e, -1.0, 1.0)
    E = np.arccos(cosE)
    if (dr @ dv) < 0.0:                                 # inbound -> past apocentre, E in (pi,2pi)
        E = 2*np.pi - E
    M_anom = E - e*np.sin(E)                            # mean anomaly since last pericentre
    return (2*np.pi - M_anom)/n                         # time forward to the next pericentre


def compute_fallback(r_sm, v_sm, r_bh, v_bh, dpos, dvel, M_smbh, m_bh, m_star, r_t,
                     disruptor="bh", norbits=400.0, norbits_bh=60.0, ntime=4000, nbins=42):
    """Return dM/dt onto the BH and the SMBH (Msun/yr) vs return time (yr), plus the
    total, the per-channel peak times, Delta t_delay, and the logarithmic slope.
    (norbits/norbits_bh are accepted for call-site compatibility but no longer used: the
    integration window is set self-consistently by the encounter time; see below.)"""
    # Cap the parcel count used for the fallback measurement: the encounter integration stores
    # a (6*N, ntime) array, so N=10000 x ntime=5000 would be multi-GB and thrash a login node.
    # A few thousand equal-mass tracers already resolve dM/dt and Delta t_delay; the movie keeps
    # the full N separately.
    NDEB_CAP = 4000
    N = dpos.shape[0]
    if N > NDEB_CAP:
        sidx = np.random.default_rng(0).choice(N, NDEB_CAP, replace=False)
        dpos, dvel = dpos[sidx], dvel[sidx]; N = NDEB_CAP
    m_parcel = m_star/N
    Rg_sm = G*M_smbh/C_AUYR**2; Rg_bh = G*m_bh/C_AUYR**2
    if disruptor == "smbh":
        SS, SB = 0.3*r_t, max(3.0*Rg_bh, 1e-7)
    else:
        SB, SS = 0.3*r_t, 3.0*Rg_sm

    # frozen-in energies at disruption -> host assignment (BH takes precedence) + window sizing
    rs = dpos - r_sm[None]; vs = dvel - v_sm[None]
    rb = dpos - r_bh[None]; vb = dvel - v_bh[None]
    eps_sm = 0.5*np.einsum("ij,ij->i", vs, vs) - G*M_smbh/np.linalg.norm(rs, axis=1)
    eps_bh = 0.5*np.einsum("ij,ij->i", vb, vb) - G*m_bh/np.linalg.norm(rb, axis=1)
    label = np.full(N, 2, int)
    if disruptor == "smbh":
        label[eps_sm < 0] = 1                      # all bound material falls onto the SMBH
    else:
        label[eps_sm < 0] = 1                      # 1 = SMBH-bound (central return)
        label[eps_bh < 0] = 0                      # 0 = BH-bound (mini-disk), takes precedence
    # --- fallback time per parcel = time of its FIRST return to pericentre about its host ---
    # The potential is time-varying only while the BH is near the SMBH (the "encounter",
    # ~ r_peri/v_rel).  We integrate SMBH+BH+ALL debris THROUGH that encounter (to t_dyn, a
    # few x its end), which captures the real energy redistribution and the fast mini-disk
    # returns directly.  Afterwards the BH has departed and the potential about each host is
    # ~static, so any parcel that has NOT yet returned is on a fixed Kepler ellipse and its
    # first-return time follows analytically (_kepler_t_to_peri) from its state at t_dyn --
    # no frozen-in Keplerian shortcut during the time-varying phase, and no need to integrate
    # thousands of static orbits.  This is both correct and bounded in cost.
    sep0 = float(np.linalg.norm(np.asarray(r_bh) - np.asarray(r_sm)))
    # Timescale for the separation to change = the local orbital time at sep0.  (Do NOT use
    # sep0/|v_rel|: for a near-marginal encounter |v_rel|->0 and that timescale diverges,
    # which balloons the integration window and hangs the solver.)  t_ref is always finite.
    t_ref = 2.0*np.pi*np.sqrt(sep0**3/(G*(M_smbh + m_bh)))
    tscan = np.logspace(np.log10(t_ref*1e-3), np.log10(t_ref*50.0), 400)
    Pbg, _ = _background(r_sm, v_sm, r_bh, v_bh, M_smbh, m_bh, tscan)
    sepf = np.linalg.norm(Pbg[1] - Pbg[0], axis=0)
    over = np.where(sepf >= 8.0*sep0)[0]               # BH has clearly departed -> potential static
    if over.size:
        t_static = tscan[over[0]]
    else:                                              # weak separation: fall back to 3x, else a few t_ref
        o3 = np.where(sepf >= 3.0*sep0)[0]
        t_static = tscan[o3[0]] if o3.size else 5.0*t_ref
    t_dyn = float(min(3.0*t_static, 20.0*t_ref))       # bounded window (hard ceiling at 20 t_ref)

    nbod = N + 2
    tg = np.logspace(np.log10(t_dyn*1e-5), np.log10(t_dyn), int(ntime)); tg[-1] = t_dyn
    y0 = np.concatenate([np.vstack([r_sm, r_bh, dpos]).ravel(),
                         np.vstack([v_sm, v_bh, dvel]).ravel()])
    sol = solve_ivp(_rhs, (0.0, t_dyn), y0, t_eval=tg, method="DOP853",
                    rtol=1e-9, atol=1e-12, args=(nbod, M_smbh, m_bh, SS*SS, SB*SB))
    if not sol.success:
        raise RuntimeError("encounter integration failed: " + sol.message)
    pos = sol.y[:3*nbod].reshape(nbod, 3, int(ntime)); vel = sol.y[3*nbod:].reshape(nbod, 3, int(ntime))
    logt = np.log(tg)
    rh_bh = np.linalg.norm(pos[2:] - pos[1][None], axis=1)     # (N,ntime) debris->BH
    rh_sm = np.linalg.norm(pos[2:] - pos[0][None], axis=1)     # (N,ntime) debris->SMBH

    tfb = np.full(N, np.nan)
    for k in range(N):
        if label[k] == 2:
            continue
        hb = 1 if label[k] == 0 else 0
        rh = rh_bh[k] if hb == 1 else rh_sm[k]
        tk = _first_return_time(logt, rh)              # a return WITHIN the integrated window?
        if np.isfinite(tk):
            tfb[k] = tk
        else:                                          # static Kepler continuation past t_dyn
            dr = pos[2+k, :, -1] - pos[hb, :, -1]
            dv = vel[2+k, :, -1] - vel[hb, :, -1]
            GM = G*(m_bh if hb == 1 else M_smbh)
            dtp = _kepler_t_to_peri(dr, dv, GM)
            if np.isfinite(dtp):
                tfb[k] = t_dyn + dtp

    fin = np.isfinite(tfb)
    if fin.sum() < 5:
        raise RuntimeError("too few returns; raise --ndeb/--ntime")

    # --- frozen-in Keplerian prediction, for the SAME parcels, to gauge the approximation ---
    # The old (reduce_delay-style) estimate assigns each parcel a fallback time from its
    # DISRUPTION-INSTANT energy about its host, assuming a fixed Kepler orbit about a static
    # point mass:  t_K = 2 pi G M_host (2|eps_host,0|)^{-3/2}.  We compute it here so the
    # figure/CSV can show measured-vs-frozen dM/dt and quantify the error per simulation.
    tK = np.full(N, np.nan)
    mb0 = (label == 0) & (eps_bh < 0)
    ms0 = (label == 1) & (eps_sm < 0)
    tK[mb0] = 2*np.pi*(G*m_bh)   * (2*np.abs(eps_bh[mb0]))**-1.5
    tK[ms0] = 2*np.pi*(G*M_smbh) * (2*np.abs(eps_sm[ms0]))**-1.5
    finK = np.isfinite(tK)

    # common log-time bins spanning BOTH the measured and frozen-in times
    allt = np.concatenate([tfb[fin], tK[finK]])
    lo, hi = np.nanpercentile(allt, [0.3, 99.7])
    bins = np.logspace(np.log10(lo), np.log10(hi), nbins+1)
    cen = np.sqrt(bins[:-1]*bins[1:]); dtb = np.diff(bins)
    mdot  = lambda s: m_parcel*np.histogram(tfb[s & fin],  bins=bins)[0]/dtb   # measured
    mdotK = lambda s: m_parcel*np.histogram(tK[s & finK],  bins=bins)[0]/dtb   # frozen-in Kepler
    md_bh, md_sm, md_tot = mdot(label == 0), mdot(label == 1), mdot(np.ones(N, bool))
    mK_bh, mK_sm, mK_tot = mdotK(label == 0), mdotK(label == 1), mdotK(np.ones(N, bool))
    pk = lambda md: (cen[np.argmax(md)] if (md > 0).any() else np.nan)
    tpk_bh, tpk_sm = pk(md_bh), pk(md_sm)
    tpkK_bh, tpkK_sm = pk(mK_bh), pk(mK_sm)
    with np.errstate(divide="ignore", invalid="ignore"):
        nsl = np.gradient(np.log(np.where(md_tot > 0, md_tot, np.nan)), np.log(cen))

    # per-parcel accuracy of the frozen-in approximation: rho = t_measured / t_Kepler,
    # reported SEPARATELY per channel -- the BH mini-disk is softening-sensitive and drives
    # rho->0, whereas the SMBH (central / second-flare) channel is the physically robust one.
    def _ratio(sel):
        m = sel & fin & finK
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = tfb[m]/tK[m]
        rr = rr[np.isfinite(rr) & (rr > 0)]
        return tuple(np.percentile(rr, [50, 16, 84])) if rr.size else (np.nan, np.nan, np.nan)
    rBH  = _ratio(label == 0)
    rSM  = _ratio(label == 1)
    rAll = _ratio(label < 2)
    dd   = (tpk_sm  - tpk_bh)  if np.isfinite(tpk_bh)  and np.isfinite(tpk_sm)  else np.nan
    ddK  = (tpkK_sm - tpkK_bh) if np.isfinite(tpkK_bh) and np.isfinite(tpkK_sm) else np.nan

    return dict(t=cen, dMdt_BH=md_bh, dMdt_SMBH=md_sm, dMdt_tot=md_tot, n_slope=nsl,
                dMdt_BH_K=mK_bh, dMdt_SMBH_K=mK_sm, dMdt_tot_K=mK_tot,            # frozen-in dM/dt
                t_peak_BH=tpk_bh, t_peak_SMBH=tpk_sm, Dt_delay=dd,
                t_peak_BH_K=tpkK_bh, t_peak_SMBH_K=tpkK_sm, Dt_delay_K=ddK,       # frozen-in peaks/delay
                ratio_BH=rBH, ratio_SMBH=rSM, ratio_all=rAll,                     # t_meas/t_Kepler per channel
                frac_BH=np.mean(label == 0), frac_SMBH=np.mean(label == 1),
                frac_unbound=np.mean(label == 2), frac_censored=float(np.mean(~fin & (label < 2))),
                disruptor=disruptor, N=N)


def energy_evolution(r_sm, v_sm, r_bh, v_bh, dpos, dvel, M_smbh, m_bh, r_t,
                     disruptor="bh", tmax=1.0, nsnap=6, ndeb_max=3000, seed=0,
                     rtol=1e-9, atol=1e-12):
    """Integrate ALL debris together and return the specific-energy distributions wrt
    BH1 (eps1) and the SMBH/BH2 (eps2) at `nsnap` snapshot times in [0, tmax].  Because
    the potential is time-varying, these distributions evolve; we sample them so the
    morphing is visible.  (Debris are subsampled to ndeb_max for a representative PDF.)"""
    N = dpos.shape[0]
    if N > ndeb_max:
        s = np.random.default_rng(seed).choice(N, ndeb_max, replace=False)
        dpos, dvel = dpos[s], dvel[s]; N = ndeb_max
    Rg_sm = G*M_smbh/C_AUYR**2; Rg_bh = G*m_bh/C_AUYR**2
    SS, SB = (0.3*r_t, max(3.0*Rg_bh, 1e-7)) if disruptor == "smbh" else (3.0*Rg_sm, 0.3*r_t)
    snaps = np.concatenate([[0.0], np.logspace(np.log10(tmax*1e-3), np.log10(tmax), nsnap-1)])
    snaps[-1] = tmax

    # The specific energies evolve ONLY while the potential is time-varying, i.e. while the
    # BH and SMBH are close.  Once their separation has grown well beyond its initial
    # (pericentre) value the potential is ~static and each parcel's energy wrt its host is
    # conserved -- the distributions freeze.  So we integrate the debris dynamics only
    # through that encounter epoch (t_dyn) and hold the distributions fixed afterwards;
    # this keeps the cost bounded no matter how large tmax is (the fast mini-disk parcels
    # would otherwise need thousands of orbits to reach a late tmax).
    sep0 = np.linalg.norm(np.asarray(r_bh) - np.asarray(r_sm))
    tf = np.logspace(np.log10(max(tmax, 1e-6)*1e-4), np.log10(max(tmax, 1e-6)), 400)
    Pf, _ = _background(r_sm, v_sm, r_bh, v_bh, M_smbh, m_bh, tf)
    sepf = np.linalg.norm(Pf[1] - Pf[0], axis=0)
    over = np.where(sepf >= 10.0*sep0)[0]
    t_static = tf[over[0]] if over.size else tmax
    t_dyn = float(min(tmax, max(2.0*t_static, snaps[1])))

    teval = np.unique(np.append(snaps[snaps <= t_dyn], t_dyn))
    nbod = N + 2
    y0 = np.concatenate([np.vstack([r_sm, r_bh, dpos]).ravel(),
                         np.vstack([v_sm, v_bh, dvel]).ravel()])
    sol = solve_ivp(_rhs, (0.0, t_dyn), y0, t_eval=teval, method="DOP853",
                    rtol=rtol, atol=atol, args=(nbod, M_smbh, m_bh, SS*SS, SB*SB))
    if not sol.success:
        raise RuntimeError("energy-evolution integration failed: " + sol.message)
    pos = sol.y[:3*nbod].reshape(nbod, 3, len(teval)); vel = sol.y[3*nbod:].reshape(nbod, 3, len(teval))

    def _energies(k):
        rd, vd = pos[2:, :, k], vel[2:, :, k]
        e1 = 0.5*np.sum((vd - vel[1, :, k])**2, 1) - G*m_bh/np.linalg.norm(rd - pos[1, :, k], axis=1)
        e2 = 0.5*np.sum((vd - vel[0, :, k])**2, 1) - G*M_smbh/np.linalg.norm(rd - pos[0, :, k], axis=1)
        return e1, e2
    E1, E2 = [], []
    for ts in snaps:
        k = int(np.argmin(np.abs(teval - min(ts, t_dyn))))   # late snaps reuse the frozen t_dyn state
        e1, e2 = _energies(k); E1.append(e1); E2.append(e2)
    return snaps, E1, E2


def plot_energy_evolution(snaps, E1, E2, path, title="", nbins=60):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10", "DejaVu Serif"],
                         "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
                         "axes.unicode_minus": False, "xtick.direction": "in", "ytick.direction": "in"})
    allE = np.concatenate([np.abs(np.concatenate(E1)), np.abs(np.concatenate(E2))])
    allE = allE[np.isfinite(allE) & (allE > 0)]
    lt = np.percentile(allE, 2); amax = np.percentile(allE, 99.9)        # symlog linthresh + range
    pos = np.logspace(np.log10(lt), np.log10(amax), nbins//2)
    bins = np.concatenate([-pos[::-1], [0.0], pos])
    ns = len(snaps); nc = min(3, ns); nr = int(np.ceil(ns/nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.2*nc, 3.1*nr), squeeze=False)
    for k in range(ns):
        ax = axs.flat[k]
        ax.hist(E1[k], bins=bins, histtype="step", color="#0072B2", lw=1.6,
                density=True, label=r"$\varepsilon_1$ (wrt BH1)")
        ax.hist(E2[k], bins=bins, histtype="step", color="#D55E00", lw=1.6,
                density=True, label=r"$\varepsilon_2$ (wrt SMBH)")
        ax.axvline(0.0, color="0.6", lw=0.8, ls="--")
        ax.set_xscale("symlog", linthresh=lt); ax.set_yscale("log")
        ax.text(0.04, 0.95, rf"$t={snaps[k]:.3g}$ yr", transform=ax.transAxes,
                ha="left", va="top", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=0.9))
        ax.set_xlabel(r"specific energy  [AU$^2$ yr$^{-2}$]")
        if k % nc == 0: ax.set_ylabel(r"$dN/d\varepsilon$ (pdf)")
        if k == 0: ax.legend(frameon=False, fontsize=8, loc="upper right")
    for k in range(ns, nr*nc): axs.flat[k].axis("off")
    if title: fig.suptitle(title, y=1.0)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    return path


def plot_fallback(res, path, title=""):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["cmr10", "DejaVu Serif"],
                         "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
                         "axes.unicode_minus": False,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    c = res["t"]; fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    # draw the total FIRST (thin, pale, behind) so the two channel curves sit on top of it
    if res["disruptor"] != "smbh":
        gt = res["dMdt_tot"] > 0
        ax[0].plot(c[gt], res["dMdt_tot"][gt], "-", lw=3.2, color="0.78", label="total", zorder=1)
    series = [(res["dMdt_SMBH"], "onto SMBH (central)", "#D55E00")]
    if res["disruptor"] != "smbh":
        series = [(res["dMdt_BH"], "onto BH (companion)", "#0072B2")] + series
    for md, lab, col in series:
        g = md > 0
        ax[0].plot(c[g], md[g], "o-", ms=3, lw=1.5, color=col, label=lab, zorder=3)
    for tp, col in [(res["t_peak_BH"], "#0072B2"), (res["t_peak_SMBH"], "#D55E00")]:
        if np.isfinite(tp): ax[0].axvline(tp, color=col, ls=":", lw=1, zorder=2)
    # reference power laws on the dM/dt panel, anchored at the central (SMBH) peak so the
    # decline is compared against the late-time tail; span whatever range follows the peak.
    mdr = res["dMdt_SMBH"]; gr = mdr > 0
    if gr.sum() > 3 and np.isfinite(res["t_peak_SMBH"]):
        t_anch = res["t_peak_SMBH"]; ia = int(np.argmin(np.abs(c - t_anch)))
        y_anch = mdr[ia] if mdr[ia] > 0 else mdr[gr][0]
        tt = c[c >= t_anch]
        for p, ls, lab in [(-5.0/3.0, "--", r"$t^{-5/3}$"), (-9.0/4.0, (0, (1, 1)), r"$t^{-9/4}$")]:
            if tt.size > 2:
                yy = y_anch*(tt/t_anch)**p
                ax[0].plot(tt, yy, ls=ls, color="k", lw=1.1, alpha=0.8, zorder=4)
                ax[0].text(tt[-1], yy[-1], " "+lab, fontsize=8, va="center")
    # frozen-in Keplerian prediction (reduce_delay-style) overlaid as thin DASHED curves in
    # the same channel colours, so the measured-vs-approximation difference is read directly.
    kseries = [(res.get("dMdt_SMBH_K"), "#D55E00")]
    if res["disruptor"] != "smbh":
        kseries = [(res.get("dMdt_BH_K"), "#0072B2")] + kseries
    for mk, col in kseries:
        if mk is None:
            continue
        g = mk > 0
        ax[0].plot(c[g], mk[g], ls="--", lw=1.3, color=col, alpha=0.9, zorder=2)
    for tp, col in [(res.get("t_peak_BH_K"), "#0072B2"), (res.get("t_peak_SMBH_K"), "#D55E00")]:
        if tp is not None and np.isfinite(tp): ax[0].axvline(tp, color=col, ls=(0, (1, 1)), lw=1, alpha=0.7, zorder=1)
    # legend handles: solid = measured, dashed = frozen-in Keplerian
    from matplotlib.lines import Line2D
    extra = [Line2D([0], [0], color="0.3", ls="-", lw=1.5, label="measured (direct)"),
             Line2D([0], [0], color="0.3", ls="--", lw=1.3, label="frozen-in Kepler")]
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("return time [yr]"); ax[0].set_ylabel(r"$dM/dt\;[M_\odot\,{\rm yr}^{-1}]$")
    h0, l0 = ax[0].get_legend_handles_labels()
    ax[0].legend(h0 + extra, l0 + [e.get_label() for e in extra], frameon=False, fontsize=8, ncol=1)
    dd, ddK = res["Dt_delay"], res.get("Dt_delay_K", np.nan)
    ttl = title
    if np.isfinite(dd):
        ttl += rf"  $\Delta t_{{\rm delay}}$: {dd:.3g} yr (direct)"
        if np.isfinite(ddK): ttl += rf" vs {ddK:.3g} (frozen)"
    ax[0].set_title(ttl, fontsize=10)
    # annotate the per-parcel accuracy of the frozen-in approximation (per channel)
    def _rtxt(lab, tup, col):
        m, lo16, hi84 = tup
        if not np.isfinite(m):
            return None
        return (col, rf"{lab}: $t_{{\rm meas}}/t_{{\rm Kep}}$ = {m:.2f}$^{{+{hi84-m:.2f}}}_{{-{m-lo16:.2f}}}$")
    lines = [_rtxt("SMBH", res.get("ratio_SMBH", (np.nan,)*3), "#D55E00")]
    if res["disruptor"] != "smbh":
        lines = [_rtxt("BH", res.get("ratio_BH", (np.nan,)*3), "#0072B2")] + lines
    yoff = 0.03
    for item in lines:
        if item is None:
            continue
        col, txt = item
        ax[0].text(0.03, yoff, txt, transform=ax[0].transAxes, ha="left", va="bottom",
                   fontsize=8, color=col,
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.8", alpha=0.9))
        yoff += 0.075
    g = np.isfinite(res["n_slope"])
    ax[1].plot(c[g], res["n_slope"][g], "k-")
    for y, lab in [(-5/3, r"$-5/3$"), (-9/4, r"$-9/4$")]:
        ax[1].axhline(y, ls="--", lw=1, color="0.5"); ax[1].text(c[g][-1], y, " "+lab, va="center", fontsize=9)
    ax[1].set_xscale("log"); ax[1].set_xlabel("return time [yr]")
    ax[1].set_ylabel(r"$n=d\ln\dot M/d\ln t$"); ax[1].set_ylim(-3.2, 0.6)
    fig.tight_layout(); fig.savefig(path, dpi=140, bbox_inches="tight"); plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw"); ap.add_argument("out_prefix")
    ap.add_argument("--Msmbh", type=float, default=1e6); ap.add_argument("--mbh", type=float, default=10.0)
    ap.add_argument("--mstar", type=float, default=1.0); ap.add_argument("--Rstar", type=float, default=1.0)
    ap.add_argument("--ndeb", type=int, default=10000); ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--disruptor", choices=["bh", "smbh"], default="bh")
    ap.add_argument("--norbits", type=float, default=400.0, help="(accepted for compatibility; unused)")
    ap.add_argument("--norbits-bh", type=float, default=60.0, help="(accepted for compatibility; unused)")
    ap.add_argument("--ntime", type=int, default=4000, help="time samples across the encounter integration")
    ap.add_argument("--nsnap", type=int, default=6, help="energy-distribution snapshots in time")
    ap.add_argument("--snap-tmax", type=float, default=0.0, dest="snap_tmax",
                    help="max time for energy snapshots [yr]; 0 = auto (2x central-fallback peak)")
    ap.add_argument("--snap-ndeb", type=int, default=3000, dest="snap_ndeb")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    tr = parse_traj.parse(a.raw)
    pos = tr["pos"]; vel = tr["vel"]*KMS_TO_AUYR
    r_sm, r_bh, r_st = pos[-1, 0], pos[-1, 1], pos[-1, 2]
    v_sm, v_bh, v_st = vel[-1, 0], vel[-1, 1], vel[-1, 2]
    R_star = a.Rstar*RSUN_AU
    if a.disruptor == "smbh":
        r_t = R_star*(a.Msmbh/a.mstar)**(1.0/3.0); r_disr, v_disr = r_sm, v_sm
    else:
        r_t = R_star*(a.mbh/a.mstar)**(1.0/3.0); r_disr, v_disr = r_bh, v_bh
    dpos, dvel = df.debris_ic(r_st, v_st, r_disr, v_disr, a.mstar, R_star, a.ndeb, rng)
    res = compute_fallback(r_sm, v_sm, r_bh, v_bh, dpos, dvel, a.Msmbh, a.mbh, a.mstar, r_t,
                           disruptor=a.disruptor, norbits=a.norbits, norbits_bh=a.norbits_bh, ntime=a.ntime)
    hdr = "t_yr,dMdt_BH,dMdt_SMBH,dMdt_tot,dMdt_BH_Kfrozen,dMdt_SMBH_Kfrozen,dMdt_tot_Kfrozen"
    np.savetxt(a.out_prefix + ".fallback.csv",
               np.column_stack([res["t"], res["dMdt_BH"], res["dMdt_SMBH"], res["dMdt_tot"],
                                res["dMdt_BH_K"], res["dMdt_SMBH_K"], res["dMdt_tot_K"]]),
               delimiter=",", header=hdr, comments="", fmt="%.6g")
    plot_fallback(res, a.out_prefix + ".fallback.png", title=f"{a.disruptor}-TDE")
    print(f"[{a.disruptor}] fractions BH={res['frac_BH']:.2f} SMBH={res['frac_SMBH']:.2f} "
          f"unb={res['frac_unbound']:.2f} censored={res['frac_censored']:.2f}")
    print(f"  t_peak_BH={res['t_peak_BH']:.3e} yr  t_peak_SMBH={res['t_peak_SMBH']:.3e} yr  "
          f"Dt_delay={res['Dt_delay']:.3e} yr  (direct)")
    print(f"  t_peak_BH={res['t_peak_BH_K']:.3e} yr  t_peak_SMBH={res['t_peak_SMBH_K']:.3e} yr  "
          f"Dt_delay={res['Dt_delay_K']:.3e} yr  (frozen-in Kepler)")
    rB, rS = res["ratio_BH"], res["ratio_SMBH"]
    print(f"  accuracy t_meas/t_Kepler  BH={rB[0]:.3f} [16-84%: {rB[1]:.3f},{rB[2]:.3f}]  "
          f"SMBH={rS[0]:.3f} [16-84%: {rS[1]:.3f},{rS[2]:.3f}]")
    print(f"  wrote {a.out_prefix}.fallback.csv (now incl. frozen-in columns) and .fallback.png")

    # ---- scalar comparison: direct (measured) vs frozen-in Kepler, one row, machine-readable ----
    rA = res["ratio_all"]
    chdr = ("disruptor,t_peak_BH_direct,t_peak_BH_frozen,t_peak_SMBH_direct,t_peak_SMBH_frozen,"
            "Dt_delay_direct,Dt_delay_frozen,Dt_delay_ratio,"
            "ratio_BH_med,ratio_BH_16,ratio_BH_84,ratio_SMBH_med,ratio_SMBH_16,ratio_SMBH_84,"
            "ratio_all_med,ratio_all_16,ratio_all_84,frac_BH,frac_SMBH,frac_unbound")
    dratio = (res["Dt_delay"]/res["Dt_delay_K"]) if np.isfinite(res["Dt_delay"]) and \
             np.isfinite(res["Dt_delay_K"]) and res["Dt_delay_K"] != 0 else float("nan")
    with open(a.out_prefix + ".compare.csv", "w") as f:
        f.write(chdr + "\n")
        f.write(f"{a.disruptor},{res['t_peak_BH']:.6g},{res['t_peak_BH_K']:.6g},"
                f"{res['t_peak_SMBH']:.6g},{res['t_peak_SMBH_K']:.6g},"
                f"{res['Dt_delay']:.6g},{res['Dt_delay_K']:.6g},{dratio:.6g},"
                f"{rB[0]:.6g},{rB[1]:.6g},{rB[2]:.6g},{rS[0]:.6g},{rS[1]:.6g},{rS[2]:.6g},"
                f"{rA[0]:.6g},{rA[1]:.6g},{rA[2]:.6g},"
                f"{res['frac_BH']:.4g},{res['frac_SMBH']:.4g},{res['frac_unbound']:.4g}\n")
    print(f"  wrote {a.out_prefix}.compare.csv  (direct vs frozen-in Dt_delay = "
          f"{res['Dt_delay']:.3g} vs {res['Dt_delay_K']:.3g} yr)")

    # ---- debris specific-energy distributions (eps1, eps2) vs time (time-varying potential) ----
    tmax = a.snap_tmax if a.snap_tmax > 0 else (2.0*res["t_peak_SMBH"] if np.isfinite(res["t_peak_SMBH"]) else 1.0)
    try:
        snaps, E1, E2 = energy_evolution(r_sm, v_sm, r_bh, v_bh, dpos, dvel, a.Msmbh, a.mbh, r_t,
                                         disruptor=a.disruptor, tmax=tmax, nsnap=a.nsnap, ndeb_max=a.snap_ndeb)
        plot_energy_evolution(snaps, E1, E2, a.out_prefix + ".energy.png",
                              title=f"{a.disruptor}-TDE: debris specific-energy distributions vs time")
        np.savez_compressed(a.out_prefix + ".energy.npz", snaps=snaps,
                            eps1=np.array(E1), eps2=np.array(E2))
        print(f"  wrote {a.out_prefix}.energy.png and .energy.npz ({a.nsnap} snapshots to t={tmax:.3g} yr)")
    except Exception as e:
        print(f"  energy-evolution skipped: {e}")


if __name__ == "__main__":
    main()
