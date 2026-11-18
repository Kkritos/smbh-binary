#!/usr/bin/env python3
"""TDE debris post-step driven by the REAL smbh_binary_PN disruption state.

   python3 debris_follow.py stage1_bss_full.txt scenes/mutde.csv \
           [--Msmbh 1e6 --mbh 10 --mstar 1 --Rstar 1 --ndeb 10000 --norbits 10]
           [--disruptor bh|smbh]

Two disruption channels, same ballistic machinery:

  --disruptor bh   (micro-TDE, TDE:star2_by_bh1)
      The star is torn apart by its *stellar-BH companion*.  Tidal radius
      r_t = R* (m_bh/m*)^(1/3) ~ 0.01 AU.  Debris orbit/escape the SMBH; the bound/
      ejected split is set mostly by where the BH sits in its SMBH orbit at disruption.
      Main panel = SMBH-centred fallback; inset = the receding BH + its mini-disk.

  --disruptor smbh (SMBH TDE, TDE:star2_by_SMBH)
      The star is torn apart by the *SMBH itself*.  Tidal radius
      r_t = R* (M_smbh/m*)^(1/3) ~ 0.46 AU (~46x larger).  The disruptor and the body
      the debris orbit are the SAME object, so the split is the textbook frozen-in
      spread dE ~ G M_smbh R* / r_t^2 -> a near half-bound / half-ejected stream (the
      canonical t^-5/3 fallback).  Main panel = SMBH-centred stream; inset = the inner
      launch region around the SMBH.  The stellar BH is a spectator point.

In both cases the executable classifies the event but cannot emit ~10^4 particles, so
the disruption + fallback is a ballistic post-step seeded from the executable's EXACT
final 3-body state:
  * star modelled as an n=3 polytrope; `ndeb` tracers drawn from its density profile.
  * each tracer velocity = star bulk velocity  +  solid-body SPIN-UP to breakup,
    Omega_break = sqrt(G m*/R*^3), about the orbital angular momentum wrt the DISRUPTOR
    (Evans & Kochanek 1989).  The frozen-in energy spread arises self-consistently.
  * tracers move as test particles in the time-dependent BH + SMBH field; coloured by
    specific energy wrt the SMBH (blue bound, red ejected).

Outputs:
  <out>.csv          canonical SMBH/BH/star trajectory (approach + fallback)
  <out>.debris.npz   debris positions over the fallback frames + per-tracer energy
  <out>.events.json  render config (SMBH-centred main + inset)
"""
import argparse, json, csv, numpy as np
import parse_traj

G = 4.0 * np.pi**2                       # AU, yr, Msun
C_AUYR = 63241.077                       # speed of light [AU/yr]
RSUN_AU = 695700.0 / 1.495978707e8       # 0.004649 AU
KMS_TO_AUYR = 1.0 / 4.740572             # 1 km/s in AU/yr


# ---- n=3 polytrope (Lane-Emden) + sampling ----------------------------------
def lane_emden_n3(h=2e-3, xi_max=10.0):
    """Integrate theta'' + (2/xi)theta' + theta^3 = 0, theta(0)=1, theta'(0)=0 (RK4).
    Returns xi grid and theta up to the first zero (the polytrope surface)."""
    def deriv(xi, th, dth):
        return dth, -(2.0 / xi) * dth - max(th, 0.0) ** 3
    xi, th, dth = 1e-6, 1.0, 0.0
    X, TH = [xi], [th]
    while th > 0.0 and xi < xi_max:
        k1 = deriv(xi, th, dth)
        k2 = deriv(xi + h/2, th + h/2*k1[0], dth + h/2*k1[1])
        k3 = deriv(xi + h/2, th + h/2*k2[0], dth + h/2*k2[1])
        k4 = deriv(xi + h,   th + h*k3[0],   dth + h*k3[1])
        th  += h/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        dth += h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        xi  += h
        X.append(xi); TH.append(max(th, 0.0))
    return np.array(X), np.array(TH)


def sample_polytrope_n3(N, R_star, rng):
    """Draw N positions (N,3) from an n=3 polytrope of radius R_star."""
    X, TH = lane_emden_n3()
    xi1 = X[-1]
    w = X**2 * TH**3                      # mass pdf in xi: 4 pi xi^2 rho, rho ~ theta^3
    cdf = np.cumsum(w); cdf /= cdf[-1]
    r = R_star * np.interp(rng.random(N), cdf, X) / xi1
    cost = rng.uniform(-1, 1, N); phi = rng.uniform(0, 2*np.pi, N)
    sint = np.sqrt(np.clip(1 - cost*cost, 0, 1))
    dirs = np.column_stack([sint*np.cos(phi), sint*np.sin(phi), cost])
    return dirs * r[:, None]


def debris_ic(r_star, v_star, r_disr, v_disr, m_star, R_star, N, rng):
    """Positions/velocities (N,3) of debris at disruption: polytrope positions + star
    bulk velocity + solid-body spin-up to breakup about the orbital angular momentum of
    the star wrt the DISRUPTOR (the BH for a micro-TDE, the SMBH for an SMBH TDE)."""
    dr = sample_polytrope_n3(N, R_star, rng)                 # (N,3) about the star
    pos = r_star[None, :] + dr
    Lhat = np.cross(r_star - r_disr, v_star - v_disr)
    Lhat = Lhat / (np.linalg.norm(Lhat) + 1e-30)
    omega_break = np.sqrt(G * m_star / R_star**3)            # breakup ang. velocity
    v_spin = omega_break * np.cross(np.broadcast_to(Lhat, dr.shape), dr)
    vel = v_star[None, :] + v_spin
    return pos, vel


# ---- ballistic fallback integration (SMBH+BH 2-body + debris test particles) --
def follow(r_sm, v_sm, r_bh, v_bh, dpos, dvel, M_smbh, m_bh, r_t,
           t_follow, n_frames=320, dt=5e-5, disruptor="bh", settle_orbits=40.0):
    """Integrate the SMBH+BH 2-body exactly and the debris as test particles in their
    combined, time-dependent field (KDK leapfrog).  Softening + accretion-freezing are
    keyed to the DISRUPTOR: the disrupted-star debris form a bound disk around whichever
    body did the disrupting, and a tracer that leaves that body's tidal sphere and later
    falls back deep is 'accreted' (co-moves) so stiff bound orbits stay cheap."""
    Rg_sm = G * M_smbh / C_AUYR**2
    Rg_bh = G * m_bh   / C_AUYR**2
    if disruptor == "smbh":
        SS = 0.3 * r_t                 # soft disk around the SMBH (the disruptor)
        SB = max(3.0 * Rg_bh, 1e-7)    # the stellar BH is nearly a point here
        R_ACC_SM, R_ACC_BH = 1.2 * r_t, 0.25
        gate_sm = True                 # SMBH accretion gated by 'has_left'
    else:
        SB = 0.3 * r_t                 # soft mini-disk around the BH (the disruptor)
        SS = 3.0 * Rg_sm
        R_ACC_BH, R_ACC_SM = 1.2 * r_t, 0.25
        gate_sm = False                # BH accretion gated by 'has_left'
    N = len(dpos)
    p_sm, v_sm = r_sm.copy(), v_sm.copy()
    p_bh, v_bh = r_bh.copy(), v_bh.copy()
    dp, dv = dpos.copy(), dvel.copy()
    frozen = np.zeros(N, bool); acc_of = np.zeros(N, int)
    off = np.zeros((N, 3)); has_left = np.zeros(N, bool)

    def acc_massive(ps, pb):
        d = pb - ps; r3 = (d @ d)**1.5
        return G*m_bh*d/r3, -G*M_smbh*d/r3
    def acc_deb(dq, ps, pb):
        rs = ps - dq; rb = pb - dq
        rs3 = (np.einsum('ij,ij->i', rs, rs) + SS*SS)**1.5
        rb3 = (np.einsum('ij,ij->i', rb, rb) + SB*SB)**1.5
        return G*M_smbh*rs/rs3[:, None] + G*m_bh*rb/rb3[:, None]

    # Two-phase timestep so a LONG follow (debris returning to the SMBH) stays cheap:
    #   fine dt while the disruptor's tightly-bound debris form/return (the mini-disk),
    #   then coarse dt for the long, wide SMBH fallback.  At the transition the bound
    #   inner disk is 'settled' (frozen, co-moving) so the coarse step can't fling it.
    m_disr = m_bh if disruptor == "bh" else M_smbh
    P_disr = 2*np.pi*np.sqrt(r_t**3 / (G*m_disr))
    t_fine = min(t_follow, settle_orbits*P_disr)   # fine-step until the mini-disk forms; then settle+hold
    dt_fine = dt
    NCOARSE = 20000
    dt_coarse = max(dt_fine, (t_follow - t_fine)/NCOARSE) if t_follow > t_fine else dt_fine
    # frame times: mild early bias so the disruption shows; the rest over the fallback
    frame_times = t_follow * (np.linspace(0.0, 1.0, n_frames)**1.3)

    SM, BH, DEB, TT = [], [], [], []
    def _record(tc):
        SM.append(p_sm.copy()); BH.append(p_bh.copy()); DEB.append(dp.copy()); TT.append(tc)

    a_sm, a_bh = acc_massive(p_sm, p_bh); a_d = acc_deb(dp, p_sm, p_bh)
    t = 0.0; fi = 0; settled = False
    max_steps = int(t_fine/dt_fine) + NCOARSE + 16
    for _ in range(max_steps):
        while fi < n_frames and t >= frame_times[fi] - 1e-15:
            _record(t); fi += 1
        if fi >= n_frames:
            break
        if (not settled) and t >= t_fine:        # settle the inner (disk) debris once
            settled = True
            pd = p_bh if disruptor == "bh" else p_sm
            vd = v_bh if disruptor == "bh" else v_sm
            rb = dp - pd; vb = dv - vd
            eps_b = 0.5*np.einsum('ij,ij->i', vb, vb) - G*m_disr/np.linalg.norm(rb, axis=1)
            disk = (~frozen) & (eps_b < 0) & (np.linalg.norm(rb, axis=1) < 5.0*r_t)
            off[disk] = dp[disk] - pd; acc_of[disk] = (0 if disruptor == "bh" else 1); frozen |= disk
        dtn = dt_fine if t < t_fine else dt_coarse
        fr = ~frozen
        v_sm += 0.5*dtn*a_sm; v_bh += 0.5*dtn*a_bh; dv[fr] += 0.5*dtn*a_d[fr]
        p_sm += dtn*v_sm; p_bh += dtn*v_bh; dp[fr] += dtn*dv[fr]
        fb = frozen & (acc_of == 0); fs = frozen & (acc_of == 1)
        dp[fb] = p_bh + off[fb]; dp[fs] = p_sm + off[fs]
        a_sm, a_bh = acc_massive(p_sm, p_bh); a_d = acc_deb(dp, p_sm, p_bh)
        v_sm += 0.5*dtn*a_sm; v_bh += 0.5*dtn*a_bh; dv[fr] += 0.5*dtn*a_d[fr]
        dbh = np.linalg.norm(dp - p_bh[None, :], axis=1)
        dsm = np.linalg.norm(dp - p_sm[None, :], axis=1)
        has_left |= (dsm > 2.5*r_t) if gate_sm else (dbh > 2.5*r_t)
        gbh = has_left if not gate_sm else np.ones(N, bool)   # BH accretion gate
        gsm = has_left if gate_sm else np.ones(N, bool)       # SMBH accretion gate
        nb = (~frozen) & gbh & (dbh < R_ACC_BH)
        ns = (~frozen) & gsm & (dsm < R_ACC_SM)
        off[nb] = dp[nb] - p_bh; acc_of[nb] = 0; frozen |= nb
        off[ns] = dp[ns] - p_sm; acc_of[ns] = 1; frozen |= ns
        t += dtn
    while fi < n_frames:                        # pad any remaining frames
        _record(t); fi += 1

    SM = np.array(SM); BH = np.array(BH); DEB = np.array(DEB); TT = np.array(TT)
    # per-tracer specific energy wrt the SMBH (for colouring: <0 bound, >0 ejected)
    v_now = dv.copy()
    v_now[frozen & (acc_of == 0)] = v_bh
    v_now[frozen & (acc_of == 1)] = v_sm
    rrel = dp - p_sm; vrel = v_now - v_sm
    eps = 0.5*np.einsum('ij,ij->i', vrel, vrel) - G*M_smbh/np.linalg.norm(rrel, axis=1)
    # specific energy wrt the BH -> which debris are bound to / returning onto the stellar BH
    rrel_b = dp - p_bh; vrel_b = v_now - v_bh
    eps_bh = 0.5*np.einsum('ij,ij->i', vrel_b, vrel_b) - G*m_bh/np.linalg.norm(rrel_b, axis=1)
    n_bound_bh = int((eps_bh < 0).sum())
    n_acc_bh   = int((frozen & (acc_of == 0)).sum())      # debris that fell back onto the BH
    n_acc_sm   = int((frozen & (acc_of == 1)).sum())      # debris that fell back onto the SMBH
    return dict(sm=SM, bh=BH, deb=DEB, t=TT, eps=eps, frozen=frozen, acc_of=acc_of,
                n_bound_bh=n_bound_bh, n_acc_bh=n_acc_bh, n_acc_sm=n_acc_sm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw"); ap.add_argument("out_csv")
    ap.add_argument("--Msmbh", type=float, default=1e6)
    ap.add_argument("--mbh", type=float, default=10.0)
    ap.add_argument("--mstar", type=float, default=1.0)
    ap.add_argument("--Rstar", type=float, default=1.0, help="stellar radius [Rsun]")
    ap.add_argument("--ndeb", type=int, default=10000)
    ap.add_argument("--norbits", type=float, default=10.0, help="SMBH/fallback orbits to follow")
    ap.add_argument("--nframes", type=int, default=420, help="fallback frames (smoothness)")
    ap.add_argument("--max1", type=int, default=240, help="approach frames (uniform time)")
    ap.add_argument("--disruptor", choices=["bh", "smbh"], default="bh",
                    help="bh = micro-TDE (companion BH); smbh = star disrupted by the SMBH")
    ap.add_argument("--settle-orbits", type=float, default=40.0, dest="settle_orbits",
                    help="freeze+hold the BH1 mini-disk in the inset after this many disruptor orbits")
    a = ap.parse_args()
    rng = np.random.default_rng(7)

    # ---- stage-1 (approach to disruption) from the real trajectory, in 3D --------
    tr = parse_traj.parse(a.raw)
    t_raw = tr["t"]; pos_raw = tr["pos"]; vel_raw = tr["vel"] * KMS_TO_AUYR    # AU, AU/yr
    # disruption state = exact final snapshot
    r_sm, r_bh, r_st = pos_raw[-1, 0], pos_raw[-1, 1], pos_raw[-1, 2]
    v_sm, v_bh, v_st = vel_raw[-1, 0], vel_raw[-1, 1], vel_raw[-1, 2]
    # binary separation over the approach, to size the (micro-TDE) inset
    sep = np.linalg.norm(pos_raw[:, 1] - pos_raw[:, 2], axis=1)
    sep_max = np.nanmax(sep[np.isfinite(sep)]) if np.isfinite(sep).any() else 0.0
    # smooth approach for display
    t1, pos = parse_traj.resample_uniform(t_raw, pos_raw, a.max1)

    R_star = a.Rstar * RSUN_AU
    if a.disruptor == "smbh":
        r_t = R_star * (a.Msmbh / a.mstar)**(1.0/3.0)
        r_disr, v_disr = r_sm, v_sm
    else:
        r_t = R_star * (a.mbh / a.mstar)**(1.0/3.0)
        r_disr, v_disr = r_bh, v_bh

    # ---- immediate-disruption guard: synthesize a short straight-line lead-in ----
    if len(t1) < 30:
        npad = 40
        speed = max(np.linalg.norm(v_st - v_disr), 1e-6)
        padT = (20.0 * r_t) / speed                          # travel ~20 r_t before disruption
        tb = np.linspace(-padT, 0.0, npad)
        pos = np.empty((npad, 3, 3))
        for bi, (rr, vv) in enumerate([(r_sm, v_sm), (r_bh, v_bh), (r_st, v_st)]):
            pos[:, bi, :] = rr[None, :] + vv[None, :] * tb[:, None]
        t1 = tb - tb[0]

    # ---- debris ICs -------------------------------------------------------------
    dpos, dvel = debris_ic(r_st, v_st, r_disr, v_disr, a.mstar, R_star, a.ndeb, rng)

    # Follow time from the disruption kinematics.
    if a.disruptor == "smbh":
        # classic frozen-in spread across the stellar radius in the SMBH potential:
        # dE ~ G M_smbh R* / r_t^2  ->  most-bound a_mb = r_t^2 / (2 R*), period P_mb.
        dE = G * a.Msmbh * R_star / r_t**2
        a_mb = r_t**2 / (2.0 * R_star)
    else:
        # micro-TDE: the binary's bulk motion decides who is bound to the SMBH;
        # dE ~ v_BH,SMBH * v_disrupt + 0.5 v_disrupt^2, v_disrupt = sqrt(G m_bh / r_t).
        v_rel = np.linalg.norm(v_bh - v_sm)
        dv_disrupt = np.sqrt(G * a.mbh / r_t)
        dE = v_rel * dv_disrupt + 0.5 * dv_disrupt**2
        a_mb = G * a.Msmbh / (2.0 * dE)
    P_mb = 2*np.pi*np.sqrt(a_mb**3 / (G*a.Msmbh))
    t_follow = a.norbits * P_mb

    F = follow(r_sm, v_sm, r_bh, v_bh, dpos, dvel, a.Msmbh, a.mbh, r_t,
               t_follow, n_frames=a.nframes, disruptor=a.disruptor, settle_orbits=a.settle_orbits)
    n_bound = int((F["eps"] < 0).sum())
    print(f"[{a.disruptor}] r_t={r_t:.5f} AU  a_mb={a_mb:.3g} AU  P_mb={P_mb:.4g} yr  "
          f"t_follow={t_follow:.4g} yr  N={a.ndeb}  bound(SMBH)={n_bound}  ejected(SMBH)={a.ndeb-n_bound}")
    if a.disruptor == "bh":
        print(f"           returning to the stellar BH: bound(BH)={F['n_bound_bh']}  "
              f"accreted(BH)={F['n_acc_bh']}  accreted(SMBH)={F['n_acc_sm']}")

    # ---- write canonical CSV (SMBH/BH/star; star NaN during fallback) -----------
    N1 = len(t1); SM, BH, DEB, TT = F["sm"], F["bh"], F["deb"], F["t"]
    N2 = len(TT)
    cols = ["t", "smbh_x", "smbh_y", "bh_x", "bh_y", "star_x", "star_y"]
    with open(a.out_csv, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for i in range(N1):
            w.writerow([f"{v:.6g}" for v in
                        [t1[i], pos[i,0,0], pos[i,0,1], pos[i,1,0], pos[i,1,1],
                         pos[i,2,0], pos[i,2,1]]])
        for i in range(N2):
            w.writerow([f"{t1[-1]+TT[i]:.6g}", f"{SM[i,0]:.6g}", f"{SM[i,1]:.6g}",
                        f"{BH[i,0]:.6g}", f"{BH[i,1]:.6g}", "nan", "nan"])

    # ---- debris positions (xy) over the fallback frames + colouring energy ------
    out_stem = a.out_csv.rsplit(".csv", 1)[0]
    np.savez_compressed(out_stem + ".debris.npz",
                        xy=DEB[:, :, :2].astype(np.float32),     # (N2, Ndeb, 2)
                        energy=F["eps"].astype(np.float32),      # (Ndeb,)
                        start_frame=np.int64(N1))

    # ---- render config ----------------------------------------------------------
    if a.disruptor == "smbh":
        # main = SMBH-centred stream; inset = the inner launch region around the SMBH
        inset_hw = max(0.15 * a_mb, 10.0 * r_t)
        meta = {"center": "smbh", "zoom_on": ["debris"], "zoom_pct": 75,
                "inset": {"center": "smbh", "show": ["smbh", "debris"], "hw": float(inset_hw),
                          "label": "SMBH + debris (zoom)"},
                "debris": {"npz": out_stem.split("/")[-1] + ".debris.npz"},
                "notes": f"SMBH TDE: n=3 polytrope, {a.ndeb} tracers; r_t={r_t:.4f} AU"}
    else:
        inset_hw = max(1.2*sep_max, 10.0*r_t)       # resolve the binary + the mini-disk
        meta = {"center": "smbh", "zoom_on": ["bh", "debris"], "zoom_pct": 65,
                "inset": {"center": "bh", "show": ["bh", "star", "debris"], "hw": float(inset_hw),
                          "label": "stellar BH + debris (zoom)"},
                "debris": {"npz": out_stem.split("/")[-1] + ".debris.npz"},
                "notes": f"micro-TDE: n=3 polytrope, {a.ndeb} tracers; r_t={r_t:.4f} AU"}
    json.dump(meta, open(out_stem + ".events.json", "w"), indent=1)
    print(f"wrote {a.out_csv}  ({N1}+{N2} frames)  and {out_stem}.debris.npz")


if __name__ == "__main__":
    main()
