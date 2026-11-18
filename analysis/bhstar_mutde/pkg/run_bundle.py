#!/usr/bin/env python3
# =============================================================================
#  run_bundle.py  —  run the debris tracer on mu-TDE events taken straight from
#                    the SeaWulf scattering bundle (their saved disruption-time
#                    phase space are the initial conditions).
#
#  Usage:
#     python3 run_bundle.py BUNDLE.h5 [--rows 0 1 2 | --nfirst K | --all]
#                           [--ndeb 1000] [--nframes 300] [--poly 3]
#                           [--movie] [--outdir debris_out]
#
#  For each selected mu-TDE row it:
#     * reads the disruption state (read_bundle_row),
#     * seeds N polytrope debris with the frozen-in tidal spin,
#     * integrates the 2 BHs + debris (vectorized Newtonian),
#     * saves the FULL trajectories to <outdir>/debris_row<idx>.h5,
#     * writes an energy/fate figure, and (optionally) a movie for the first row,
#     * prints the bound-BH1 / bound-BH2 / unbound fractions and accuracy checks.
#
#  The vectorized Newtonian integrator runs anywhere (numpy/scipy); no SpaceHub
#  or cluster needed.  Use the SpaceHub C++ driver only for the consistency check.
# =============================================================================
import os, sys, argparse, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import debris_tracer as dt
try:
    import style as _sty; _sty.use_paper_style(base=12)
except Exception:
    pass

def run_one(bundle, idx, args, make_movie=False):
    st = dt.read_bundle_row(bundle, row_index=idx)
    r, v, omega, info = dt.make_debris_ic(st, args.ndeb, n=args.poly, spin=True, seed=1)
    # integration time: several most-bound-debris orbits about BH1, and ~1 BH1-BH2 orbit
    a_bound = st["m_BH1"]/(2.0*info["Deps_scale"])
    T_bound = 2*np.pi*np.sqrt(abs(a_bound)**3/st["m_BH1"])
    sep12   = np.linalg.norm(st["r_BH1"]-st["r_BH2"])
    T_pair  = 2*np.pi*np.sqrt(sep12**3/(st["M2"]+st["m_BH1"]))
    t_end   = max(10*T_bound, 1.0*T_pair)
    res = dt.integrate_debris(st, r, v, t_end, n_frames=args.nframes)
    os.makedirs(args.outdir, exist_ok=True)
    dt.save_h5(res, os.path.join(args.outdir, f"debris_row{idx}.h5"))
    cl  = dt.classify(res, frame=-1); cl0 = dt.classify(res, frame=0)
    print(f"  row {idx}: M2={st['M2']:.2g} m1={st['m_BH1']:.0f}  "
          f"BH1={cl['frac']['BH1']*100:.1f}% BH2={cl['frac']['BH2']*100:.1f}% "
          f"unb={cl['frac']['unbound']*100:.1f}%  |dE/E|={res['dE_bh'].max():.1e}  "
          f"t_end={t_end*dt.YR_PER_T:.3g} yr")
    # --- energy + fate figure ---
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    w = info["Deps_scale"]
    ax[0].hist(cl0["eps_BH1"], bins=60, range=(-3*w, 3*w), color="steelblue", edgecolor="none")
    ax[0].axvline(0, color="k", ls="--", lw=1); ax[0].axvspan(-w, w, color="orange", alpha=0.12)
    ax[0].set_xlabel(r"frozen-in $\varepsilon_{\rm BH1}$ [$V_1^2$]"); ax[0].set_ylabel(r"$dN/d\varepsilon$")
    ax[0].set_title(f"fallback energy (row {idx})")
    deb = res["traj_deb"][:, -1, :3]; bh = res["traj_bh"][:, -1, :3]
    col = np.array(["#d62728", "#1f77b4", "#7f7f7f"])
    for lab, nm in [(0,"BH1"),(1,"BH2"),(2,"unb")]:
        m = cl["label"] == lab
        ax[1].scatter(deb[m,0], deb[m,1], s=4, c=col[lab], alpha=0.6, label=f"{nm} ({m.sum()})")
    ax[1].scatter(*bh[0,:2], marker="o", s=40, c="k"); ax[1].scatter(*bh[1,:2], marker="o", s=170, c="k")
    ax[1].set_xlabel("x [AU]"); ax[1].set_ylabel("y [AU]"); ax[1].set_aspect("equal","datalim")
    ax[1].legend(fontsize=7); ax[1].set_title(f"final fate (t={res['t'][-1]*dt.YR_PER_T:.2g} yr)")
    fig.tight_layout(); fig.savefig(os.path.join(args.outdir, f"debris_row{idx}.png"), dpi=140); plt.close(fig)
    # publication figure: BH1 orbit + representative debris trajectories (+ energy inset)
    dt.paper_figure(res, os.path.join(args.outdir, f"debris_paper_row{idx}.png"), n_curves=48)

    if make_movie:
        from matplotlib import animation
        tb, td = res["traj_bh"], res["traj_deb"]; span = 60*st.get("r_t_star", 0.01)
        fm, axm = plt.subplots(figsize=(6.4,6.0))
        sc = axm.scatter([], [], s=3, c=[], cmap="coolwarm", vmin=-1, vmax=1)
        p0, = axm.plot([], [], "k*", ms=13)
        axm.set_xlim(-span,span); axm.set_ylim(-span,span); axm.set_aspect("equal")
        axm.set_xlabel("x$-$x$_{BH1}$ [AU]"); axm.set_ylabel("y$-$y$_{BH1}$ [AU]"); tt = axm.set_title("")
        def fr(i):
            r1=tb[0,i,:3]; v1=tb[0,i,3:]; rel=td[:,i,:3]-r1; dv=td[:,i,3:]-v1
            e1=0.5*np.einsum("ij,ij->i",dv,dv)-st["m_BH1"]/np.linalg.norm(td[:,i,:3]-r1,axis=1)
            sc.set_offsets(rel[:,:2]); sc.set_array(np.tanh(e1/(np.std(e1)+1e-30))); p0.set_data([0],[0])
            tt.set_text(f"t={res['t'][i]*dt.YR_PER_T*dt.SEC_PER_YR/86400:.2f} d"); return sc,p0,tt
        an = animation.FuncAnimation(fm, fr, frames=args.nframes, interval=40, blit=True)
        mp4 = os.path.join(args.outdir, f"debris_row{idx}.mp4")
        try: an.save(mp4, writer=animation.FFMpegWriter(fps=25, bitrate=2400)); print(f"    movie -> {mp4}")
        except Exception: an.save(mp4.replace(".mp4",".gif"), writer=animation.PillowWriter(fps=20))
        plt.close(fm)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle"); ap.add_argument("--rows", type=int, nargs="+")
    ap.add_argument("--nfirst", type=int); ap.add_argument("--all", action="store_true")
    ap.add_argument("--ndeb", type=int, default=1000); ap.add_argument("--nframes", type=int, default=300)
    ap.add_argument("--poly", type=float, default=3.0); ap.add_argument("--movie", action="store_true")
    ap.add_argument("--outdir", default="debris_out")
    args = ap.parse_args()
    mu = dt.mutde_indices(args.bundle)
    print(f"{len(mu)} mu-TDE rows in {args.bundle}")
    if args.all:        sel = list(range(len(mu)))
    elif args.nfirst:   sel = list(range(min(args.nfirst, len(mu))))
    elif args.rows:     sel = args.rows
    else:               sel = [0]
    for j, idx in enumerate(sel):
        run_one(args.bundle, idx, args, make_movie=(args.movie and j == 0))

if __name__ == "__main__":
    main()
