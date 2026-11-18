#!/usr/bin/env python3
"""Write a smbh_binary_PN RESTART file from the final snapshot of a stage-1 run.
Used to continue the integration (e.g. resolve the EMRI after the star is ejected):
the restart starts from stage-1's exact final state, so the seam is continuous.

   python3 make_restart.py bss_..._full.txt restart.txt [--Msmbh 1e6 --mbh 10 --mstar 1]

Restart line format (from --help):
   mass_Msun radius_Rsun type x_AU y_AU z_AU vx_kms vy_kms vz_kms
order = body0(SMBH) body1(BH) body2(star).  Velocities are taken from the last two
snapshots (AU/yr -> km/s), which is robust to the dump's own velocity units."""
import argparse, numpy as np
import parse_traj

AU_YR_TO_KMS = 4.740572
def two_rg_Rsun(M):        # 2 r_g in Rsun (matches the sweep's auto radii)
    return 4.245e-6 * M

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw"); ap.add_argument("out")
    ap.add_argument("--Msmbh", type=float, default=1e6)
    ap.add_argument("--mbh", type=float, default=10.0)
    ap.add_argument("--mstar", type=float, default=1.0)
    a = ap.parse_args()
    tr = parse_traj.parse(a.raw); pos = tr["pos"]; t = tr["t"]
    p = pos[-1]                                          # (3,3) AU  (last snapshot)
    if tr.get("vel") is not None:                        # exact state from the dump
        vel = tr["vel"][-1]                              # (3,3) km/s
    else:                                                # fallback: finite-difference
        dt = t[-1] - t[-2]
        vel = (pos[-1] - pos[-2]) / dt * AU_YR_TO_KMS    # (3,3) km/s
    masses = [a.Msmbh, a.mbh, a.mstar]
    radii  = [two_rg_Rsun(a.Msmbh), two_rg_Rsun(a.mbh), 1.0]
    types  = ["bh", "bh", "star"]
    with open(a.out, "w") as f:
        f.write("# restart: body0=SMBH body1=BH body2=star\n")
        f.write("# mass_Msun radius_Rsun type x_AU y_AU z_AU vx_kms vy_kms vz_kms\n")
        for i in range(3):
            f.write(f"{masses[i]:.10g} {radii[i]:.6g} {types[i]} "
                    f"{p[i,0]:.10g} {p[i,1]:.10g} {p[i,2]:.10g} "
                    f"{vel[i,0]:.10g} {vel[i,1]:.10g} {vel[i,2]:.10g}\n")
    print(f"wrote {a.out} from final snapshot (t={t[-1]:.6g})")
    # orbital period about the SMBH of each companion (for sizing the restart integration)
    G = 4*np.pi**2
    def period(i, m):
        r = p[i]-p[0]; vrel = (vel[i]-vel[0])/AU_YR_TO_KMS   # AU, AU/yr
        mu = G*(a.Msmbh+m); rr = np.linalg.norm(r)
        sma = 1.0/(2.0/rr - (vrel@vrel)/mu)
        return sma, 2*np.pi*np.sqrt(abs(sma)**3/mu)
    sma_bh, P_bh = period(1, a.mbh); sma_st, P_st = period(2, a.mstar)
    print(f"SMA_AU={sma_bh:.6g}"); print(f"PERIOD_YR={P_bh:.6g}")     # BH (back-compat)
    print(f"PERIOD_BH_YR={P_bh:.6g}"); print(f"PERIOD_STAR_YR={P_st:.6g}")

if __name__ == "__main__":
    main()
