#!/usr/bin/env python3
"""EJECTION handoff: read a stage-1 canonical trajectory, take the BH's exact state
(r, v) relative to the SMBH at the final frame, and emit the osculating 2-body
orbital elements for the stage-2 EMRI run (BH + SMBH, full PN incl. 2.5PN in SpaceHub).

   python3 handoff.py scenes/ejection_stage1.csv [--Msmbh 1e6 --mbh 10]

Prints a, e, argument of pericenter, true anomaly and writes handoff.json. Feed these
to your 2-body SpaceHub run; starting stage 2 from this exact state makes the N-body ->
EMRI seam continuous (position, velocity, a, e all match)."""
import json, argparse, numpy as np
import scene_lib
from nbody import osculating

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage1_csv")
    ap.add_argument("--Msmbh", type=float, default=1e6)
    ap.add_argument("--mbh", type=float, default=10.0)
    ap.add_argument("-o", "--out", default="handoff.json")
    a = ap.parse_args()
    s = scene_lib.load(a.stage1_csv); B = s["bodies"]; t = s["t"]
    # BH relative to SMBH; velocity by finite difference of the last two frames
    rbh = B["bh"][-1] - B["smbh"][-1]
    rbh_prev = B["bh"][-2] - B["smbh"][-2]
    dt = t[-1] - t[-2]
    vbh = (rbh - rbh_prev) / dt
    osc = osculating(rbh, vbh, a.Msmbh + a.mbh)
    res = dict(a_AU=osc["a"], e=osc["e"], rp_AU=osc["rp"],
               omega_rad=osc["omega"], nu_rad=osc["nu"],
               r_xy=list(rbh), v_xy=list(vbh), Msmbh=a.Msmbh, mbh=a.mbh)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"handoff: a={osc['a']:.3f} AU  e={osc['e']:.5f}  rp={osc['rp']:.4f} AU  "
          f"omega={np.degrees(osc['omega']):.1f} deg  nu={np.degrees(osc['nu']):.1f} deg")
    print(f"BH state (SMBH frame): r={rbh} AU  v={vbh} AU/yr")
    print(f"-> wrote {a.out}. Use (a,e,omega,nu) as the stage-2 BH+SMBH orbit; the full")
    print("   2.5PN SpaceHub run then inspirals to merger, continuous across the seam.")

if __name__ == "__main__":
    main()
