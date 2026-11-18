#!/usr/bin/env python3
# =============================================================================
#  check_bundle.py  --  one-shot sanity check of smbh_bhstar_bundle.h5
#  -----------------------------------------------------------------------------
#  Run this FIRST on the new bundle.  It reports:
#    * n_rows and the full /results column list,
#    * whether the disruption-time PHASE-SPACE columns the debris pipeline needs
#      are present (x_smbh_AU ... vz_b2_kms, m1, m2, R2, v_inf),
#    * the mu-TDE count and overall fraction,
#    * the rp_max/r_t convention stored in /grid_meta (should be 10 for the v2 run),
#      cross-checked against scales.RP_MAX_OVER_RT.
#
#  If the phase-space columns are MISSING, the scattering campaign was written
#  without the final Cartesian state dump; the debris notebooks (debris_fig1,
#  run_bundle, reduce_delay) then cannot run and the campaign must be re-dumped
#  with full final-state output (the cross-section / branching / loss-cone
#  products do NOT need the phase space and still work).
#
#  Usage:  python3 check_bundle.py smbh_bhstar_bundle.h5
# =============================================================================
import sys, os
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scales as sc

PHASE = (tuple(f"{c}_{p}_AU" for p in ("smbh", "b1", "b2") for c in "xyz")
         + tuple(f"v{c}_{p}_kms" for p in ("smbh", "b1", "b2") for c in "xyz")
         + ("m1", "m2", "R2", "v_inf"))
FITCOLS = ("M_smbh", "a_bin", "e_bin", "beta", "v_inf", "outcome")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "smbh_bhstar_bundle.h5"
    with h5py.File(path, "r") as f:
        R = f["results"]
        cols = sorted(R.keys())
        n = R[cols[0]].shape[0]
        print(f"=== {path} ===")
        print(f"/results : {n:,} rows, {len(cols)} columns")
        print("  columns:", ", ".join(cols))

        print("\n-- analysis columns (cross sections / branching / loss cone) --")
        miss_fit = [c for c in FITCOLS if c not in R]
        print("  present" if not miss_fit else f"  MISSING: {miss_fit}")

        print("\n-- debris phase-space columns (debris_fig1 / run_bundle / reduce_delay) --")
        miss_ph = [c for c in PHASE if c not in R]
        if not miss_ph:
            print("  present -> debris pipeline can run on this bundle")
        else:
            print(f"  MISSING {len(miss_ph)}: {miss_ph}")
            print("  -> debris notebooks cannot run; re-dump the campaign with full final state.")

        if "outcome" in R:
            oc = R["outcome"].asstr()[:]
            nmu = int(np.sum(oc == "TDE:star2_by_bh1"))
            print(f"\nmu-TDEs: {nmu:,}  (f_mu = {nmu/n:.4f})")

        if "grid_meta" in f:
            ratio, ok = sc.check_rpmax_convention(path)
            print(f"\ngrid_meta rp_max/r_t (tidal cells): {ratio:.3f}  "
                  f"(scales.RP_MAX_OVER_RT={sc.RP_MAX_OVER_RT}, {'OK' if ok else 'MISMATCH'})")
            gm = sc.read_grid_meta(path)
            if "b_max_AU" in gm:
                print(f"grid_meta b_max_AU: {np.nanmin(gm['b_max_AU']):.3g} .. "
                      f"{np.nanmax(gm['b_max_AU']):.3g} AU")
        else:
            print("\n(no /grid_meta table in bundle; rp_max convention taken from scales.RP_MAX_OVER_RT)")


if __name__ == "__main__":
    main()
