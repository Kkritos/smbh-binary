#!/usr/bin/env python3
# =============================================================================
#  make_test_bundle.py  --  a small, SCHEMA-CORRECT synthetic bundle.
#  -----------------------------------------------------------------------------
#  Writes smbh_bhstar_bundle.h5-shaped data (a few hundred cells x realisations)
#  so the whole v2 pipeline -- reduce_cellsummary, reduce_delay, and all five
#  notebooks -- can be exercised end-to-end WITHOUT the real 9 GB bundle.
#
#  It is a faithful stand-in, not a physics model:
#    * impact parameter sampled b^2 ~ U[0, b_max^2]  (the real estimator),
#    * penetration beta = r_t/r_p reconstructed from the sampled b^2,
#    * outcomes assigned by a KNOWN closed-form threshold beta_crit(rho,e) with
#      coefficients CF_TRUE, so analyze_mutde_cross_sections should ~recover CF_TRUE,
#    * the rp_max = 10 r_t (v2) convention, written into /grid_meta,
#    * full disruption-time phase space (x_smbh_AU ... vz_b2_kms) for the mu-TDE
#      rows, built from debris_tracer.build_parabolic_state.
#
#  Usage:  python3 make_test_bundle.py [out.h5] [--n-per-cell 200] [--seed 0]
# =============================================================================
import sys, os, argparse
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scales as sc
import debris_tracer as dt

# a known beta_crit(rho,e) law the refit should recover:  log10 bc = CF.[1,L,L^2,e,e^2]
CF_TRUE = np.array([-0.26, 1.47, -0.21, 0.18, -0.66])
OUTCOMES_STRONG = ["EXCHANGE_bh1", "EXCHANGE_star2", "BOUND_TRIPLE"]

STR = h5py.string_dtype("utf-8")


def beta_crit_true(a, e):
    L = np.log10(a * (1.0 - e) / sc.R_TDE)
    return 10.0**(CF_TRUE[0] + CF_TRUE[1]*L + CF_TRUE[2]*L**2 + CF_TRUE[3]*e + CF_TRUE[4]*e**2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default="smbh_bhstar_bundle.h5")
    ap.add_argument("--n-per-cell", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    # grid matched to the REAL campaign node structure so every hardcoded notebook
    # slice (M half-decades, a=0.3, e=0.5/0.9, v/vorb=0.01 & 3.0) resolves to a node.
    M_NODES = np.logspace(1, 10, 19)                                  # 10^1 .. 10^10, 19 half-decades
    A_NODES = np.logspace(np.log10(0.03), np.log10(30.0), 13)         # 0.03 .. 30 AU, 13 nodes (incl 0.3)
    E_NODES = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99])
    FR_NODES = sc.FRAC_NODES                                          # the 7 canonical v_inf/v_orb nodes
    Npc = args.n_per_cell

    res = {k: [] for k in
           ["M_smbh", "m1", "m2", "R2", "a_bin", "e_bin", "beta", "v_inf", "outcome"]}
    phase_cols = (tuple(f"{c}_{p}_AU" for p in ("smbh", "b1", "b2") for c in "xyz")
                  + tuple(f"v{c}_{p}_kms" for p in ("smbh", "b1", "b2") for c in "xyz"))
    for k in phase_cols:
        res[k] = []

    gm = {k: [] for k in ["M_smbh", "a_bin_AU", "frac", "v_orb_kms", "v_inf_kms",
                          "r_t_bin_AU", "r_g_smbh_AU", "rp_max_AU", "b_max_AU", "N"]}

    def zeros_phase():
        return {k: 0.0 for k in phase_cols}

    for M in M_NODES:
        for a in A_NODES:
            for fr in FR_NODES:
                r_t, r_foc, rp_max, b2max = sc.cell_scales(M, a, fr)
                v_orb = sc.V1 * np.sqrt(sc.M_BIN / a); v_inf = fr * v_orb
                gm["M_smbh"].append(M); gm["a_bin_AU"].append(a); gm["frac"].append(fr)
                gm["v_orb_kms"].append(v_orb); gm["v_inf_kms"].append(v_inf)
                gm["r_t_bin_AU"].append(r_t)
                gm["r_g_smbh_AU"].append((sc.V1/sc.C_KMS)**2 * (M + sc.M_BIN))
                gm["rp_max_AU"].append(rp_max); gm["b_max_AU"].append(np.sqrt(b2max))
                for e in E_NODES:
                    bc = beta_crit_true(a, e)
                    # b^2 ~ U[0, b_max^2]  ->  r_p  ->  beta
                    u = rng.random(Npc)
                    b2 = u * b2max
                    r_p = (-r_foc + np.sqrt(r_foc**2 + 4.0*b2)) / 2.0
                    beta = r_t / r_p
                    for bk in beta:
                        res["M_smbh"].append(M); res["m1"].append(10.0); res["m2"].append(1.0)
                        res["R2"].append(1.0); res["a_bin"].append(a); res["e_bin"].append(e)
                        res["beta"].append(bk); res["v_inf"].append(v_inf)
                        # outcome ladder
                        soft = fr >= 1.0
                        if bk < 0.3:
                            oc = "IONIZATION" if (soft and rng.random() < 0.3) else "FLYBY"
                        elif bk < bc:
                            oc = OUTCOMES_STRONG[rng.integers(3)]
                        else:
                            if M >= 1e9 and rng.random() < 0.5:
                                oc = rng.choice(["TDE:star2_by_SMBH", "PLUNGE:bh1_into_SMBH",
                                                 "PLUNGE:star2_into_SMBH"])
                            else:
                                oc = "TDE:star2_by_bh1"
                        res["outcome"].append(oc)
                        # phase space: only mu-TDE rows carry a real disruption state
                        if oc == "TDE:star2_by_bh1":
                            st = dt.build_parabolic_state(M2=M, a_bin=a,
                                                          beta_pen=float(np.clip(bk, 1.0, 50.0)))
                            for pre, rk, vk in [("smbh", "r_BH2", "v_BH2"),
                                                ("b1", "r_BH1", "v_BH1"),
                                                ("b2", "r_star", "v_star")]:
                                r = st[rk]; v = st[vk] * sc.V1
                                res[f"x_{pre}_AU"].append(r[0]); res[f"y_{pre}_AU"].append(r[1])
                                res[f"z_{pre}_AU"].append(r[2])
                                res[f"vx_{pre}_kms"].append(v[0]); res[f"vy_{pre}_kms"].append(v[1])
                                res[f"vz_{pre}_kms"].append(v[2])
                        else:
                            z = zeros_phase()
                            for k in phase_cols:
                                res[k].append(z[k])

    gm["N"] = [Npc * len(E_NODES)] * len(gm["M_smbh"])

    # ecc_meta: survival ceiling per a_bin
    ecc = {"a_bin_AU": [], "e_max_survival": [], "n_e_nodes_used": []}
    for a in A_NODES:
        emax = float(np.clip(1.0 - sc.R_TDE / a, 0.6, 0.995))
        ecc["a_bin_AU"].append(a); ecc["e_max_survival"].append(emax)
        ecc["n_e_nodes_used"].append(len(E_NODES))

    with h5py.File(args.out, "w") as f:
        g = f.create_group("results")
        for k, v in res.items():
            if k == "outcome":
                g.create_dataset(k, data=np.array(v, dtype=STR))
            else:
                g.create_dataset(k, data=np.asarray(v, float), compression="gzip")
        g.attrs["n_rows"] = len(res["M_smbh"])
        gg = f.create_group("grid_meta")
        for k, v in gm.items():
            gg.create_dataset(k, data=np.asarray(v, float))
        ge = f.create_group("ecc_meta")
        for k, v in ecc.items():
            ge.create_dataset(k, data=np.asarray(v, float))

    n = len(res["M_smbh"]); nmu = int(np.sum(np.array(res["outcome"]) == "TDE:star2_by_bh1"))
    print(f"wrote {args.out}: {n:,} encounters, {nmu:,} mu-TDEs "
          f"({len(gm['M_smbh'])} grid_meta rows); CF_TRUE={CF_TRUE.tolist()}")


if __name__ == "__main__":
    main()
