#!/usr/bin/env python3
# =============================================================================
#  reduce_cellsummary.py  --  bundle  ->  mutde_cellsummary.csv
#  -----------------------------------------------------------------------------
#  Collapse the 27.6M-row scattering bundle into a tiny per-(M,a,e,v) cell table
#  that the cross-section and loss-cone notebooks consume.  This is the ONLY heavy
#  step on the full bundle; run it ONCE and the notebooks then work on the ~20k-row
#  CSV.  Reads the bundle in chunks so it fits in modest RAM (a laptop is fine).
#
#  Output columns (exactly what analyze_mutde_cross_sections.ipynb expects):
#     M_smbh, a_bin, e_bin, frac, N_tot, N_mu, beta_med
#     - frac    : v_inf/v_orb snapped to the 7 canonical nodes
#     - N_tot   : encounters in the cell
#     - N_mu    : encounters ending in TDE:star2_by_bh1 (the mu-TDE)
#     - beta_med: median penetration beta over the cell's mu-TDE events
#                 (diagnostic only; not used in the f_mu fit)
#
#  Usage:
#     python3 reduce_cellsummary.py smbh_bhstar_bundle.h5 [--out mutde_cellsummary.csv]
#                                   [--chunk 2000000] [--mu TDE:star2_by_bh1]
# =============================================================================
import sys, os, argparse, csv
from collections import defaultdict
import numpy as np, h5py

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scales as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bundle")
    ap.add_argument("--out", default="mutde_cellsummary.csv")
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--mu", default="TDE:star2_by_bh1")
    args = ap.parse_args()

    with h5py.File(args.bundle, "r") as f:
        R = f["results"]
        need = ["M_smbh", "a_bin", "e_bin", "beta", "v_inf", "outcome"]
        missing = [c for c in need if c not in R]
        if missing:
            raise KeyError(f"bundle /results missing {missing}; have {sorted(R.keys())}")
        n = R["M_smbh"].shape[0]
        print(f"{args.bundle}: {n:,} encounters; reducing in chunks of {args.chunk:,}")

        n_tot = defaultdict(int)             # (logM, a, e, logfrac) -> count
        n_mu  = defaultdict(int)
        beta_mu = defaultdict(list)          # mu-TDE betas per cell (small)

        for s in range(0, n, args.chunk):
            t = min(s + args.chunk, n)
            M  = R["M_smbh"][s:t]; a = R["a_bin"][s:t]; e = R["e_bin"][s:t]
            vinf = R["v_inf"][s:t]; beta = R["beta"][s:t]
            oc = R["outcome"].asstr()[s:t]
            frac = sc.frac_from_vinf(vinf, a)
            keys = list(zip(np.round(np.log10(M), 4), np.round(a, 5),
                            np.round(e, 4), np.round(np.log10(frac), 4)))
            ismu = (oc == args.mu)
            for k in keys:
                n_tot[k] += 1
            for k, b, mu in zip(keys, beta, ismu):
                if mu:
                    n_mu[k] += 1
                    beta_mu[k].append(b)
            print(f"  rows {s:,}..{t:,}  cells so far: {len(n_tot):,}")

    # --- recover the (M,a,e,frac) physical values from the rounded keys ------
    rows = []
    for k, Nt in sorted(n_tot.items()):
        lM, a, e, lfr = k
        M = 10.0**lM; frac = 10.0**lfr
        Nmu = n_mu.get(k, 0)
        bmed = float(np.median(beta_mu[k])) if beta_mu.get(k) else np.nan
        rows.append((M, a, e, frac, Nt, Nmu, bmed))

    with open(args.out, "w", newline="") as g:
        w = csv.writer(g)
        w.writerow(["M_smbh", "a_bin", "e_bin", "frac", "N_tot", "N_mu", "beta_med"])
        for M, a, e, frac, Nt, Nmu, bmed in rows:
            w.writerow([f"{M:.6g}", f"{a:.6g}", f"{e:.6g}", f"{frac:.6g}",
                        Nt, Nmu, f"{bmed:.6g}"])

    Nt_tot = sum(n_tot.values()); Nmu_tot = sum(n_mu.values())
    print(f"\nwrote {args.out}: {len(rows):,} cells, "
          f"{Nt_tot:,} encounters, {Nmu_tot:,} mu-TDEs (overall f_mu={Nmu_tot/Nt_tot:.4f})")


if __name__ == "__main__":
    main()
