#!/usr/bin/env python3
import os, sys, numpy as np, pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else \
      os.path.join(os.environ.get("SMBH_SCRATCH", "/path/to/scratch"),
                   "smbh_bhstar_bmaxtest", "smbh_bhstar_results_test.csv")
RP_TIDAL_FAC  = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
RP_PLUNGE_FAC = 2.0
V1, C_KMS, M_BIN = 29.7847, 299792.458, 11.0

d = pd.read_csv(CSV)
for c in ["M_smbh", "a_bin", "v_inf", "beta"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
d = d.dropna(subset=["M_smbh", "a_bin", "v_inf", "beta"])

M, a, vinf, beta = d.M_smbh.values, d.a_bin.values, d.v_inf.values, d.beta.values
r_t    = a * (M / M_BIN) ** (1/3)
r_g    = (V1*V1 / C_KMS**2) * (M + M_BIN)
r_foc  = 2*V1*V1 * (M + M_BIN) / vinf**2
rp_max = np.maximum(RP_TIDAL_FAC * r_t, RP_PLUNGE_FAC * 4.0 * r_g)   # = max(20 r_t, 8 r_g)
b_max2 = rp_max**2 + r_foc * rp_max
rp     = r_t / beta
bb     = np.sqrt(np.clip((rp**2 + r_foc*rp) / b_max2, 0, 1))

print(f"file: {CSV}")
print(f"rows: {len(d)}   RP_TIDAL_FAC={RP_TIDAL_FAC:g}\n")
print(f"{'channel':16s} {'events':>8}  {'>0.9':>7}  {'>0.8':>7}")
for lbl, s in [("muTDE","TDE:star2_by_bh1"), ("bound_triple","BOUND_TRIPLE"),
               ("plunge_star","PLUNGE:star2_into_SMBH"), ("TDE_central","TDE:star2_by_SMBH"),
               ("plunge_BH","PLUNGE:bh1_into_SMBH"), ("exchange_bh","EXCHANGE_bh1"),
               ("exchange_star","EXCHANGE_star2"), ("flyby","FLYBY")]:
    bx = bb[(d.outcome == s).values]
    if bx.size:
        print(f"{lbl:16s} {bx.size:8d}  {100*(bx>0.9).mean():6.1f}%  {100*(bx>0.8).mean():6.1f}%")
