#!/usr/bin/env python3
# Emits one line per smbh_binary_* run for the Yu & Lai (2024) Fig.5-6 verification.
# Each line is the full argument string for ./smbh_binary.
#   beta in {0.5,0.6,0.7,0.8,1.0,3.0,5.0,10.0} x alpha 3..15
#   x 3 inclinations (30,60,90 deg, stratified) x N_PER_INCL realizations
#   = 8 * 13 * 3 * 1000 = 312,000 runs
# Sampling: Omega ~ U[0,2pi), nu ~ U[0,2pi), omega=0 (e=0), seed 42.

import math
import random

M_SMBH = 1.0e6
M_STAR = 1.0
R_STAR = 1.0
E_BIN  = 0.0
RTOL   = 1.0e-9

N_PER_INCL = 1000
INCLS_DEG  = [30.0, 60.0, 90.0]
SEED       = 42

BETAS  = [0.5, 0.6, 0.7, 0.8, 1.0, 3.0, 5.0, 10.0]
ALPHAS = list(range(3, 16))          # 3..15 inclusive

R_SUN_AU = 0.00465047
V_FAC    = 29.7847
R_col_AU = (2.0 * R_STAR) * R_SUN_AU


def emit():
    rng = random.Random(SEED)
    for beta in BETAS:
        for alpha in ALPHAS:
            a_bin = alpha * R_col_AU
            v_orb = math.sqrt(2.0 * M_STAR / a_bin) * V_FAC
            v_inf = 0.1 * v_orb
            for i_deg in INCLS_DEG:
                incl = math.radians(i_deg)
                for _ in range(N_PER_INCL):
                    Omega = rng.uniform(0.0, 2.0 * math.pi)
                    omega = 0.0
                    nu    = rng.uniform(0.0, 2.0 * math.pi)
                    print(
                        f"--M_smbh {M_SMBH:.6e} "
                        f"--m1 {M_STAR} --type1 star --R1 {R_STAR} "
                        f"--m2 {M_STAR} --type2 star --R2 {R_STAR} "
                        f"--a_bin {a_bin:.10g} --e_bin {E_BIN} "
                        f"--v_inf {v_inf:.10g} --beta {beta} "
                        f"--incl {incl:.10g} --Omega {Omega:.10g} "
                        f"--omega {omega} --nu {nu:.10g} "
                        f"--mode classify --rtol {RTOL:.1e} --report 0"
                    )


if __name__ == "__main__":
    emit()
