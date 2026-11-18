# reduce_bundle.py  —  run on SeaWulf, outputs mutde_cellsummary.csv
# b_max recomputed analytically per cell (same formula as the sweep); no grid_meta join.
import h5py, numpy as np, csv

RP_TIDAL_FAC  = 10.0          # <-- MUST match BHSTAR_RP_TIDAL_FAC used for the run
RP_PLUNGE_FAC = 2.0
V1, C_KMS, M_BIN = 29.7847, 299792.458, 11.0

f = h5py.File('smbh_bhstar_bundle.h5', 'r'); R = f['results']
M    = R['M_smbh'][:]; a = R['a_bin'][:]; e = R['e_bin'][:]
vinf = R['v_inf'][:];  beta = R['beta'][:]
oc   = R['outcome'][:].astype(str)
frac = vinf / (V1 * np.sqrt(M_BIN / a))            # v_inf / v_orb

mu  = (oc == 'TDE:star2_by_bh1')                   # micro-TDE
key = np.stack([M, a, e, np.round(frac, 6)], 1)
uniq, inv = np.unique(key, axis=0, return_inverse=True)
inv  = inv.ravel()                                 # guard vs numpy>=2.0 shape
Ntot = np.bincount(inv, minlength=len(uniq))
Nmu  = np.bincount(inv, weights=mu.astype(float), minlength=len(uniq))

def b_max_of(Mc, ac, vc):
    r_t    = ac * (Mc / M_BIN) ** (1.0 / 3.0)
    r_g    = (V1 * V1 / C_KMS**2) * (Mc + M_BIN)
    r_foc  = 2.0 * V1 * V1 * (Mc + M_BIN) / (vc * vc)
    rp_max = max(RP_TIDAL_FAC * r_t, RP_PLUNGE_FAC * 4.0 * r_g)
    return np.sqrt(rp_max * rp_max + r_foc * rp_max)

with open('mutde_cellsummary.csv', 'w', newline='') as g:
    w = csv.writer(g)
    w.writerow(['M_smbh','a_bin','e_bin','frac','N_tot','N_mu','b_max_AU',
                'beta_min','beta_med','beta_max'])
    for i, (Mc, ac, ec, fc) in enumerate(uniq):
        m  = (inv == i)
        bb = beta[m]
        vc = np.median(vinf[m])                    # cell v_inf (constant within a cell)
        w.writerow([Mc, ac, ec, fc, int(Ntot[i]), int(Nmu[i]),
                    b_max_of(Mc, ac, vc),
                    bb.min(), np.median(bb), bb.max()])
print('wrote mutde_cellsummary.csv  cells=', len(uniq), ' total sims=', Ntot.sum())
