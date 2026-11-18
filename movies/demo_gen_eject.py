"""Generate the ejection+capture+EMRI canonical trajectory (PREVIEW).
Stage 1: real N-body 3-body slingshot (SMBH+BH+star), star ejected, BH captured.
Stage 2: GW inspiral of the captured BH about the SMBH, continued from the exact
handoff state (osculating elements) -> rendered as a shrinking, precessing orbit.
PREVIEW ONLY: production uses SpaceHub direct 2.5PN; here the inspiral e-a track is
Peters(1964) and precession is 1PN, time-compressed."""
import json, numpy as np
from nbody import (make_ic, integrate, classify, osculating, peters_inspiral,
                   kepler_solve, ellipse_point, G, C_AUYR)

M_smbh, m_bh, m_star = 1e6, 10.0, 1.0
a_bin = 0.1; m_b = m_bh + m_star
r_t_bin = a_bin*(M_smbh/m_b)**(1/3)
rp_out = 0.6*r_t_bin; e_out = 1.0
r0 = 12*rp_out
nu0 = -np.arccos(np.clip(((rp_out*(1+e_out))/r0 - 1)/e_out, -1, 1))
phase, sense = np.radians(30.0), +1

# ---------- stage 1: N-body ----------
y0, m = make_ic(M_smbh, m_bh, m_star, a_bin, rp_out, e_out, nu0, phase, bin_sense=sense)
t1 = 0.055
sol = integrate(y0, m, (0.0, t1), 160)
Y = sol.y; tt = sol.t
smbh = Y[0:2].T; bh = Y[2:4].T; star = Y[4:6].T
vsm = Y[6:8].T; vbh = Y[8:10].T; vst = Y[10:12].T
# SMBH-centric
bh_c = bh - smbh; star_c = star - smbh
# ejection frame index: first time star is unbound from SMBH
ej = None
for i in range(len(tt)):
    dr = star[i]-smbh[i]; dv = vst[i]-vsm[i]
    E = 0.5*(dv@dv) - G*M_smbh/np.linalg.norm(dr)
    if E > 0 and np.linalg.norm(dr) > 3*rp_out:
        ej = i; break

c = classify(sol, m)
print("end-of-stage1:", {k: (round(v['E'],1), v['bound']) for k,v in c.items()})

# ---------- handoff: osculating elements of BH about SMBH ----------
rbh = (bh[-1]-smbh[-1]); vrel = (vbh[-1]-vsm[-1])
Mtot = M_smbh + m_bh
osc = osculating(rbh, vrel, Mtot)
a0, e0, om0, nu_h = osc["a"], osc["e"], osc["omega"], osc["nu"]
print(f"handoff: a0={a0:.1f} AU  e0={e0:.4f}  rp={osc['rp']:.3f} AU")

# ---------- stage 2: GW inspiral rendered as N shrinking orbits ----------
Rg = G*M_smbh/C_AUYR**2
a_merge = 8*Rg                      # stop ~ ISCO scale
T, A, E = peters_inspiral(a0, e0, m_bh, M_smbh, a_merge)
t_insp_yr = T[-1]
print(f"Peters inspiral time ~ {t_insp_yr:.3e} yr  ({t_insp_yr/1e6:.2f} Myr)")

# smooth reconstruction: continuous a(s),e(s) over N_orb orbits, physical motion
# (equal mean-anomaly steps -> BH lingers at apocenter, whips through pericenter),
# with 1PN pericenter precession. Continuous at the seam: starts at (a0,e0,nu_h).
N_orb, ppo = 13, 40
s = np.linspace(0.0, N_orb, N_orb*ppo)
a_s = a0 * (a_merge/a0)**(s/N_orb)                     # smooth geometric shrink
e_s = np.interp(a_s[::-1], A[::-1], E[::-1])[::-1]
t_s = np.interp(a_s[::-1], A[::-1], T[::-1])[::-1]
# mean anomaly continuous from the handoff true anomaly
E_h = 2*np.arctan2(np.sqrt(1-e0)*np.sin(nu_h/2), np.sqrt(1+e0)*np.cos(nu_h/2))
M_h = E_h - e0*np.sin(E_h)
M_s = M_h + 2*np.pi*s
# accumulate precession: dom/ds = 6 pi G Mtot / (c^2 a (1-e^2))
dom = 6*np.pi*G*Mtot/(C_AUYR**2 * a_s*(1-e_s**2))
om_s = om0 + np.concatenate([[0], np.cumsum(0.5*(dom[1:]+dom[:-1])*np.diff(s))])
xs2 = np.empty((len(s), 2))
for i in range(len(s)):
    Ei = kepler_solve(((M_s[i]+np.pi) % (2*np.pi)) - np.pi, e_s[i])
    nui = 2*np.arctan2(np.sqrt(1+e_s[i])*np.sin(Ei/2), np.sqrt(1-e_s[i])*np.cos(Ei/2))
    xs2[i] = ellipse_point(a_s[i], e_s[i], om_s[i], nui, Mtot)
ts2 = t_s

# ---------- assemble canonical CSV ----------
nan = np.nan
N1 = len(tt); N2 = len(xs2)
rows = []
for i in range(N1):
    rows.append((tt[i], 0.0, 0.0, bh_c[i,0], bh_c[i,1], star_c[i,0], star_c[i,1]))
for i in range(N2):
    rows.append((ts2[i], 0.0, 0.0, xs2[i,0], xs2[i,1], nan, nan))
import csv
out = "scenes/ejection.csv"
import os; os.makedirs("scenes", exist_ok=True)
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["t","smbh_x","smbh_y","bh_x","bh_y","star_x","star_y"])
    for r in rows: w.writerow([f"{v:.6g}" if np.isfinite(v) else "nan" for v in r])

events = []
if ej is not None: events.append({"frame": ej, "label": "star ejected  (HVS)"})
events.append({"frame": N1, "label": "BH bound → GW inspiral"})
events.append({"frame": N1+N2-1, "label": "merger"})
meta = {"center": "smbh", "zoom_on": ["bh"], "events": events,
        "inset": {"center": "com", "com": {"bh": 10, "star": 1}, "show": ["bh", "star"], "hw": 0.22,
                  "label": "binary (zoom)"},
        "notes": f"a0={a0:.0f}AU e0={e0:.3f}, inspiral~{t_insp_yr/1e6:.1f}Myr (compressed)"}
json.dump(meta, open("scenes/ejection.events.json","w"), indent=1)
print(f"wrote {out}  ({N1}+{N2}={N1+N2} frames)")
