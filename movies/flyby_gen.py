"""Generate the flyby canonical trajectory (PREVIEW): a distant SMBH passage that
perturbs but does NOT disrupt the BH-star binary. SMBH-centric + binary inset."""
import json, os, csv, numpy as np
from nbody import make_ic, integrate, classify, G

M_smbh, m_bh, m_star = 1e6, 10.0, 1.0
a_bin, e_bin = 0.3, 0.2; m_b = m_bh+m_star
r_t_bin = a_bin*(M_smbh/m_b)**(1/3)
rp_out = 2.0*r_t_bin           # beta ~ 0.5 : survives, mildly perturbed
e_out = 1.0; r0 = 6*rp_out
nu0 = -np.arccos(np.clip(((rp_out*(1+e_out))/r0-1)/e_out,-1,1))
phase = np.radians(40.0)

y0, m = make_ic(M_smbh, m_bh, m_star, a_bin, rp_out, e_out, nu0, phase, e_bin=e_bin)
sol = integrate(y0, m, (0.0, 0.38), 260, rtol=1e-11, atol=1e-13)
c = classify(sol, m)
print("binary bound at end:", c["binary"]["bound"], " sep=%.3f AU"%c["binary"]["sep"])
Y = sol.y; tt = sol.t
smbh=Y[0:2].T; bh=Y[2:4].T; star=Y[4:6].T
# SMBH-centric
bh_c = bh-smbh; star_c = star-smbh

os.makedirs("scenes", exist_ok=True)
with open("scenes/flyby.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["t","smbh_x","smbh_y","bh_x","bh_y","star_x","star_y"])
    for i in range(len(tt)):
        w.writerow([f"{v:.6g}" for v in (tt[i],0.0,0.0,bh_c[i,0],bh_c[i,1],star_c[i,0],star_c[i,1])])
# pericenter frame (closest CoM approach) for the event label
d = np.linalg.norm(bh_c, axis=1); pf = int(np.argmin(d))
meta={"center":"smbh","zoom_on":["bh"],
      "events":[{"frame":pf,"label":"pericenter"}],
      "inset":{"center":"com","com":{"bh":10,"star":1},"show":["bh","star"],"label":"binary (zoom)"},
      "notes":f"a_bin={a_bin} e_bin={e_bin}, rp_out={rp_out:.1f}AU (beta~0.5)"}
json.dump(meta, open("scenes/flyby.events.json","w"), indent=1)
print(f"wrote scenes/flyby.csv ({len(tt)} frames), pericenter frame {pf}")
