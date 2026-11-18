"""micro-TDE + debris (PREVIEW), two-scale with debris FATES.
The SMBH encounter drives the star to r_t around its BH companion -> disruption.
Debris are followed in the SMBH+BH field (softened to tame the stiff inner orbits)
and classified by fate: bound to the stellar BH (fallback), captured by the SMBH,
or ejected. Main panel = SMBH-centric (the three populations spreading); inset =
BH-centric zoom (the disruption + fallback onto the stellar BH)."""
import json, os, csv, numpy as np
from scipy.integrate import solve_ivp
from nbody import make_ic, integrate, G, RSUN_AU, spawn_debris

M_smbh, m_bh, m_star = 1e6, 10.0, 1.0
R_star = 1.0*RSUN_AU
r_t = R_star*(m_bh/m_star)**(1/3)
a_bin, e_bin = 0.1, 0.85; m_b = m_bh+m_star
r_t_bin = a_bin*(M_smbh/m_b)**(1/3)
# BOUND outer orbit: a nuclear-cluster binary on an eccentric orbit about the SMBH,
# disrupted at its (penetrating) pericenter -> debris mostly stay bound to the SMBH.
rp_out = 0.5*r_t_bin; e_out = 0.95
a_out = rp_out/(1-e_out)
nu0 = np.pi                              # start at apocenter, fall in to pericenter

def disrupt_event(t, y, m, soft):
    pos = y[:6].reshape(3,2); return np.linalg.norm(pos[2]-pos[1]) - r_t
disrupt_event.terminal = True; disrupt_event.direction = -1

chosen = None
for phase in np.linspace(0, 2*np.pi, 32, endpoint=False):
    y0, m = make_ic(M_smbh, m_bh, m_star, a_bin, rp_out, e_out, nu0, phase, e_bin=e_bin)
    sol = integrate(y0, m, (0.0, 0.15), 120, events=disrupt_event, rtol=1e-10, atol=1e-13)
    if sol.t_events[0].size > 0: chosen = (phase, sol.t_events[0][0]); break
if chosen is None: raise SystemExit("no disruption found")
phase, t_dis = chosen
print(f"micro-TDE at t={t_dis:.4f} yr, phase={np.degrees(phase):.0f} deg")

# stage 1 to disruption
y0, m = make_ic(M_smbh, m_bh, m_star, a_bin, rp_out, e_out, nu0, phase, e_bin=e_bin)
sol = integrate(y0, m, (0.0, t_dis*0.999), 170, rtol=1e-11, atol=1e-13)
Y=sol.y; tt=sol.t
smbh=Y[0:2].T; bh=Y[2:4].T; star=Y[4:6].T
vsm=Y[6:8].T; vbh=Y[8:10].T; vst=Y[10:12].T

# spawn debris
dpos, dvel = spawn_debris(star[-1], vst[-1], bh[-1], vbh[-1], m_bh, n=160, spread=0.5)
Nd = len(dpos)

# stage 2: velocity-Verlet leapfrog for SMBH+BH (2-body) + debris test particles,
# lightly softened, with ACCRETION FREEZING to remove the stiff inner orbits: a debris
# reaching the BH (<2 r_t) or deep near the SMBH is frozen onto (co-moves with) that
# body ("accreted"). Lets us follow several SMBH orbits cheaply without tiny timesteps.
Rg = G*M_smbh/(63241.0**2)
SB = 0.3*r_t; SS = 3.0*Rg                      # softening (keeps BH binding tight)
R_ACC_BH = 1.2*r_t; R_ACC_SM = 0.25            # accretion radii (AU)
t2 = 0.9; dt = 3e-5; nsteps = int(t2/dt); n2 = 340
p_sm, v_sm = smbh[-1].copy(), vsm[-1].copy()
p_bh, v_bh = bh[-1].copy(), vbh[-1].copy()
dp = dpos.copy(); dv = dvel.copy()
frozen = np.zeros(Nd, bool); acc_of = np.zeros(Nd, int); off = np.zeros((Nd, 2)); has_left = np.zeros(Nd, bool)

def acc_massive(ps, pb):
    d = pb-ps; r3 = (d@d)**1.5; return G*m_bh*d/r3, -G*M_smbh*d/r3
def acc_deb(dq, ps, pb):
    rs = ps[None,:]-dq; rb = pb[None,:]-dq
    rs3 = (np.sum(rs*rs,1)+SS*SS)**1.5; rb3 = (np.sum(rb*rb,1)+SB*SB)**1.5
    return G*M_smbh*rs/rs3[:,None] + G*m_bh*rb/rb3[:,None]

rec = set(np.linspace(0, nsteps, n2).astype(int).tolist())
smbh2=[]; bh2=[]; deb2=[]; trec=[]
a_sm, a_bh = acc_massive(p_sm, p_bh); a_d = acc_deb(dp, p_sm, p_bh)
for step in range(nsteps+1):
    if step in rec:
        smbh2.append(p_sm.copy()); bh2.append(p_bh.copy()); deb2.append(dp.copy()); trec.append(step*dt)
    fr = ~frozen
    v_sm += 0.5*dt*a_sm; v_bh += 0.5*dt*a_bh; dv[fr] += 0.5*dt*a_d[fr]
    p_sm += dt*v_sm; p_bh += dt*v_bh; dp[fr] += dt*dv[fr]
    fb = frozen & (acc_of==0); fs = frozen & (acc_of==1)
    dp[fb] = p_bh + off[fb]; dp[fs] = p_sm + off[fs]
    a_sm, a_bh = acc_massive(p_sm, p_bh); a_d = acc_deb(dp, p_sm, p_bh)
    v_sm += 0.5*dt*a_sm; v_bh += 0.5*dt*a_bh; dv[fr] += 0.5*dt*a_d[fr]
    dbh = np.linalg.norm(dp-p_bh[None,:], axis=1); dsm = np.linalg.norm(dp-p_sm[None,:], axis=1)
    has_left |= (dbh > 2.5*r_t)
    nb = (~frozen)&has_left&(dbh<R_ACC_BH); ns = (~frozen)&(dsm<R_ACC_SM)
    off[nb] = dp[nb]-p_bh; acc_of[nb] = 0; frozen |= nb
    off[ns] = dp[ns]-p_sm; acc_of[ns] = 1; frozen |= ns
smbh2=np.array(smbh2); bh2=np.array(bh2); deb2=np.array(deb2); trec=np.array(trec); N2=len(trec)

# classify fate (accreted -> by accretor; else by final binding energy)
a_bh_orbit=np.linalg.norm(p_bh-p_sm); r_hill=a_bh_orbit*(m_bh/(3*M_smbh))**(1/3)
fate=[]
for k in range(Nd):
    if frozen[k]:
        fate.append("debB" if acc_of[k]==0 else "debS"); continue
    dbh=np.linalg.norm(dp[k]-p_bh); dsm=np.linalg.norm(dp[k]-p_sm)
    Ebh=0.5*np.sum((dv[k]-v_bh)**2)-G*m_bh/(dbh+SB)
    Esm=0.5*np.sum((dv[k]-v_sm)**2)-G*M_smbh/(dsm+SS)
    if Ebh<0 and dbh<r_hill:   fate.append("debB")
    elif Esm<0:                fate.append("debS")
    else:                      fate.append("debE")
import collections
print("fates:", dict(collections.Counter(fate)), " r_hill=%.4f AU  frozen=%d"%(r_hill, frozen.sum()))

# ---------- canonical (debris named by fate) ----------
N1=len(tt)
cols=["t","smbh_x","smbh_y","bh_x","bh_y","star_x","star_y"]+\
     [f"{fate[k]}{k:03d}_{ax}" for k in range(Nd) for ax in ("x","y")]
os.makedirs("scenes",exist_ok=True)
with open("scenes/mutde.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(cols)
    for i in range(N1):
        row=[tt[i],smbh[i,0],smbh[i,1],bh[i,0],bh[i,1],star[i,0],star[i,1]]+["nan"]*(2*Nd)
        w.writerow([f"{v:.6g}" if not isinstance(v,str) and np.isfinite(v) else (v if isinstance(v,str) else "nan") for v in row])
    for i in range(N2):
        base=[t_dis+trec[i],smbh2[i,0],smbh2[i,1],bh2[i,0],bh2[i,1],"nan","nan"]
        d=[]
        for k in range(Nd): d+=[f"{deb2[i,k,0]:.6g}",f"{deb2[i,k,1]:.6g}"]
        w.writerow([f"{base[0]:.6g}"]+[f"{v:.6g}" if not isinstance(v,str) else v for v in base[1:]]+d)

meta={"center":"smbh","zoom_on":["bh","star","deb"],"zoom_pct":82,
      "events":[{"frame":N1-1,"label":"micro-TDE: star disrupted by BH"}],
      "inset":{"center":"com","com":{"bh":10,"star":1},"show":["bh","star","deb"],"hw":0.25,"label":"stellar BH (zoom)"},
      "notes":f"debris fates bound-BH/SMBH/ejected; r_t={r_t:.4f}AU"}
json.dump(meta,open("scenes/mutde.events.json","w"),indent=1)
print(f"wrote scenes/mutde.csv ({N1}+{N2} frames, {Nd} debris)")
