#!/usr/bin/env python3
"""reduce_delay.py -- survey the two-peak fallback delay Dt_delay for EVERY micro-TDE in a bundle,
WITHOUT integrating orbits.  For each mu-TDE we seed N polytrope debris at the disruption state,
take the frozen-in energies about BH1 and BH2, and read the two peak fallback times off the
high-energy edges of dM/deps:   t_peak = 2pi M_c (2|eps|_edge)^{-3/2}.

    Dt_delay   = t_peak,BH2 - t_peak,BH1      [yr]
    ratio      = t_peak,BH2 / t_peak,BH1

Run (on the full bundle, e.g. SeaWulf):
    module load python/3.12.0
    python3 reduce_delay.py smbh_bhstar_bundle.h5            # -> delay_per_event.csv + delay_cellsummary.csv
Needs debris_tracer.py on the path.
"""
import sys, os, time, functools, numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import debris_tracer as dt

BUNDLE  = sys.argv[1] if len(sys.argv) > 1 else "smbh_bhstar_bundle.h5"
N       = int(sys.argv[2]) if len(sys.argv) > 2 else 500      # debris per event (edge needs only ~hundreds)
EDGE_PCT= 90.0                                                 # "most-bound edge" percentile of |eps|
POLY_N  = 3.0
V1, YR, RSUN_AU = dt.V1_KMS, dt.YR_PER_T, dt.RSUN_AU
FRAC_NODES = np.array([1e-3,1e-2,1e-1,3e-1,1.0,3.0,10.0]); M_BIN=11.0

dt.lane_emden = functools.lru_cache(maxsize=8)(dt.lane_emden)   # cache the Lane-Emden solve

# ---- load everything once ----
POS=lambda p:(f"x_{p}_AU",f"y_{p}_AU",f"z_{p}_AU"); VEL=lambda p:(f"vx_{p}_kms",f"vy_{p}_kms",f"vz_{p}_kms")
with h5py.File(BUNDLE,"r") as f:
    R=f["results"]; oc=R["outcome"].asstr()[:]
    mu=np.where(oc=="TDE:star2_by_bh1")[0]
    cols={}
    for k in (POS("smbh")+VEL("smbh")+POS("b1")+VEL("b1")+POS("b2")+VEL("b2")
              +("M_smbh","m1","m2","R2","a_bin","e_bin")):
        cols[k]=R[k][:]
    cols["v_inf"]=R["v_inf"][:] if "v_inf" in R else np.full(len(oc), np.nan)  # may be absent in a test bundle
print(f"{len(mu)} micro-TDE events in {BUNDLE};  seeding N={N} each")

def row_state(i):
    g=lambda keys: np.array([cols[k][i] for k in keys], float)
    return dict(M2=float(cols["M_smbh"][i]), m_BH1=float(cols["m1"][i]), m_star=float(cols["m2"][i]),
                R_star=float(cols["R2"][i])*RSUN_AU,
                r_BH2=g(POS("smbh")), v_BH2=g(VEL("smbh"))/V1,
                r_BH1=g(POS("b1")),   v_BH1=g(VEL("b1"))/V1,
                r_star=g(POS("b2")),  v_star=g(VEL("b2"))/V1,
                r_t_star=np.linalg.norm(g(POS("b2"))-g(POS("b1"))))

def peaks(st):
    r,v,om,info = dt.make_debris_ic(st, N, n=POLY_N, spin=True, seed=1)
    d1=r-st["r_BH1"]; eps1=0.5*np.einsum("ij,ij->i",v-st["v_BH1"],v-st["v_BH1"])-st["m_BH1"]/np.linalg.norm(d1,axis=1)
    d2=r-st["r_BH2"]; eps2=0.5*np.einsum("ij,ij->i",v-st["v_BH2"],v-st["v_BH2"])-st["M2"]/np.linalg.norm(d2,axis=1)
    b1=eps1<0; b2=(~b1)&(eps2<0)
    out=dict(f_BH1=b1.mean(), f_BH2=b2.mean(), f_unb=((~b1)&(eps2>=0)).mean(), beta=info["beta_pen"])
    if b1.sum()>=5 and b2.sum()>=5:
        e1=np.percentile(np.abs(eps1[b1]),EDGE_PCT); e2=np.percentile(np.abs(eps2[b2]),EDGE_PCT)
        t1=2*np.pi*st["m_BH1"]*(2*e1)**-1.5*YR; t2=2*np.pi*st["M2"]*(2*e2)**-1.5*YR
        out.update(t_peak_BH1=t1, t_peak_BH2=t2, Dt_delay=t2-t1, ratio=t2/t1)
    else:
        out.update(t_peak_BH1=np.nan, t_peak_BH2=np.nan, Dt_delay=np.nan, ratio=np.nan)
    return out

T0=time.time(); rows=[]
for j,i in enumerate(mu):
    st=row_state(i); p=peaks(st)
    if np.isfinite(cols["v_inf"][i]):
        vorb=V1*np.sqrt(M_BIN/cols["a_bin"][i]); frac=FRAC_NODES[np.argmin(np.abs(np.log10(cols["v_inf"][i]/vorb)-np.log10(FRAC_NODES)))]
    else:
        frac=np.nan
    rows.append((cols["M_smbh"][i],cols["a_bin"][i],cols["e_bin"][i],frac,cols["v_inf"][i],
                 p["t_peak_BH1"],p["t_peak_BH2"],p["Dt_delay"],p["ratio"],p["f_BH1"],p["f_BH2"],p["f_unb"]))
    if (j+1)%2000==0: print(f"  {j+1}/{len(mu)}  ({(time.time()-T0)/(j+1)*1e3:.1f} ms/event)")
rows=np.array(rows)
hdr="M_smbh,a_bin,e_bin,frac,v_inf,t_peak_BH1,t_peak_BH2,Dt_delay,ratio,f_BH1,f_BH2,f_unb"
np.savetxt("delay_per_event.csv", rows, delimiter=",", header=hdr, comments="", fmt="%.6g")

# ---- per-cell summary (median over each M,a,e,frac cell) ----
import csv
key=lambda r:(round(np.log10(r[0]),4),round(r[1],5),round(r[2],4),round(np.log10(r[3]),4))
from collections import defaultdict
cells=defaultdict(list)
for r in rows: cells[key(r)].append(r)
with open("delay_cellsummary.csv","w",newline="") as g:
    w=csv.writer(g); w.writerow(["M_smbh","a_bin","e_bin","frac","N","t_peak_BH1_med","t_peak_BH2_med",
                                 "Dt_delay_med","Dt_delay_16","Dt_delay_84","ratio_med"])
    for k,v in cells.items():
        v=np.array(v); dd=v[:,7][np.isfinite(v[:,7])]
        if dd.size==0: continue
        w.writerow([v[0,0],v[0,1],v[0,2],v[0,3],len(v),
                    np.nanmedian(v[:,5]),np.nanmedian(v[:,6]),
                    np.median(dd),np.percentile(dd,16),np.percentile(dd,84),np.nanmedian(v[:,8])])
dd=rows[:,7][np.isfinite(rows[:,7])]
print(f"done in {time.time()-T0:.0f}s.  Dt_delay median={np.median(dd):.3e} yr, 16-84=[{np.percentile(dd,16):.2e},{np.percentile(dd,84):.2e}]")
print(f"wrote delay_per_event.csv ({len(rows)} rows) and delay_cellsummary.csv")
