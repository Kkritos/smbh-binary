#!/usr/bin/env python3
"""Concatenate stage-1 (N-body) and stage-2 (EMRI) canonical CSVs into one movie,
star set to NaN through stage 2. The seam is continuous because stage 2 started from
stage 1's exact final BH state (see handoff.py).

   python3 concat_stages.py stage1.csv stage2.csv out.csv
"""
import argparse, csv, json, numpy as np
import scene_lib

ap = argparse.ArgumentParser()
ap.add_argument("stage1"); ap.add_argument("stage2"); ap.add_argument("out")
ap.add_argument("--drop", default="star", choices=["star", "bh"],
                help="the EJECTED body to blank during stage 2 (so the inset can hide)")
A = ap.parse_args()
s1 = scene_lib.load(A.stage1); s2 = scene_lib.load(A.stage2); out = A.out
B1, B2 = s1["bodies"], s2["bodies"]; t1, t2 = s1["t"], s2["t"]
cols = ["t","smbh_x","smbh_y","bh_x","bh_y","star_x","star_y"]
rows = []
for i in range(len(t1)):
    rows.append([t1[i], *B1["smbh"][i], *B1["bh"][i], *B1["star"][i]])
t0 = t1[-1]
nan2 = (np.nan, np.nan)
for i in range(len(t2)):
    bh = nan2 if A.drop == "bh" else B2["bh"][i]
    st = nan2 if A.drop == "star" else B2["star"][i]
    rows.append([t0+(t2[i]-t2[0]), *B2["smbh"][i], *bh, *st])
with open(out,"w",newline="") as f:
    w = csv.writer(f); w.writerow(cols)
    for r in rows: w.writerow([f"{v:.6g}" if np.isfinite(v) else "nan" for v in r])
keep = "star" if A.drop == "bh" else "bh"           # the captured body, for the zoom
json.dump({"center":"smbh","zoom_on":[keep]},
          open(out.rsplit(".csv",1)[0]+".events.json","w"), indent=1)
print(f"wrote {out}  ({len(t1)}+{len(t2)} frames)")
