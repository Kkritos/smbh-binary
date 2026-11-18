#!/usr/bin/env python3
"""Stitch canonical trajectory CSVs and resample them in WARPED time so playback
accelerates: dense (slow) early frames resolve the binary disruption, then time
speeds up to show the bodies flying apart.  Used by the IONIZATION scene
(keep integrating past the breakup + time-lapse the separation).

   python3 timewarp.py out.csv in1.csv [in2.csv ...] [--nframes 340 --accel 2.6]

Each input is a canonical CSV (t,smbh_x,smbh_y,bh_x,bh_y,star_x,star_y).  Successive
files are assumed continuous (each starts from the previous one's final state, as a
--restart does) and are stitched with a cumulative time offset.  Frame times are

   t_k = t0 + (tN - t0) * (k/(N-1))^accel          (accel > 1 => slow start, fast finish)

and every column is interpolated there with shape-preserving PCHIP (SciPy if present).
NaN bodies stay NaN.  The caller usually writes the events.json afterwards.
"""
import argparse, csv, os, numpy as np


def load_canon(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    data = np.array([[float(x) for x in r] for r in rows[1:]], dtype=float)
    return head, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("inp", nargs="+")
    ap.add_argument("--nframes", type=int, default=340)
    ap.add_argument("--accel", type=float, default=2.6,
                    help=">1 slows the start (disruption) and speeds up the finish (fly-apart)")
    a = ap.parse_args()

    head, d0 = load_canon(a.inp[0])
    ti = head.index("t")
    datas = [d0]
    toff = d0[-1, ti]
    for p in a.inp[1:]:
        h, d = load_canon(p)
        if h != head:
            raise SystemExit(f"column mismatch between {a.inp[0]} and {p}")
        d = d.copy()
        d[:, ti] = d[:, ti] - d[0, ti] + toff          # continue after the previous file
        toff = d[-1, ti]
        datas.append(d)
    data = np.vstack(datas)

    order = np.argsort(data[:, ti], kind="stable"); data = data[order]
    keep = np.concatenate(([True], np.diff(data[:, ti]) > 0)); data = data[keep]
    t = data[:, ti]
    if t.size < 2:
        raise SystemExit("timewarp: <2 distinct time samples")
    t0, tN = t[0], t[-1]
    u = np.linspace(0.0, 1.0, a.nframes) ** a.accel
    tq = t0 + (tN - t0) * u

    try:
        from scipy.interpolate import PchipInterpolator
        def ip(x, y, xq): return PchipInterpolator(x, y, extrapolate=False)(xq)
    except Exception:
        def ip(x, y, xq): return np.interp(xq, x, y, left=np.nan, right=np.nan)

    out = np.full((a.nframes, len(head)), np.nan)
    out[:, ti] = tq
    for j in range(len(head)):
        if j == ti:
            continue
        y = data[:, j]; good = np.isfinite(y)
        if good.sum() >= 2:
            out[:, j] = ip(t[good], y[good], tq)

    with open(a.out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(head)
        for r in out:
            w.writerow([f"{v:.6g}" if np.isfinite(v) else "nan" for v in r])
    print(f"wrote {a.out}  ({t.size} src rows -> {a.nframes} warped frames, accel={a.accel})")


if __name__ == "__main__":
    main()
