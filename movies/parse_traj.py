#!/usr/bin/env python3
"""Parse a smbh_binary_PN trajectory file and convert it to the canonical CSV the
animator reads.  Also used by make_restart.py.

   python3 parse_traj.py bss_..._full.txt canonical.csv [--plane xy]

WHERE THE TRAJECTORY COMES FROM
-------------------------------
`smbh_binary_PN --mode full` writes the trajectory to a FILE (not stdout) via
SpaceHub's DefaultWriter.  The file is named
    bss_<outer>_<type1><m1>_<type2><m2>_i<incl>_full.txt
and lands in the run's working directory.  (stdout only carries the human summary
and the one-line classify CSV; `--report 0|1` toggles only that stdout summary.)

FILE FORMAT (SpaceHub SizeParticles dump)
-----------------------------------------
A header line:
    time,id,mass,radius,px,py,pz,vx,vy,vz
then, per integrator output, ONE LINE PER BODY (ids 0,1,2 = SMBH, BH, star):
    <time>,<id>,<mass>,<radius>,<px>,<py>,<pz>,<vx>,<vy>,<vz>
so each time step is three consecutive rows.

UNITS  (SpaceHub internal, from macros.hpp: AU = 1, Ms = 1, year = 2*pi)
    px,py,pz : already AU        -> no conversion
    time     : internal          -> years = time / (2*pi)
    vx,vy,vz : internal speed     -> km/s  = v * (1 AU/yr = 4.7406 km/s) / (1/2pi)
               i.e. 1 internal speed unit = sqrt(G Msun / AU) = 29.7846 km/s
"""
import sys, argparse, numpy as np

PLANES = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}

TWO_PI       = 2.0 * np.pi
_KM_PER_AU   = 149597870.7
_SEC_PER_YR  = 365.25636042 * 86400.0
# 1 internal speed unit in km/s  (= AU per (internal time), internal time: yr = 2*pi)
V_INT_TO_KMS = _KM_PER_AU * TWO_PI / _SEC_PER_YR     # ~= 29.7846
AU_PER_YR_TO_KMS = 4.740572                          # 1 AU/yr in km/s (for finite-diff)

# Column indices in the SpaceHub dump line (0-based).
C_TIME, C_ID, C_PX, C_PZ, C_VX, C_VZ = 0, 1, 4, 6, 7, 9


def _has_header(path):
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            # numeric first field => no header; 'time'/letters => header
            try:
                float(s.split(",")[0].split()[0])
                return False
            except ValueError:
                return True
    return False


def load_rows(path):
    """Load the numeric table (comma- or whitespace-separated), skipping a header."""
    skip = 1 if _has_header(path) else 0
    # SpaceHub writes comma-separated; be tolerant of whitespace too.
    try:
        data = np.loadtxt(path, delimiter=",", skiprows=skip, comments="#")
    except Exception:
        data = np.loadtxt(path, skiprows=skip, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return data


def parse(path):
    """-> dict(t[years], pos[N,nb,3] AU, vel[N,nb,3] km/s).  Bodies 0=SMBH,1=BH,2=star.

    Pivots the long (one-row-per-body) dump into per-step arrays by grouping on id,
    which is robust to how many steps were written."""
    data = load_rows(path)
    if data.ndim != 2 or data.shape[0] < 2:
        raise SystemExit(f"could not read a numeric trajectory table from {path}")
    if data.shape[1] < 10:
        raise SystemExit(f"{path}: {data.shape[1]} columns; expected >=10 "
                         f"(time,id,mass,radius,px,py,pz,vx,vy,vz). First row:\n{data[0]}")

    ids = np.rint(data[:, C_ID]).astype(int)
    uids = np.unique(ids)
    nb = uids.size
    # rows for each body, in the order written (= time order)
    per = [data[ids == b] for b in sorted(uids)]
    ns = min(len(r) for r in per)             # guard against a truncated final step
    per = [r[:ns] for r in per]

    t = per[0][:, C_TIME] / TWO_PI                            # -> years
    pos = np.stack([r[:, C_PX:C_PX + 3] for r in per], axis=1)      # (ns, nb, 3) AU
    vel = np.stack([r[:, C_VX:C_VX + 3] for r in per], axis=1) * V_INT_TO_KMS  # km/s
    return dict(t=t, pos=pos, vel=vel)


def resample_uniform(t, pos, N):
    """Resample to N frames spaced evenly in the ANGLE SWEPT around the SMBH (body 0) by
    the orbiting body. Uniform-time sampling under-resolves the fast pericenter of an
    eccentric orbit (looks stroboscopic / jumps); sampling by swept angle puts frames
    where the body turns fastest, so pericenter AND apocenter are both smooth. Falls back
    to arc length (then time) when nothing orbits. Shape-preserving PCHIP if SciPy is
    present. NaN bodies (e.g. an ejected companion) stay NaN outside their span."""
    t = np.asarray(t, float)
    if pos.shape[0] < 2:
        return t, pos
    # AR-chain dumps can repeat or momentarily reverse the time stamp; sort + dedupe.
    order = np.argsort(t, kind="stable"); t = t[order]; pos = pos[order]
    keep = np.concatenate(([True], np.diff(t) > 0)); t = t[keep]; pos = pos[keep]
    if pos.shape[0] < 2:
        return t, pos
    nb = pos.shape[1]

    # ---- monotonic sampling parameter ----
    param = None; best_span = 0.5                       # need >~0.5 rad of winding to "orbit"
    c = pos[:, 0, :2]                                   # SMBH (x, y)
    for b in range(1, nb):
        rel = pos[:, b, :2] - c
        if not np.isfinite(rel).all():
            continue
        ang = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
        span = abs(ang[-1] - ang[0])
        if span > best_span:
            best_span = span
            dth = np.abs(np.diff(ang))                  # angle swept per step
            arc = np.linalg.norm(np.diff(pos[:, b, :2], axis=0), axis=1)
            w = dth + 0.05 * (dth.max() or 1.0) * (arc / (arc.max() + 1e-30))
            param = np.concatenate(([0.0], np.cumsum(w)))
    if param is None:                                   # nothing orbits -> arc length
        step = np.zeros(pos.shape[0])
        for b in range(nb):
            dd = np.linalg.norm(np.diff(pos[:, b, :], axis=0), axis=1)
            step[1:] = np.maximum(step[1:], np.nan_to_num(dd))
        param = np.cumsum(step)
    if param[-1] <= 0:
        param = t - t[0]
        if param[-1] <= 0:
            return t, pos

    keep2 = np.concatenate(([True], np.diff(param) > 0))
    param = param[keep2]; t = t[keep2]; pos = pos[keep2]
    if len(param) < 2:
        return t, pos
    pu = np.linspace(param[0], param[-1], N)
    try:
        from scipy.interpolate import PchipInterpolator
        def ip(x, y, xq): return PchipInterpolator(x, y, extrapolate=False)(xq)
    except Exception:
        def ip(x, y, xq): return np.interp(xq, x, y, left=np.nan, right=np.nan)
    tu = ip(param, t, pu)
    out = np.full((N,) + pos.shape[1:], np.nan)
    for b in range(nb):
        for d in range(pos.shape[2]):
            y = pos[:, b, d]; good = np.isfinite(y)
            if good.sum() >= 2:
                yq = ip(param[good], y[good], pu)
                yq[(pu < param[good][0]) | (pu > param[good][-1])] = np.nan
                out[:, b, d] = yq
    return tu, out


def to_canonical(path, out, plane="xy", max_frames=450):
    tr = parse(path); ix, iy = PLANES[plane]
    t, pos = resample_uniform(tr["t"], tr["pos"], max_frames)    # uniform time -> smooth
    cols = ["t", "smbh_x", "smbh_y", "bh_x", "bh_y", "star_x", "star_y"]
    M = np.column_stack([t,
        pos[:, 0, ix], pos[:, 0, iy], pos[:, 1, ix], pos[:, 1, iy],
        pos[:, 2, ix], pos[:, 2, iy]])
    np.savetxt(out, M, delimiter=",", header=",".join(cols), comments="", fmt="%.6g")
    print(f"wrote {out}  ({M.shape[0]} frames, plane={plane})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("raw"); ap.add_argument("out")
    ap.add_argument("--plane", default="xy", choices=list(PLANES))
    a = ap.parse_args()
    to_canonical(a.raw, a.out, a.plane)
