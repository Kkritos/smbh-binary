"""Load a canonical trajectory CSV and render a looping GIF (+ mp4).

Canonical CSV (one row per frame):   t,<body>_x,<body>_y,...
Body names drive styling by role: 'smbh', 'bh', 'star'.  Missing bodies are NaN.
A sidecar <name>.events.json may hold {"center","zoom_on","zoom_pct","inset","debris"}.

Debris (micro-TDE) are carried in a compact sidecar <name>.debris.npz referenced by
meta["debris"]["npz"]: arrays xy (Nfollow, Ndeb, 2), energy (Ndeb,), start_frame (int).
They are drawn as one scatter coloured by specific energy wrt the SMBH, with a colorbar.

Design: fixed encounter-type title (top-center), a white scale bar + label in BOTH the
main window and the inset, the inset auto-placed in the emptiest corner, NO transient
event text.
"""
import json, os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from matplotlib.colors import Normalize

# ---- role styling (dark "space" theme; colorblind-safe accents) -------------
ROLE = {
    "smbh": dict(color="#000000", edge="#D55E00", size=320, z=5),
    "bh":   dict(color="#0b1021", edge="#56B4E9", size=70,  z=4),
    "star": dict(color="#F0E442", edge="#E69F00", size=55,  z=4),
    "deb":  dict(color="#E69F00", edge=None,      size=4,   z=3),  # demo CSV debris
}
TRAIL = {"bh": "#56B4E9", "star": "#E69F00", "smbh": None}
CORNERS = ("TL", "TR", "BL", "BR")


def role_of(name):
    if name.startswith("deb"):
        return "deb"
    return name if name in ROLE else "star"


def load(csv_path):
    import csv as _csv
    with open(csv_path) as f:
        rows = list(_csv.reader(f))
    head = rows[0]; data = np.array(rows[1:], dtype=float)
    t = data[:, head.index("t")]
    bodies = {}
    for j, col in enumerate(head):
        if col.endswith("_x"):
            name = col[:-2]
            xy = np.column_stack([data[:, j], data[:, head.index(name+"_y")]])
            bodies[name] = xy
    try:
        meta = json.load(open(csv_path.rsplit(".csv", 1)[0] + ".events.json"))
    except OSError:
        meta = {}
    debris = None
    dcfg = meta.get("debris")
    if dcfg and dcfg.get("npz"):
        npz = os.path.join(os.path.dirname(csv_path), dcfg["npz"])
        try:
            z = np.load(npz)
            debris = dict(xy=z["xy"], energy=z["energy"], start=int(z["start_frame"]),
                          cbar_label=dcfg.get("cbar_label", "specific energy wrt SMBH"))
        except OSError:
            debris = None
    return dict(t=t, bodies=bodies, meta=meta, debris=debris)


def _ffmpeg_path():
    """Locate an ffmpeg binary: system first, else the pip-bundled imageio-ffmpeg."""
    import shutil
    p = shutil.which("ffmpeg")
    if p:
        return p
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()      # bundled static binary
    except Exception:
        return None


def _nice_len(x):
    """Round x down to a 1/2/5 x 10^k 'nice' number for a scale bar."""
    if x <= 0 or not np.isfinite(x):
        return 1.0
    k = np.floor(np.log10(x)); base = x / 10**k
    nice = 1 if base < 2 else (2 if base < 5 else 5)
    return nice * 10**k


def _fmt_len(L):
    if L >= 1:   return f"{L:g} AU"
    if L >= 1/23455: return f"{L:.2g} AU"          # still AU
    return f"{L*23455:.2g} R$_\\odot$"             # (not normally hit)


class _ScaleBar:
    """A white horizontal scale bar + label drawn in axes-fraction coords."""
    def __init__(self, ax, x0, y0, ha="left", fs=9):
        self.ax, self.x0, self.y0, self.ha, self.fs = ax, x0, y0, ha, fs
        self.line = Line2D([], [], color="white", lw=2.0, transform=ax.transAxes,
                           solid_capstyle="butt", zorder=12)
        ax.add_line(self.line)
        self.txt = ax.text(x0, y0+0.028, "", transform=ax.transAxes, color="white",
                           fontsize=fs, ha=ha, va="bottom", zorder=12)

    def update(self, half_width):
        L = _nice_len(0.5 * half_width)                 # ~25% of the full view
        frac = min(0.45, L / (2.0 * half_width))
        if self.ha == "left":
            xa, xb, xt = self.x0, self.x0+frac, self.x0+frac/2
        else:
            xa, xb, xt = self.x0-frac, self.x0, self.x0-frac/2
        self.line.set_data([xa, xb], [self.y0, self.y0])
        self.txt.set_position((xt, self.y0+0.028)); self.txt.set_ha("center")
        self.txt.set_text(_fmt_len(L))


def _occupancy(centers, hws, tracks):
    """Fraction of tracked-body samples landing in each corner of the main view."""
    cnt = {c: 0 for c in CORNERS}; tot = 0
    for i in range(len(hws)):
        c = centers[i]; h = hws[i]
        for xy in tracks:
            p = xy[i]
            if not np.isfinite(p[0]):
                continue
            u = (p[0]-c[0])/h; v = (p[1]-c[1])/h
            if abs(u) > 1 or abs(v) > 1:
                continue
            key = ("T" if v >= 0 else "B") + ("R" if u >= 0 else "L")
            cnt[key] += 1; tot += 1
    if tot == 0:
        return {c: 0.0 for c in CORNERS}
    return {c: cnt[c]/tot for c in CORNERS}


def _inset_boxes(pb):
    """Four corner inset rectangles that sit INSIDE the plot box pb = [x0,y0,w,h]."""
    x0, y0, w, h = pb
    iw = 0.33 if w >= 0.99 else 0.30
    m = 0.012
    return {"TL": [x0+m, y0+h-m-iw, iw, iw], "TR": [x0+w-m-iw, y0+h-m-iw, iw, iw],
            "BL": [x0+m, y0+m, iw, iw],       "BR": [x0+w-m-iw, y0+m, iw, iw]}


def render(scene, out_gif, title="", fps=24, trail=40, figsize=6.2,
           center="smbh", zoom_on=("bh",), hw_floor=1.0, hw_pad=1.18,
           fontsize=12, inset=False, inset_center="pair", inset_show=None,
           inset_hw=None, inset_label="binary (zoom)", inset_com=None,
           zoom_pct=100, mp4=True):
    t = scene["t"]; B = scene["bodies"]; meta = scene.get("meta", {})
    debris = scene.get("debris")
    center = meta.get("center", center)
    zoom_on = tuple(meta.get("zoom_on", zoom_on))
    zoom_pct = meta.get("zoom_pct", zoom_pct)
    nf = len(t)
    icfg = meta.get("inset")
    if icfg:
        inset = True
        inset_center = icfg.get("center", inset_center)
        inset_show = icfg.get("show", inset_show)
        inset_hw = icfg.get("hw", inset_hw)
        inset_label = icfg.get("label", inset_label)
        inset_com = icfg.get("com", inset_com)

    def ctr(i):
        return B[center][i] if center in B else np.array([0.0, 0.0])

    # ---- main-view half-width: smoothed, instant-expand / slow-shrink -----------
    zoom_names = [nm for nm in zoom_on if nm in B and nm != "debris"]
    use_deb_zoom = ("debris" in zoom_on) and (debris is not None)
    centers = np.array([ctr(i) for i in range(nf)])
    raw = np.full(nf, hw_floor)
    for i in range(nf):
        ds = []
        for nm in zoom_names:                            # point bodies (e.g. the BH)
            d = np.linalg.norm(B[nm][i] - centers[i])
            if np.isfinite(d): ds.append(d)
        if use_deb_zoom and debris is not None:
            j = i - debris["start"]
            if 0 <= j < debris["xy"].shape[0]:
                dd = np.linalg.norm(debris["xy"][j] - centers[i], axis=1)
                ds.extend(dd[np.isfinite(dd)].tolist())
        # one percentile over the POOLED distances: tracks the BH in during the approach,
        # frames the debris bulk during fallback, and ignores the escaping BH / far ejecta
        raw[i] = max(hw_floor, np.percentile(ds, zoom_pct)) if ds else hw_floor
    raw *= hw_pad
    hw = np.copy(raw)
    for i in range(1, nf):
        hw[i] = max(raw[i], 0.99*hw[i-1])

    # ---- occupancy of the frame (drives inset corner AND colorbar side) ---------
    tracks = [B[nm] for nm in zoom_names]
    if use_deb_zoom and debris is not None:
        pad = np.full((nf, 2), np.nan)
        pad[debris["start"]:debris["start"]+debris["xy"].shape[0]] = \
            np.nanmean(debris["xy"], axis=1)
        tracks.append(pad)
    # fixed layout (per request): inset lower-LEFT, title top-center, time lower-CENTER,
    # scale bar lower-RIGHT — regardless of where the orbits/debris happen to be.
    inset_corner = "BL"

    # ---- figure -----------------------------------------------------------------
    # debris scenes reserve a right margin for the colorbar, so it is OUTSIDE the data
    # region and can never overlap the orbits/debris (which can recede either way).
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})
    PB = [0.0, 0.09, 0.82, 0.82] if debris is not None else [0.0, 0.0, 1.0, 1.0]
    INSET = _inset_boxes(PB)
    fig = plt.figure(figsize=(figsize, figsize), facecolor="black")
    ax = fig.add_axes(PB); ax.set_facecolor("black"); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])

    deb_csv = [nm for nm in B if nm.startswith("deb")]      # demo path (CSV debris)
    scat = {}; trails = {}
    for nm, xy in B.items():
        if nm in deb_csv:
            continue
        r = role_of(nm); st = ROLE[r]
        scat[nm] = ax.scatter([], [], s=st["size"], c=st["color"],
                              edgecolors=st["edge"], linewidths=1.4, zorder=st["z"])
        tl = TRAIL.get(r)
        if tl:
            (trails[nm],) = ax.plot([], [], "-", color=tl, lw=1.4, alpha=0.5,
                                    zorder=st["z"]-0.5)
    scat_debcsv = ax.scatter([], [], s=4, c=ROLE["deb"]["color"], alpha=0.6,
                             edgecolors="none", zorder=3) if deb_csv else None
    glow_smbh = ax.scatter([], [], s=1600, c="#D55E00", alpha=0.10,
                           edgecolors="none", zorder=1)

    # ---- debris scatter + colorbar ---------------------------------------------
    deb_scat = None; deb_colors = None
    if debris is not None:
        e = debris["energy"].astype(float)
        emax = np.percentile(np.abs(e[np.isfinite(e)]), 95) or 1.0   # normalization scale
        en = np.clip(e/emax, -1.0, 1.0)              # O(1): -1 bound ... +1 ejected
        norm = Normalize(vmin=-1.0, vmax=1.0)
        try:
            cmap = matplotlib.colormaps["coolwarm"]       # mpl >= 3.6
        except (AttributeError, KeyError):
            cmap = cm.get_cmap("coolwarm")                # older mpl
        deb_colors = cmap(norm(en))
        deb_scat = ax.scatter([], [], s=3.0, edgecolors="none", zorder=3)
        cax = fig.add_axes([0.855, 0.30, 0.02, 0.40])    # right margin, outside the data box
        sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
        cb.ax.yaxis.set_ticks_position("left")           # ticks/label toward the plot
        cb.ax.yaxis.set_label_position("left")
        cb.ax.tick_params(colors="0.85", labelsize=fontsize-3)
        cb.outline.set_edgecolor("0.4")
        cb.set_label(r"$\epsilon/\epsilon_{\max}$", color="0.9", fontsize=fontsize+2, labelpad=4)

    # ---- overlay text: fixed positions -----------------------------------------
    # title top-center (inset is lower-left, so the top is always clear), time
    # lower-center, scale bar lower-right.
    ax.text(0.5, 0.965, title, transform=ax.transAxes, color="white",
            fontsize=fontsize+1, ha="center", va="top", zorder=12)
    txt_t = ax.text(0.5, 0.045, "", transform=ax.transAxes, color="0.8",
                    fontsize=fontsize-1, ha="center", va="bottom", zorder=12)  # lower-center
    sbar = _ScaleBar(ax, 0.97, 0.055, ha="right", fs=fontsize-2)               # lower-right

    # ---- inset ------------------------------------------------------------------
    axi = None
    if inset:
        ishow = list(inset_show) if inset_show is not None else ["bh", "star"]
        ipt = [nm for nm in ishow if nm in B]
        show_deb_inset = ("debris" in ishow) and (debris is not None)
        if inset_hw is not None:
            ihw = inset_hw
        elif len(ipt) >= 2:                              # auto-size to the binary
            seps = np.linalg.norm(B[ipt[0]] - B[ipt[1]], axis=1)
            fs = seps[np.isfinite(seps)]
            # robust to a body escaping (use the typical separation, not the max)
            ihw = 1.2*np.percentile(fs, 75) if fs.size else hw_floor
        else:
            ihw = hw_floor
        axi = fig.add_axes(INSET[inset_corner]); axi.set_facecolor("#05070f")
        axi.set_xticks([]); axi.set_yticks([]); axi.set_aspect("equal")
        for sp in axi.spines.values(): sp.set_color("0.35")
        axi.text(0.5, 0.965, inset_label, transform=axi.transAxes, color="0.82",
                 ha="center", va="top", fontsize=fontsize-3, zorder=13)   # inside (never clipped)
        iscat = {nm: axi.scatter([], [], s=ROLE[role_of(nm)]["size"]*0.9,
                                 c=ROLE[role_of(nm)]["color"], edgecolors=ROLE[role_of(nm)]["edge"],
                                 linewidths=1.1, zorder=4) for nm in ipt}
        itrail = {nm: axi.plot([], [], "-", color=TRAIL.get(role_of(nm), "0.5"),
                               lw=1.0, alpha=0.5)[0] for nm in ipt if TRAIL.get(role_of(nm))}
        idebscat = axi.scatter([], [], s=5.0, edgecolors="none", zorder=3) if show_deb_inset else None
        isbar = _ScaleBar(axi, 0.06, 0.07, ha="left", fs=fontsize-4)

        def inset_center_at(i):
            if inset_center == "com":
                mm = inset_com or {nm: 1.0 for nm in ipt}
                num = np.zeros(2); tot = 0.0
                for nm, w in mm.items():
                    if nm in B and np.isfinite(B[nm][i][0]): num += w*B[nm][i]; tot += w
                if tot > 0: return num/tot
            if inset_center == "pair" and len(ipt) >= 2 \
               and np.isfinite(B[ipt[0]][i][0]) and np.isfinite(B[ipt[1]][i][0]):
                return 0.5*(B[ipt[0]][i] + B[ipt[1]][i])
            if inset_center in B and np.isfinite(B[inset_center][i][0]):
                return B[inset_center][i]
            for nm in ipt:
                if np.isfinite(B[nm][i][0]): return B[nm][i]
            return np.array([0.0, 0.0])

    def deb_at(i):
        if debris is None: return None
        j = i - debris["start"]
        if 0 <= j < debris["xy"].shape[0]:
            return debris["xy"][j]
        return None

    def frame(i):
        c = ctr(i)
        ax.set_xlim(c[0]-hw[i], c[0]+hw[i]); ax.set_ylim(c[1]-hw[i], c[1]+hw[i])
        for nm, xy in B.items():
            if nm in deb_csv:
                continue
            p = xy[i]
            if np.isfinite(p[0]):
                scat[nm].set_offsets([p])
                if nm in trails:
                    lo = max(0, i-trail); seg = xy[lo:i+1]
                    trails[nm].set_data(seg[:, 0], seg[:, 1])
            else:
                scat[nm].set_offsets(np.empty((0, 2)))
                if nm in trails: trails[nm].set_data([], [])
        if scat_debcsv is not None:
            pts = [B[nm][i] for nm in deb_csv if np.isfinite(B[nm][i][0])]
            scat_debcsv.set_offsets(np.array(pts) if pts else np.empty((0, 2)))
        if deb_scat is not None:
            d = deb_at(i)
            if d is not None:
                m = np.isfinite(d[:, 0])
                deb_scat.set_offsets(d[m]); deb_scat.set_facecolors(deb_colors[m])
            else:
                deb_scat.set_offsets(np.empty((0, 2)))
        if "smbh" in B and np.isfinite(B["smbh"][i][0]):
            glow_smbh.set_offsets([B["smbh"][i]])
        txt_t.set_text(f"t = {t[i]:.4g} yr")
        sbar.update(hw[i])
        if axi is not None:
            # inset is meaningful only while there's a binary to zoom (or debris);
            # it auto-hides once the star is gone (e.g. after ejection).
            npt_fin = sum(np.isfinite(B[nm][i][0]) for nm in ipt)
            deb_here = (idebscat is not None) and (deb_at(i) is not None)
            active = deb_here or (npt_fin >= 2)
            axi.set_visible(active)
            if not active:
                return []
            mc = inset_center_at(i)
            axi.set_xlim(mc[0]-ihw, mc[0]+ihw); axi.set_ylim(mc[1]-ihw, mc[1]+ihw)
            for nm in ipt:
                p = B[nm][i]
                iscat[nm].set_offsets([p] if np.isfinite(p[0]) else np.empty((0, 2)))
                if nm in itrail:
                    lo = max(0, i-trail); seg = B[nm][lo:i+1]
                    itrail[nm].set_data(seg[:, 0], seg[:, 1])
            if idebscat is not None:
                d = deb_at(i)
                if d is not None:
                    m = np.isfinite(d[:, 0])
                    idebscat.set_offsets(d[m]); idebscat.set_facecolors(deb_colors[m])
                else:
                    idebscat.set_offsets(np.empty((0, 2)))
            isbar.update(ihw)
        return []

    anim = FuncAnimation(fig, frame, frames=nf, interval=1000/fps, blit=False)
    anim.save(out_gif, writer=PillowWriter(fps=fps))
    out_mp4 = None
    if mp4:
        fp = _ffmpeg_path()
        if fp:
            out_mp4 = out_gif.rsplit(".gif", 1)[0] + ".mp4"
            try:
                from matplotlib.animation import FFMpegWriter
                matplotlib.rcParams["animation.ffmpeg_path"] = fp
                anim.save(out_mp4, writer=FFMpegWriter(fps=fps, bitrate=2400))
            except Exception as ex:
                print(f"  (mp4 skipped: {ex})"); out_mp4 = None
        else:
            print("  (mp4 skipped: no ffmpeg; `pip install imageio-ffmpeg` "
                  "or `module load ffmpeg` to enable)")
    plt.close(fig)
    return out_gif, out_mp4
