#!/usr/bin/env python3
# =============================================================================
#  style.py  --  shared plotting conventions for the whole mu-TDE analysis.
#  -----------------------------------------------------------------------------
#  Matches the conventions used in BOTH the debris-pipeline Figure 1 and the
#  smbh_movies render (scene_lib.py), so every figure and movie in the paper
#  reads as one set:
#    * Computer-Modern serif (cmr10 / mathtext 'cm'),
#    * ticks pointing IN on all four sides, with minor ticks everywhere
#      (including colorbars),
#    * no grid, one fixed box aspect,
#    * colour-blind-safe Okabe-Ito line palette + Wong, cividis for heatmaps,
#    * inside-axes panel titles in a white bbox,
#    * log axes: a major tick every decade, labels only every 2 decades.
# =============================================================================
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# numpy 1.x / 2.x trapezoid
trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

PANEL_ASPECT = 0.82

# Okabe-Ito (lines/markers) and the Wong palette alias (same 8 hues)
OI = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
          blue="#0072B2", vermillion="#D55E00", purple="#CC79A7", grey="#8C8C8C",
          gold="#B8860B")
WONG = [OI["black"], OI["orange"], OI["skyblue"], OI["green"],
        OI["blue"], OI["vermillion"], OI["purple"]]

CIVIDIS = mpl.colormaps["cividis"].copy()
CIVIDIS.set_bad("#e9e9e9")        # empty cells -> light grey


def use_paper_style(base=13):
    """Install the shared rcParams.  `base` scales the font sizes."""
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 150,
        "font.size": base, "axes.titlesize": base + 1.5, "axes.labelsize": base + 2.5,
        "legend.fontsize": base - 1.5, "xtick.labelsize": base, "ytick.labelsize": base,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 7.0, "ytick.major.size": 7.0,
        "xtick.minor.size": 4.0, "ytick.minor.size": 4.0,
        "xtick.major.width": 0.9, "ytick.major.width": 0.9,
        "xtick.minor.width": 0.7, "ytick.minor.width": 0.7,
        "axes.grid": False, "axes.linewidth": 0.9,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "font.family": "serif", "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
    })


def style_axes(ax, logx=False, logy=False, aspect=PANEL_ASPECT):
    ax.set_box_aspect(aspect)
    if logx: ax.set_xscale("log")
    if logy: ax.set_yscale("log")
    return ax


# ---- log-axis tick helpers: major every dex, labels every 2 dex -------------
def _label_even_pow(val, pos):
    if val <= 0: return ""
    e = np.log10(val)
    if abs(e - round(e)) > 1e-3: return ""
    e = int(round(e)); return rf"$10^{{{e}}}$" if e % 2 == 0 else ""

def log_axis_dex(axis):
    axis.set_major_locator(mticker.LogLocator(base=10, numticks=500))
    axis.set_major_formatter(mticker.FuncFormatter(_label_even_pow))
    axis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=500))
    axis.set_minor_formatter(mticker.NullFormatter())

def log_minor(axis):
    axis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=500))
    axis.set_minor_formatter(mticker.NullFormatter())


def whiten_ticks(ax):
    """For dark heatmaps: white tick marks, dark spines (ticks read on any cmap)."""
    ax.tick_params(which="both", color="white")
    for sp in ax.spines.values(): sp.set_edgecolor("0.2")

def whiten_cbar(cb):
    cb.ax.tick_params(which="both", color="white", direction="in")
    cb.outline.set_edgecolor("0.2")
    cb.ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())


def inside_title(ax, text, loc="upper left", fontsize=11, pad=0.03):
    """Panel title INSIDE the axes in a white bbox (matches debris Fig 1 / movies)."""
    xy = {"upper left": (pad, 1 - pad, "left", "top"),
          "upper right": (1 - pad, 1 - pad, "right", "top"),
          "lower left": (pad, pad, "left", "bottom"),
          "lower right": (1 - pad, pad, "right", "bottom")}[loc]
    ax.text(xy[0], xy[1], text, transform=ax.transAxes, ha=xy[2], va=xy[3],
            fontsize=fontsize, color="black",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9))


def edges(centers, log=True):
    """Cell edges from centres (log-spaced by default) for pcolormesh."""
    c = np.log10(np.asarray(centers, float)) if log else np.asarray(centers, float)
    e = np.empty(c.size + 1); e[1:-1] = 0.5 * (c[:-1] + c[1:])
    e[0] = c[0] - (c[1] - c[0]) / 2; e[-1] = c[-1] + (c[-1] - c[-2]) / 2
    return 10**e if log else e
