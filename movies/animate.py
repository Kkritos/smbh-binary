#!/usr/bin/env python3
"""Render a looping GIF (+ mp4) from a canonical trajectory CSV.

    python3 animate.py scenes/ejection.csv -o scenes/ejection.gif --title "..." \
            --center smbh --zoom bh --fps 24 --trail 48 [--inset] [--no-mp4]

Styling/zoom/inset/debris come mostly from a sidecar <name>.events.json; CLI flags
override. 'deb'/'debris' in --zoom tracks the debris cloud (micro-TDE).
"""
import argparse, scene_lib

def main():
    p = argparse.ArgumentParser(description="animate a canonical trajectory CSV -> GIF (+mp4)")
    p.add_argument("csv")
    p.add_argument("-o", "--out", default=None, help="output .gif (default: alongside the csv)")
    p.add_argument("--title", default="")
    p.add_argument("--center", default=None, help="body to center on (e.g. smbh, bh)")
    p.add_argument("--zoom", default=None, help="comma list of bodies the zoom tracks; 'debris' = debris cloud")
    p.add_argument("--fps", type=int, default=24)
    p.add_argument("--trail", type=int, default=40)
    p.add_argument("--figsize", type=float, default=6.2)
    p.add_argument("--hw-floor", type=float, default=1.0, help="minimum view half-width (AU)")
    p.add_argument("--inset", action="store_true", help="add a zoom inset")
    p.add_argument("--no-mp4", action="store_true", help="skip the mp4 export")
    a = p.parse_args()

    out = a.out or (a.csv.rsplit(".csv", 1)[0] + ".gif")
    scene = scene_lib.load(a.csv)
    kw = dict(title=a.title, fps=a.fps, trail=a.trail, figsize=a.figsize,
              hw_floor=a.hw_floor, inset=a.inset, mp4=not a.no_mp4)
    if a.center: kw["center"] = a.center
    if a.zoom:   kw["zoom_on"] = tuple(a.zoom.split(","))
    gif, mp4 = scene_lib.render(scene, out, **kw)
    print("wrote", gif + (f" and {mp4}" if mp4 else ""))

if __name__ == "__main__":
    main()
