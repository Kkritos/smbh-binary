#!/usr/bin/env python3
"""Tile a set of scene GIFs into ONE grid montage video (gif + mp4).

   python3 montage.py scenes/survey/*.gif -o scenes/survey/montage.gif \
           [--cols N --cell 240 --frames 160 --fps 20]

Reads GIF frames with PIL (no extra deps), loops shorter clips to a common length,
labels each cell by its filename stem, tiles into a grid, and writes a gif (+ mp4 when
ffmpeg is available).  Good as a single overview slide of the whole gallery.
"""
import argparse, glob, os, numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def read_gif(path, cell):
    im = Image.open(path)
    try:
        n = im.n_frames
    except Exception:
        n = 1
    frames = []
    for i in range(n):
        im.seek(i)
        frames.append(np.asarray(im.convert("RGB").resize((cell, cell), Image.BILINEAR)))
    return np.array(frames)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gifs", nargs="+")
    ap.add_argument("-o", "--out", default="montage.gif")
    ap.add_argument("--cols", type=int, default=0, help="0 = auto (~square)")
    ap.add_argument("--cell", type=int, default=240)
    ap.add_argument("--frames", type=int, default=160)
    ap.add_argument("--fps", type=int, default=20)
    a = ap.parse_args()

    paths = []
    for g in a.gifs:
        paths += sorted(glob.glob(g))
    paths = [p for p in paths if "montage" not in os.path.basename(p)]
    seen = set(); paths = [p for p in paths if not (p in seen or seen.add(p))]
    if not paths:
        raise SystemExit("no gifs matched")
    labels = [os.path.basename(p).rsplit(".gif", 1)[0] for p in paths]
    clips = [read_gif(p, a.cell) for p in paths]

    n = len(clips)
    cols = a.cols or int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    cell, NF = a.cell, a.frames
    H, W = rows * cell, cols * cell

    fig = plt.figure(figsize=(W/100, H/100), dpi=100, facecolor="black")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    imart = ax.imshow(np.zeros((H, W, 3), dtype=np.uint8))
    for idx, lab in enumerate(labels):            # static cell labels (on a readable strip)
        r, c = idx // cols, idx % cols
        ax.text(c*cell + cell*0.5, r*cell + 12, lab, color="white",
                fontsize=8, ha="center", va="top", family="monospace", zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc="black", ec="none", alpha=0.6))

    def frame(k):
        canv = np.zeros((H, W, 3), dtype=np.uint8)
        for idx, clip in enumerate(clips):
            r, c = idx // cols, idx % cols
            f = clip[int(k/NF * len(clip)) % len(clip)]
            canv[r*cell:(r+1)*cell, c*cell:(c+1)*cell] = f
        imart.set_array(canv)
        return [imart]

    anim = FuncAnimation(fig, frame, frames=NF, interval=1000/a.fps, blit=False)
    anim.save(a.out, writer=PillowWriter(fps=a.fps))
    try:
        import scene_lib
        fp = scene_lib._ffmpeg_path()
        if fp:
            from matplotlib.animation import FFMpegWriter
            matplotlib.rcParams["animation.ffmpeg_path"] = fp
            anim.save(a.out.rsplit(".gif", 1)[0] + ".mp4",
                      writer=FFMpegWriter(fps=a.fps, bitrate=3500))
    except Exception as ex:
        print(f"  (montage mp4 skipped: {ex})")
    plt.close(fig)
    print(f"wrote {a.out}  ({n} clips, {rows}x{cols} grid)")


if __name__ == "__main__":
    main()
