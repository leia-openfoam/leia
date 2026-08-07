#!/usr/bin/env python3
"""Volume-fraction + velocity-glyph animation of a reconstructed 2D droplet case.

Renders one frame per written time directory (alpha.water as the background,
the alpha = 0.5 contour, and velocity glyphs sub-sampled and re-normalised per
frame so the pattern stays readable while max|U| grows by decades), then
assembles an mp4 with ffmpeg.

Uniform-hex 2D cases only (the frame is built from the (N, N) cell array).

Usage (from the case directory or with an explicit path):
    python3 workflow/scripts/make_droplet_animation.py <case> <N> <out.mp4> \
        [--sub 5] [--fps 5]
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_scalar(path, n):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(", txt)
    if not m:
        u = re.search(r"internalField\s+uniform\s+([0-9eE+.-]+)", txt)
        return np.full(n, float(u.group(1)))
    cnt = int(m.group(1))
    start = txt.index("(", m.end() - 1) + 1
    return np.array([float(v) for v in txt[start:].split(")")[0].split()[:cnt]])


def read_vector(path, n):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n(\d+)\s*\n\(", txt)
    if not m:
        u = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", txt)
        return np.tile([float(v) for v in u.group(1).split()], (n, 1))
    cnt = int(m.group(1))
    tuples = re.findall(r"\(([^)]+)\)", txt[m.end():])[:cnt]
    return np.array([[float(x) for x in t.split()] for t in tuples])


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("N", type=int)
    ap.add_argument("out")
    ap.add_argument("--sub", type=int, default=5, help="glyph sub-sampling")
    ap.add_argument("--fps", type=int, default=5)
    ap.add_argument("--alpha", default="alpha.water")
    ap.add_argument("--L", type=float, default=0.01, help="domain edge [m]")
    a = ap.parse_args(argv)

    N, nc = a.N, a.N*a.N
    times = sorted(
        (d for d in os.listdir(a.case)
         if re.fullmatch(r"[0-9]+(\.[0-9eE+-]+)?", d)
         and os.path.isdir(os.path.join(a.case, d))),
        key=float)
    if not times:
        print(f"[anim] no time directories in {a.case}"); return 1

    x = (np.arange(N) + 0.5)*a.L/N
    X, Y = np.meshgrid(x, x)
    s = np.s_[a.sub//2::a.sub, a.sub//2::a.sub]

    with tempfile.TemporaryDirectory() as tmp:
        for fi, t in enumerate(times):
            td = os.path.join(a.case, t)
            al = read_scalar(os.path.join(td, a.alpha), nc).reshape(N, N)
            U = read_vector(os.path.join(td, "U"), nc)
            ux, uy = U[:, 0].reshape(N, N), U[:, 1].reshape(N, N)
            umag = np.hypot(ux, uy)
            scale = max(umag.max(), 1e-30)

            fig, ax = plt.subplots(figsize=(6.4, 6.4))
            ax.pcolormesh(X, Y, al, cmap="Blues", vmin=0, vmax=1,
                          shading="auto", rasterized=True)
            ax.contour(X, Y, al, levels=[0.5], colors="#D55E00", linewidths=1.2)
            ax.quiver(X[s], Y[s], ux[s]/scale, uy[s]/scale, umag[s]/scale,
                      cmap="viridis", scale=25, width=2.2e-3, clim=(0, 1))
            ax.set_aspect("equal")
            ax.set_xlim(0, a.L); ax.set_ylim(0, a.L)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"t = {float(t):.3f} s    max|U| = {umag.max():.2e} m/s"
                         "   (glyphs normalized per frame)", fontsize=10)
            fig.tight_layout()
            fig.savefig(os.path.join(tmp, f"frame_{fi:04d}.png"), dpi=110)
            plt.close(fig)
        print(f"[anim] rendered {len(times)} frames")

        subprocess.run(
            ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps),
             "-i", os.path.join(tmp, "frame_%04d.png"),
             "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "22", a.out],
            check=True)
    print(f"[anim] wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
