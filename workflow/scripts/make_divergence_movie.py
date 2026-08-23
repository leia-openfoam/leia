#!/usr/bin/env python3
"""Look at the parasitic current: where on the interface does it actually live?

Every diagnostic in this campaign so far has been a NORM -- max|U|, an L2 of the
delivered remainder, a per-step gain. A norm cannot distinguish "the whole
interface is slightly wrong" from "four cells are very wrong", and those have
different causes. This renders the field itself.

The question it answers: is the growing mode LOCKED TO THE MESH? On a Cartesian
mesh the angle between the interface normal and the nearest face normal sweeps
0-90 degrees four times around a circle, so a mesh-locked mode should put max|U|
at a small set of preferred angles and keep it there. A mode that is not
mesh-locked will let the maximum wander.

Outputs, into <case>/figures/:
  divergence_frames.png    montage of |U| with the psi = 0 contour and the
                           location of max|U| ringed, over the run
  divergence_maxloc.png    the polar angle and radius of the argmax of |U|
                           against time, with max|U| itself on a log axis
  divergence_maxloc.csv    the same as data (time, max|U|, x, y, theta, r)

Colour follows the job: |U| is magnitude, so one hue light-to-dark (Blues) and
never a rainbow. The two overlays are the only categorical marks and they carry
shape as well as hue -- the interface is a line, the argmax is a ring -- so
identity never rests on colour alone. Axis furniture is recessive; text is neutral
ink, never the overlay colour.

Usage:
  make_divergence_movie.py --case studies/divergenceMovie2D/N150 \
      [--radius 1e-3] [--center 0.003,0.003] [--frames 9]

Assumes a single uniform blockMesh hex block, one cell thick, written by a SERIAL
run (so fields are at <case>/<time>/, needing no reconstruction). The cell
ordering is x fastest then y, which the loader VERIFIES against the analytic
signed distance rather than assuming.
"""

import argparse
import csv
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

INK = "#1f2933"          # primary text
INK_MUTED = "#7b8794"    # axes, grid, secondary labels
INTERFACE = "#E8590C"    # overlay 1: the psi = 0 contour (a LINE)
ARGMAX = "#7048E8"       # overlay 2: the location of max|U| (a RING)


def read_field(path, ncells=None):
    """Parse an OpenFOAM ascii volField internalField into a numpy array.

    A field that is spatially constant is written as `internalField uniform ...`
    -- which happens for U at t = 0 on this case -- so a uniform entry is
    broadcast to the full cell count rather than returned as a bare 3-vector.
    """
    with open(path) as fh:
        txt = fh.read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", txt)
    if not m:
        raise RuntimeError("no internalField in %s" % path)
    if m.group(1) == "uniform":
        v = re.search(r"internalField\s+uniform\s+(\(?[^;]+?\)?)\s*;", txt).group(1)
        nums = [float(x) for x in re.findall(r"[-+0-9.eE]+", v)]
        if ncells is None:
            return np.array(nums), True
        comp = np.array(nums)
        if comp.size == 1:
            return np.full(ncells, comp[0]), True
        return np.tile(comp, (ncells, 1)), True
    body = txt[m.end():]
    n = int(re.search(r"(\d+)", body).group(1))
    start = body.index("(", body.index(str(n)))
    # vector entries look like "(a b c)"; scalars are bare numbers
    if body[start + 1:start + 40].lstrip().startswith("("):
        vals = re.findall(r"\(([^)]*)\)", body[start + 1:])
        arr = np.array([[float(t) for t in v.split()] for v in vals[:n]])
    else:
        chunk = body[start + 1:]
        end = chunk.index(")")
        arr = np.array([float(t) for t in chunk[:end].split()][:n])
    if len(arr) != n:
        raise RuntimeError("expected %d entries, parsed %d in %s" % (n, len(arr), path))
    return arr, False


def times(case):
    out = []
    for d in os.listdir(case):
        p = os.path.join(case, d)
        if not os.path.isdir(p):
            continue
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isfile(os.path.join(p, "U")) and os.path.isfile(os.path.join(p, "psi")):
            out.append((t, p))
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--center", default=None, help="x,y of the droplet centre")
    ap.add_argument("--frames", type=int, default=9)
    a = ap.parse_args()

    ts = times(a.case)
    if len(ts) < 3:
        raise SystemExit("need at least 3 written times in %s (found %d)"
                         % (a.case, len(ts)))

    psi0, _ = read_field(os.path.join(ts[0][1], "psi"))
    ncells = len(psi0)
    n = int(round(math.sqrt(ncells)))
    if n * n != ncells:
        raise SystemExit("%d cells is not a square 2D block" % ncells)

    # Domain length and centre: recover L from the case tokens if available,
    # else assume the psi field's own extent. Verify the assumed cell ordering
    # against the analytic signed distance before trusting any of it.
    L = None
    pj = os.path.join(a.case, "case_params.json")
    if os.path.isfile(pj):
        import json
        tok = json.load(open(pj))["tokens"]
        L = float(tok.get("DOMAIN_LENGTH") or 0.01)
    if L is None:
        L = 0.01
    h = L / n
    cx, cy = (L / 2, L / 2)
    if a.center:
        cx, cy = [float(v) for v in a.center.split(",")]

    xs = (np.arange(n) + 0.5) * h
    X, Y = np.meshgrid(xs, xs)                     # x fastest, then y
    psi_exact = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2) - a.radius
    psi_grid = psi0.reshape(n, n)
    err = np.abs(psi_grid - psi_exact).max() / a.radius
    print("cell-ordering check: max|psi_file - psi_exact|/R = %.3e" % err)
    if err > 0.25:
        raise SystemExit("cell ordering or geometry assumption is wrong "
                         "(relative mismatch %.3e); refusing to plot" % err)

    # ---- track the argmax of |U| over every written time -------------------
    rows = []
    for t, p in ts:
        U, _ = read_field(os.path.join(p, "U"), ncells)
        mag = np.sqrt((U ** 2).sum(axis=1)).reshape(n, n)
        j, i = np.unravel_index(np.argmax(mag), mag.shape)
        x, y = xs[i], xs[j]
        th = math.degrees(math.atan2(y - cy, x - cx)) % 360.0
        rows.append(dict(time=t, maxMagU=float(mag.max()), x=float(x), y=float(y),
                         theta=th, r=float(math.hypot(x - cx, y - cy) / a.radius)))

    figdir = os.path.join(a.case, "figures")
    os.makedirs(figdir, exist_ok=True)
    with open(os.path.join(figdir, "divergence_maxloc.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- montage of |U| with the interface and the argmax -------------------
    idx = np.unique(np.linspace(0, len(ts) - 1, a.frames).astype(int))
    ncol = min(3, len(idx))
    nrow = int(math.ceil(len(idx) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 3.7 * nrow),
                             squeeze=False)
    for k, ax in enumerate(axes.ravel()):
        if k >= len(idx):
            ax.axis("off")
            continue
        t, p = ts[idx[k]]
        U, _ = read_field(os.path.join(p, "U"), ncells)
        ps, _ = read_field(os.path.join(p, "psi"), ncells)
        mag = np.sqrt((U ** 2).sum(axis=1)).reshape(n, n)
        vmax = max(mag.max(), 1e-30)
        im = ax.imshow(mag, origin="lower", cmap="Blues",
                       extent=[0, L * 1e3, 0, L * 1e3],
                       norm=LogNorm(vmin=vmax / 1e4, vmax=vmax))
        ax.contour(X * 1e3, Y * 1e3, ps.reshape(n, n), levels=[0.0],
                   colors=INTERFACE, linewidths=1.4)
        j, i = np.unravel_index(np.argmax(mag), mag.shape)
        ax.plot(xs[i] * 1e3, xs[j] * 1e3, "o", mfc="none", mec=ARGMAX, mew=1.8,
                ms=11)
        ax.set_title("t = %.4g ms   max|U| = %.3g m/s" % (t * 1e3, vmax),
                     fontsize=9, color=INK)
        # zoom to the droplet plus a margin: the far field is featureless
        pad = 2.2 * a.radius * 1e3
        ax.set_xlim(cx * 1e3 - pad, cx * 1e3 + pad)
        ax.set_ylim(cy * 1e3 - pad, cy * 1e3 + pad)
        ax.tick_params(colors=INK_MUTED, labelsize=7, length=2)
        for s in ax.spines.values():
            s.set_color(INK_MUTED)
            s.set_linewidth(0.6)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(colors=INK_MUTED, labelsize=6, length=2)
        cb.outline.set_edgecolor(INK_MUTED)
        cb.outline.set_linewidth(0.6)
    fig.suptitle("Velocity magnitude, interface (line) and location of max|U| (ring)",
                 fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    f1 = os.path.join(figdir, "divergence_frames.png")
    fig.savefig(f1, dpi=140)
    plt.close(fig)

    # ---- where the maximum sits, against time ------------------------------
    tms = [r["time"] * 1e3 for r in rows]
    fig, axs = plt.subplots(3, 1, figsize=(8.2, 8.0), sharex=True)
    axs[0].semilogy(tms, [r["maxMagU"] for r in rows], color=INK, lw=1.6)
    axs[0].set_ylabel("max|U|  [m/s]", color=INK, fontsize=9)
    axs[1].scatter(tms, [r["theta"] for r in rows], s=9, color=ARGMAX)
    for g in (0, 45, 90, 135, 180, 225, 270, 315, 360):
        axs[1].axhline(g, color=INK_MUTED, lw=0.5, ls=":", alpha=0.7)
    axs[1].set_ylim(0, 360)
    axs[1].set_yticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
    axs[1].set_ylabel("angle of max|U|  [deg]", color=INK, fontsize=9)
    axs[2].scatter(tms, [r["r"] for r in rows], s=9, color=INTERFACE)
    axs[2].axhline(1.0, color=INK_MUTED, lw=0.8)
    axs[2].set_ylabel("radius of max|U|  [R]", color=INK, fontsize=9)
    axs[2].set_xlabel("time  [ms]", color=INK, fontsize=9)
    for ax in axs:
        ax.grid(True, color=INK_MUTED, alpha=0.18, lw=0.5)
        ax.tick_params(colors=INK_MUTED, labelsize=8)
        for s in ax.spines.values():
            s.set_color(INK_MUTED)
            s.set_linewidth(0.6)
    axs[1].set_title("dotted lines mark the mesh-locked angles (multiples of 45 deg)",
                     fontsize=8, color=INK_MUTED)
    fig.suptitle("Does the growing mode sit at preferred angles?", fontsize=11,
                 color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    f2 = os.path.join(figdir, "divergence_maxloc.png")
    fig.savefig(f2, dpi=140)
    plt.close(fig)

    # ---- a plain-numbers summary, since a picture is not a measurement -----
    th = np.array([r["theta"] for r in rows])
    late = th[len(th) // 2:]
    d45 = np.minimum(late % 45.0, 45.0 - (late % 45.0))
    print("\nargmax angle over the SECOND HALF of the run (%d frames):" % len(late))
    print("  distance to the nearest multiple of 45 deg: mean %.2f, median %.2f"
          % (d45.mean(), np.median(d45)))
    print("  (uniform-random would average 11.25 deg; mesh-locked tends to 0)")
    print("  distinct angles within 5 deg tolerance: %d"
          % len(np.unique(np.round(late / 5.0))))
    print("  radius of argmax, second half: mean %.3f R, min %.3f, max %.3f"
          % (np.mean([r["r"] for r in rows[len(rows) // 2:]]),
             min(r["r"] for r in rows[len(rows) // 2:]),
             max(r["r"] for r in rows[len(rows) // 2:])))
    print("\nwrote %s\n      %s" % (f1, f2))


if __name__ == "__main__":
    main()
