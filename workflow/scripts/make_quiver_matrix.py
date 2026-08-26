#!/usr/bin/env python3
"""Velocity GLYPHS over the volume fraction: what the |U| heatmap cannot show.

A magnitude heatmap renders a stagnation point and a hole in the data the same
way -- as a light spot. This figure settles which is which: arrows show the
DIRECTION of the parasitic current everywhere, the volume fraction alpha is the
background (so any phase where none should exist is immediately visible), and
the psi = 0 contour stays as the interface mark.

Encoding, stated once and printed on the figure:
  * arrow COLOUR  = |U| on ONE log scale shared by every panel (single hue,
    clipped to the dark half of Blues so every arrow reads on the light
    background) -- magnitude is comparable across panels;
  * arrow LENGTH  = |U| normalised by the PANEL's own maximum -- geometry is
    visible in every panel regardless of how many decades the run spans;
  * background    = alpha.water in light greys (white = air, grey = water);
  * orange line   = psi = 0.
A stagnation point of a recirculation now looks like arrows CIRCLING a spot;
a data defect (phase island, source) would look like alpha off {0,1} or arrows
DIVERGING from a point.

Usage:
  make_quiver_matrix.py --study studies/parasiticMatrix2D [--radius 1e-3]
      [--times 2.5,7.5,12.5,18.75,25] [--detail N150] [--out <dir>]

--detail adds a second figure: FULL-DOMAIN panels for one case (denser arrows),
for structures that live outside the droplet zoom. Same serial single-block
assumptions as make_parasitic_matrix.py; cell ordering verified against the
analytic circle before anything is drawn.
"""

import argparse
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap

INK = "#1f2933"
INK_MUTED = "#7b8794"
INTERFACE = "#E8590C"

# single hue, dark half only: readable on the near-white alpha background
ARROWS = LinearSegmentedColormap.from_list(
    "bluesDark", plt.get_cmap("Blues")(np.linspace(0.38, 1.0, 256)))


def read_field(path, ncells=None):
    txt = open(path).read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", txt)
    if not m:
        raise RuntimeError("no internalField in %s" % path)
    if m.group(1) == "uniform":
        comp = np.array([float(x) for x in re.findall(
            r"[-+0-9.eE]+", re.search(
                r"internalField\s+uniform\s+(\(?[^;]+?\)?)\s*;", txt).group(1))])
        if ncells is None:
            return comp
        return np.full(ncells, comp[0]) if comp.size == 1 else np.tile(comp, (ncells, 1))
    body = txt[m.end():]
    n = int(re.search(r"(\d+)", body).group(1))
    start = body.index("(", body.index(str(n)))
    if body[start + 1:start + 40].lstrip().startswith("("):
        vals = re.findall(r"\(([^)]*)\)", body[start + 1:])
        return np.array([[float(t) for t in v.split()] for v in vals[:n]])
    chunk = body[start + 1:]
    return np.array([float(t) for t in chunk[:chunk.index(")")].split()][:n])


def load_case(case, radius):
    import json
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    N = int(float(tok["N_CELLS"]))
    L = float(tok.get("DOMAIN_LENGTH") or 0.01)
    times = {}
    for d in os.listdir(case):
        p = os.path.join(case, d)
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, "U")):
            times[t] = p
    return radius / (L / N), L, N, times


def draw_panel(ax, tdir, n, L, radius, norm, sub, zoom):
    """One panel: alpha background, psi contour, direction glyphs."""
    h = L / n
    xs = (np.arange(n) + 0.5) * h
    X, Y = np.meshgrid(xs, xs)
    U = read_field(os.path.join(tdir, "U"), n * n).reshape(n, n, 3)
    ps = read_field(os.path.join(tdir, "psi"), n * n).reshape(n, n)
    al = read_field(os.path.join(tdir, "alpha.water"), n * n).reshape(n, n)
    mag = np.sqrt(U[..., 0] ** 2 + U[..., 1] ** 2)
    pmax = max(mag.max(), 1e-30)

    ax.imshow(al, origin="lower", cmap="Greys", vmin=0, vmax=4.5,
              extent=[0, L * 1e3, 0, L * 1e3])          # alpha=1 -> light grey
    ax.contour(X * 1e3, Y * 1e3, ps, levels=[0.0],
               colors=INTERFACE, linewidths=1.1)

    s = slice(sub // 2, None, sub)
    Xs, Ys = X[s, s] * 1e3, Y[s, s] * 1e3
    Us, Vs, Ms = U[s, s, 0], U[s, s, 1], mag[s, s]
    keep = Ms > 1e-4 * pmax
    with np.errstate(invalid="ignore", divide="ignore"):
        ux, uy = np.where(keep, Us / Ms, 0), np.where(keep, Vs / Ms, 0)
    rel = np.where(keep, Ms / pmax, 0)                  # per-panel length
    glyph = 0.92 * sub * h * 1e3                        # longest arrow [mm]
    ax.quiver(Xs[keep], Ys[keep], (ux * rel * glyph)[keep],
              (uy * rel * glyph)[keep], Ms[keep],
              cmap=ARROWS, norm=norm, angles="xy", scale_units="xy", scale=1.0,
              width=0.004, headwidth=3.2, headlength=3.6, headaxislength=3.2,
              minlength=0.5)
    if zoom:
        pad = 1.9 * radius * 1e3
        ax.set_xlim(L * 1e3 / 2 - pad, L * 1e3 / 2 + pad)
        ax.set_ylim(L * 1e3 / 2 - pad, L * 1e3 / 2 + pad)
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(INK_MUTED)
        sp.set_linewidth(0.6)
    ax.text(0.03, 0.03, "max %.2e" % pmax, transform=ax.transAxes, fontsize=6.5,
            color=INK, bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.4))
    return pmax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True)
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--times", default="2.5,7.5,12.5,18.75,25")
    ap.add_argument("--detail", default=None, help="case dir name for full-domain detail")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    want = [float(v) * 1e-3 for v in a.times.split(",")]
    cases = []
    for d in sorted(os.listdir(a.study)):
        p = os.path.join(a.study, d)
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, "case_params.json")):
            cases.append((d,) + load_case(p, a.radius))
    cases = [c for c in cases if len(c[4]) >= 2]
    cases.sort(key=lambda c: c[1])

    # shared colour scale from the global max over all panels drawn
    gmax = 0.0
    sel = {}
    for name, RoH, L, n, times in cases:
        avail = sorted(times)
        for w in want:
            t = min(avail, key=lambda s: abs(s - w))
            if abs(t - w) > 0.05 * max(want):
                continue
            U = read_field(os.path.join(times[t], "U"), n * n)
            gmax = max(gmax, float(np.sqrt((U ** 2).sum(axis=1)).max()))
            sel[(name, w)] = (t, times[t])
    norm = LogNorm(vmin=gmax / 1e5, vmax=gmax)

    nrow, ncol = len(cases), len(want)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 1.5, 2.75 * nrow),
                             squeeze=False)
    for r, (name, RoH, L, n, times) in enumerate(cases):
        for c, w in enumerate(want):
            ax = axes[r][c]
            if (name, w) not in sel:
                ax.axis("off")
                continue
            t, tdir = sel[(name, w)]
            draw_panel(ax, tdir, n, L, a.radius, norm, max(1, n // 26), True)
            if r == 0:
                ax.set_title("t = %.2f ms" % (t * 1e3), fontsize=10, color=INK)
            if c == 0:
                ax.set_ylabel("R/h = %.1f" % RoH, fontsize=10, color=INK)
    fig.suptitle("Direction glyphs of the parasitic current over the volume fraction\n"
                 "arrow colour = |U|, one log scale for all panels; arrow length = |U| "
                 "normalised per panel; grey = water, white = air; orange = $\\psi$ = 0",
                 fontsize=11, color=INK)
    sm = plt.cm.ScalarMappable(cmap=ARROWS, norm=norm)
    fig.tight_layout(rect=[0, 0.02, 0.92, 0.94])
    cax = fig.add_axes([0.935, 0.10, 0.016, 0.76])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("|U| [m/s]", fontsize=8, color=INK_MUTED)
    cb.ax.tick_params(colors=INK_MUTED, labelsize=7.5, length=2)
    cb.outline.set_edgecolor(INK_MUTED)

    out = a.out or os.path.join(a.study, "figures")
    os.makedirs(out, exist_ok=True)
    f1 = os.path.join(out, "quiver_matrix.png")
    fig.savefig(f1, dpi=150)
    plt.close(fig)
    print("wrote", f1)

    if a.detail:
        name = a.detail
        match = [c for c in cases if c[0] == name]
        if match:
            _, RoH, L, n, times = match[0]
            dt_want = [w for w in want if (name, w) in sel][-3:]
            fig, axes = plt.subplots(1, len(dt_want),
                                     figsize=(5.4 * len(dt_want) + 1.2, 5.6),
                                     squeeze=False)
            for c, w in enumerate(dt_want):
                t, tdir = sel[(name, w)]
                draw_panel(axes[0][c], tdir, n, L, a.radius, norm,
                           max(1, n // 40), False)
                axes[0][c].set_title("t = %.2f ms" % (t * 1e3), fontsize=11,
                                     color=INK)
            fig.suptitle("FULL DOMAIN, R/h = %.1f -- same encoding" % RoH,
                         fontsize=12, color=INK)
            sm = plt.cm.ScalarMappable(cmap=ARROWS, norm=norm)
            fig.tight_layout(rect=[0, 0.02, 0.93, 0.93])
            cax = fig.add_axes([0.945, 0.12, 0.014, 0.72])
            cb = fig.colorbar(sm, cax=cax)
            cb.set_label("|U| [m/s]", fontsize=8, color=INK_MUTED)
            cb.ax.tick_params(colors=INK_MUTED, labelsize=7.5, length=2)
            cb.outline.set_edgecolor(INK_MUTED)
            f2 = os.path.join(out, "quiver_detail_%s.png" % name)
            fig.savefig(f2, dpi=150)
            plt.close(fig)
            print("wrote", f2)


if __name__ == "__main__":
    main()
