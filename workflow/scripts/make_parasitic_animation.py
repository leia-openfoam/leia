#!/usr/bin/env python3
"""Animate the parasitic current from a dense-output serial run.

Two panels per frame, both with the psi = 0 contour in orange:
  left  -- direction glyphs over the volume fraction (white air, grey water),
           arrow colour = |U| on ONE log scale shared across ALL frames (so
           growth is visible as arrows darkening over time), arrow length = |U|
           normalised per frame (so the flow geometry stays readable while the
           amplitude spans decades); zoomed to the droplet.
  right -- log |U| heatmap of the FULL domain on the same shared scale, so the
           global recirculation and its stagnation points are visible.
The title carries t and the current max|U|; a progress bar along the bottom
shows position in the run.

Frames are rendered in parallel to PNGs and encoded with ffmpeg (h264, yuv420p).

Usage:
  make_parasitic_animation.py --case studies/divergenceMovie2D/N150 \
      [--radius 1e-3] [--fps 12] [--out <file.mp4>]
"""

import argparse
import math
import os
import re
import subprocess
import tempfile
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap

INK = "#1f2933"
INK_MUTED = "#7b8794"
INTERFACE = "#E8590C"
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


CFG = {}


def render_frame(job):
    k, t, tdir = job
    n = CFG["n"]; L = CFG["L"]; R = CFG["radius"]
    norm = CFG["norm"]; tmax = CFG["tmax"]; nfr = CFG["nframes"]
    h = L / n
    xs = (np.arange(n) + 0.5) * h
    X, Y = np.meshgrid(xs, xs)
    U = read_field(os.path.join(tdir, "U"), n * n).reshape(n, n, 3)
    # Translating cases: the interesting field is the DISTURBANCE U - U_ref, not U.
    # With U_ref = (0.05, 0, 0) the free stream is 0.05 m/s while the parasitic current
    # is ~1e-2, so plotting |U| shows a uniform blue rectangle and nothing else.
    U = (U - np.asarray(CFG["uref"])) / CFG["vscale"]
    ps = read_field(os.path.join(tdir, "psi"), n * n).reshape(n, n)
    al = read_field(os.path.join(tdir, "alpha.water"), n * n).reshape(n, n)
    mag = np.sqrt(U[..., 0] ** 2 + U[..., 1] ** 2)
    pmax = max(mag.max(), 1e-30)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 5.6))

    # ---- left: glyphs over alpha, zoomed --------------------------------
    axL.imshow(al, origin="lower", cmap="Greys", vmin=0, vmax=4.5,
               extent=[0, L * 1e3, 0, L * 1e3])
    axL.contour(X * 1e3, Y * 1e3, ps, levels=[0.0], colors=INTERFACE,
                linewidths=1.2)
    sub = max(1, n // 30)
    s = slice(sub // 2, None, sub)
    Xs, Ys, Us, Vs, Ms = X[s, s] * 1e3, Y[s, s] * 1e3, U[s, s, 0], U[s, s, 1], mag[s, s]
    keep = Ms > 1e-4 * pmax
    with np.errstate(invalid="ignore", divide="ignore"):
        ux, uy = np.where(keep, Us / Ms, 0), np.where(keep, Vs / Ms, 0)
    rel = np.where(keep, Ms / pmax, 0)
    glyph = 0.92 * sub * h * 1e3
    axL.quiver(Xs[keep], Ys[keep], (ux * rel * glyph)[keep],
               (uy * rel * glyph)[keep], Ms[keep], cmap=ARROWS, norm=norm,
               angles="xy", scale_units="xy", scale=1.0, width=0.0038,
               headwidth=3.2, headlength=3.6, headaxislength=3.2, minlength=0.5)
    pad = 1.9 * R * 1e3
    # Follow the droplet: a translating case leaves a box centred on L/2 within a few ms.
    # The centroid is alpha-weighted, and falls back to the domain centre if the phase
    # has left the domain (alpha sums to ~0), so a frame can never blow up the render.
    aw = al.sum()
    if aw > 1e-9:
        cx, cy = (al * X).sum() / aw * 1e3, (al * Y).sum() / aw * 1e3
    else:
        cx = cy = L * 1e3 / 2
    axL.set_xlim(cx - pad, cx + pad)
    axL.set_ylim(cy - pad, cy + pad)
    axL.set_title("direction glyphs over $\\alpha$  (length: per-frame, "
                  "colour: shared)", fontsize=9, color=INK)

    # ---- right: |U| heatmap, full domain, shared scale -------------------
    im = axR.imshow(mag, origin="lower", cmap="Blues", norm=norm,
                    extent=[0, L * 1e3, 0, L * 1e3])
    axR.contour(X * 1e3, Y * 1e3, ps, levels=[0.0], colors=INTERFACE,
                linewidths=1.0)
    axR.set_title("%s full domain, shared log scale" % CFG["qlabel"],
                  fontsize=10, color=INK)
    # Explicit colorbar axes with RESERVED canvas room. The first cut used
    # fig.colorbar(..., ax=axR) with subplots_adjust(right=0.99), which pushed
    # the tick labels off the canvas -- invisible in the rendered frames, worst
    # after the GIF downscale. Fonts are sized to stay legible at 880 px.
    cax = fig.add_axes([0.915, 0.13, 0.018, 0.72])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("%s%s" % (CFG["qlabel"], CFG["qunit"]), fontsize=11, color=INK)
    cb.ax.tick_params(colors=INK, labelsize=10, length=3)
    cb.outline.set_edgecolor(INK_MUTED)

    for ax in (axL, axR):
        ax.set_xticks([]), ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(INK_MUTED)
            sp.set_linewidth(0.6)

    fig.suptitle("%s      t = %6.2f ms      max %s = %s"
                 % (CFG["label"], t * 1e3, CFG["qlabel"],
                    ("%.3f" % pmax) if CFG["vscale"] != 1.0 else ("%.3e m/s" % pmax)),
                 fontsize=12, color=INK, family="monospace")
    # progress bar
    axp = fig.add_axes([0.075, 0.035, 0.85, 0.012])
    axp.axhspan(0, 1, 0, t / tmax, color=INTERFACE, alpha=0.85)
    axp.set_xlim(0, 1), axp.set_xticks([]), axp.set_yticks([])
    for sp in axp.spines.values():
        sp.set_color(INK_MUTED)
        sp.set_linewidth(0.5)

    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.03, right=0.90,
                        wspace=0.06)
    out = os.path.join(CFG["framedir"], "f%04d.png" % k)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--uref", type=float, nargs=3, default=[0.0, 0.0, 0.0],
                    help="reference velocity subtracted before plotting; for a "
                         "translating droplet pass the translation velocity so the "
                         "plotted field is the disturbance U - U_ref")
    ap.add_argument("--label", default="Stationary droplet, R/h = 25, filter off",
                    help="suptitle text")
    ap.add_argument("--normalize", action="store_true",
                    help="divide the plotted magnitude by |uref|, so the colour bar is "
                         "the DIMENSIONLESS |U-U0|/U0 rather than m/s -- 0.2 reads as "
                         "'20 percent of the translation speed', which 1e-2 m/s does not")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import json
    tok = json.load(open(os.path.join(a.case, "case_params.json")))["tokens"]
    n = int(float(tok["N_CELLS"]))
    L = float(tok.get("DOMAIN_LENGTH") or 0.01)

    frames = []
    for d in os.listdir(a.case):
        p = os.path.join(a.case, d)
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, "U")) \
           and os.path.isfile(os.path.join(p, "alpha.water")):
            frames.append((t, p))
    frames.sort()
    if len(frames) < 10:
        raise SystemExit("only %d frames in %s" % (len(frames), a.case))

    # global colour scale over all frames
    uref = np.asarray(a.uref)
    vscale = float(np.linalg.norm(uref)) if (a.normalize and np.linalg.norm(uref) > 0) else 1.0
    if not uref.any():
        qlabel, qunit = "|U|", " [m/s]"
    elif vscale != 1.0:
        qlabel, qunit = "|U-U0|/U0", ""
    else:
        qlabel, qunit = "|U-U0|", " [m/s]"
    gmax = 0.0
    for t, p in frames[::5] + [frames[-1]]:
        U = (read_field(os.path.join(p, "U"), n * n) - uref) / vscale
        gmax = max(gmax, float(np.sqrt((U[:, 0] ** 2 + U[:, 1] ** 2)).max()))
    print("frames: %d   global max %s = %.3e" % (len(frames), qlabel, gmax))

    framedir = tempfile.mkdtemp(prefix="anim_")
    CFG.update(n=n, L=L, radius=a.radius, framedir=framedir,
               norm=LogNorm(vmin=gmax / 1e5, vmax=gmax),
               tmax=frames[-1][0], nframes=len(frames),
               uref=uref, qlabel=qlabel, qunit=qunit, label=a.label, vscale=vscale)

    jobs = [(k, t, p) for k, (t, p) in enumerate(frames)]
    with Pool(a.jobs) as pool:
        for i, _ in enumerate(pool.imap_unordered(render_frame, jobs)):
            if (i + 1) % 25 == 0:
                print("  rendered %d/%d" % (i + 1, len(jobs)))

    out = a.out or os.path.join(a.case, "figures", "parasitic.mp4")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps),
         "-i", os.path.join(framedir, "f%04d.png"),
         "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
         "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "22", out],
        check=True)
    print("wrote", out, "(%.1f MB)" % (os.path.getsize(out) / 1e6))

    # GIF is THE deliverable format (renders inline and downloads without a
    # player); the mp4 above stays for decks. Palette-optimised, 880 px wide.
    gif = os.path.splitext(out)[0] + ".gif"
    pal = os.path.join(framedir, "pal.png")
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", "-i", out,
         "-vf", "fps=10,scale=880:-1:flags=lanczos,palettegen", pal],
        check=True)
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", "-i", out, "-i", pal,
         "-lavfi", "fps=10,scale=880:-1:flags=lanczos [x]; "
         "[x][1:v] paletteuse=dither=bayer:bayer_scale=4", gif],
        check=True)
    print("wrote", gif, "(%.1f MB)" % (os.path.getsize(gif) / 1e6))


if __name__ == "__main__":
    main()
