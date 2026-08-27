#!/usr/bin/env python3
"""Parasitic-current panels for an UNSTRUCTURED (polyhedral, 3D) droplet run.

The 2D animator (make_parasitic_animation.py) reshapes the internal field to
n x n, so it only ever worked on a Cartesian hex mesh. This script is the
mesh-agnostic twin: it reads CELL CENTRES from the case (written once by
`postProcess -func writeCellCentres`, mesh is static) and treats every field as
a POINT CLOUD, so hex, polyhedral and any decomposition are handled the same
way. Parallel cases are read straight out of processor*/ -- no reconstructPar.

Per frame, two panels, matching the 2D figures:
  left  -- direction glyphs over the volume fraction in a z-slab through the
           droplet centre; arrow colour = |U| on ONE log scale shared across
           all frames (growth shows as darkening), arrow length normalised per
           frame (geometry stays readable while the amplitude spans decades).
  right -- log |U| over the full slab on the same shared scale, so the outer
           recirculation and its stagnation points stay visible.
Plus a max|U|(t) trace with a cursor, read from the solver's per-step CSV.

Usage:
  make_poly_droplet_panels.py --case studies/<study>/<arm> [--radius 1e-3]
      [--slab 1.0] [--times all|t1,t2,...] [--gif] [--out-dir DIR]
"""

import argparse
import glob
import os
import re
import subprocess

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


def read_field(path):
    """internalField of an OpenFOAM ascii field -> (n,) or (n,3) array."""
    txt = open(path).read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", txt)
    if not m:
        raise RuntimeError("no internalField in %s" % path)
    if m.group(1) == "uniform":
        comp = [float(x) for x in re.findall(
            r"[-+0-9.eE]+",
            re.search(r"internalField\s+uniform\s+(\(?[^;]+?\)?)\s*;",
                      txt).group(1))]
        return np.array(comp)          # caller broadcasts
    body = txt[m.end():]
    n = int(re.search(r"(\d+)", body).group(1))
    start = body.index("(", body.index(str(n)))
    if body[start + 1:start + 40].lstrip().startswith("("):
        vals = re.findall(r"\(([^)]*)\)", body[start + 1:])
        return np.array([[float(t) for t in v.split()] for v in vals[:n]])
    chunk = body[start + 1:]
    return np.array([float(t) for t in chunk[:chunk.index(")")].split()][:n])


def rank_dirs(case):
    """processor*/ in rank order, or [case] for a serial case."""
    procs = sorted(glob.glob(os.path.join(case, "processor*")),
                   key=lambda p: int(re.search(r"(\d+)$", p).group(1)))
    return procs if procs else [case]


def _open_mesh_file(meshdir, name):
    """polyMesh file, transparently gz."""
    p = os.path.join(meshdir, name)
    if os.path.isfile(p):
        return open(p, "rb").read().decode("utf-8", "replace")
    if os.path.isfile(p + ".gz"):
        import gzip
        return gzip.open(p + ".gz", "rb").read().decode("utf-8", "replace")
    raise SystemExit(f"missing {p}")


def _read_labels(txt):
    """Leading count + flat list of labels from a polyMesh owner/neighbour."""
    body = txt[txt.index("// * * *"):] if "// * * *" in txt else txt
    m = re.search(r"^\s*(\d+)\s*\n\s*\(", body, re.M)
    start = body.index("(", m.start())
    end = body.index(")", start)
    return np.fromstring(body[start + 1:end], dtype=np.int64, sep=" ")


def cell_centres_from_mesh(root):
    """Approximate cell centres straight from constant/polyMesh.

    Face centre = mean of its points, cell centre = mean of its face centres.
    Exact OpenFOAM centroids use an area/volume-weighted tet decomposition;
    for slab selection and plotting this approximation is ample, and it makes
    the script independent of OpenFOAM (no mpirun, no reconstructPar), so a
    case pulled from the cluster -- which has no processor*/system -- reads
    exactly like a local one.
    """
    meshdir = os.path.join(root, "constant", "polyMesh")
    ptxt = _open_mesh_file(meshdir, "points")
    pm = re.search(r"^\s*(\d+)\s*\n\s*\(", ptxt, re.M)
    pstart = ptxt.index("(", pm.start())
    pts = np.array([[float(v) for v in s.split()]
                    for s in re.findall(r"\(([^()]*)\)", ptxt[pstart + 1:])],
                   dtype=float)
    ftxt = _open_mesh_file(meshdir, "faces")
    fm = re.search(r"^\s*(\d+)\s*\n\s*\(", ftxt, re.M)
    fstart = ftxt.index("(", fm.start())
    faces = [np.fromstring(s, dtype=np.int64, sep=" ")
             for s in re.findall(r"\(([^()]*)\)", ftxt[fstart + 1:])]
    fc = np.array([pts[f].mean(axis=0) for f in faces])
    own = _read_labels(_open_mesh_file(meshdir, "owner"))
    nei = _read_labels(_open_mesh_file(meshdir, "neighbour"))
    ncells = int(max(own.max(), nei.max() if nei.size else 0)) + 1
    acc = np.zeros((ncells, 3)); cnt = np.zeros(ncells)
    np.add.at(acc, own, fc[:len(own)]); np.add.at(cnt, own, 1.0)
    if nei.size:
        np.add.at(acc, nei, fc[:len(nei)]); np.add.at(cnt, nei, 1.0)
    return acc / np.maximum(cnt, 1)[:, None]


def ensure_cell_centres(case, preamble):
    """Cache cell centres per rank as 0/C-derived .npy (static mesh -> once)."""
    for r in rank_dirs(case):
        npy = os.path.join(r, "cellCentres.npy")
        if os.path.isfile(npy):
            continue
        cpath = os.path.join(r, "0", "C")
        C = read_field(cpath) if os.path.isfile(cpath) \
            else cell_centres_from_mesh(r)
        np.save(npy, C)


def field_times(case):
    """Write times present in every rank, sorted, excluding 0."""
    roots = rank_dirs(case)
    sets = []
    for r in roots:
        ts = set()
        for p in glob.glob(os.path.join(r, "*")):
            b = os.path.basename(p)
            if re.fullmatch(r"[0-9]+(\.[0-9]+)?(e[-+]?[0-9]+)?", b) \
               and os.path.isfile(os.path.join(p, "U")):
                ts.add(b)
        sets.append(ts)
    common = set.intersection(*sets) if sets else set()
    return sorted(common, key=float)


def load_time(case, t, alpha_name):
    """Concatenate C, U, alpha over ranks -> (x, U, alpha) point cloud."""
    X, Uv, Al = [], [], []
    for r in rank_dirs(case):
        c = np.load(os.path.join(r, "cellCentres.npy"))
        u = read_field(os.path.join(r, t, "U"))
        ap = os.path.join(r, t, alpha_name)
        a = read_field(ap) if os.path.isfile(ap) else np.zeros(len(c))
        if u.ndim == 1:                       # uniform
            u = np.tile(u, (len(c), 1))
        if np.ndim(a) == 0 or (a.ndim == 1 and a.size == 1):
            a = np.full(len(c), float(np.ravel(a)[0]))
        X.append(c); Uv.append(u); Al.append(a)
    return np.vstack(X), np.vstack(Uv), np.concatenate(Al)


def read_trace(case):
    """(t, max|U|) from the solver's per-step CSV."""
    for name in ("leiaSemiLagrangianLevelSetTwoPhaseFoam.csv",
                 "leiaLevelSetTwoPhaseFoam.csv"):
        p = os.path.join(case, name)
        if os.path.isfile(p):
            import csv
            rows = [r for r in csv.DictReader(open(p))
                    if r.get("maxMagU") not in (None, "")]
            return (np.array([float(r["TIME"]) for r in rows]),
                    np.array([float(r["maxMagU"]) for r in rows]))
    return np.array([]), np.array([])


def render(case, t, cfg, out_png):
    X, U, A = load_time(case, t, cfg["alpha_name"])
    zc = cfg["zc"]
    keep = np.abs(X[:, 2] - zc) < cfg["slab"]
    x, y = X[keep, 0], X[keep, 1]
    ux, uy = U[keep, 0], U[keep, 1]
    a = A[keep]
    mag = np.sqrt(U[keep, 0] ** 2 + U[keep, 1] ** 2 + U[keep, 2] ** 2)
    norm = cfg["norm"]
    R, xc, yc = cfg["radius"], cfg["xc"], cfg["yc"]
    mm = 1e3

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.4, 6.2))
    fig.patch.set_facecolor("white")

    # ---- left: glyphs over alpha, zoomed to the droplet ------------------
    zoom = 1.9 * R
    z = (np.abs(x - xc) < zoom) & (np.abs(y - yc) < zoom)
    axL.scatter((x[z] - xc) * mm, (y[z] - yc) * mm, c=a[z], cmap="Greys",
                vmin=0, vmax=1, s=cfg["dot"], marker="s", linewidths=0)
    pmax = max(mag[z].max(), 1e-30) if z.any() else 1.0
    rel = mag[z] / pmax
    glyph = cfg["glyph"]
    axL.quiver((x[z] - xc) * mm, (y[z] - yc) * mm,
               ux[z] / np.maximum(mag[z], 1e-30) * rel * glyph,
               uy[z] / np.maximum(mag[z], 1e-30) * rel * glyph,
               np.maximum(mag[z], norm.vmin), cmap=ARROWS, norm=norm,
               angles="xy", scale_units="xy", scale=1, width=0.0038,
               headwidth=3.4, headlength=4.0)
    circ = plt.Circle((0, 0), R * mm, fill=False, color=INTERFACE, lw=1.8)
    axL.add_patch(circ)
    axL.set_xlim(-zoom * mm, zoom * mm); axL.set_ylim(-zoom * mm, zoom * mm)
    axL.set_aspect("equal")
    axL.set_xlabel("x - x$_c$  [mm]", color=INK); axL.set_ylabel("y - y$_c$  [mm]", color=INK)
    axL.set_title("direction glyphs over $\\alpha$  (length: per-frame, "
                  "colour: shared log |U|)", fontsize=10, color=INK)

    # ---- right: |U| over the full slab -----------------------------------
    sc = axR.scatter((x - xc) * mm, (y - yc) * mm,
                     c=np.maximum(mag, norm.vmin), cmap="magma_r", norm=norm,
                     s=cfg["dot_full"], marker="s", linewidths=0)
    axR.add_patch(plt.Circle((0, 0), R * mm, fill=False, color=INTERFACE, lw=1.4))
    half = cfg["half"] * mm
    axR.set_xlim(-half, half); axR.set_ylim(-half, half); axR.set_aspect("equal")
    axR.set_xlabel("x - x$_c$  [mm]", color=INK)
    axR.set_title("|U| over the full domain slab", fontsize=10, color=INK)
    cb = fig.colorbar(sc, ax=axR, fraction=0.046, pad=0.03)
    cb.set_label("|U|  [m/s]", color=INK)

    tv = float(t)
    mx = mag.max()
    fig.suptitle(f"{cfg['label']}    t = {tv*1e3:.2f} ms    "
                 f"max|U| (slab) = {mx:.3e} m/s",
                 fontsize=12.5, color=INK, y=0.985)

    # ---- trace strip ------------------------------------------------------
    tt, uu = cfg["trace"]
    if tt.size:
        axp = fig.add_axes([0.09, 0.035, 0.84, 0.075])
        axp.semilogy(tt * 1e3, np.maximum(uu, 1e-12), color=INK, lw=1.2)
        axp.axvline(tv * 1e3, color=INTERFACE, lw=1.6)
        axp.set_xlim(0, max(tt.max(), tv) * 1e3)
        axp.set_ylabel("max|U|", fontsize=8, color=INK_MUTED)
        axp.set_xlabel("t [ms]", fontsize=8, color=INK_MUTED)
        axp.tick_params(labelsize=7, colors=INK_MUTED)
        for s in axp.spines.values():
            s.set_color(INK_MUTED)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.17, wspace=0.16)
    fig.savefig(out_png, dpi=cfg["dpi"])
    plt.close(fig)
    return mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--slab", type=float, default=1.0,
                    help="half-thickness of the z-slab, in cell sizes")
    ap.add_argument("--alpha-name", default="alpha.water")
    ap.add_argument("--times", default="all")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--label", default=None)
    ap.add_argument("--gif", action="store_true")
    ap.add_argument("--fps", type=int, default=3)
    ap.add_argument("--dpi", type=int, default=115)
    ap.add_argument("--preamble",
                    default="source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc")
    a = ap.parse_args()

    case = os.path.abspath(a.case)
    out_dir = a.out_dir or os.path.join(case, "figures")
    os.makedirs(out_dir, exist_ok=True)
    ensure_cell_centres(case, a.preamble)

    times = field_times(case)
    if a.times != "all":
        want = set(a.times.split(","))
        times = [t for t in times if t in want]
    if not times:
        raise SystemExit(f"no field times with U in {case}")

    # geometry + one shared colour scale over ALL frames (growth must be
    # visible as darkening, so the scale cannot be per-frame)
    X0, _, _ = load_time(case, times[0], a.alpha_name)
    xc, yc, zc = X0.mean(axis=0)
    half = 0.5 * (X0[:, 0].max() - X0[:, 0].min())
    ncell = len(X0)
    h = (2 * half) / max(round((ncell) ** (1.0 / 3.0)), 1)
    lo, hi = np.inf, 0.0
    for t in times:
        _, U, _ = load_time(case, t, a.alpha_name)
        m = np.linalg.norm(U, axis=1)
        m = m[m > 0]
        if m.size:
            lo = min(lo, np.percentile(m, 5)); hi = max(hi, m.max())
    lo = max(lo, hi * 1e-6) if np.isfinite(lo) else hi * 1e-6

    cfg = dict(radius=a.radius, slab=a.slab * h, xc=xc, yc=yc, zc=zc,
               half=half, alpha_name=a.alpha_name, dpi=a.dpi,
               norm=LogNorm(vmin=max(lo, 1e-12), vmax=max(hi, 1e-11)),
               trace=read_trace(case), dot=max(4.0, 900.0 * h / (2 * half)),
               dot_full=max(1.5, 320.0 * h / (2 * half)),
               glyph=0.9 * h * 1e3,
               label=a.label or os.path.basename(case))

    pngs = []
    for k, t in enumerate(times):
        p = os.path.join(out_dir, f"frame_{k:03d}_t{float(t):.5f}.png")
        mx = render(case, t, cfg, p)
        pngs.append(p)
        print(f"  t = {float(t)*1e3:7.3f} ms   max|U| = {mx:.4e}   -> {p}")

    if a.gif and len(pngs) > 1:
        gif = os.path.join(out_dir, "parasitic.gif")
        # GIF, not mp4: mp4 downloads fail for this user (memory: animations-as-gif)
        subprocess.run(["bash", "-lc",
                        f"ffmpeg -y -framerate {a.fps} -pattern_type glob -i "
                        f"'{out_dir}/frame_*.png' -vf "
                        f"'scale=1100:-1:flags=lanczos,split[s0][s1];"
                        f"[s0]palettegen[p];[s1][p]paletteuse' {gif}"],
                       check=False)
        print("  gif:", gif)


if __name__ == "__main__":
    main()
