#!/usr/bin/env python3
"""Render the PLIC surfaces of a DAVOF study side by side.

    python3 workflow/scripts/render_davof_interfaces.py studies/davof/sphereNormal3D \
        [--theme davof] [--alpha-source exactSphere] [--models davof plicRDF gradAlpha] \
        [--max-resolutions 3] [--field eNormal] [--elev 22 --azim -55]

Every case of the study written by leiaTestDavofNormal carries
postProcessing/davofInterface/<time>/plic.<model>.vtk: one polygon per interface
cell (the DAVOF plane, or the PLIC plane of a geometricVoF scheme run on the
identical alpha), with the per-cell error fields eNormal (|n_c - n_exact|) and
ePos (distance [m] of the polygon centroid from the exact surface) attached as
CELL_DATA. This script draws a grid of panels, rows = models, columns = the
coarsest --max-resolutions resolutions of one alpha source, the polygons
coloured by --field on one logarithmic scale per figure, and writes
figures/davof_plic_<study>.{pdf,png} into the theme data (paths.figs_dir).

Self-contained: the ASCII legacy polydata is parsed here, and the polygons are
projected orthographically (painter's algorithm) onto 2D axes, so neither a
VTK library nor mplot3d is needed (mplot3d is broken on hosts that mix an apt
and a pip matplotlib).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

MODEL_LABEL = {"davof": "DAVOF", "plicRDF": "plicRDF", "gradAlpha": "gradAlpha",
               "isoAlpha": "isoAlpha"}
FIELD_LABEL = {"eNormal": r"$|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
               "ePos": "distance of the polygon centroid from the surface [m]"}


def read_legacy_vtk_polydata(path):
    """ASCII legacy VTK polydata -> (points (n,3), [index arrays], {name: array})."""
    with open(path) as fh:
        tok = fh.read().split()
    i, n = 0, len(tok)
    points, polys, fields, ncell = None, [], {}, 0
    while i < n:
        t = tok[i]
        if t == "POINTS":
            npts = int(tok[i + 1]); i += 3
            points = np.array(tok[i:i + 3*npts], dtype=float).reshape(npts, 3); i += 3*npts
        elif t == "POLYGONS":
            npoly, size = int(tok[i + 1]), int(tok[i + 2]); i += 3
            data = np.array(tok[i:i + size], dtype=int); i += size
            j = 0
            for _ in range(npoly):
                k = data[j]; polys.append(data[j + 1:j + 1 + k]); j += 1 + k
        elif t == "CELL_DATA":
            ncell = int(tok[i + 1]); i += 2
        elif t == "SCALARS":
            name = tok[i + 1]; i += 4            # SCALARS name type ncomp
            if i < n and tok[i] == "LOOKUP_TABLE":
                i += 2
            fields[name] = np.array(tok[i:i + ncell], dtype=float); i += ncell
        else:
            i += 1
    return points, polys, fields


def view_matrix(elev_deg, azim_deg):
    """Rows: screen x, screen y, depth (towards the viewer) of an orthographic
    camera at elevation elev and azimuth azim (matplotlib's convention)."""
    el, az = np.radians(elev_deg), np.radians(azim_deg)
    view = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])
    up = np.array([0.0, 0.0, 1.0])
    right = np.cross(up, view); right /= np.linalg.norm(right)
    sup = np.cross(view, right)
    return np.vstack([right, sup, view])


def collect(study_dir, alpha_source, models, max_res):
    """[(N, {model: vtk path})] for the coarsest max_res resolutions, and the
    alpha sources seen."""
    cases, sources = {}, set()
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        d = os.path.dirname(meta)
        tok = json.load(open(meta)).get("tokens", {})
        sources.add(tok.get("DAVOF_ALPHA_SOURCE"))
        if alpha_source and tok.get("DAVOF_ALPHA_SOURCE") != alpha_source:
            continue
        try:
            N = int(float(tok.get("N_CELLS")))
        except (TypeError, ValueError):
            continue
        dirs = sorted(glob.glob(os.path.join(d, "postProcessing", "davofInterface", "*")))
        if not dirs:
            continue
        files = {m: os.path.join(dirs[-1], f"plic.{m}.vtk") for m in models}
        files = {m: f for m, f in files.items() if os.path.isfile(f)}
        if files:
            cases[N] = files
    return [(N, cases[N]) for N in sorted(cases)[:max_res]], sources


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="davof")
    ap.add_argument("--alpha-source", default=None,
                    help="DAVOF_ALPHA_SOURCE to draw (default: exactSphere if present, "
                         "else linearInterpolant, else planePhaseIndicator)")
    ap.add_argument("--models", nargs="+", default=["davof", "plicRDF", "gradAlpha"])
    ap.add_argument("--max-resolutions", type=int, default=3)
    ap.add_argument("--field", default="eNormal", choices=["eNormal", "ePos"])
    ap.add_argument("--elev", type=float, default=None, dest="elev_set",
                    help="camera elevation [deg] (default 22, or head-on for a planar surface)")
    ap.add_argument("--azim", type=float, default=None, dest="azim_set",
                    help="camera azimuth [deg] (default -55, or head-on for a planar surface)")
    a = ap.parse_args(argv)
    a.elev = 22.0 if a.elev_set is None else a.elev_set
    a.azim = -55.0 if a.azim_set is None else a.azim_set

    study = os.path.basename(os.path.normpath(a.study_dir))
    cases, sources = collect(a.study_dir, a.alpha_source, a.models, a.max_resolutions)
    if not cases and a.alpha_source is None:
        for pref in ("exactSphere", "linearInterpolant", "planePhaseIndicator"):
            if pref in sources:
                cases, _ = collect(a.study_dir, pref, a.models, a.max_resolutions)
                a.alpha_source = pref
                if cases:
                    break
    if not cases:
        print(f"[render] no plic.<model>.vtk files under {a.study_dir} "
              f"(alpha sources seen: {sorted(s for s in sources if s)}); nothing drawn")
        return 0
    if a.alpha_source is None:
        a.alpha_source = sorted(s for s in sources if s)[0]

    # Load everything first: one colour scale per figure.
    data, vmin, vmax = {}, np.inf, -np.inf
    for N, files in cases:
        for m, f in files.items():
            pts, polys, flds = read_legacy_vtk_polydata(f)
            e = flds.get(a.field, np.zeros(len(polys)))
            data[(N, m)] = (pts, polys, e)
            pos = e[e > 0]
            if pos.size:
                vmin, vmax = min(vmin, pos.min()), max(vmax, pos.max())
    if not np.isfinite(vmin):
        vmin, vmax = 1e-16, 1.0
    vmin = max(vmin, 1e-16)
    if vmax <= vmin:
        vmax = vmin*10
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("Blues")

    # The camera: --elev/--azim, unless the surface is (nearly) planar, in which
    # case it is viewed head-on along its mean normal (a plane seen obliquely
    # is a sliver). Newell normals of the coarsest DAVOF (or first) surface.
    elev, azim = a.elev, a.azim
    if not any(x is not None for x in (a.elev_set, a.azim_set)):
        N0, files0 = cases[0]
        m0 = "davof" if (N0, "davof") in data else next(m for m in a.models if (N0, m) in data)
        pts, polys, _ = data[(N0, m0)]
        mean_n, total = np.zeros(3), 0.0
        for p in polys:
            v = pts[p]
            nrm = np.zeros(3)
            for k in range(len(v)):
                p1, p2 = v[k], v[(k + 1) % len(v)]
                nrm += np.cross(p1, p2)
            area = np.linalg.norm(nrm)
            if area > 0:
                mean_n += nrm; total += area
        if total > 0 and np.linalg.norm(mean_n)/total > 0.9:
            n_hat = mean_n/np.linalg.norm(mean_n)
            elev = float(np.degrees(np.arcsin(np.clip(n_hat[2], -1, 1))))
            azim = float(np.degrees(np.arctan2(n_hat[1], n_hat[0])))
            print(f"[render] planar surface: viewing along its normal (elev {elev:.1f}, azim {azim:.1f})")
    M = view_matrix(elev, azim)

    models = [m for m in a.models if any((N, m) in data for N, _ in cases)]
    nrow, ncol = len(models), len(cases)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1*ncol + 1.3, 3.0*nrow), squeeze=False)
    for r, m in enumerate(models):
        for c, (N, _) in enumerate(cases):
            ax = axes[r][c]
            ax.set_axis_off()
            if (N, m) not in data:
                continue
            pts, polys, e = data[(N, m)]
            proj = pts @ M.T                    # (n, 3): screen x, screen y, depth
            depth = np.array([proj[p, 2].mean() for p in polys])
            order = np.argsort(depth)           # far first (painter's algorithm)
            verts = [proj[polys[k], :2] for k in order]
            colors = cmap(norm(np.maximum(e[order], vmin)))
            ax.add_collection(PolyCollection(verts, facecolors=colors, edgecolors="k",
                                             linewidths=0.08))
            lo, hi = proj[:, :2].min(axis=0), proj[:, :2].max(axis=0)
            span = (hi - lo).max()*1.04
            mid = 0.5*(lo + hi)
            ax.set_xlim(mid[0] - span/2, mid[0] + span/2)
            ax.set_ylim(mid[1] - span/2, mid[1] + span/2)
            ax.set_aspect("equal")
            ax.set_title(f"{MODEL_LABEL.get(m, m)}, N = {N} ({len(polys)} cells)",
                         fontsize=9)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.7, pad=0.02)
    cbar.set_label(FIELD_LABEL.get(a.field, a.field), fontsize=9)
    fig.suptitle(f"PLIC surfaces on the identical alpha ({a.alpha_source}), "
                 f"coloured by the per-cell {a.field}", fontsize=10)
    figs = paths.figs_dir(a.theme)
    os.makedirs(figs, exist_ok=True)
    for ext in ("pdf", "png"):
        out = os.path.join(figs, f"davof_plic_{study}.{ext}")
        fig.savefig(out, dpi=170, bbox_inches="tight")
        print(f"[render] wrote {out}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
