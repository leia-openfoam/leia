#!/usr/bin/env python3
"""Stationary-droplet FIELD montages across mesh resolution (2D).

Reads the solved cases of the ``dropletFields2D`` study (one case per resolution,
N = 32 / 64 / 128) and renders two paper/slide figures, coarse -> fine left to
right, with the mesh drawn on every panel so the resolution is visible:

  sl_droplet_velocity_field.png
      volume-fraction background (alpha.water: drop vs ambient) + velocity glyphs
      (quiver, sub-sampled, per-panel normalised) at a representative time -- the
      parasitic-current pattern and its amplitude (annotated max|u|) vs resolution.

  sl_droplet_curvature_field.png
      the reconstructed mean curvature kappa in the snGrad(alpha) band, on a fixed
      colour scale centred on the exact 1/R, at the same time -- shows the curvature
      is a clean ~1/R ring on the coarse meshes and roughens on the finer one.

Fields are read with VTK's vtkOpenFOAMReader (same robust reader as
make_isosurface_fig.py); cell-centred data + cell centres are pulled to numpy and
drawn with matplotlib (Agg, no GL context needed). The target time is the last
reconstructed write time by default (the study writes t=0 and one representative
t>0), overridable with --time.

Usage:  python3 make_droplet_fields_fig.py <study_dir> [--time T]
"""
import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import paths

# Domain / droplet geometry of cases/stationaryDroplet2D (SI units, metres).
BOX = 0.01
CENTRE = (0.005, 0.005)
R = 1.0e-3
KAPPA_EXACT = 1.0 / R              # 2D mean curvature of the circle [1/m]
ZOOM_HALF = 2.2 * R                # half-width of the view box around the drop


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def study_cases(study_dir):
    """[(N, case_dir)] for the study, coarse -> fine; skip cases stuck at t=0."""
    rows = []
    for cp in sorted(glob.glob(os.path.join(study_dir, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp))
        n = _f(d.get("tokens", {}).get("N_CELLS"))
        cdir = os.path.dirname(cp)
        if n:
            rows.append((int(n), cdir))
    rows.sort(key=lambda r: r[0])
    return rows


def _internal(mb, array):
    """First unstructured-grid block carrying the named CELL array."""
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        blk = it.GetCurrentDataObject()
        if blk is not None and blk.IsA("vtkUnstructuredGrid") \
           and blk.GetCellData().GetArray(array) is not None:
            return blk
        it.GoToNextItem()
    return None


def read_fields(case_dir, want_time=None, kappa_name="kappa"):
    """Return dict {x,y,alpha,U,kappa,time} of cell-centred data at want_time
    (nearest reconstructed time; default the last). ``kappa_name`` selects the
    curvature field to read (``kappa`` from the two-phase run, or ``kappaDiv``
    from the standalone clean-SDF curvature test)."""
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(foam)
    reader.CreateCellToPointOff()          # keep native cell-centred data
    reader.SetSkipZeroTime(0)
    reader.Update()
    reader.EnableAllCellArrays()
    reader.Update()
    tv = reader.GetTimeValues()
    tvals = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())] if tv else []
    if not tvals:
        return None
    t = (min(tvals, key=lambda x: abs(x - want_time)) if want_time is not None
         else max(tvals))
    reader.UpdateTimeStep(t)
    mb = reader.GetOutput()
    internal = _internal(mb, "alpha.water")
    if internal is None:
        return None
    cc = vtk.vtkCellCenters(); cc.SetInputData(internal); cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    cd = internal.GetCellData()

    def arr(name):
        a = cd.GetArray(name)
        return vtk_to_numpy(a) if a is not None else None

    U = arr("U")
    return {
        "x": pts[:, 0], "y": pts[:, 1],
        "alpha": arr("alpha.water"),
        "psi": arr("psi"),                 # level set: psi=0 is the true interface
        "Ux": U[:, 0] if U is not None else None,
        "Uy": U[:, 1] if U is not None else None,
        "kappa": arr(kappa_name),
        "time": t,
    }


def _draw_mesh(ax, N, x0, x1, y0, y1):
    """Uniform blockMesh grid lines within the view box (light grey)."""
    h = BOX / N
    xs = np.arange(0.0, BOX + 0.5 * h, h)
    ys = xs
    for xv in xs[(xs >= x0) & (xs <= x1)]:
        ax.plot([xv, xv], [y0, y1], color="0.75", lw=0.3, zorder=0)
    for yv in ys[(ys >= y0) & (ys <= y1)]:
        ax.plot([x0, x1], [yv, yv], color="0.75", lw=0.3, zorder=0)


def _view_box():
    cx, cy = CENTRE
    return cx - ZOOM_HALF, cx + ZOOM_HALF, cy - ZOOM_HALF, cy + ZOOM_HALF


def velocity_figure(cases, out):
    x0, x1, y0, y1 = _view_box()
    n = len(cases)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.3), squeeze=False)
    tref = None
    for j, (N, cdir) in enumerate(cases):
        ax = axes[0][j]
        f = read_fields(cdir)
        if f is None:
            ax.text(0.5, 0.5, "n/a", ha="center", va="center"); continue
        tref = f["time"]
        # alpha background: fill the drop (alpha>1/2) a solid light blue, ambient white.
        # A single band (not a 0..1 ramp) reads cleanly for the near-sharp indicator.
        ax.tricontourf(f["x"], f["y"], f["alpha"], levels=[0.5, 2.0],
                       colors=["#bcd4ee"], zorder=1)
        _draw_mesh(ax, N, x0, x1, y0, y1)
        # interface outline = the level-set zero contour psi=0 (the true interface
        # location; smoother than the alpha=1/2 phase-indicator step).
        if f.get("psi") is not None:
            ax.tricontour(f["x"], f["y"], f["psi"], levels=[0.0],
                          colors="k", linewidths=1.0, zorder=3)
        # velocity glyphs, sub-sampled to ~18 arrows across the box, normalised so
        # the PATTERN shows regardless of amplitude (annotated separately).
        m = (f["x"] >= x0) & (f["x"] <= x1) & (f["y"] >= y0) & (f["y"] <= y1)
        xs, ys = f["x"][m], f["y"][m]
        us, vs = f["Ux"][m], f["Uy"][m]
        mag = np.hypot(us, vs)
        umax = float(mag.max()) if mag.size else 0.0
        # stride: aim for ~20 glyphs per axis inside the box
        span_cells = int(round(2 * ZOOM_HALF / (BOX / N)))
        stride = max(1, span_cells // 20)
        # sub-sample on the structured index is awkward on scattered data; instead
        # thin by a coarse spatial grid: keep the cell nearest each grid node.
        xi = np.linspace(x0, x1, 20); yi = np.linspace(y0, y1, 20)
        keep = []
        used = set()
        for gx in xi:
            for gy in yi:
                k = int(np.argmin((xs - gx) ** 2 + (ys - gy) ** 2))
                if k not in used and mag[k] > 0:
                    used.add(k); keep.append(k)
        keep = np.array(keep, dtype=int) if keep else np.array([], dtype=int)
        if keep.size and umax > 0:
            ax.quiver(xs[keep], ys[keep], us[keep] / umax, vs[keep] / umax,
                      color="#B22222", scale=22, width=0.005,
                      pivot="mid", zorder=4)
        ax.set_title(rf"$N={N}$   $\max|\mathbf{{u}}|={umax:.1e}$ m/s", fontsize=11)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ttl = rf"Parasitic current over volume fraction ($t={tref:g}$ s)" if tref else \
        "Parasitic current over volume fraction"
    fig.suptitle(ttl, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"[fields] wrote {out}")


def curvature_figure(cases, out, kappa_name="kappa", clean=False):
    x0, x1, y0, y1 = _view_box()
    n = len(cases)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 4.6), squeeze=False)
    # fixed diverging scale centred on the exact 1/R so panels are comparable.
    norm = TwoSlopeNorm(vmin=0.0, vcenter=KAPPA_EXACT, vmax=2.0 * KAPPA_EXACT)
    tref = None
    sc = None
    for j, (N, cdir) in enumerate(cases):
        ax = axes[0][j]
        f = read_fields(cdir, kappa_name=kappa_name)
        if f is None or f["kappa"] is None:
            ax.text(0.5, 0.5, "n/a", ha="center", va="center"); continue
        tref = f["time"]
        k = f["kappa"]
        band = np.abs(k) > 1e-9                     # kappa is nonzero only in the band
        _draw_mesh(ax, N, x0, x1, y0, y1)
        if band.any():
            sc = ax.scatter(f["x"][band], f["y"][band], c=k[band], cmap="RdBu_r",
                            norm=norm, s=(140.0 / N) ** 2 * 6, marker="s", zorder=3)
        # level-set zero contour psi=0 = the true interface, drawn ON TOP of the
        # curvature markers so its location is clear.
        if f.get("psi") is not None:
            ax.tricontour(f["x"], f["y"], f["psi"], levels=[0.0],
                          colors="k", linewidths=1.1, zorder=5)
        ax.set_title(rf"$N={N}$", fontsize=11)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if sc is not None:
        cb = fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.8, pad=0.02)
        cb.set_label(r"$\kappa$ [1/m]  (exact $1/R=1000$)")
    if clean:
        ttl = "Reconstructed curvature on the initialised interface (signed distance)"
    else:
        ttl = rf"Reconstructed curvature in the band ($t={tref:g}$ s)" if tref else \
            "Reconstructed curvature in the band"
    fig.suptitle(ttl, fontsize=13)
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"[fields] wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study_dir")
    ap.add_argument("--time", type=float, default=None,
                    help="target time (default: last reconstructed write)")
    args = ap.parse_args()
    cases = study_cases(args.study_dir)
    if not cases:
        print(f"[fields] no cases in {args.study_dir}; skip"); return 1
    figs = paths.figs_dir("semi-lagrangian-level-set")
    velocity_figure(cases, os.path.join(figs, "sl_droplet_velocity_field.png"))

    # Curvature QUALITY is a clean-signed-distance property, so render it from the
    # standalone leiaTestMeanCurvature study (curvatureDroplet2D, field kappaDiv on
    # the initialised interface) -- there the curvature CONVERGES with refinement
    # (N=128 is the smoothest). Rendering it from the two-phase run at t>0 instead
    # would conflate resolution with run-time drift (the fine mesh degrades first),
    # misleadingly showing N=128 as the noisiest. Fall back to the run-time kappa
    # only if the curvature study is absent.
    curv_study = os.path.join(os.path.dirname(os.path.normpath(args.study_dir)),
                              "curvatureDroplet2D")
    curv_cases = ([c for c in study_cases(curv_study) if c[0] in (32, 64, 128)]
                  if os.path.isdir(curv_study) else [])
    out_kappa = os.path.join(figs, "sl_droplet_curvature_field.png")
    if curv_cases:
        curvature_figure(curv_cases, out_kappa, kappa_name="kappaDiv", clean=True)
    else:
        print("[fields] curvatureDroplet2D not found; curvature from run-time kappa")
        curvature_figure(cases, out_kappa)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
