#!/usr/bin/env python3
"""WHERE does the stationary-droplet parasitic current / curvature blow-up start?

Reads ONE solved two-phase case (uniform 2D blockMesh, N x N x 1) at a set of
write times bracketing the divergence onset and localizes the runaway in space:

  * |u|            -- parasitic velocity magnitude (which phase, which angle?)
  * kappa - 1/R    -- reconstructed-curvature error in the snGrad(alpha) band
  * |grad psi|     -- signed-distance quality (1 = clean SDF; collapse = drift),
                      computed from psi by centred finite differences on the grid

For each requested time it renders a montage row (|u|, curvature error, |grad psi|)
zoomed on the interface with the psi=0 contour, and an ANGULAR-PROFILE panel: the
interface-band cells binned by polar angle theta about the drop centre, so we can
read off whether the blow-up seeds at the poles / equator / 45-deg diagonals
(grid-aligned vs grid-diagonal), and a signed split showing the light (air,
psi>0 outside) vs heavy (water, psi<0 inside) side.

Self-contained (VTK OpenFOAM reader + numpy + matplotlib Agg); no leia imports.

Usage:
  python3 make_droplet_localization.py <case_dir> --times 0.03 0.04 0.05 0.06 \
      [--out sl_droplet_localization.png] [--box 0.01 --cx 0.005 --cy 0.005 --R 1e-3]
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def _internal(mb, array):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        blk = it.GetCurrentDataObject()
        if blk is not None and blk.IsA("vtkUnstructuredGrid") \
           and blk.GetCellData().GetArray(array) is not None:
            return blk
        it.GoToNextItem()
    return None


def read_case(case_dir):
    """Open the case once; return (reader, sorted list of time values)."""
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(foam)
    reader.CreateCellToPointOff()
    reader.SetSkipZeroTime(0)
    reader.Update()
    reader.EnableAllCellArrays()
    reader.Update()
    tv = reader.GetTimeValues()
    tvals = sorted(tv.GetValue(i) for i in range(tv.GetNumberOfTuples())) if tv else []
    return reader, tvals


def fields_at(reader, t):
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
        "alpha": arr("alpha.water"), "psi": arr("psi"),
        "Ux": U[:, 0] if U is not None else None,
        "Uy": U[:, 1] if U is not None else None,
        "kappa": arr("kappa"),
    }


def to_grid(f, N, box):
    """Bin scattered cell-centred data onto the uniform N x N grid (y row, x col)."""
    h = box / N
    ix = np.clip(np.round((f["x"] - 0.5 * h) / h).astype(int), 0, N - 1)
    iy = np.clip(np.round((f["y"] - 0.5 * h) / h).astype(int), 0, N - 1)

    def grid(v):
        g = np.full((N, N), np.nan)
        g[iy, ix] = v
        return g

    G = {k: grid(f[k]) for k in ("alpha", "psi", "Ux", "Uy", "kappa")}
    G["h"] = h
    G["Umag"] = np.hypot(G["Ux"], G["Uy"])
    # |grad psi| from psi via centred differences on the uniform grid
    gy, gx = np.gradient(G["psi"], h, h)
    G["gradPsiMag"] = np.hypot(gx, gy)
    xs = (np.arange(N) + 0.5) * h
    G["X"], G["Y"] = np.meshgrid(xs, xs)
    return G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--times", type=float, nargs="+", required=True)
    ap.add_argument("--out", default="sl_droplet_localization.png")
    ap.add_argument("--box", type=float, default=0.01)
    ap.add_argument("--cx", type=float, default=0.005)
    ap.add_argument("--cy", type=float, default=0.005)
    ap.add_argument("--R", type=float, default=1.0e-3)
    ap.add_argument("--zoom", type=float, default=2.2, help="view half-width in units of R")
    args = ap.parse_args()

    reader, tvals = read_case(args.case_dir)
    if not tvals:
        print("[loc] no time values"); return 1
    # snap requested times to nearest available
    times = [min(tvals, key=lambda x: abs(x - tt)) for tt in args.times]
    kap_exact = 1.0 / args.R
    zh = args.zoom * args.R
    x0, x1 = args.cx - zh, args.cx + zh
    y0, y1 = args.cy - zh, args.cy + zh

    nrow = len(times)
    fig, axes = plt.subplots(nrow, 4, figsize=(17, 4.0 * nrow), squeeze=False)
    for r, t in enumerate(times):
        f = fields_at(reader, t)
        if f is None or f["Ux"] is None:
            for c in range(4):
                axes[r][c].text(0.5, 0.5, "n/a", ha="center")
            continue
        # infer N from cell count (2D square)
        N = int(round(np.sqrt(f["x"].size)))
        G = to_grid(f, N, args.box)
        ext = [0, args.box, 0, args.box]

        # --- panel 0: |u| (log) ---
        ax = axes[r][0]
        um = np.nan_to_num(G["Umag"])
        vmax = max(um.max(), 1e-12)
        vmin = max(vmax * 1e-4, 1e-12)
        im = ax.imshow(np.clip(um, vmin, None), origin="lower", extent=ext,
                       cmap="inferno", norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.contour(G["X"], G["Y"], G["psi"], levels=[0], colors="cyan", linewidths=1.0)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label(r"$|\mathbf{u}|$ [m/s]")
        ax.set_title(rf"$t={t:g}$  |u|  (max={um.max():.2e})")

        # --- panel 1: curvature error in band ---
        ax = axes[r][1]
        kerr = np.where(np.abs(G["kappa"]) > 1e-9, G["kappa"] - kap_exact, np.nan)
        amax = np.nanmax(np.abs(kerr)) if np.isfinite(kerr).any() else kap_exact
        norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
        im = ax.imshow(kerr, origin="lower", extent=ext, cmap="RdBu_r", norm=norm)
        ax.contour(G["X"], G["Y"], G["psi"], levels=[0], colors="k", linewidths=1.0)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label(r"$\kappa - 1/R$ [1/m]")
        ax.set_title(rf"curvature error (|max|={np.nanmax(np.abs(kerr)):.0f})")

        # --- panel 2: |grad psi| ---
        ax = axes[r][2]
        gpm = G["gradPsiMag"]
        norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=2.0)
        im = ax.imshow(np.clip(gpm, 0, 2), origin="lower", extent=ext,
                       cmap="PuOr", norm=norm)
        ax.contour(G["X"], G["Y"], G["psi"], levels=[0], colors="k", linewidths=1.0)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label(r"$|\nabla\psi|$")
        ax.set_title(r"$|\nabla\psi|$ (1 = clean SDF)")

        for c in range(3):
            axes[r][c].set_xlim(x0, x1); axes[r][c].set_ylim(y0, y1)
            axes[r][c].set_aspect("equal")

        # --- panel 3: angular profiles over the interface band ---
        ax = axes[r][3]
        psi = G["psi"]; h = G["h"]
        band = np.abs(psi) < 2.0 * h
        th = np.degrees(np.arctan2(G["Y"] - args.cy, G["X"] - args.cx))
        thb = th[band]
        order = np.argsort(thb)
        umb = G["Umag"][band][order]
        keb = np.where(np.abs(G["kappa"][band]) > 1e-9,
                       np.abs(G["kappa"][band] - kap_exact), np.nan)[order]
        thb = thb[order]
        ax.plot(thb, umb / max(umb.max(), 1e-30), ".", ms=3, color="#B22222",
                label=r"$|\mathbf{u}|$ (norm.)")
        ax.plot(thb, keb / max(np.nanmax(keb), 1e-30), ".", ms=3, color="#0072B2",
                label=r"$|\kappa-1/R|$ (norm.)")
        for a in (-90, 0, 90):
            ax.axvline(a, color="0.8", lw=0.6)
        ax.set_xlabel(r"interface angle $\theta$ [deg]  (0=+x pole, 90=+y pole)")
        ax.set_ylabel("normalised")
        ax.set_title("angular profile (band)")
        ax.legend(fontsize=8, frameon=False)
        ax.set_xlim(-180, 180)

    fig.suptitle(f"Divergence localization: {os.path.basename(os.path.normpath(args.case_dir))}",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(args.out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[loc] wrote {args.out}  (times: {', '.join(f'{t:g}' for t in times)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
