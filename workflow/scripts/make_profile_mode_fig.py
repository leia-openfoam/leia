#!/usr/bin/env python3
"""SMOKING-GUN test of the profile-mode aliasing theory of the parasitic runaway.

Theory (verified on synthetic fields): the anchored 3x3 quadratic WLS fit that
evaluates the symbolic curvature kappa = (tr(H)|g|^2 - g.H.g)/|g|^3 is exactly
invariant to SMOOTH normal-profile changes of psi (any error confined to H_nn
cancels identically), but a GRID-SCALE (2h-wavelength) oscillation of the profile
-- the level-set compression wave pumped by the parasitic strain -- aliases into
the tangential Hessian components and corrupts kappa with sign-indefinite O(1/h)
errors, worst where the interface normal is DIAGONAL to the grid (incommensurate
sampling), while the grid-aligned poles cancel (Nyquist -> pure H_nn).

This script tests that mechanism on the REAL fields of a diverging run:
  (1) per band cell, the 2h normal-mode amplitude of psi:
        A2h = |psi(x+h n) - 2 psi(x) + psi(x - h n)| / h
      (second difference along the local unit normal n = grad psi/|grad psi|,
      sampled by bilinear interpolation; for a smooth SDF A2h ~ h/R << 1,
      for a 2h profile mode A2h ~ 2*eps).
  (2) the curvature APPROXIMATION error |kappa - 1/r| (vs the local isosurface
      curvature; kappa read from the run).
  (3) an OFFLINE REFIT: for the top-|error| cells, feed the MEASURED 3x3 psi
      stencil into a replica of the solver's anchored WLS fit and verify it
      reproduces the written kappa (mechanism, not just correlation).

Outputs a 4-panel figure per requested time: A2h map, kappa-error map, scatter
A2h vs |kappa err| coloured by interface angle (pole vs diagonal), and the
refit-vs-written kappa parity plot. Prints correlation + refit statistics.

Usage:
  python3 make_profile_mode_fig.py <case_dir> --times 0.06 0.09 0.1 [--out fig.png]
"""
import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import vtk
from vtk.util.numpy_support import vtk_to_numpy

BOX, R, CX, CY = 0.01, 1.0e-3, 0.005, 0.005


def _internal(mb, a):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        b = it.GetCurrentDataObject()
        if b is not None and b.IsA("vtkUnstructuredGrid") \
           and b.GetCellData().GetArray(a) is not None:
            return b
        it.GoToNextItem()
    return None


def open_case(case_dir):
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    r = vtk.vtkOpenFOAMReader(); r.SetFileName(foam)
    r.CreateCellToPointOff(); r.SetSkipZeroTime(0); r.Update()
    r.EnableAllCellArrays(); r.Update()
    tv = r.GetTimeValues()
    ts = sorted(tv.GetValue(i) for i in range(tv.GetNumberOfTuples())) if tv else []
    return r, ts


def grids_at(reader, t, N, h):
    reader.UpdateTimeStep(t)
    ig = _internal(reader.GetOutput(), "psi")
    cc = vtk.vtkCellCenters(); cc.SetInputData(ig); cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    ix = np.clip(np.round((pts[:, 0] - 0.5*h)/h).astype(int), 0, N - 1)
    iy = np.clip(np.round((pts[:, 1] - 0.5*h)/h).astype(int), 0, N - 1)

    def G(name):
        a = ig.GetCellData().GetArray(name)
        g = np.full((N, N), np.nan)
        if a is not None:
            g[iy, ix] = vtk_to_numpy(a)
        return g

    return G("psi"), G("kappa")


def bilinear(F, x, y, h, N):
    """Bilinear sample of cell-centred grid F at physical point (x, y)."""
    fx = x/h - 0.5; fy = y/h - 0.5
    i0 = int(np.floor(fx)); j0 = int(np.floor(fy))
    if i0 < 0 or j0 < 0 or i0 + 1 >= N or j0 + 1 >= N:
        return np.nan
    tx = fx - i0; ty = fy - j0
    return ((1-tx)*(1-ty)*F[j0, i0] + tx*(1-ty)*F[j0, i0+1]
            + (1-tx)*ty*F[j0+1, i0] + tx*ty*F[j0+1, i0+1])


# ---- offline replica of the solver's anchored WLS fit (3x3 CPC, IDW 1/|d|) ----
def basis(dx, dy):
    return np.array([dx, dy, 0.5*dx*dx, dx*dy, 0.5*dy*dy])  # gx gy Hxx Hxy Hyy


def refit_kappa(psi, i, j, h):
    """Anchored quadratic WLS on the MEASURED 3x3 stencil -> symbolic kappa."""
    A, W, ds = [], [], []
    pc = psi[j, i]
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            v = psi[j+dj, i+di]
            if not np.isfinite(v):
                return np.nan
            A.append(basis(di*h, dj*h))
            W.append(1.0/math.hypot(di*h, dj*h))
            ds.append(v - pc)
    A = np.array(A); W = np.array(W); ds = np.array(ds)
    WA = A*W[:, None]; col = np.linalg.norm(WA, axis=0)
    cf, *_ = np.linalg.lstsq(WA/col, ds*W, rcond=None); cf = cf/col
    g = np.array([cf[0], cf[1]]); H = np.array([[cf[2], cf[3]], [cf[3], cf[4]]])
    gm = np.linalg.norm(g)
    if gm < 1e-30:
        return np.nan
    return (np.trace(H)*gm*gm - g @ (H @ g))/gm**3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--times", type=float, nargs="+", required=True)
    ap.add_argument("--out", default="sl_profile_mode.png")
    args = ap.parse_args()

    reader, tvals = open_case(args.case_dir)
    if not tvals:
        print("[profile] no time values"); return 1
    times = [min(tvals, key=lambda x: abs(x - tt)) for tt in args.times]

    nrow = len(times)
    fig, axes = plt.subplots(nrow, 4, figsize=(18, 4.2*nrow), squeeze=False)
    zh = 1.6e-3
    x0, x1, y0, y1 = CX - zh, CX + zh, CY - zh, CY + zh

    for r, t in enumerate(times):
        # infer N from the first grid read (assume square 2D)
        reader.UpdateTimeStep(t)
        ig = _internal(reader.GetOutput(), "psi")
        ncell = ig.GetNumberOfCells()
        N = int(round(math.sqrt(ncell)))
        h = BOX/N
        psi, kap = grids_at(reader, t, N, h)

        xs = (np.arange(N) + 0.5)*h
        X, Y = np.meshgrid(xs, xs)
        rr = np.hypot(X - CX, Y - CY)
        gy, gx = np.gradient(psi, h, h)
        gm = np.hypot(gx, gy)
        band = (np.abs(kap) > 1e-9) & (gm > 0.1)

        # (1) 2h normal-mode amplitude: second difference along the unit normal
        A2h = np.full((N, N), np.nan)
        js, is_ = np.where(band)
        for j, i in zip(js, is_):
            nx, ny = gx[j, i]/gm[j, i], gy[j, i]/gm[j, i]
            pp = bilinear(psi, X[j, i] + h*nx, Y[j, i] + h*ny, h, N)
            pm = bilinear(psi, X[j, i] - h*nx, Y[j, i] - h*ny, h, N)
            if np.isfinite(pp) and np.isfinite(pm):
                A2h[j, i] = abs(pp - 2*psi[j, i] + pm)/h

        # (2) curvature approximation error vs the local isosurface curvature
        kerr = np.where(band, np.abs(kap - 1.0/rr), np.nan)
        th = np.degrees(np.arctan2(Y - CY, X - CX))
        # distance (in deg) to the nearest grid-aligned pole (0/90/180/270)
        dpole = np.minimum.reduce([np.abs(((th - c + 180) % 360) - 180)
                                   for c in (0, 90, 180, -90)])

        ext = [0, BOX, 0, BOX]
        ax = axes[r][0]
        im = ax.imshow(A2h, origin="lower", extent=ext, cmap="viridis",
                       norm=LogNorm(vmin=max(np.nanmin(A2h), 1e-4),
                                    vmax=max(np.nanmax(A2h), 1e-3)))
        ax.contour(X, Y, psi, levels=[0], colors="w", linewidths=0.8)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label(r"$A_{2h}$ (2h profile mode)")
        ax.set_title(rf"$t={t:g}$: 2h normal-mode amplitude")

        ax = axes[r][1]
        im = ax.imshow(kerr, origin="lower", extent=ext, cmap="inferno",
                       norm=LogNorm(vmin=10, vmax=max(np.nanmax(kerr), 100)))
        ax.contour(X, Y, psi, levels=[0], colors="c", linewidths=0.8)
        fig.colorbar(im, ax=ax, shrink=0.8).set_label(r"$|\kappa - 1/r|$ [1/m]")
        ax.set_title("curvature approximation error")
        for c in range(2):
            axes[r][c].set_xlim(x0, x1); axes[r][c].set_ylim(y0, y1)
            axes[r][c].set_aspect("equal")
            axes[r][c].set_xticks([]); axes[r][c].set_yticks([])

        # (3) scatter: A2h vs kappa error, coloured by pole distance
        ax = axes[r][2]
        m = band & np.isfinite(A2h) & np.isfinite(kerr) & (kerr > 0) & (A2h > 0)
        sc = ax.scatter(A2h[m], kerr[m], c=dpole[m], cmap="coolwarm", s=8,
                        vmin=0, vmax=45)
        fig.colorbar(sc, ax=ax, shrink=0.8).set_label(
            r"angle to nearest pole [deg] (45 = diagonal)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$A_{2h}$"); ax.set_ylabel(r"$|\kappa-1/r|$ [1/m]")
        # correlation on logs
        la, lk = np.log(A2h[m]), np.log(kerr[m])
        cc_all = np.corrcoef(la, lk)[0, 1] if m.sum() > 3 else np.nan
        diag = m & (dpole > 25)
        cc_diag = (np.corrcoef(np.log(A2h[diag]), np.log(kerr[diag]))[0, 1]
                   if diag.sum() > 3 else np.nan)
        ax.set_title(f"corr(log-log): all={cc_all:.2f}, diagonal={cc_diag:.2f}")

        # (4) offline refit parity for the top-40 |kappa err| cells
        ax = axes[r][3]
        order = np.argsort(np.where(m, -kerr, np.inf), axis=None)[:40]
        wr, rf = [], []
        for k in order:
            j, i = np.unravel_index(k, kerr.shape)
            if not m[j, i]:
                continue
            kr = refit_kappa(psi, i, j, h)
            if np.isfinite(kr):
                wr.append(kap[j, i]); rf.append(kr)
        wr = np.array(wr); rf = np.array(rf)
        if wr.size:
            ax.plot(wr, rf, "o", ms=4, color="#B22222", alpha=0.7)
            lim = max(np.abs(np.concatenate([wr, rf])))*1.05
            ax.plot([-lim, lim], [-lim, lim], "k-", lw=0.7)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            rel = np.abs(rf - wr)/np.maximum(np.abs(wr), 1e-30)
            frac_sign = np.mean(np.sign(rf) == np.sign(wr))
            ax.set_title(f"refit parity: med rel diff {np.median(rel)*100:.1f}%, "
                         f"sign match {frac_sign*100:.0f}%")
            print(f"[profile] t={t:g}: corr all={cc_all:.3f} diag={cc_diag:.3f}; "
                  f"refit med rel diff={np.median(rel)*100:.2f}% "
                  f"sign match={frac_sign*100:.0f}% (n={wr.size}); "
                  f"max A2h={np.nanmax(A2h):.3f} (smooth-SDF ~ {h/R:.3f})")
        ax.set_xlabel(r"$\kappa$ written by the solver")
        ax.set_ylabel(r"$\kappa$ offline anchored-WLS refit")

    fig.suptitle("Profile-mode aliasing test: 2h normal-mode amplitude vs curvature "
                 "error (theory: correlated, diagonal-dominant, refit-reproducible)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[profile] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
