#!/usr/bin/env python3
"""Angular structure of the CLEAN-SDF reconstructed-curvature error (pole diagnosis).

Reads the leiaTestMeanCurvature study (curvatureDroplet2D: one case per resolution,
curvature evaluated on the INITIALISED circle signed-distance field -- no run-time
drift) and localizes the curvature error by interface polar angle theta, so we can
see WHERE around the circle the approximation is worst (grid-aligned poles theta =
0/90/180/270 vs grid-diagonals 45/135/...).

For each resolution it plots, for the band cells (|kappa|>0), the SIGNED relative
error (kappa - 1/R)/(1/R) vs theta for the no-extension cell-centre curvature
(kappaNoExt) -- the variant the solver now uses -- with kappaDiv (foot-point
extension) and kappaLap (tr H) overlaid for reference, coloured by the signed
normal distance psi (inside/outside). Prints a pole-vs-diagonal error summary.

Usage: python3 make_curvature_pole.py <curvatureDroplet2D_study_dir> [--out fig.png]
"""
import argparse, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

R = 1.0e-3
CX, CY = 0.005, 0.005
KAP = 1.0 / R
BOX = 0.01


def _internal(mb, array):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        blk = it.GetCurrentDataObject()
        if blk is not None and blk.IsA("vtkUnstructuredGrid") \
           and blk.GetCellData().GetArray(array) is not None:
            return blk
        it.GoToNextItem()
    return None


def read(case_dir, fields):
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    r = vtk.vtkOpenFOAMReader(); r.SetFileName(foam)
    r.CreateCellToPointOff(); r.SetSkipZeroTime(0); r.Update()
    r.EnableAllCellArrays(); r.Update()
    tv = r.GetTimeValues()
    tvals = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())] if tv else []
    if not tvals:
        return None
    r.UpdateTimeStep(max(tvals)); mb = r.GetOutput()
    key = next((f for f in fields if _internal(mb, f) is not None), None)
    ig = _internal(mb, key) if key else _internal(mb, "psi")
    if ig is None:
        return None
    cc = vtk.vtkCellCenters(); cc.SetInputData(ig); cc.Update()
    pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
    cd = ig.GetCellData()
    avail = [cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())]

    def arr(n):
        a = cd.GetArray(n)
        return vtk_to_numpy(a) if a is not None else None

    return {"x": pts[:, 0], "y": pts[:, 1], "avail": avail,
            "psi": arr("psi"),
            "kappaNoExt": arr("kappaNoExt"), "kappaDiv": arr("kappaDiv"),
            "kappaLap": arr("kappaLap")}


def cases(study):
    out = []
    for cp in sorted(glob.glob(os.path.join(study, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp))
        n = d.get("tokens", {}).get("N_CELLS")
        if n:
            out.append((int(float(n)), os.path.dirname(cp)))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--out", default="sl_curvature_pole.png")
    args = ap.parse_args()
    cs = [c for c in cases(args.study) if c[0] in (32, 64, 128, 256)]
    if not cs:
        print("no cases"); return 1

    fig, axes = plt.subplots(len(cs), 2, figsize=(13, 3.4 * len(cs)), squeeze=False)
    # Error metric is the PURE approximation error kappa_fit - kappa_exact, where
    # kappa_exact = 1/r is the curvature of the level set THROUGH that cell (a
    # circle of radius r about the centre) -- NOT 1/R. Comparing to 1/R would
    # conflate the genuine geometric 1/r variation across the band with the
    # discretization error. Pole/diagonal summary is on the NEAR-INTERFACE layer
    # (|psi| < 0.75 h, the cells whose kappa the CSF force actually interpolates).
    print(f"{'N':>5} {'field':>10} {'poleErr%':>9} {'diagErr%':>9} {'maxErr%':>9} {'argmax_theta':>12}")
    for r, (N, cdir) in enumerate(cs):
        f = read(cdir, ["kappaNoExt"])
        if f is None:
            axes[r][0].text(0.5, 0.5, "n/a", ha="center"); continue
        h = BOX / N
        rr = np.hypot(f["x"] - CX, f["y"] - CY)
        kap_exact = 1.0 / rr
        th = np.degrees(np.arctan2(f["y"] - CY, f["x"] - CX))
        near0 = np.abs(f["psi"]) < 0.75 * h                     # CSF-relevant layer
        for name, col in [("kappaNoExt", "#0072B2"), ("kappaDiv", "#D55E00"),
                          ("kappaLap", "#009E73")]:
            k = f.get(name)
            if k is None:
                continue
            band = (np.abs(k) > 1e-9) & near0
            relerr = (k[band] - kap_exact[band]) / kap_exact[band] * 100.0
            thb = th[band]
            o = np.argsort(thb)
            axes[r][0].plot(thb[o], relerr[o], ".", ms=3, color=col, alpha=0.7,
                            label=name)
            if name == "kappaNoExt":
                # pole bins (within 15 deg of 0/90/180/-90) vs diagonal bins (45..)
                def binmask(centres, tol=15):
                    m = np.zeros_like(thb, dtype=bool)
                    for c0 in centres:
                        d = np.abs(((thb - c0 + 180) % 360) - 180)
                        m |= d < tol
                    return m
                pole = binmask([0, 90, 180, -90])
                diag = binmask([45, 135, -45, -135])
                ap_ = np.mean(np.abs(relerr[pole])) if pole.any() else np.nan
                ad_ = np.mean(np.abs(relerr[diag])) if diag.any() else np.nan
                im = np.argmax(np.abs(relerr)) if relerr.size else 0
                print(f"{N:>5} {'noExt':>10} {ap_:>9.1f} {ad_:>9.1f} "
                      f"{(np.abs(relerr).max() if relerr.size else 0):>9.1f} "
                      f"{(thb[im] if relerr.size else 0):>12.1f}")
        axes[r][0].axhline(0, color="k", lw=0.5)
        for a in (-90, 0, 90):
            axes[r][0].axvline(a, color="0.8", lw=0.6)
        axes[r][0].set_title(rf"$N={N}$: approximation error at the interface layer")
        axes[r][0].set_xlabel(r"$\theta$ [deg] (0,$\pm$90,180 = grid-aligned poles)")
        axes[r][0].set_ylabel(r"$(\kappa-1/r)/(1/r)$ [%]")
        axes[r][0].set_xlim(-180, 180)
        axes[r][0].legend(fontsize=8, frameon=False, ncol=3)

        # spatial map of the no-extension approximation error (band vs local 1/r)
        ax = axes[r][1]
        k = f["kappaNoExt"]; band = np.abs(k) > 1e-9
        relerr = (k[band] - kap_exact[band]) / kap_exact[band] * 100.0
        amax = np.percentile(np.abs(relerr), 98) if relerr.size else 50
        sc = ax.scatter(f["x"][band], f["y"][band], c=relerr, cmap="RdBu_r",
                        vmin=-amax, vmax=amax, s=(160.0 / N) ** 2 * 5, marker="s")
        if f.get("psi") is not None:
            from matplotlib.tri import Triangulation
            try:
                ax.tricontour(f["x"], f["y"], f["psi"], levels=[0], colors="k",
                              linewidths=0.8)
            except Exception:
                pass
        fig.colorbar(sc, ax=ax, shrink=0.8).set_label(r"$(\kappa-1/r)/(1/r)$ [%]")
        ax.set_title(rf"$N={N}$: no-extension error (band)")
        ax.set_xlim(CX - 1.5e-3, CX + 1.5e-3); ax.set_ylim(CY - 1.5e-3, CY + 1.5e-3)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("Clean-SDF reconstructed curvature error around the circle "
                 "(no normal extension)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(args.out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[pole] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
