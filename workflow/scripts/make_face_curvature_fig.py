#!/usr/bin/env python3
"""Face-centered curvature convergence figure + order table (faceCurvatureDroplet2D).

Reads each case's leiaTestFaceCurvature.csv (tidy: one row per MODEL x FOOT_POINT,
|Sf|-weighted error norms of the CSF-applied face curvature kappa_f vs the exact
1/R over the ACTIVE |snGrad(alpha)| > 0 faces) plus case_params.json for N_CELLS,
and emits into the method-comparison theme's data source:

    data/figures/face_curvature_convergence.png   log-log active-face L2 vs h
    data/tables/face_curvature_orders.csv          every model, fitted orders
    data/tables/face_curvature_orders.tex          booktabs body (\\input-able)

The figure carries the seven decision-relevant deliveries; models that duplicate
a plotted line to within measurement noise (trHessian ~ quadraticCellCentre,
fvmTraceGradGradPsi ~ fvmDivGradPsi, quadraticNewtonFoot/quadraticClosestPoint ~
quadraticCellCentre, isoConormal/isoParabola diagnostics) are reported in the
order table only. Dashed = the raw CSF delivery, solid = the same model
re-referenced to the interface with the stabilized foot point; the constant
interface-mean diagnostic is drawn as a gray reference.

Usage (from repo root):
    python3 workflow/scripts/make_face_curvature_fig.py studies/faceCurvatureDroplet2D
"""
import csv
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

THEME = "method-comparison"
CSV_NAME = "leiaTestFaceCurvature.csv"

# Entity -> fixed hue (validated categorical order; color follows the model,
# solid/dashed carries the foot-point axis).
COLORS = {
    "quadraticCellCentre":     "#0072B2",
    "fvmDivGradPsi":           "#D55E00",
    "fvmDivGradAlpha":         "#009E73",
    "stableFootPoint":         "#E69F00",
    "footPointHeightFunction": "#56B4E9",
    "connectedInterface":      "#CC79A7",
    "kangQuadratic":           "#882255",
    "cutCellInverse":          "#117733",
    "cellMeanInverse":         "#44AA99",
}
MARKERS = {
    "quadraticCellCentre":     "o",
    "fvmDivGradPsi":           "s",
    "fvmDivGradAlpha":         "D",
    "stableFootPoint":         "^",
    "footPointHeightFunction": "v",
    "connectedInterface":      "P",
    "kangQuadratic":           "X",
    "cutCellInverse":          "*",
    "cellMeanInverse":         "h",
}
LABELS = {
    "quadraticCellCentre":     "quadratic cell centre",
    "trHessian":               "tr(H) Laplacian",
    "fvmDivGradPsi":           r"FVM div($\nabla\psi/|\nabla\psi|$)",
    "fvmTraceGradGradPsi":     r"FVM tr(grad $n_\psi$)",
    "fvmDivGradAlpha":         r"FVM $-$div($\nabla\alpha/|\nabla\alpha|$)",
    "kangQuadratic":           "Kang face interpolation",
    "quadraticNewtonFoot":     "Newton-foot extension",
    "quadraticClosestPoint":   "closest-point extension",
    "footPointHeightFunction": "foot-point height function",
    "connectedInterface":      "connected interface (face-native)",
    "interfaceMean":           "interface mean (constant)",
    "isoConormal":             "iso-geometric conormal",
    "isoParabola":             "iso-geometric parabola",
    "stableFootPoint":         "stabilized foot point (per-side variant)",
    "scalarInverse2D":         "2D scalar inverse (control, no Gaussian term)",
    "cutCellInverse":          "one inverted value per cut cell",
    "cellMeanInverse":         "cut-cell mean of per-face inversions",
}
# (entity, foot_point) pairs drawn in the figure; the rest go to the table.
PLOT_PAIRED = ["quadraticCellCentre", "fvmDivGradPsi", "fvmDivGradAlpha"]
PLOT_SINGLE = [("stableFootPoint", 1), ("footPointHeightFunction", 0),
               ("connectedInterface", 0), ("kangQuadratic", 0),
               ("cutCellInverse", 1), ("cellMeanInverse", 1)]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _fit(h, y):
    m = [(a, b) for a, b in zip(h, y) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None
    lh, ly = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, ly, 1)
    r = ly - np.polyval(p, lh)
    ss = float((r**2).sum()); st = float(((ly - ly.mean())**2).sum())
    return float(p[0]), (1.0 - ss/st if st > 0 else 1.0)


def _series(study_dir):
    """{(model, fp): {"N": [...], "h": [...], "L1": [...], "L2": [...],
    "Linf": [...], "L2w": [...]}} sorted by N; plus kappaExact."""
    recs = {}
    kappa_exact = None
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        cdir = os.path.dirname(meta)
        cpath = os.path.join(cdir, CSV_NAME)
        if not os.path.isfile(cpath) or os.path.getsize(cpath) == 0:
            continue
        with open(meta) as fh:
            n = _f(json.load(fh).get("tokens", {}).get("N_CELLS"))
        if not n:
            continue
        with open(cpath, newline="") as fh:
            for row in csv.DictReader(fh):
                # A truncated row from an interrupted run must not kill the
                # whole report: skip anything without the load-bearing fields.
                model, fp = row.get("MODEL"), _f(row.get("FOOT_POINT"))
                h = _f(row.get("DELTA_X"))
                if not model or fp is None or not h:
                    continue
                key = (model, int(fp))
                s = recs.setdefault(key, {"N": [], "h": [], "L1": [], "L2": [],
                                          "Linf": [], "L2w": []})
                s["N"].append(int(n))
                s["h"].append(h)
                s["L1"].append(_f(row.get("E_L1")))
                s["L2"].append(_f(row.get("E_L2")))
                s["Linf"].append(_f(row.get("E_LINF")))
                s["L2w"].append(_f(row.get("E_L2_FORCEW")))
                kappa_exact = _f(row.get("KAPPA_EXACT")) or kappa_exact
    for s in recs.values():
        order = np.argsort(s["N"])
        for k in s:
            s[k] = [s[k][i] for i in order]
    return recs, kappa_exact


def main(argv):
    if not argv:
        print("usage: make_face_curvature_fig.py <study_dir>"); return 1
    recs, kappa_exact = _series(argv[0])
    if not recs:
        print(f"[facecurv] no completed cases under {argv[0]}"); return 1

    # 3D sphere gate -> its own artifact names (auto from the study name).
    # Artifact suffix per GATE, so the circle, sphere and varying-curvature
    # gates write side by side instead of overwriting one another.
    _base = os.path.basename(os.path.normpath(argv[0])).lower()
    if "3d" in _base or "sphere" in _base:
        suffix = "_3d"
    elif "ellipse" in _base:
        suffix = "_ellipse"
    else:
        suffix = ""

    figs, tables = paths.figs_dir(THEME), paths.tables_dir(THEME)
    orders = {k: _fit(s["h"], s["L2"]) for k, s in recs.items()}

    # --- figure --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.8, 5.4))

    def _plot(model, fp, ls, label_suffix=""):
        s = recs.get((model, fp))
        if not s:
            return
        p, _ = orders[(model, fp)]
        lab = LABELS[model] + label_suffix
        if p is not None:
            lab += rf"  ($\propto h^{{{p:.2f}}}$)"
        ax.loglog(s["h"], s["L2"], ls, color=COLORS[model],
                  marker=MARKERS[model], ms=5, lw=1.6, label=lab)

    for model in PLOT_PAIRED:
        _plot(model, 0, "--")
        _plot(model, 1, "-", " + stabilized foot point")
    for model, fp in PLOT_SINGLE:
        _plot(model, fp, "-")

    # The K-less control: identical to the corrected line in 2D (skip the
    # duplicate); in 3D it exposes the first-order failure of the 2D formula.
    ctrl, corr = recs.get(("scalarInverse2D", 1)), recs.get(("quadraticCellCentre", 1))
    if ctrl and corr and any(
            abs(a - b) > 0.01*max(abs(b), 1e-30)
            for a, b in zip(ctrl["L2"], corr["L2"])):
        p, _ = orders[("scalarInverse2D", 1)]
        lab = LABELS["scalarInverse2D"]
        if p is not None:
            lab += rf"  ($\propto h^{{{p:.2f}}}$)"
        ax.loglog(ctrl["h"], ctrl["L2"], "--", color="#7f7f7f", marker="x",
                  ms=6, lw=1.3, label=lab)

    # Constant interface-mean diagnostic: a control, not a competing series.
    s = recs.get(("interfaceMean", 0))
    if s:
        p, _ = orders[("interfaceMean", 0)]
        lab = LABELS["interfaceMean"]
        if p is not None:
            lab += rf"  ($\propto h^{{{p:.2f}}}$)"
        ax.loglog(s["h"], s["L2"], ":", color="#7f7f7f", marker=".", ms=6,
                  lw=1.3, label=lab)

    # h^1 / h^2 order guides anchored to the plotted data range.
    hs = sorted({h for s in recs.values() for h in s["h"] if h})
    if hs and kappa_exact:
        h0, h1 = hs[0], hs[-1]
        y1 = 0.12*kappa_exact
        ax.loglog([h0, h1], [y1*h0/h1, y1], color="#7f7f7f", ls="-", lw=0.8,
                  alpha=0.7)
        ax.annotate(r"$h^1$", xy=(h0*1.05, 1.15*y1*h0/h1), color="#7f7f7f",
                    fontsize=9)
        y2 = 0.015*kappa_exact
        ax.loglog([h0, h1], [y2*(h0/h1)**2, y2], color="#7f7f7f", ls="-",
                  lw=0.8, alpha=0.7)
        ax.annotate(r"$h^2$", xy=(h0*1.05, 1.15*y2*(h0/h1)**2), color="#7f7f7f",
                    fontsize=9)

    ax.set_xlabel(r"$h$ [m]")
    ax.set_ylabel(r"active-face $L_2$ error of $\kappa_f$ [m$^{-1}$]")
    geom = "sphere" if suffix else "circle"
    ax.set_title("Face-centered curvature on the CSF force support\n"
                 rf"(static {geom}, exact SDF; faces with $|\mathrm{{snGrad}}\,\alpha| > 0$)",
                 fontsize=11)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(frameon=False, fontsize=8, loc="center left",
              bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    ppath = os.path.join(figs, f"face_curvature_convergence{suffix}.png")
    fig.savefig(ppath, dpi=200); plt.close(fig)
    print(f"[facecurv] wrote {ppath}")

    # --- order table: EVERY model x foot-point row ---------------------------
    def _key(k):
        model, fp = k
        return (list(LABELS).index(model) if model in LABELS else 99, fp)

    rows = []
    for k in sorted(recs, key=_key):
        model, fp = k
        s = recs[k]
        p2, r2 = orders[k]
        pinf, _ = _fit(s["h"], s["Linf"])
        rows.append({
            "model": model, "footPoint": fp,
            "L2_finest": s["L2"][-1], "Linf_finest": s["Linf"][-1],
            # the |snGrad(alpha)||Sf|-weighted norm: how the error enters the
            # integrated force flux G_sigma,f.
            "L2_forceWeighted_finest": s["L2w"][-1],
            "order_L2": p2, "R2": r2, "order_Linf": pinf,
            "N_finest": s["N"][-1],
        })

    cpath = os.path.join(tables, f"face_curvature_orders{suffix}.csv")
    with open(cpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"[facecurv] wrote {cpath} ({len(rows)} rows)")

    tpath = os.path.join(tables, f"face_curvature_orders{suffix}.tex")
    with open(tpath, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_face_curvature_fig.py"
                 " -- do not edit.\n")
        fh.write("\\begin{tabular}{llcccc}\n\\toprule\n")
        fh.write("face curvature & foot point & $L_2$ ($N{=}%d$) & $L_\\infty$"
                 " ($N{=}%d$) & $p(L_2)$ & $R^2$ \\\\\n\\midrule\n"
                 % (rows[0]["N_finest"], rows[0]["N_finest"]))
        for r in rows:
            p = f"{r['order_L2']:.2f}" if r["order_L2"] is not None else "--"
            q = f"{r['R2']:.3f}" if r["R2"] is not None else "--"
            fh.write(f"{LABELS.get(r['model'], r['model'])} & "
                     f"{'yes' if r['footPoint'] else 'no'} & "
                     f"{r['L2_finest']:.3e} & {r['Linf_finest']:.3e} & "
                     f"{p} & {q} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[facecurv] wrote {tpath}")

    for r in rows:
        print(f"[facecurv] {r['model']}{' +fp' if r['footPoint'] else '':>4}: "
              f"p(L2) = {r['order_L2']}, L2@{r['N_finest']} = {r['L2_finest']:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
