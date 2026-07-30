#!/usr/bin/env python3
"""Mean-curvature accuracy table + log-log figure (curvatureDroplet2D study).

Reads each case of the curvature study (leiaTestMeanCurvature.csv = band-restricted
error norms of the reconstructed curvature vs the exact 1/R, plus case_params.json
for N_CELLS) and emits into the semi-Lagrangian theme's data source:

    data/tables/curvature_error.csv    one row per resolution
    data/tables/curvature_error.tex    booktabs body (\\input by the article)
    data/figures/sl_curvature_error.png  log-log L1/L2 error vs h (div n and tr H)

The error is measured in the snGrad(alpha) band (the cells that carry the surface-
tension force). Two estimators are reported: div = div(grad psi/|grad psi|) (the
method the solver uses) and lap = tr(H) = Laplacian(psi) (the signed-distance
simplification). The convergence order is the least-squares slope of log(L2) vs log h.

Usage (from repo root):
    python3 workflow/scripts/make_curvature_table.py studies/curvatureDroplet2D
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

BOX_L = 0.01          # domain edge length [m] -> h = BOX_L / N
THEME = "semi-lagrangian-level-set"
CSV_NAME = "leiaTestMeanCurvature.csv"


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _final_row(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows[-1] if rows else {}


def _fit(h, y):
    m = [(a, b) for a, b in zip(h, y) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None
    lh, ly = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, ly, 1)
    r = ly - np.polyval(p, lh)
    ss = float((r**2).sum()); st = float(((ly - ly.mean())**2).sum())
    return float(p[0]), (1.0 - ss/st if st > 0 else 1.0)


def _rows(study_dir):
    recs = []
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        cdir = os.path.dirname(meta)
        cpath = os.path.join(cdir, CSV_NAME)
        if not os.path.isfile(cpath) or os.path.getsize(cpath) == 0:
            continue
        with open(meta) as fh:
            n = _f(json.load(fh).get("tokens", {}).get("N_CELLS"))
        row = _final_row(cpath)
        if not n or not row:
            continue
        ke = _f(row.get("KAPPA_EXACT")) or 1.0
        recs.append({
            "N": int(n), "h": BOX_L/n, "kappaExact": ke,
            "nBand": int(_f(row.get("N_BAND")) or 0),
            "L1div": _f(row.get("E_L1_DIV")), "L2div": _f(row.get("E_L2_DIV")),
            "Linfdiv": _f(row.get("E_LINF_DIV")),
            "L1lap": _f(row.get("E_L1_LAP")), "L2lap": _f(row.get("E_L2_LAP")),
            # div(n) evaluated at the cell centre, WITHOUT the normal extension:
            "L1noext": _f(row.get("E_L1_NOEXT")), "L2noext": _f(row.get("E_L2_NOEXT")),
            # face-average of the interpolated normal-extended kappa (force-effective):
            "L1faceavg": _f(row.get("E_L1_FACEAVG")), "L2faceavg": _f(row.get("E_L2_FACEAVG")),
            # foot-point height-function curvature:
            "L1foot": _f(row.get("E_L1_FOOT")), "L2foot": _f(row.get("E_L2_FOOT")),
            # iso-agnostic geometric curvature: conormal integral + local parabola fit:
            "L1iso": _f(row.get("E_L1_ISO")), "L2iso": _f(row.get("E_L2_ISO")),
            "L1parab": _f(row.get("E_L1_PARAB")), "L2parab": _f(row.get("E_L2_PARAB")),
        })
    recs.sort(key=lambda r: r["N"])
    return recs


def main(argv):
    if not argv:
        print("usage: make_curvature_table.py <study_dir>"); return 1
    recs = _rows(argv[0])
    if not recs:
        print(f"[curv] no completed cases under {argv[0]}"); return 1

    h = [r["h"] for r in recs]
    o2d, r2d = _fit(h, [r["L2div"] for r in recs])       # with normal extension
    o1d, _ = _fit(h, [r["L1div"] for r in recs])
    o2l, _ = _fit(h, [r["L2lap"] for r in recs])         # tr(H)
    o2n, r2n = _fit(h, [r["L2noext"] for r in recs])     # without normal extension
    o2f, _ = _fit(h, [r["L2faceavg"] for r in recs])     # face-averaged extended
    o2ft, _ = _fit(h, [r["L2foot"] for r in recs])       # foot-point height function
    o2i, r2i = _fit(h, [r["L2iso"] for r in recs])       # iso-agnostic conormal integral
    o2pb, r2pb = _fit(h, [r["L2parab"] for r in recs])   # iso-agnostic local parabola
    ke = recs[0]["kappaExact"]

    figs, tables = paths.figs_dir(THEME), paths.tables_dir(THEME)

    # --- CSV -----------------------------------------------------------------
    cpath = os.path.join(tables, "curvature_error.csv")
    cols = ["N", "h", "nBand", "kappaExact", "L1div", "L2div", "Linfdiv", "L1lap", "L2lap"]
    with open(cpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows({c: r.get(c, "") for c in cols} for r in recs)
    print(f"[curv] wrote {cpath} ({len(recs)} rows)")

    # --- LaTeX booktabs body -------------------------------------------------
    tpath = os.path.join(tables, "curvature_error.tex")
    with open(tpath, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_curvature_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{rccccc}\n\\toprule\n")
        fh.write("$N$ & $h$ [m] & band cells & $\\|\\kappa-\\kappa_e\\|_{1}$ & "
                 "$\\|\\kappa-\\kappa_e\\|_{2}$ & rel.\\ $L_2$ \\\\\n\\midrule\n")
        for r in recs:
            rel = r["L2div"]/r["kappaExact"] if r["L2div"] is not None else 0
            fh.write(f"{r['N']} & {r['h']:.3e} & {r['nBand']} & {r['L1div']:.3e} & "
                     f"{r['L2div']:.3e} & {100*rel:.1f}\\% \\\\\n")
        fh.write("\\midrule\n")
        order_txt = (f"$p = {o2d:.2f}$ ($R^2={r2d:.3f}$)" if (o2d is not None and len(recs) >= 3)
                     else (f"$p\\approx{o2d:.1f}$" if o2d is not None else "--"))
        fh.write("\\multicolumn{6}{l}{\\footnotesize Least-squares order of the $L_2$ curvature "
                 "error in $h$: " + order_txt + f"; $\\kappa_e=\\sigma$-exact $={ke:.0f}$ m$^{{-1}}$.}} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[curv] wrote {tpath}")

    # --- figure: L2 curvature error, WITH vs WITHOUT the normal extension, plus
    #     the face-averaged (force-effective) curvature. -------------------------
    def _lab(name, o):
        return f"{name}" + (rf" ($\propto h^{{{o:.2f}}}$)" if o is not None else "")
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.loglog(h, [r["L2div"] for r in recs], "o-", color="#0072B2",
              label=_lab("with normal extension", o2d))
    ax.loglog(h, [r["L2noext"] for r in recs], "s--", color="#D55E00",
              label=_lab("without extension (cell centre)", o2n))
    ax.loglog(h, [r["L2faceavg"] for r in recs], "D-.", color="#009E73",
              label=_lab("face-averaged (force-effective)", o2f))
    ax.loglog(h, [r["L2lap"] for r in recs], "^:", color="gray",
              label=r"$\mathrm{tr}(H)=\Delta\psi$")
    if any(r.get("L2iso") for r in recs):
        ax.loglog(h, [r["L2iso"] for r in recs], "v-", color="#CC79A7",
                  label=_lab("iso-agnostic conormal integral", o2i))
    if any(r.get("L2parab") for r in recs):
        ax.loglog(h, [r["L2parab"] for r in recs], "P--", color="#E69F00",
                  label=_lab("iso-agnostic local parabola", o2pb))
    ax.set_xlabel(r"$h$ [m]")
    ax.set_ylabel(r"band $L_2$ curvature error [m$^{-1}$]")
    ax.set_title("Mean-curvature error in the snGrad$(\\alpha)$ band")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    ppath = os.path.join(figs, "sl_curvature_error.png")
    fig.savefig(ppath, dpi=200); plt.close(fig)
    print(f"[curv] wrote {ppath}")
    print(f"[curv] L2 order: with-ext p={o2d}; no-ext p={o2n}; face-avg p={o2f}; "
          f"tr(H) p={o2l}; foot p={o2ft}; iso-conormal p={o2i} (R2={r2i}); "
          f"iso-parabola p={o2pb} (R2={r2pb})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
