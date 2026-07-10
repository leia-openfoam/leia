#!/usr/bin/env python3
"""Comparison figure: quadraticWeightedLeastSquares vs defectCorrectedIDW departure-point
reconstruction on the three SL cases (2D reversed vortex, 3D shear, 3D LeVeque
deformation). Reads the curated <study>_errors.csv of the focused idwCompare*
studies and plots shape error + |grad psi|-1 error vs h (log-log), one line per
(reconstruction, CFL). Writes doc/slides/figures/sl_idw_vs_quadratic.png.

Shows the finding: quadraticWeightedLeastSquares (value least-squares fit) converges and keeps
|grad psi| bounded, while defectCorrectedIDW (gradient-based Taylor correction) is
2nd-order accurate for a STATIC field but unstable under iterated SL transport --
its |grad psi| blows up and the shape error does not converge.

    python3 workflow/scripts/make_idw_compare_fig.py
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIGDIR = os.path.join(REPO, "doc", "slides", "figures")

CASES = [
    ("idwCompare2D",                                                 "2D reversed vortex (T=2)"),
    ("idwCompare3Dshear",                                            "3D shear (T=3)"),
    (["idwCompare3Ddeformation", "idwCompare3Ddeformation128"],      "3D deformation (T=3)"),
]

# Historical short names -> spelled-out OpenFOAM class names, so CSVs written before the
# rename (and the memory-heavy 128^3 run launched with the old token) merge into one series.
_RECON_ALIAS = {
    "quadraticWLSQ": "quadraticWeightedLeastSquares",
    "bandQuadraticWLSQ": "bandQuadraticWeightedLeastSquares",
}
RECON = [
    ("quadraticWeightedLeastSquares",     "quadraticWeightedLeastSquares",     "o", "-"),
    ("defectCorrectedIDW","defectCorrectedIDW","s", "--"),
]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rows(study):
    studies = [study] if isinstance(study, str) else list(study)
    out = []
    for s in studies:
        p = os.path.join(REPO, "studies", s, s + "_errors.csv")
        if not os.path.isfile(p):
            continue
        for r in csv.DictReader(open(p)):
            rec = r.get("reconstruction")
            r["reconstruction"] = _RECON_ALIAS.get(rec, rec)
            out.append(r)
    return out


def series(rs, rec, cfl, col):
    pts = sorted(
        (_f(r["h"]), _f(r[col])) for r in rs
        if r.get("reconstruction") == rec and r.get("cfl") == cfl
        and _f(r.get(col)) not in (None,) and _f(r["h"])
    )
    pts = [(h, e) for h, e in pts if h and e is not None and e > 0]
    return zip(*pts) if pts else ([], [])


def main():
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    any_data = False
    for j, (study, title) in enumerate(CASES):
        rs = rows(study)
        if rs:
            any_data = True
        for irow, col, ylab in [(0, "shapeError", "geometric shape error"),
                                (1, "gradientError", r"$\||\nabla\psi|-1\|_{L_2}$")]:
            ax = axes[irow][j]
            for rec, lab, mk, ls in RECON:
                for cfl in ["0.5", "1.0"]:
                    h, e = series(rs, rec, cfl, col)
                    if not h:
                        continue
                    ax.loglog(h, e, marker=mk, ls=ls, ms=5,
                              label=f"{lab} CFL {cfl}")
            if irow == 0:
                ax.set_title(title)
            ax.set_xlabel(r"$h=1/N$"); ax.set_ylabel(ylab)
            ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle("Departure-point reconstruction: quadraticWeightedLeastSquares vs defectCorrectedIDW "
                 "(top: shape error; bottom: |grad psi|-1)", fontsize=13)
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, "sl_idw_vs_quadratic.png")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[idwcmp] wrote {out}")
    return 0 if any_data else 1


if __name__ == "__main__":
    sys.exit(main())
