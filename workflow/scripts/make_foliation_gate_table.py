#!/usr/bin/env python3
"""Curate the FOLIATION face-curvature gate (config/faceCurvatureEllipsoidPsi2D.yaml).

The study runs the SAME 2:1 ellipse interface with psi initialized two ways:
  signedDistanceEllipse : TRUE signed distance -- parallel foliation, D == 0,
                          so the delivered-curvature error is the FIT error.
  implicitEllipsoid     : the quadratic form -- beta = |grad psi| varies by
                          a/b = 2 ALONG the interface. psi is EXACTLY quadratic,
                          so the fit is exact to machine precision and the
                          error is the parallel-surface inverse's FOLIATION
                          BIAS (plan sec. 11.2, first order in the offset d)
                          alone.
Each arm therefore isolates one term of the delivered kappa_f error budget.

Reads the per-case leiaTestFaceCurvature.csv rows, fits the L2 order over the
last three ladder levels (N=64 is a trend point, per the case header), and
writes docs/method-comparison/method-comparison-article/data/tables/
face_curvature_orders_foliation.csv (+ .tex twin).

Usage: make_foliation_gate_table.py <study_dir> <out_dir>
"""

import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 1
    study, outdir = argv

    # (surface, model, footPoint) -> {N: (E_L2, E_LINF)}
    data = defaultdict(dict)
    for d in sorted(glob.glob(os.path.join(study, "ellipseDroplet2D_0*"))):
        fv = open(os.path.join(d, "system", "fvSolution")).read()
        surf = re.search(
            r"implicitSurface\s*\{[^}]*type\s+(\w+)", fv, re.S
        ).group(1)
        n = int(re.search(
            r"n_cells\s+(\d+)",
            open(os.path.join(d, "system", "blockMeshDict")).read()
        ).group(1))
        p = os.path.join(d, "leiaTestFaceCurvature.csv")
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p)):
            key = (surf, r["MODEL"], r["FOOT_POINT"])
            data[key][n] = (float(r["E_L2"]), float(r["E_LINF"]))

    rows = []
    for (surf, model, fp), byN in sorted(data.items()):
        ns = sorted(byN)
        if len(ns) < 3:
            continue
        # Order fit over the last three levels (skip the N=64 trend point).
        fit = ns[-3:]
        e = [byN[n][0] for n in fit]
        if min(e) <= 0:
            continue
        num = math.log(e[0]/e[-1])
        den = math.log(fit[-1]/fit[0])
        rows.append(dict(
            surface=surf, model=model, footPoint=fp,
            N_finest=ns[-1],
            L2_finest=byN[ns[-1]][0], Linf_finest=byN[ns[-1]][1],
            order_L2=num/den,
        ))

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "face_curvature_orders_foliation.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[foliation-gate] wrote {out} ({len(rows)} rows)")

    tex = os.path.join(outdir, "face_curvature_orders_foliation.tex")
    with open(tex, "w") as f:
        f.write("\\begin{tabular}{l l c r r r}\n\\hline\n"
                "$\\psi$ & model & foot & $E_{L_2}$ ($N{=}512$) & "
                "$E_{L_\\infty}$ & order \\\\\n\\hline\n")
        for r in rows:
            f.write(
                f"{r['surface']} & {r['model']} & {r['footPoint']} & "
                f"{r['L2_finest']:.3g} & {r['Linf_finest']:.3g} & "
                f"{r['order_L2']:.2f} \\\\\n"
            )
        f.write("\\hline\n\\end{tabular}\n")
    print(f"[foliation-gate] wrote {tex}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
