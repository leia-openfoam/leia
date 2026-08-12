#!/usr/bin/env python3
"""Curvature-delivery GAIN table (leiaTestCurvatureNoiseGain).

The face-curvature convergence study answers "how accurate is kappa_f on an
exact signed distance field". In the coupled stationary-droplet solve psi is
never exact, and MEASURED (docs/plan-curvature-stabilization.md sec. 10) that
accuracy does not predict coupled stability: over the three deliveries with
blow-up times at N=128 the two orders are INVERTED -- the most accurate
delivery blows up first.

The second number is the gain,

    G = || d kappa_f ||_2 / || d psi ||_inf      [1/m^2]

on the CSF force support, i.e. how far the delivered face curvature moves per
unit perturbation of psi. Reported here as G h^2, the amplification relative to
the h^-2 that any curvature operator must pay, which makes it comparable across
resolutions (it is constant to ~3% over the ladder) and across deliveries.

Reads studies/<gate>/*/leiaTestCurvatureNoiseGain.csv (+ case_params.json for
N_CELLS) and writes into the method-comparison theme's data source:

    data/tables/curvature_gain.csv    one row per (model, N, eps)
    data/tables/curvature_gain.tex    booktabs body, the finest eps-linear row
                                      per model: accuracy AND gain side by side

Usage (from repo root):
    python3 workflow/scripts/make_curvature_gain_table.py studies/faceCurvatureDroplet2D
"""
import csv
import glob
import json
import os
import sys

import paths

THEME = "method-comparison"
CSV_NAME = "leiaTestCurvatureNoiseGain.csv"

LABELS = {
    "arithmetic":     "arithmetic (interpolated cell curvature)",
    "perFaceInverse": "per-face parallel-surface inverse",
    "cutCellInverse": "one inverted value per cut cell",
    "cellMeanInverse": "cut-cell mean of per-face inversions",
}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main(argv):
    if not argv:
        print("usage: make_curvature_gain_table.py <study_dir>")
        return 1
    study = argv[0]
    # Artifact suffix per GATE, so the circle, sphere and varying-curvature
    # gates write side by side instead of overwriting one another.
    _base = os.path.basename(os.path.normpath(study)).lower()
    if "3d" in _base or "sphere" in _base:
        suffix = "_3d"
    elif "ellipse" in _base:
        suffix = "_ellipse"
    else:
        suffix = ""

    rows = []
    for meta in sorted(glob.glob(os.path.join(study, "*", "case_params.json"))):
        cpath = os.path.join(os.path.dirname(meta), CSV_NAME)
        if not os.path.isfile(cpath) or os.path.getsize(cpath) == 0:
            continue
        with open(meta) as fh:
            n = _f(json.load(fh).get("tokens", {}).get("N_CELLS"))
        if not n:
            continue
        with open(cpath, newline="") as fh:
            for r in csv.DictReader(fh):
                model, h = r.get("MODEL"), _f(r.get("DELTA_X"))
                gd = _f(r.get("GAIN_DIMLESS"))
                if not model or not h or gd is None:
                    continue
                rows.append({
                    "MODEL": model, "N_CELLS": int(n), "DELTA_X": h,
                    "EPS": _f(r.get("EPS")), "N_SEEDS": r.get("N_SEEDS"),
                    "E_L2": _f(r.get("E_L2")), "GAIN_L2": _f(r.get("GAIN_L2")),
                    "GAIN_LINF": _f(r.get("GAIN_LINF")), "GAIN_DIMLESS": gd,
                })
    if not rows:
        print(f"[curvgain] no {CSV_NAME} under {study}")
        return 1

    tables = paths.tables_dir(THEME)
    out_csv = os.path.join(tables, f"curvature_gain{suffix}.csv")
    cols = ["MODEL", "N_CELLS", "DELTA_X", "EPS", "N_SEEDS",
            "E_L2", "GAIN_L2", "GAIN_LINF", "GAIN_DIMLESS"]
    rows.sort(key=lambda r: (r["MODEL"], r["N_CELLS"], r["EPS"]))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    # The .tex body: the finest mesh, at the middle amplitude (the response is
    # linear in eps -- the CSV carries every amplitude so that can be checked).
    nMax = max(r["N_CELLS"] for r in rows)
    eps_vals = sorted({r["EPS"] for r in rows if r["EPS"]})
    eps_mid = eps_vals[len(eps_vals)//2] if eps_vals else None
    body = [r for r in rows if r["N_CELLS"] == nMax and r["EPS"] == eps_mid]
    body.sort(key=lambda r: r["GAIN_DIMLESS"])

    out_tex = os.path.join(tables, f"curvature_gain{suffix}.tex")
    with open(out_tex, "w") as fh:
        for r in body:
            fh.write(
                f"{LABELS.get(r['MODEL'], r['MODEL'])} & "
                f"{r['E_L2']:.3g} & {r['GAIN_L2']:.4g} & "
                f"{r['GAIN_DIMLESS']:.3f} \\\\\n"
            )

    print(f"[curvgain] {len(rows)} rows -> {out_csv}")
    for r in body:
        print(f"[curvgain] N={r['N_CELLS']} {r['MODEL']:<16}"
              f" E_L2 = {r['E_L2']:.4g} 1/m,"
              f" gain = {r['GAIN_DIMLESS']:.3f} x a second difference")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
