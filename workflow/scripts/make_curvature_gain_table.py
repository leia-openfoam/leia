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
    # The shipped production curvature (CURVATURE_EXTENSION cellCentreInverse), 2026-09-29.
    "cellCentreInverse":    "cell-centre inverse (production, K-aware)",
    "cellCentreInverseNoK": "cell-centre inverse without the Gaussian term (control)",
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
    # A study that varies the psi surface (the ellipsoid gate: signed distance and implicit)
    # gets one table per surface, suffix + "_<surface>"; the two are never mixed (2026-09-29).
    surfaces = set()
    for meta in glob.glob(os.path.join(study, "*", "case_params.json")):
        with open(meta) as fh:
            surfaces.add(str(json.load(fh).get("tokens", {}).get("PSI_SURFACE", "")))
    if len(surfaces) <= 1:
        return _curate(study, None)
    rc = 0
    for s in sorted(surfaces):
        rc = max(rc, _curate(study, s))
    return rc


def _curate(study, surface):
    # Artifact suffix per GATE, so the circle, sphere and varying-curvature
    # gates write side by side instead of overwriting one another.
    # "ellipsoid" first: faceCurvatureEllipsoid3D contains "3d" and faceCurvatureEllipsoidPsi2D
    # no "ellipse", so they used to overwrite the sphere's and the circle's tables (2026-09-29);
    # a polyhedral twin gets its own suffix instead of overwriting the hex study's.
    _base = os.path.basename(os.path.normpath(study)).lower()
    if "ellipsoid" in _base:
        suffix = "_ellipsoid3d" if "3d" in _base else "_ellipsoidPsi"
    elif "3d" in _base or "sphere" in _base:
        suffix = "_3d"
    elif "ellipse" in _base:
        suffix = "_ellipse"
    else:
        suffix = ""
    if "poly" in _base and suffix:
        suffix += "_poly"
    if surface is not None:
        suffix += "_" + surface

    rows = []
    for meta in sorted(glob.glob(os.path.join(study, "*", "case_params.json"))):
        cpath = os.path.join(os.path.dirname(meta), CSV_NAME)
        if not os.path.isfile(cpath) or os.path.getsize(cpath) == 0:
            continue
        with open(meta) as fh:
            _tokens = json.load(fh).get("tokens", {})
        if surface is not None and str(_tokens.get("PSI_SURFACE", "")) != surface:
            continue
        n = _f(_tokens.get("N_CELLS"))
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
                    # The part of the response the pressure projection cannot
                    # remove: sigma*kappa_f enters the momentum equation only
                    # through face-to-face DIFFERENCES, so the mean and smooth
                    # parts of d kappa_f are absorbed exactly and GAIN_DIMLESS
                    # over-counts them. Measured: these two do NOT separate the
                    # deliveries either (0.799/0.800/0.806 for arithmetic /
                    # per-face inverse / footEval, whose coupled blow-up times
                    # are 3x apart) -- what differs by 112x is the delivered
                    # field's own roughness on a clean psi, not its response.
                    "GAIN_DRIVER_ACROSS_DIMLESS":
                        _f(r.get("GAIN_DRIVER_ACROSS_DIMLESS")),
                    "GAIN_DRIVER_ALONG_DIMLESS":
                        _f(r.get("GAIN_DRIVER_ALONG_DIMLESS")),
                })
    if not rows:
        print(f"[curvgain] no {CSV_NAME} under {study}")
        return 1

    tables = paths.tables_dir(THEME)
    out_csv = os.path.join(tables, f"curvature_gain{suffix}.csv")
    cols = ["MODEL", "N_CELLS", "DELTA_X", "EPS", "N_SEEDS",
            "E_L2", "GAIN_L2", "GAIN_LINF", "GAIN_DIMLESS",
            "GAIN_DRIVER_ACROSS_DIMLESS", "GAIN_DRIVER_ALONG_DIMLESS"]
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
