#!/usr/bin/env python3
"""Consolidate the semi-Lagrangian convergence studies into publication CSVs.

Reads the curated ``<study>_errors.csv`` of the uniform-refinement SL studies
(2D vortex, 3D shear/deformation, hex + polyhedral) and emits two tidy CSVs the
paper and slides consume directly:

    doc/paper/sl_convergence.csv
        one row per resolution:
        case, mesh, reconstruction, cfl, h, maxCellSize,
        gradientError, shapeError, volumeError
    doc/paper/sl_convergence_orders.csv
        one row per (case, mesh, reconstruction, cfl) with the least-squares
        convergence order (slope of log error vs log h) of each error metric:
        case, mesh, reconstruction, cfl, nLevels, hMin, hMax,
        gradientOrder, shapeOrder, volumeOrder

Least-squares order only (the uniform refinement ladder makes the slope clean);
no per-triplet Richardson / GCI. Missing studies are skipped with a note.

Usage (from repo root):
    python3 workflow/scripts/make_convergence_table.py
"""
import csv
import os
import sys

import numpy as np

import paths  # thematic docs layout (single source of truth for output dirs)

REPO = paths.REPO
# Curated convergence CSVs belong to the semi-Lagrangian theme's data/tables.
OUTDIR = paths.tables_dir("semi-lagrangian-level-set")

# study dir -> (case label, mesh label) for the publication tables.
STUDIES = [
    ("uncachedConv2Dvortex",         "2Dvortex",      "hex"),
    ("uncachedConv3Dshear",          "3Dshear",       "hex"),
    ("uncachedConv3Ddeformation",    "3Ddeformation", "hex"),
    ("uncachedConv3DshearPoly",      "3Dshear",       "poly"),
    ("uncachedConv3DdeformationPoly", "3Ddeformation", "poly"),
]

ERRORS = ["gradientError", "shapeError", "volumeError"]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _order(pts):
    """Least-squares order (slope of log e vs log h) over (h, e) pairs."""
    pts = [(h, e) for h, e in pts if h and e and h > 0 and e > 0]
    if len(pts) < 2:
        return None
    h, e = zip(*pts)
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def main():
    conv_rows, order_rows = [], []
    for study, case, mesh in STUDIES:
        path = os.path.join(REPO, "studies", study, f"{study}_errors.csv")
        if not os.path.isfile(path):
            print(f"[convtable] no {path}; skip {case}/{mesh}")
            continue
        rows = _rows(path)
        # tidy per-resolution rows
        for r in rows:
            conv_rows.append({
                "case": case, "mesh": mesh,
                "reconstruction": r.get("reconstruction", ""),
                "cfl": r.get("cfl", ""),
                "h": r.get("h", ""),
                "maxCellSize": r.get("maxCellSize", ""),
                "gradientError": r.get("gradientError", ""),
                "shapeError": r.get("shapeError", ""),
                "volumeError": r.get("volumeError", ""),
            })
        # per (reconstruction, cfl) group -> LSQ order of each error metric
        groups = {}
        for r in rows:
            key = (r.get("reconstruction", ""), r.get("cfl", ""))
            groups.setdefault(key, []).append(r)
        for (rec, cfl), grp in sorted(groups.items()):
            hs = [_f(r.get("h")) for r in grp]
            hs = [h for h in hs if h]
            row = {"case": case, "mesh": mesh, "reconstruction": rec, "cfl": cfl,
                   "nLevels": len(grp),
                   "hMin": f"{min(hs):.6g}" if hs else "",
                   "hMax": f"{max(hs):.6g}" if hs else ""}
            for col, name in zip(ERRORS, ("gradientOrder", "shapeOrder", "volumeOrder")):
                p = _order([(_f(r.get("h")), _f(r.get(col))) for r in grp])
                row[name] = f"{p:.3f}" if p is not None else ""
            order_rows.append(row)
        print(f"[convtable] {study}: {len(rows)} rows, {len(groups)} group(s)")

    if not conv_rows:
        print("[convtable] no studies found; nothing written")
        return 1

    os.makedirs(OUTDIR, exist_ok=True)
    conv_path = os.path.join(OUTDIR, "sl_convergence.csv")
    with open(conv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(conv_rows[0].keys()))
        w.writeheader(); w.writerows(conv_rows)
    order_path = os.path.join(OUTDIR, "sl_convergence_orders.csv")
    with open(order_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(order_rows[0].keys()))
        w.writeheader(); w.writerows(order_rows)
    print(f"[convtable] wrote {conv_path} ({len(conv_rows)} rows)")
    print(f"[convtable] wrote {order_path} ({len(order_rows)} rows)")

    # Also emit a booktabs LaTeX tabular so the article table stays a single source
    # (regenerated here; \input'd by the .tex). Body only -- the article supplies the
    # surrounding table float + caption.
    pretty = {"2Dvortex": "2D reversed vortex", "3Dshear": "3D shear",
              "3Ddeformation": "3D deformation"}
    tex_path = os.path.join(OUTDIR, "convergence_orders.tex")
    with open(tex_path, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_convergence_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{llccccc}\n\\toprule\n")
        fh.write("Case & Mesh & CFL & Levels & shape & volume & "
                 "$\\bigl\\|\\,|\\nabla\\psi|-1\\bigr\\|$ \\\\\n\\midrule\n")
        for r in order_rows:
            fh.write(" & ".join([
                pretty.get(r["case"], r["case"]), r["mesh"], str(r["cfl"]),
                str(r["nLevels"]),
                r["shapeOrder"] or "--", r["volumeOrder"] or "--",
                r["gradientOrder"] or "--"]) + " \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[convtable] wrote {tex_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
