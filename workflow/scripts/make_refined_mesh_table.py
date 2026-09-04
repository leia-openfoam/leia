#!/usr/bin/env python3
"""Curate the statically refined meshes: one row per arm, from the drivers' own records.

Reads, per case directory of the named studies: case_params.json (tokens, mesh kind,
commit), refinement.csv (one row per refinement pass: dilations, 2:1 growth, cell counts),
refinedBand.csv (the band check on the FINAL mesh: complete fine layers at the worst
point, interface cells outside the fine region, fractions per level, checkMesh shape and
quality counts) and, when the arm ran the constant-curvature well-balanced force, the
maximum over time of the mean and L2 parasitic velocity from the solver CSV. Nothing is
grepped out of a solver log.

The uniform-equivalent cell count is N_CELLS^3 for hex (the twin's count) and, for a
polyhedral arm, the twin uniform study's count read from its owner-file header.

Usage:  python3 workflow/scripts/make_refined_mesh_table.py [--root .] [--studies A,B,...]
            [--out docs/.../refined_mesh_stats.csv]
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr  # noqa: E402

DEFAULT_STUDIES = ("stationaryDroplet3DrefinedWB", "stationaryDroplet3Drefined",
                   "stationaryDroplet3DrefinedL2", "stationaryDroplet3DrefinedBall",
                   "polyDroplet3Drefined_r13p8")
POLY_TWINS = {"polyDroplet3Drefined_r13p8": "polyDroplet3D_r13p8"}
SOLVER_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"


def _f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def _one_row(path):
    if not os.path.isfile(path):
        return {}
    rows = list(csv.DictReader(open(path)))
    return rows[0] if rows else {}


def _rows(path):
    return list(csv.DictReader(open(path))) if os.path.isfile(path) else []


def _max_over_time(path, col):
    best = None
    for r in _rows(path):
        v = _f(r.get(col))
        if v is not None and (best is None or v > best):
            best = v
    return best


def _twin_cells(root, study):
    twin = POLY_TWINS.get(study)
    if not twin:
        return None
    for d in sorted(glob.glob(os.path.join(root, "studies", twin, "*_0*"))):
        try:
            return lr.n_cells(d)
        except (OSError, RuntimeError):
            continue
    return None


def curate(root, studies):
    out = []
    for study in studies:
        for d in sorted(glob.glob(os.path.join(root, "studies", study, "*_[0-9]*"))):
            cp = os.path.join(d, "case_params.json")
            if not os.path.isfile(cp):
                continue
            meta = json.load(open(cp))
            tok = meta.get("tokens", {})
            band = _one_row(os.path.join(d, "refinedBand.csv"))
            passes = _rows(os.path.join(d, "refinement.csv"))
            if not band or not passes:
                print(f"[refined_mesh_table] {d}: no refinement.csv/refinedBand.csv yet, skipped")
                continue
            n = _f(tok.get("N_CELLS"))
            L = _f(tok.get("DOMAIN_LENGTH"))
            R = _f(tok.get("DROPLET_RADIUS"), 1e-3)
            mesh = meta.get("mesh")
            n_cells = int(_f(band.get("nCells"), 0))
            n_uni = int(round(n ** 3)) if (mesh == "hexRefined" and n) else _twin_cells(root, study)
            wb = tok.get("SURFACE_TENSION_FORCE") == "constantCurvatureSurfaceTension"
            csvp = os.path.join(d, SOLVER_CSV)
            row = {
                "study": study, "arm": os.path.basename(d), "mesh": mesh,
                "source": tok.get("REFINE_SOURCE", ""), "levels": tok.get("REFINE_LEVELS", ""),
                "bandCells": tok.get("REFINE_BAND_CELLS", ""),
                "N_fine": int(n) if n else "", "R_over_h": (R * n / L) if (n and L) else "",
                "hFine": band.get("hFine", ""),
                "nCells": n_cells, "nCellsUniform": n_uni if n_uni else "",
                "cellRatio": (n_uni / n_cells) if (n_uni and n_cells) else "",
                "dilationsPerPass": " ".join(p.get("dilations", "") for p in passes),
                "growth2to1PerPass": " ".join(str(p.get("growth2to1", "")) for p in passes),
                "fineLayersPredictedPerPass": " ".join(
                    f"{_f(p.get('fineLayersPredicted'), float('nan')):.2f}" for p in passes),
                "fineLayersWorst": band.get("fineLayersWorst", ""),
                "nBand": band.get("nBand", ""),
                "nBandOutsideFine": band.get("nBandOutsideFine", ""),
                "fracFine": band.get("fracFine", ""), "fracPerLevel": band.get("fracPerLevel", ""),
                "nHexahedra": band.get("nHexahedra", ""), "nPolyhedra": band.get("nPolyhedra", ""),
                "maxNonOrtho": band.get("maxNonOrtho", ""), "maxSkewness": band.get("maxSkewness", ""),
                "pinRelError": band.get("pinRelError", ""), "bandStatus": band.get("status", ""),
                "wellBalanced": "yes" if wb else "no",
                "wbMaxMeanMagUPrime": _max_over_time(csvp, "meanMagUPrime") if wb else "",
                "wbMaxL2MagUPrime": _max_over_time(csvp, "l2MagUPrime") if wb else "",
                "gitCommit": meta.get("gitCommit", ""),
            }
            out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--studies", default=",".join(DEFAULT_STUDIES))
    ap.add_argument("--out", default="docs/semi-lagrangian-level-set/sl-level-set-article/"
                                     "data/tables/refined_mesh_stats.csv")
    a = ap.parse_args()
    rows = curate(a.root, [s for s in a.studies.split(",") if s])
    if not rows:
        raise SystemExit("[refined_mesh_table] nothing to curate")
    os.makedirs(os.path.dirname(os.path.join(a.root, a.out)), exist_ok=True)
    lr.write_csv(os.path.join(a.root, a.out), rows)
    hdr = ("study", "arm", "source", "levels", "N_fine", "nCells", "nCellsUniform", "cellRatio",
           "dilationsPerPass", "fineLayersWorst", "nBandOutsideFine", "nPolyhedra", "maxNonOrtho",
           "bandStatus", "wbMaxL2MagUPrime")
    print("  ".join(hdr))
    for r in rows:
        print("  ".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in hdr))
    print(f"[refined_mesh_table] wrote {a.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
