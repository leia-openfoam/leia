#!/usr/bin/env python3
"""3D semi-Lagrangian convergence figures for the reveal deck / paper.

For each 3D case (shear, deformation) this overlays the **hexahedral** and
**polyhedral** (cfMesh) convergence of the production ``uncachedQuadraticWeightedLeastSquares``
reconstruction on the same axes, reading the curated ``<study>_errors.csv`` of
the hex + poly studies and writing one triptych per case into
``doc/slides/figures/``:

    sl_3Dshear_convergence.png
    sl_3Ddeformation_convergence.png

Each triptych: gradient / shape / volume error vs h (log-log), one line per
(mesh, CFL) -- hex solid, poly dashed -- an O(h^2) guide, and the least-squares
measured order in each legend entry. h = 1/N for hex, h = maxCellSize for poly
(both are the mesh characteristic length, so the two lines are directly
comparable). Resilient: plots whatever (h, CFL) rows exist, so a missing study
or a failed resolution just leaves a shorter line (or drops that mesh).

Usage (from repo root):
    python3 workflow/scripts/make_sl_3d_fig.py
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import paths  # thematic docs layout (single source of truth for output dirs)

REPO = paths.REPO
# 3D SL convergence figures belong to the semi-Lagrangian theme's data/figures.
FIGDIR = paths.figs_dir("semi-lagrangian-level-set")

SERIES = [
    ("gradientError", r"$\||\nabla\psi|-1\|_{L_2}$"),
    ("shapeError",    "geometric shape error"),
    ("volumeError",   "volume-conservation error"),
]

# One figure per case; each overlays the hex + poly study of that case.
# (2D is hexahedral only; the 3D cases overlay hex + polyhedral.)
CASES = [
    {"out": "sl_2Dvortex_convergence.png",
     "title": r"Semi-Lagrangian $\mathtt{quadraticWeightedLeastSquares}$ — 2D reversed vortex ($T=2$)",
     "meshes": [("hex", "uncachedConv2Dvortex", "-", "o")]},
    {"out": "sl_3Dshear_convergence.png",
     "title": r"Semi-Lagrangian $\mathtt{quadraticWeightedLeastSquares}$ — 3D shear ($[0,1]^2\times[0,2]$)",
     "meshes": [("hex",  "uncachedConv3Dshear",      "-",  "o"),
                ("poly", "uncachedConv3DshearPoly",  "--", "s")]},
    {"out": "sl_3Ddeformation_convergence.png",
     "title": r"Semi-Lagrangian $\mathtt{quadraticWeightedLeastSquares}$ — 3D deformation (unit cube)",
     "meshes": [("hex",  "uncachedConv3Ddeformation",     "-",  "o"),
                ("poly", "uncachedConv3DdeformationPoly", "--", "s")]},
]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _order(h, e):
    """Least-squares convergence order (slope of log e vs log h)."""
    if len(h) < 2:
        return None
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def _load(study):
    """(rows, sorted CFL list) for a study, or ([], []) if it has no errors.csv."""
    path = os.path.join(REPO, "studies", study, f"{study}_errors.csv")
    if not os.path.isfile(path):
        print(f"[sl3d] no {path}; skipping that mesh")
        return [], []
    rows = _rows(path)
    cfls = sorted({_f(r["cfl"]) for r in rows if _f(r.get("cfl")) is not None})
    return rows, cfls


def make(case):
    # Load both meshes up front; skip the whole figure only if neither exists.
    loaded = [(mesh, style, mk, *_load(study))
              for mesh, study, style, mk in case["meshes"]]
    if not any(rows for *_, rows, _cfls in loaded):
        print(f"[sl3d] no studies for {case['out']}; skip")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, (col, ttl) in zip(axes, SERIES):
        allh, alle = [], []
        for mesh, style, mk, rows, cfls in loaded:
            for cfl in cfls:
                pts = sorted(
                    (_f(r["h"]), _f(r[col])) for r in rows
                    if _f(r.get("cfl")) == cfl and _f(r.get(col)) not in (None, 0.0)
                )
                pts = [(h, e) for h, e in pts if h and e]
                if not pts:
                    continue
                h, e = zip(*pts)
                p = _order(h, e)
                lbl = f"{mesh} CFL {cfl:g}" + (f"  (order {p:.2f})" if p is not None else "")
                ax.loglog(h, e, linestyle=style, marker=mk, ms=5, label=lbl)
                allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0 * (hr / h0) ** 2, "k:", lw=0.9, alpha=0.6, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$ (= $1/N$ hex, maxCellSize poly)"); ax.set_title(ttl)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(case["title"], fontsize=14)
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, case["out"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=130); plt.close(fig)
    nrows = sum(len(rows) for *_, rows, _cfls in loaded)
    print(f"[sl3d] wrote {out}  ({nrows} rows across "
          f"{sum(1 for *_, rows, _cfls in loaded if rows)} mesh type(s))")
    return out


if __name__ == "__main__":
    any_ok = False
    for case in CASES:
        if make(case):
            any_ok = True
    sys.exit(0 if any_ok else 1)
