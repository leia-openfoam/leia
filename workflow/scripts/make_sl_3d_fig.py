#!/usr/bin/env python3
"""3D semi-Lagrangian convergence figures for the reveal deck.

Reads the curated ``<study>_errors.csv`` of the two 3D SL studies
(config/3DshearSL.yaml, config/3DdeformationSL.yaml) and writes one convergence
triptych per case into ``doc/slides/figures/``:

    sl_3Dshear_convergence.png
    sl_3Ddeformation_convergence.png

Each triptych: gradient / shape / volume error vs h (log-log), ONE line per CFL,
an O(h^2) guide, and the least-squares measured order in each legend entry. Only
the production ``quadraticWeightedLeastSquares`` reconstruction is present in these studies.
Resilient: plots whatever (h, CFL) rows exist, so a missing/failed resolution
(e.g. an N=128 that ran out of memory) just leaves a shorter line + a note.

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

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIGDIR = os.path.join(REPO, "doc", "slides", "figures")

SERIES = [
    ("gradientError", r"$\||\nabla\psi|-1\|_{L_2}$"),
    ("shapeError",    "geometric shape error"),
    ("volumeError",   "volume-conservation error"),
]

CASES = [
    ("3DshearSL",       "sl_3Dshear_convergence.png",
     r"Semi-Lagrangian $\mathtt{quadraticWeightedLeastSquares}$ — 3D shear ($[0,1]^2\times[0,2]$)"),
    ("3DdeformationSL", "sl_3Ddeformation_convergence.png",
     r"Semi-Lagrangian $\mathtt{quadraticWeightedLeastSquares}$ — 3D deformation (unit cube)"),
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


def make(study, out_name, title):
    errors = os.path.join(REPO, "studies", study, f"{study}_errors.csv")
    if not os.path.isfile(errors):
        print(f"[sl3d] no {errors}; skip {out_name}")
        return None
    rows = _rows(errors)
    cfls = sorted({_f(r["cfl"]) for r in rows if _f(r.get("cfl")) is not None})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, (col, ttl) in zip(axes, SERIES):
        allh, alle = [], []
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
            lbl = f"CFL {cfl:g}" + (f"  (order {p:.2f})" if p is not None else "")
            ax.loglog(h, e, marker="o", ms=5, label=lbl)
            allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0 * (hr / h0) ** 2, "k:", lw=0.9, alpha=0.6, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h = 1/N$"); ax.set_title(ttl)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(title, fontsize=14)
    os.makedirs(FIGDIR, exist_ok=True)
    out = os.path.join(FIGDIR, out_name)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=130); plt.close(fig)
    ns = sorted({_f(r["h"]) for r in rows if _f(r.get("h"))}, reverse=True)
    print(f"[sl3d] wrote {out}  ({len(rows)} rows, CFL={cfls}, "
          f"resolutions h={['%.4g'%x for x in ns]})")
    return out


if __name__ == "__main__":
    any_ok = False
    for study, out_name, title in CASES:
        if make(study, out_name, title):
            any_ok = True
    sys.exit(0 if any_ok else 1)
