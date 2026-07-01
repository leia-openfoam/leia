#!/usr/bin/env python3
"""Post-processing figures for a leia study, using foamlib to read the cases.

Produces, into ``<study>/figures/`` and (copied) into the reveal slide folder:
  * ``convergence_<model>.png`` -- one triptych per velocity-extension model
    (gradient / shape / volume error vs h, one line per T, O(h)/O(h^2) guides),
    with the velocity-extension model in the figure title;
  * ``interface_grid.png`` -- reversed-interface montage: initial (dashed) vs
    final (solid) alpha=0.5 contour for each (T, h) case  [2D structured cases,
    a single representative velocity-extension model].

Driven by the curated ``*_errors.csv`` (columns velocityExtension,T,h,
gradientError,shapeError,volumeError) written by aggregate.py, plus the per-case
``case_params.json`` for the montage geometry.
"""
import csv
import glob
import json
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from foamlib import FoamCase

# (column, subplot title) for the three convergence measures.
SERIES = [
    ("gradientError", r"$\||\nabla\psi|-1\|_{L_2}$"),
    ("shapeError",    "geometric shape error"),
    ("volumeError",   "volume-conservation error"),
]


def _rows(errors_csv):
    with open(errors_csv, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def convergence_triptych(rows, model, figdir):
    """One figure per velocity-extension model: gradient/shape/volume vs h,
    one line per period T, with O(h)/O(h^2) guides. Model name in the title."""
    sub = [r for r in rows if r.get("velocityExtension", "none") == model]
    if not sub:
        return None
    Ts = sorted({_f(r["T"]) for r in sub if _f(r["T"]) is not None})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, (col, title) in zip(axes, SERIES):
        allh, alle = [], []
        for T in Ts:
            pts = sorted(
                (_f(r["h"]), _f(r[col])) for r in sub
                if _f(r["T"]) == T and _f(r.get(col)) not in (None, 0.0)
            )
            pts = [(h, e) for h, e in pts if h and e]
            if not pts:
                continue
            h, e = zip(*pts)
            ax.loglog(h, e, marker="o", ms=5, label=f"T = {T:g}")
            allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.8, alpha=0.6, label=r"$O(h)$")
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.8, alpha=0.6, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$"); ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"Velocity extension: {model}", fontsize=14)
    p = os.path.join(figdir, f"convergence_{model}.png")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def _cases_2d(study_dir):
    """Discover materialized 2D cases from their case_params.json.
    Returns list of {dir, model, T, N}."""
    out = []
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        try:
            meta = json.load(open(meta_path))
        except Exception:
            continue
        if meta.get("dims") != 2:
            continue
        tok = meta.get("tokens", {})
        T, N = _f(tok.get("END_TIME")), _f(tok.get("N_CELLS"))
        if T is None or N is None:
            continue
        out.append({
            "dir": os.path.dirname(meta_path),
            "model": tok.get("VELOCITY_EXTENSION", "none"),
            "T": T, "N": int(N),
        })
    return out


def interface_montage(study_dir, figdir, model="none"):
    """Reversed-interface montage (initial dashed vs final solid alpha=0.5) over
    the (T, h) grid, for a single representative velocity-extension model.
    Domain geometry is inferred per case (square cells, x in [0,1])."""
    cases = [c for c in _cases_2d(study_dir) if c["model"] == model]
    if not cases:
        allc = _cases_2d(study_dir)
        if not allc:
            return None
        model = sorted({c["model"] for c in allc})[0]
        cases = [c for c in allc if c["model"] == model]
    grid = {(c["T"], c["N"]): c["dir"] for c in cases}
    Ts = sorted({c["T"] for c in cases})
    Ns = sorted({c["N"] for c in cases})
    fig, axes = plt.subplots(
        len(Ts), len(Ns),
        figsize=(3.0*len(Ns), 2.0*len(Ts)), squeeze=False,
    )
    any_ok = False
    for ir, T in enumerate(Ts):
        for ic, N in enumerate(Ns):
            ax = axes[ir][ic]
            ax.set_xticks([]); ax.set_yticks([])
            cdir = grid.get((T, N))
            if not cdir:
                ax.axis("off"); continue
            try:
                case = FoamCase(cdir)
                a0 = np.asarray(case[0]["alpha"].internal_field, float)
                aT = np.asarray(case[-1]["alpha"].internal_field, float)
                Nx = N; Ny = a0.size // Nx           # square cells, x in [0,1]
                ymax = Ny / Nx
                a0 = a0.reshape(Ny, Nx); aT = aT.reshape(Ny, Nx)
                xc = (np.arange(Nx)+0.5)/Nx
                yc = (np.arange(Ny)+0.5)/Ny*ymax
                X, Y = np.meshgrid(xc, yc)
                ax.contour(X, Y, aT, [0.5], colors="#D85A30", linewidths=2.2)
                ax.contour(X, Y, a0, [0.5], colors="k", linestyles="--", linewidths=1.0)
                ax.set_xlim(0, 1); ax.set_ylim(0, ymax); ax.set_aspect("equal")
                ax.set_title(f"T={T:g}, h={1.0/Nx:.3g}", fontsize=8)
                any_ok = True
            except Exception as exc:
                ax.text(0.5, 0.5, "n/a", ha="center", va="center", fontsize=8)
                ax.axis("off")
                print(f"[plots] {cdir}: {type(exc).__name__}: {exc}")
    if not any_ok:
        plt.close(fig); return None
    fig.suptitle(
        f"Reversed vortex ({model}): initial (dashed) vs final (solid) interface",
        fontsize=11)
    p = os.path.join(figdir, "interface_grid.png")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def indicator_comparison(rows, figdir):
    """Overlay convergence of the phase indicators (e.g. geometric vs
    detrixheAslam) on one triptych: gradient/shape/volume vs h, colour per T,
    line style + marker per indicator. Coinciding curves ⇒ the analytic
    (tolerance-free) indicator matches the geometric one. Also prints the max
    relative indicator disagreement per error (they integrate the same LLS plane,
    so this should be tiny)."""
    inds = sorted({r.get("phaseIndicator", "") for r in rows if r.get("phaseIndicator")})
    if len(inds) < 2:
        return None
    Ts = sorted({_f(r["T"]) for r in rows if _f(r["T"]) is not None})
    lstyles = ["-", "--", ":", "-."]
    markers = ["o", "x", "s", "^"]
    istyle = {ind: (lstyles[k % 4], markers[k % 4]) for k, ind in enumerate(inds)}
    cmap = plt.get_cmap("viridis")
    tcol = {T: cmap(i/max(1, len(Ts)-1)) for i, T in enumerate(Ts)}

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    for ax, (col, title) in zip(axes, SERIES):
        allh, alle = [], []
        for T in Ts:
            for ind in inds:
                pts = sorted(
                    (_f(r["h"]), _f(r[col])) for r in rows
                    if r.get("phaseIndicator") == ind and _f(r["T"]) == T
                    and _f(r.get(col)) not in (None, 0.0)
                )
                pts = [(h, e) for h, e in pts if h and e]
                if not pts:
                    continue
                h, e = zip(*pts)
                ls, mk = istyle[ind]
                ax.loglog(h, e, color=tcol[T], linestyle=ls, marker=mk, ms=5,
                          label=f"{ind}, T={T:g}")
                allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.7, alpha=0.5, label=r"$O(h)$")
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.7, alpha=0.5, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$"); ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle("Phase indicator: " + " vs ".join(inds)
                 + " (bulk shear vortex)", fontsize=14)
    p = os.path.join(figdir, "convergence_indicators.png")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(p, dpi=130); plt.close(fig)

    # Quantify indicator disagreement (same T,h; relative to the first indicator).
    ref = inds[0]
    for col, _ in SERIES:
        worst = 0.0
        for T in Ts:
            for r in rows:
                if r.get("phaseIndicator") != ref or _f(r["T"]) != T:
                    continue
                a = _f(r.get(col))
                b = next((_f(o.get(col)) for o in rows
                          if o.get("phaseIndicator") == inds[1]
                          and _f(o["T"]) == T and o["h"] == r["h"]), None)
                if a not in (None, 0.0) and b is not None:
                    worst = max(worst, abs(b-a)/abs(a))
        print(f"[plots] indicator max rel. diff ({inds[1]} vs {ref}) {col}: {worst:.3e}")
    return p


def main(study_dir, slides_dir):
    hits = glob.glob(os.path.join(study_dir, "*_errors.csv"))
    if not hits:
        print(f"[plots] no *_errors.csv in {study_dir}; nothing to plot")
        return []
    rows = _rows(hits[0])
    figdir = os.path.join(study_dir, "figures")
    os.makedirs(figdir, exist_ok=True)
    if slides_dir:
        os.makedirs(slides_dir, exist_ok=True)

    inds = sorted({r.get("phaseIndicator", "") for r in rows if r.get("phaseIndicator")})
    if len(inds) > 1:
        # Study varies the phase indicator -> overlay comparison.
        outputs = [p for p in [indicator_comparison(rows, figdir)] if p]
    else:
        # Study varies the velocity-extension model -> one triptych per model.
        models = sorted({r.get("velocityExtension", "none") for r in rows})
        outputs = [p for p in (convergence_triptych(rows, m, figdir) for m in models) if p]

    montage = interface_montage(study_dir, figdir, "none")
    if montage:
        outputs.append(montage)

    if slides_dir:
        for p in outputs:
            shutil.copy(p, os.path.join(slides_dir, os.path.basename(p)))
    dest = slides_dir if slides_dir else "(study only)"
    print(f"[plots] wrote {len(outputs)} figures to {figdir} and {dest}")
    return outputs


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else sys.argv[1])
