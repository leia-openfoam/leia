#!/usr/bin/env python3
"""Semi-Lagrangian convergence figures for the articles + decks.

The two semi-Lagrangian METHOD LINES are handled as independent methods, each
with its own study set, docs theme and figure prefix:

    --method quadratic   uncachedConv* studies  -> semi-lagrangian-level-set (sl_*)
    --method linear      linearConv*   studies  -> linear-semi-lagrangian-level-set (lsl_*)

For each case (2D vortex, 3D shear, 3D deformation) this overlays the
**hexahedral** and **polyhedral** (cfMesh) convergence of the method's
reconstruction on the same axes, reading the curated ``<study>_errors.csv`` of
the hex + poly studies and writing per case into the theme's ``data/figures/``
(see ``paths.figs_dir``):

    <prefix>_<case>_convergence.png
        triptych of the PRIMARY errors: gradient / shape / volume error vs h
        (log-log), one line per (mesh, CFL) -- hex solid, poly dashed -- guide
        slopes (O(h^2) for quadratic; O(h) AND O(h^2) for linear) and the
        least-squares measured order in each legend entry.
    <prefix>_<case>_convergence_diagnostics.png
        the remaining error measures, completing the ALL-ERRORS picture:
        band-restricted gradient error at t=T, gradient + band-gradient error
        at the half time t=T/2 (maximum deformation), and the half-time volume
        error.

h = 1/N for hex, h = maxCellSize for poly (both are the mesh characteristic
length, so the two lines are directly comparable). Resilient: plots whatever
(h, CFL) rows exist, so a missing study or a failed resolution just leaves a
shorter line (or drops that mesh).

Usage (from repo root):
    python3 workflow/scripts/make_sl_3d_fig.py [--method quadratic|linear]
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import paths  # thematic docs layout (single source of truth for output dirs)

REPO = paths.REPO

# Primary triptych: the headline error measures at t=T.
SERIES = [
    ("gradientError", r"$\||\nabla\psi|-1\|_{L_2}$"),
    ("shapeError",    "geometric shape error"),
    ("volumeError",   "volume-conservation error"),
]

# Diagnostics: band-restricted + half-time (maximum deformation) error measures.
DIAG_SERIES = [
    ("gradientErrorBand",     r"$\||\nabla\psi|-1\|_{L_2}$ band, $t=T$"),
    ("gradientErrorHalf",     r"$\||\nabla\psi|-1\|_{L_2}$, $t=T/2$"),
    ("gradientErrorBandHalf", r"$\||\nabla\psi|-1\|_{L_2}$ band, $t=T/2$"),
    ("volumeErrorHalf",       r"volume-conservation error, $t=T/2$"),
]

# Signed, bounded, non-monotone residuals: plotted but no fitted order quoted.
NO_ORDER = {"volumeError", "volumeErrorHalf"}

# Stable-envelope detection. The value-fitting reconstructions
# (linearWeightedLeastSquares, quadraticWeightedLeastSquares) are stable by
# construction -- they never differentiate the psi field -- so the full ladder is
# normally inside the envelope. The guard is retained as a safety net (and it is
# what truncates the field-differentiation Taylor family, whose signed-distance
# gradient runs away on the finest meshes). A resolution is inside the stable
# envelope of its (mesh, CFL) series iff, walking coarse -> fine, the gradient
# defect has not run away (< GRAD_CEILING; the elevated-but-bounded polyhedral
# baseline stays well under it) AND the shape error has not turned upward
# (<= running minimum x SHAPE_TOL). The convergence figures fit and draw only
# this prefix; any first destabilized resolution is marked separately.
GRAD_CEILING = 10.0   # |grad psi|-1 defect above this = runaway, not a mesh-noise floor
SHAPE_TOL    = 1.30   # shape error rising past 1.3x its running min = interface degrading


def _stable_prefix(rows, cfl):
    """Return (set of stable h, first-unstable h or None) for one (mesh, CFL)
    series, from the coarse->fine prefix where neither the gradient defect has
    run away nor the shape error has turned upward."""
    pts = sorted(
        ((_f(r.get("h")), _f(r.get("shapeError")), _f(r.get("gradientError")))
         for r in rows if _f(r.get("cfl")) == cfl),
        key=lambda t: (t[0] if t[0] else 0.0), reverse=True)   # coarse (large h) -> fine
    stable, limit, shape_min = set(), None, None
    for h, shape, grad in pts:
        if h is None:
            continue
        ok = True
        if grad is not None and abs(grad) > GRAD_CEILING:
            ok = False
        if shape is not None:
            if shape_min is not None and shape > shape_min * SHAPE_TOL:
                ok = False
            shape_min = shape if shape_min is None else min(shape_min, shape)
        if ok:
            stable.add(h)
        else:
            limit = h            # first destabilized resolution
            break                # destabilization is monotone: all finer are out too
    return stable, limit

# Per method line: theme, figure prefix, guide slopes, and per-case study sets.
# (2D is hexahedral only; the 3D cases overlay hex + polyhedral.)
METHODS = {
    "quadratic": {
        "theme": "semi-lagrangian-level-set",
        "prefix": "sl",
        "guides": (2,),
        "label": "UQWLSR",
        "cases": [
            {"case": "2Dvortex",
             "title": r"Semi-Lagrangian uncached quadratic weighted least-squares reconstruction (UQWLSR) — 2D reversed vortex ($T=2$)",
             "figsize": (15.0, 5.2),   # larger high-res panels for the 2D convergence figure
             "meshes": [("hex", "uncachedConv2Dvortex", "-", "o")]},
            {"case": "3Dshear",
             "title": r"Semi-Lagrangian UQWLSR reconstruction — 3D shear ($[0,1]^3$)",
             "meshes": [("hex",  "uncachedConv3Dshear",      "-",  "o"),
                        ("poly", "uncachedConv3DshearPoly",  "--", "s")]},
            {"case": "3Ddeformation",
             "title": r"Semi-Lagrangian UQWLSR reconstruction — 3D deformation (unit cube)",
             "meshes": [("hex",  "uncachedConv3Ddeformation",     "-",  "o"),
                        ("poly", "uncachedConv3DdeformationPoly", "--", "s")]},
        ],
    },
    "linear": {
        "theme": "linear-semi-lagrangian-level-set",
        "prefix": "lsl",
        "guides": (1, 2),
        # NOTE (2026-07-23): the lsl "linear line" framing is UNDER REVISION. linearWLSQ
        # FAILED SL advection (|grad psi| runaway); the linearConv3D* studies currently hold
        # the quadraticTaylor+clip campaign data. Label kept as the historical "nestedLSQ"
        # (= quadraticTaylor via RTS alias, what the configs still select) until the framing
        # decision is made in the lsl docs -- do not read this as an endorsement.
        "label": "nestedLSQ",
        "cases": [
            {"case": "2Dvortex",
             "title": r"Semi-Lagrangian least-squares linear (nestedLSQ) reconstruction — 2D reversed vortex ($T=2$)",
             "figsize": (15.0, 5.2),
             "meshes": [("hex", "linearConv2Dvortex", "-", "o")]},
            {"case": "3Dshear",
             "title": r"Semi-Lagrangian nestedLSQ reconstruction — 3D shear ($[0,1]^3$)",
             "meshes": [("hex",  "linearConv3Dshear",      "-",  "o"),
                        ("poly", "linearConv3DshearPoly",  "--", "s")]},
            {"case": "3Ddeformation",
             "title": r"Semi-Lagrangian nestedLSQ reconstruction — 3D deformation (unit cube)",
             "meshes": [("hex",  "linearConv3Ddeformation",     "-",  "o"),
                        ("poly", "linearConv3DdeformationPoly", "--", "s")]},
        ],
    },
}


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


def _guides(ax, allh, alle, guides):
    """Reference slopes anchored at the coarsest point."""
    if not allh:
        return
    h0, e0 = max(allh), max(alle)
    hr = np.array([min(allh), max(allh)])
    styles = {1: ("k--", r"$O(h)$"), 2: ("k:", r"$O(h^2)$"), 3: ("k-.", r"$O(h^3)$")}
    for p in guides:
        st, lbl = styles.get(p, ("k:", rf"$O(h^{p})$"))
        ax.loglog(hr, e0 * (hr / h0) ** p, st, lw=0.9, alpha=0.6, label=lbl)


def _panels(fig, axes, series, loaded, guides):
    """Fill one panel per error metric. Only the STABLE envelope of each
    (mesh, CFL) series is drawn as the convergence line and used for the order
    fit; the stability limit (first destabilized resolution) is marked with a
    hollow red cross so it is called out separately without blowing up the axis
    scale with the runaway values."""
    # stable prefix per (mesh, CFL) -- shared across all panels of the figure.
    masks = {}
    for mesh, style, mk, rows, cfls in loaded:
        for cfl in cfls:
            masks[(mesh, cfl)] = _stable_prefix(rows, cfl)
    for ax, (col, ttl) in zip(axes, series):
        allh, alle = [], []
        limit_pts = []
        for mesh, style, mk, rows, cfls in loaded:
            for cfl in cfls:
                stable_h, limit = masks[(mesh, cfl)]
                pts = sorted(
                    (_f(r["h"]), _f(r.get(col))) for r in rows
                    if _f(r.get("cfl")) == cfl and _f(r.get(col)) not in (None, 0.0)
                )
                # magnitudes: signed residuals (volume) plot on the log axis too
                pts = [(h, abs(e)) for h, e in pts if h and e]
                stab = [(h, e) for h, e in pts if h in stable_h]
                if not stab:
                    continue
                h, e = zip(*stab)
                p = _order(h, e) if col not in NO_ORDER else None
                lbl = f"{mesh} CFL {cfl:g}" + (f"  (order {p:.2f})" if p is not None else "")
                if limit is not None:
                    lbl += f"  [stable h≥{min(h):.1e}]"
                ax.loglog(h, e, linestyle=style, marker=mk, ms=6, label=lbl)
                allh += list(h); alle += list(e)
                # mark the stability limit: a red cross at the last stable point's
                # value, placed at the first-unstable h (the runaway value itself
                # is off-scale and deliberately not plotted).
                if limit is not None:
                    limit_pts.append((limit, e[-1]))
        for hx, ey in limit_pts:
            ax.loglog([hx], [ey], marker="x", ms=11, mew=2.4, color="#d32f2f",
                      linestyle="none",
                      label="stability limit" if hx == limit_pts[0][0] else None)
        _guides(ax, allh, alle, guides)
        ax.set_xlabel(r"$h$ (= $1/N$ hex, maxCellSize poly)", fontsize=12)
        ax.set_title(ttl, fontsize=13)
        ax.tick_params(labelsize=10)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)


def make(case, method):
    """Primary triptych + diagnostics figure for one case of one method line."""
    figdir = paths.figs_dir(method["theme"])
    prefix, guides = method["prefix"], method["guides"]
    # Load both meshes up front; skip the whole figure only if neither exists.
    loaded = [(mesh, style, mk, *_load(study))
              for mesh, study, style, mk in case["meshes"]]
    if not any(rows for *_, rows, _cfls in loaded):
        print(f"[sl3d] no studies for {case['case']} ({prefix}); skip")
        return None

    written = []
    for series, suffix, width in (
            (SERIES, "convergence", case.get("figsize", (14.0, 4.6))),
            (DIAG_SERIES, "convergence_diagnostics", (18.5, 4.6))):
        fig, axes = plt.subplots(1, len(series), figsize=width)
        _panels(fig, axes, series, loaded, guides)
        fig.suptitle(case["title"], fontsize=15)
        out = os.path.join(figdir, f"{prefix}_{case['case']}_{suffix}.png")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(out, dpi=200); plt.close(fig)
        written.append(out)
        nrows = sum(len(rows) for *_, rows, _cfls in loaded)
        print(f"[sl3d] wrote {out}  ({nrows} rows across "
              f"{sum(1 for *_, rows, _cfls in loaded if rows)} mesh type(s))")
    return written


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=sorted(METHODS), default="quadratic",
                    help="semi-Lagrangian method line (default: quadratic)")
    args = ap.parse_args()
    method = METHODS[args.method]
    any_ok = False
    for case in method["cases"]:
        if make(case, method):
            any_ok = True
    sys.exit(0 if any_ok else 1)
