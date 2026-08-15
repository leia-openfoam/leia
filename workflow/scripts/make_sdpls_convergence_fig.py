#!/usr/bin/env python3
"""Convergence panels for the BEST SDPLS setup, on the two cases where a
convergence claim is defensible.

WHICH CASES. The 2D reversed vortex and the 3D shear. The 3D LeVeque
deformation is deliberately absent: at t=T/2 its interface still carries
cell-scale structure at N=203 (the thin rim degenerates into a comb of
cell-width teeth; mean sheet thickness 5.4 cells but the rim is ~1), so its
fitted orders describe under-resolution rather than the scheme. It is a
robustness showcase, not a convergence study.

WHICH SETUP, and why each choice is the measured one:
  source          R          beta relaxes to g* = beta - a, a target set by the
                             FLOW, so it cannot converge (2D order -0.04).
  linearization   simpleImp  the two admissible linearizations agree to four
                             digits on every entry; either would do.
  cut-off         none       both cut-offs fix the far field and destroy the
                             band error. The |psi| threshold degrades 3D shear
                             (volume +0.23 -> -0.97); the topological band flattens
                             the far field to -0.02 but takes the band error to
                             ~1.5 with order -0.03, near the sourceless value.
                             nLayers 1..5 changes it by <10%, so it is not a
                             width problem: the level set ADVECTS through the
                             band, and material arrives carrying an uncorrected
                             gradient.
  CFL 0.5, dc 3   converged  measured: refining dt 4x at fixed h moves the band
                             error by 0.3-1.6%, and nDefCorr 3 -> 10 by <0.4%.

WHICH INSTANT. The two-instant protocol. On a reversed flow the interfacial
stretching integrates to zero over the period, so the t=T gradient error
compares a residual against an exactly-cancelling control and inverts the
ranking. Gradient errors are therefore read at t=T/2 (maximal accumulated
deformation, nothing cancelled yet); shape error only at t=T, where it is
defined against the initial condition; volume at both.

Usage (from repo root):
    python3 workflow/scripts/make_sdpls_convergence_fig.py
"""
import csv
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

REPO = paths.REPO
THEME = "sdpls-level-set"

# (study, case label, arm-matching predicate on the method string)
CASES = [
    ("benchVortexEulerT2", "2D reversed vortex, $T=2$"),
    ("sdplsConv3Dshear",   "3D shear, $T=3$"),
]

# (csv column, panel title, instant tag). Every error the curated table carries
# for these runs, each at the instant where it means something.
PANELS = [
    ("gradientErrorBandHalf", r"band $\||\nabla\psi|-1\|_{L_2}$",      "$t=T/2$"),
    ("gradientErrorHalf",     r"global $\||\nabla\psi|-1\|_{L_2}$",    "$t=T/2$"),
    ("volumeErrorHalf",       r"volume error",                          "$t=T/2$"),
    ("volumeError",           r"volume error",                          "$t=T$"),
    ("shapeError",            r"shape error $E_{\mathrm{geom}}$",       "$t=T$"),
    # DIAGNOSTIC, not an error: the ideal value is 1, not 0. Plotted without
    # reference slopes and without a fitted "order", because a slope on a
    # quantity that should approach a nonzero constant means nothing.
    ("minGradPsiBandHalf",    r"band $\min|\nabla\psi|$  (flattening; ideal $=1$)", "$t=T/2$"),
]

# Columns whose target is 1, not 0.
DIAGNOSTIC = {"minGradPsiBandHalf", "minGradPsiBand"}

ARMS = [
    ("euler+SDPLS:R/simpleImp", "SDPLS: R",  "#1f4e9c", "o", "-"),
    ("euler+div",               "no source", "#b03030", "s", "--"),
]


def _f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def _order(hs, es):
    pts = [(h, abs(e)) for h, e in zip(hs, es) if h and e and abs(e) > 0]
    if len(pts) < 2:
        return None
    h, e = zip(*pts)
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def _rows(study):
    p = os.path.join(REPO, "studies", study, f"{study}_errors.csv")
    if not os.path.isfile(p):
        p = os.path.join(paths.tables_dir(THEME), f"{study}_errors.csv")
    if not os.path.isfile(p):
        return None
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def make(study, title):
    rows = _rows(study)
    if not rows:
        print(f"[sdplsconv] no data for {study}; skipped")
        return None

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.6))
    fig.suptitle(f"SDPLS convergence -- {title}\n"
                 r"$R$ (strain-cancelling), simpleLinearImplicit, no cut-off, "
                 r"CFL 0.5, $n_{\mathrm{defCorr}}=3$",
                 fontsize=12.5, y=0.985)

    for ax, (col, name, inst) in zip(axes.ravel(), PANELS):
        plotted = False
        for key, label, colour, marker, ls in ARMS:
            sel = [r for r in rows if r["method"].startswith(key)]
            # collapse duplicate (h) rows: an inert axis crossed against an arm
            seen, hs, es = set(), [], []
            for r in sorted(sel, key=lambda r: -(_f(r["h"]) or 0)):
                h, e = _f(r["h"]), _f(r.get(col))
                if h is None or h in seen:
                    continue
                seen.add(h)
                if e is not None and e > 0:
                    hs.append(h); es.append(e)
            if len(hs) < 2:
                continue
            if col in DIAGNOSTIC:
                lab = f"{label}   ({es[0]:.2f} at coarsest, {es[-1]:.2f} at finest)"
            else:
                p = _order(hs, es)
                lab = f"{label}" + (f"   $p={p:+.2f}$" if p is not None else "")
            ax.loglog(hs, es, marker=marker, ls=ls, color=colour, lw=1.8,
                      ms=6, label=lab)
            plotted = True

        if plotted and col in DIAGNOSTIC:
            ax.axhline(1.0, color="0.35", ls="-", lw=1.2)
            ax.annotate("signed distance", (ax.get_xlim()[0]*1.3, 1.0),
                        color="0.35", fontsize=8.5, va="bottom")
        elif plotted:
            # first- and second-order reference slopes, anchored to the axes
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            xr = np.array([x0*1.25, x1*0.8])
            for order, style in ((1.0, ":"), (2.0, "-.")):
                yr = (y0*y1)**0.5 * (xr/xr[0])**order
                ax.loglog(xr, yr, style, color="0.55", lw=1.1)
                ax.annotate(f"$h^{{{order:.0f}}}$", (xr[-1], yr[-1]),
                            color="0.45", fontsize=8.5,
                            ha="left", va="center")
        ax.set_xlabel(r"$h$")
        ax.set_title(f"{name}   ({inst})", fontsize=10.5)
        ax.grid(True, which="both", alpha=0.25, lw=0.6)
        ax.legend(fontsize=8.5, loc="best", framealpha=0.9)

    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    out = os.path.join(paths.figs_dir(THEME), f"sdpls_convergence_{study}.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[sdplsconv] wrote {out}")
    return out


def main():
    made = [make(s, t) for s, t in CASES]
    return 0 if any(made) else 1


if __name__ == "__main__":
    sys.exit(main())
