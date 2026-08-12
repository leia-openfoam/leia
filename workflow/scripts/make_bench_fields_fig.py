#!/usr/bin/env python3
"""Field atlases for the method-decision benchmark (uniform 2D hex cases).

For every method in the benchVortex* / sdpls* studies and every field in
{alpha, psi, graderr}, one montage: rows = resolutions N (32..512), cols =
snapshots {0, T/2, T}. Rendering follows plots.py::_field_montage (foamlib
FoamCase + reshape; ||grad psi|-1| via np.gradient).

These are DEBUGGING atlases: alpha shows where the phase indicator smears, psi
where the level set deforms, and graderr where the signed-distance error lives.
A convergence order says a method got better; only these say WHERE.

CAVEAT, and it belongs in every caption: `graderr` here is ||grad psi|-1|
computed with np.gradient on the reshaped uniform grid. That is NOT the solver's
`gradPsiMetric` operator, which is what the curated CSV's gradientError columns
measure. The atlas localises the error; the CSV measures it. They will not agree
digit-for-digit and must not be presented as if they do.

Output theme follows the STUDY, matching the repo's per-theme data convention:
sdpls* studies feed docs/sdpls-level-set, everything else method-comparison
(where methodComparison.tex and the comparison deck already reference 36
existing bench_*.png by name -- do not move those).

    python3 make_bench_fields_fig.py [studies_glob ...]   (default: both groups)
"""
import method_label
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from foamlib import FoamCase

import paths

THEME = "method-comparison"
SDPLS_THEME = "sdpls-level-set"


def theme_of(study_dir):
    """Per-theme data, keyed on the study that produced the case."""
    return SDPLS_THEME if os.path.basename(study_dir.rstrip("/")).startswith("sdpls") \
        else THEME

FIELDS = {
    "alpha":   dict(cmap="Blues", title=r"$\alpha$"),
    "psi":     dict(cmap="RdBu",  title=r"$\psi$"),
    "graderr": dict(cmap="inferno", title=r"$||\nabla\psi|-1|$ ($\log_{10}$)"),
}


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def method_of(t):
    """Filename component. Delegates to the SHARED definition so it can never
    drift from the curated CSV's `method` column again -- this function used to
    ignore SOURCE_SCHEME, so a study sweeping both SDPLS linearizations wrote
    both arms to the SAME PNG and the second silently overwrote the first."""
    return method_label.method_slug(t)


def snapshots(case, T):
    idxs = []
    times = [float(f.time) for f in case]
    for target in (0.0, T/2, T):
        idxs.append(min(range(len(times)), key=lambda i: abs(times[i]-target)))
    return idxs, times


def read2d(case, idx, name, n):
    fld = np.asarray(case[idx][name].internal_field, dtype=float)
    if fld.size == 1:
        fld = np.full(n*n, float(fld))
    return fld.reshape(n, n)


def montage(cases, T, field, out):
    # cases: {N: casedir} sorted
    ns = sorted(cases)
    fig, axes = plt.subplots(len(ns), 3,
                             figsize=(3*2.9, len(ns)*2.9), squeeze=False)
    for i, n in enumerate(ns):
        case = FoamCase(cases[n])
        idxs, times = snapshots(case, T)
        h = 1.0/n
        for j, idx in enumerate(idxs):
            ax = axes[i][j]
            psi2d = read2d(case, idx, "psi", n)
            if field == "alpha":
                f2d = read2d(case, idx, "alpha", n)
                im = ax.imshow(f2d, origin="lower", extent=(0, 1, 0, 1),
                               cmap="Blues", vmin=0, vmax=1)
            elif field == "psi":
                im = ax.imshow(psi2d, origin="lower", extent=(0, 1, 0, 1),
                               cmap="RdBu")
            else:
                gy, gx = np.gradient(psi2d, h, h)
                dev = np.abs(np.sqrt(gx**2 + gy**2) - 1.0)
                im = ax.imshow(np.log10(np.maximum(dev, 1e-12)),
                               origin="lower", extent=(0, 1, 0, 1),
                               cmap="inferno", vmin=-6, vmax=1)
            x = (np.arange(n) + 0.5)*h
            X, Y = np.meshgrid(x, x)
            ax.contour(X, Y, psi2d, levels=[0], colors="#4fc3f7",
                       linewidths=0.8)
            ax.set_xticks([]); ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(f"N={n}", fontsize=9)
            if i == 0:
                ax.set_title(f"t = {times[idx]:.3g}", fontsize=9)
    fig.colorbar(im, ax=[a for r in axes for a in r], shrink=0.6, pad=0.01)
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[bench_fields] wrote {out}")


def main():
    repo = paths.REPO
    # Deliberately NOT `sdpls*`: sdplsStability is single-resolution and sweeps
    # CFL, which is not part of the group key below, so its cases would collide
    # on (method, T, N) and silently overwrite each other. Only resolution
    # ladders with a pinned CFL belong here.
    patterns = sys.argv[1:] or [
        os.path.join(repo, "studies", "benchVortex*"),
        os.path.join(repo, "studies", "sdplsBetaSweep"),
    ]
    study_dirs = [d for p in patterns for d in glob.glob(p)]
    # (theme, method, T) -> {N: case_dir}
    groups = {}
    for sd in study_dirs:
        theme = theme_of(sd)
        for d in glob.glob(os.path.join(sd, "2Dvortex_*")):
            pj = os.path.join(d, "case_params.json")
            if not os.path.isfile(pj):
                continue
            t = json.load(open(pj))["tokens"]
            key = (theme, method_of(t), t["END_TIME"])
            n = int(t["N_CELLS"])
            prev = groups.setdefault(key, {}).get(n)
            if prev is not None and prev != d:
                # Two cases share (method, T, N): some axis the key does not
                # carry (CFL, say) is being swept. Say so -- silently keeping one
                # would publish an atlas labelled for a configuration it is not.
                print(f"[bench_fields] COLLISION {key} N={n}: "
                      f"{os.path.basename(prev)} vs {os.path.basename(d)} -- "
                      f"keeping the first; the study sweeps an axis the figure "
                      f"filename does not encode")
                continue
            groups[key][n] = d
    if not groups:
        print(f"[bench_fields] no 2Dvortex cases under {patterns}; nothing to do")
        return
    for (theme, m, T), cases in sorted(groups.items()):
        figs_dir = paths.figs_dir(theme)
        for field in FIELDS:
            out = os.path.join(figs_dir, f"bench_{m}_T{T}_{field}.png")
            try:
                montage(cases, float(T), field, out)
            except Exception as e:
                print(f"[bench_fields] skip {m} T={T} {field}: {e}")


if __name__ == "__main__":
    main()
