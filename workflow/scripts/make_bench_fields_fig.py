#!/usr/bin/env python3
"""Field atlases for the method-decision benchmark (uniform 2D hex cases).

For every method in the benchVortex* studies and every field in
{alpha, psi, graderr}, one montage: rows = resolutions N (32..512), cols =
snapshots {0, T/2, T}. Rendering follows plots.py::_field_montage (foamlib
FoamCase + reshape; ||grad psi|-1| via np.gradient). PNGs land in the
sdpls-level-set theme data/figures for the deck field atlas and the article.

    python3 make_bench_fields_fig.py [studies_glob]   (default studies/benchVortex*)
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
    pattern = sys.argv[1] if len(sys.argv) > 1 else None
    repo = paths.REPO
    figs_dir = paths.figs_dir(THEME)
    study_dirs = glob.glob(pattern) if pattern else \
        glob.glob(os.path.join(repo, "studies", "benchVortex*"))
    # (method, T) -> {N: case_dir}
    groups = {}
    for sd in study_dirs:
        for d in glob.glob(os.path.join(sd, "2Dvortex_*")):
            pj = os.path.join(d, "case_params.json")
            if not os.path.isfile(pj):
                continue
            t = json.load(open(pj))["tokens"]
            key = (method_of(t), t["END_TIME"])
            groups.setdefault(key, {})[int(t["N_CELLS"])] = d
    for (m, T), cases in sorted(groups.items()):
        for field in FIELDS:
            out = os.path.join(figs_dir, f"bench_{m}_T{T}_{field}.png")
            try:
                montage(cases, float(T), field, out)
            except Exception as e:
                print(f"[bench_fields] skip {m} T={T} {field}: {e}")


if __name__ == "__main__":
    main()
