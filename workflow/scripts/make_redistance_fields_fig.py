#!/usr/bin/env python3
"""Field plots for the static idempotency gate (redistanceStatic2D).

For each SURFACE and a chosen resolution (default N=64), one montage across
the redistancer models: input psi (0/psi) contours, output psi (1/psi)
contours, and the per-event deviation |psi_out - psi_in| (log scale). Plain
OpenFOAM-ascii parsing of the uniform N x N hex fields -- no extra deps.

    python3 make_redistance_fields_fig.py <STUDY_DIR> [N]
"""
import csv
import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

THEME = "geometrically-redistanced-levelset"


def read_scalar_field(path, ncells):
    """internalField of an ascii volScalarField -> 1D numpy array."""
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?"
                  r"(\d+)\s*\n?\(", txt)
    if m:
        start = txt.index("(", m.end() - 1) + 1
        end = txt.index(")", start)
        vals = np.fromstring(txt[start:end], sep="\n")
        if vals.size != int(m.group(1)):
            vals = np.array([float(x) for x in txt[start:end].split()])
        return vals
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)", txt)
    if m:
        return np.full(ncells, float(m.group(1)))
    raise ValueError(f"no internalField in {path}")


def main():
    study_dir = sys.argv[1]
    nshow = sys.argv[2] if len(sys.argv) > 2 else "64"
    study = os.path.basename(os.path.normpath(study_dir))
    figs_dir = paths.figs_dir(THEME)

    # collect cases: (surface, redistancer) -> case dir at N = nshow
    cases = {}
    for d in sorted(glob.glob(os.path.join(study_dir, "*_[0-9]*"))):
        pj = os.path.join(d, "case_params.json")
        if not os.path.isfile(pj):
            continue
        t = json.load(open(pj))["tokens"]
        if t.get("N_CELLS") != nshow:
            continue
        cases[(t.get("SURFACE", ""), t["REDISTANCER"])] = d

    n = int(nshow)
    x = (np.arange(n) + 0.5)/n
    X, Y = np.meshgrid(x, x)

    for surface in sorted({s for s, _ in cases}):
        models = sorted({r for s, r in cases if s == surface})
        fig, axes = plt.subplots(2, len(models),
                                 figsize=(3.1*len(models), 6.2))
        for j, redist in enumerate(models):
            d = cases[(surface, redist)]
            try:
                psi0 = read_scalar_field(os.path.join(d, "0", "psi"), n*n)
                psi1 = read_scalar_field(os.path.join(d, "1", "psi"), n*n)
            except (FileNotFoundError, ValueError) as e:
                print(f"[fields] skip {redist}/{surface}: {e}")
                continue
            P0 = psi0.reshape(n, n)
            P1 = psi1.reshape(n, n)
            ax = axes[0, j]
            ax.contour(X, Y, P0, levels=[0], colors="k", linewidths=1.0)
            cs = ax.contour(X, Y, P1, levels=[0], colors="r",
                            linewidths=1.0, linestyles="--")
            ax.contour(X, Y, P1, levels=np.linspace(-0.4, 0.4, 9),
                       colors="gray", linewidths=0.4, alpha=0.6)
            ax.set_title(f"{redist}\n zero sets: in (black) / out (red)",
                         fontsize=8)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

            ax = axes[1, j]
            dev = np.abs(P1 - P0)
            im = ax.imshow(np.log10(np.maximum(dev, 1e-16)),
                           origin="lower", extent=(0, 1, 0, 1),
                           vmin=-16, vmax=0, cmap="magma")
            ax.contour(X, Y, P0, levels=[0], colors="w", linewidths=0.6)
            ax.set_title(r"$\log_{10}|\psi_{out}-\psi_{in}|$"
                         f"  max={dev.max():.1e}", fontsize=8)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=axes[1, :].tolist(), shrink=0.8, pad=0.01)
        fig.suptitle(f"{study} [{surface}] N={nshow}: one event, "
                     "deviation from exact-SDF input", fontsize=11)
        out = os.path.join(figs_dir, f"{study}_fields_{surface}.png")
        fig.savefig(out, dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"[make_redistance_fields_fig] wrote {out}")


if __name__ == "__main__":
    main()
