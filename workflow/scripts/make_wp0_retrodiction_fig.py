#!/usr/bin/env python3
"""WP0 retrodiction figure: the band mode spectrum A{2h,4h,8h}(t) against
max|U|(t), minGradPsiBand(t) and kErrL2Band(t) from a droplet metrics CSV --
the leading-indicator test of docs/plan-curvature-stabilization.md WP0.

Usage:
    python3 workflow/scripts/make_wp0_retrodiction_fig.py out.png "label=case" [...]
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_NAME = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00"]
KEYS = ("TIME", "maxMagU", "A2hL2Band", "A4hL2Band", "A8hL2Band",
        "minGradPsiBand", "kErrL2Band")


def _series(case):
    out = {k: [] for k in KEYS}
    with open(os.path.join(case, CSV_NAME), newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                vals = {k: float(row[k]) for k in KEYS}
            except (TypeError, ValueError, KeyError):
                continue
            for k, v in vals.items():
                out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    out = argv[0]

    fig, axes = plt.subplots(3, 1, figsize=(7.6, 8.4), sharex=True,
                             gridspec_kw={"height_ratios": [1.4, 1.4, 1.0]})
    axU, axA, axG = axes

    for i, spec in enumerate(argv[1:]):
        label, case = spec.rsplit("=", 1)
        s = _series(case)
        if not len(s["TIME"]):
            print(f"[wp0] {case}: no rows, skipped")
            continue
        t = s["TIME"]
        base = COLORS[i % len(COLORS)]
        axU.semilogy(t, np.clip(s["maxMagU"], None, 10.0), "-", color=base,
                     lw=1.1, label=f"{label}: max|U|")
        axA.semilogy(t, s["A2hL2Band"], "-", color="#0072B2", lw=1.1,
                     label="A2h L2")
        axA.semilogy(t, s["A4hL2Band"], "-", color="#D55E00", lw=1.1,
                     label="A4h L2")
        axA.semilogy(t, s["A8hL2Band"], "-", color="#009E73", lw=1.1,
                     label="A8h L2")
        axG.plot(t, s["minGradPsiBand"], "-", color=base, lw=1.1,
                 label=f"{label}: min|∇ψ| band")
        g2 = axG.twinx()
        g2.semilogy(t, s["kErrL2Band"], ":", color="#7f7f7f", lw=1.1)
        g2.set_ylabel(r"band $L_2$ $|\kappa-\kappa_e|$ [1/m] (dotted)",
                      fontsize=8, color="#7f7f7f")
        print(f"[wp0] {label}: {len(t)} steps to t={t[-1]:.4g}")

    axU.set_ylabel(r"max$_c$ $|U|$ [m/s] (capped 10)")
    axA.set_ylabel(r"band mode amplitude [m]")
    axG.set_ylabel(r"min$|\nabla\psi|$ band")
    axG.set_xlabel(r"$t$ [s]")
    for ax in (axU, axA, axG):
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    axU.set_title("WP0 retrodiction: velocity, band mode spectrum, profile",
                  fontsize=11)
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[wp0] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
