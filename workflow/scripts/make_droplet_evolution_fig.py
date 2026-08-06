#!/usr/bin/env python3
"""Parasitic-current EVOLUTION figure: max|U|(t) (+ the band curvature-error
diagnostic) from the per-step droplet metrics CSV of one or more cases.

Unlike make_droplet_table.py (settled values vs h), this plots the full
time histories -- the observable of the stationaryDropletStableFoot gate:
where the currents originate, whether they settle, and how the settled level
scales under refinement.

Usage:
    python3 workflow/scripts/make_droplet_evolution_fig.py out.png \
        "label=path/to/case" "label2=path/to/case2" [...]

Each case dir must contain leiaSemiLagrangianLevelSetTwoPhaseFoam.csv
(TIME, maxMagU, ..., kErrL2Band, ...). Labels render verbatim in the legend.
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_NAME = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"

# Fixed entity order (validated categorical palette; assign by CLI order --
# strongest hue contrast between the first two series, they carry the story).
COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7"]


def _series(case):
    t, umax, kerr = [], [], []
    with open(os.path.join(case, CSV_NAME), newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                # Parse the full row BEFORE appending anything, so a truncated
                # trailing line of a live run cannot desynchronise the arrays.
                ti = float(row["TIME"])
                ui = float(row["maxMagU"])
                ki = float(row.get("kErrL2Band") or "nan")
            except (TypeError, ValueError, KeyError):
                continue   # truncated row of a live/interrupted run
            t.append(ti); umax.append(ui); kerr.append(ki)
    return np.array(t), np.array(umax), np.array(kerr)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    out = argv[0]

    fig, (ax, axk) = plt.subplots(
        2, 1, figsize=(7.6, 6.6), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]})

    # Cap the velocity axis at a physical ceiling: a sigma-driven droplet with
    # |U| above O(1) m/s has already diverged, and the ~1e70 m/s FPE tails
    # would otherwise compress the decades where the mechanism lives. A
    # blow-up is marked where the trace first leaves the cap.
    CAP = 10.0
    for i, spec in enumerate(argv[1:]):
        # Split at the LAST '=' -- labels may contain '=' ("N=256").
        label, case = spec.rsplit("=", 1)
        t, umax, kerr = _series(case)
        if not len(t):
            print(f"[evo] {case}: no rows, skipped")
            continue
        c = COLORS[i % len(COLORS)]
        ax.semilogy(t, umax, "-", color=c, lw=1.1, label=label)
        if (umax > CAP).any():
            tb = t[int(np.argmax(umax > CAP))]
            ax.axvline(tb, color=c, ls=":", lw=1.0, alpha=0.8)
            ax.annotate(f"blow-up $t$={tb:.4f}", xy=(tb, 0.3*CAP),
                        rotation=90, ha="right", va="top", fontsize=8, color=c)
        if np.isfinite(kerr).any():
            axk.semilogy(t, kerr, "-", color=c, lw=1.1)
        print(f"[evo] {label}: {len(t)} steps to t={t[-1]:.4g}, "
              f"max|U| last={umax[-1]:.3e}, peak={umax.max():.3e}")
    ax.set_ylim(top=CAP)

    ax.set_ylabel(r"max$_c$ $|U|$ [m/s]")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.set_title("Stationary droplet: parasitic-current evolution", fontsize=11)
    axk.set_ylabel(r"band $L_2$ $|\kappa-\kappa_e|$ [1/m]")
    axk.set_xlabel(r"$t$ [s]")
    axk.grid(True, which="both", ls=":", alpha=0.5)
    fig.align_ylabels((ax, axk))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[evo] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
