#!/usr/bin/env python3
"""WP8.1 figure: where the delivered face curvature varies during a coupled run.

Plots, from the per-step droplet metrics of one or more cases:
  (top)    max|U|(t) -- the parasitic current
  (middle) the relative variation rate of the delivered kappa_f ACROSS the CSF
           support (face-pair separations aligned with the interface normal)
           and ALONG the interface, both in 1/m
  (bottom) their ratio, and the angle between the quadratic-fit normal used for
           the classification and the linear-plane normal (right axis)

The middle panel answers the WP8 ordering question: a delivery that attaches
one curvature per interface element and extends it along the normal sets the
ACROSS curve to zero by construction, so it can only help if that curve
dominates -- and keeps dominating once the flow perturbs psi.

Usage:
    python3 workflow/scripts/make_driver_split_fig.py out.png "label=case" [...]
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_NAME = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
KEYS = ("TIME", "maxMagU", "driverAcrossSupportL2", "driverAlongInterfaceL2",
        "driverNormalDisagreementMeanDeg", "minGradPsiBand")
COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00"]
UCAP = 10.0     # velocity axis ceiling: above this the run has diverged


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

    fig, (axU, axD, axR) = plt.subplots(
        3, 1, figsize=(7.8, 8.8), sharex=True,
        gridspec_kw={"height_ratios": [1.2, 1.6, 1.0]})

    for i, spec in enumerate(argv[1:]):
        label, case = spec.rsplit("=", 1)
        s = _series(case)
        if not len(s["TIME"]):
            print(f"[split] {case}: no rows, skipped")
            continue
        t = s["TIME"]
        c = COLORS[i % len(COLORS)]

        axU.semilogy(t, np.clip(s["maxMagU"], None, UCAP), "-", color=c, lw=1.1,
                     label=label)
        axD.semilogy(t, s["driverAcrossSupportL2"], "-", color=c, lw=1.3,
                     label=f"{label}: across support")
        axD.semilogy(t, s["driverAlongInterfaceL2"], "--", color=c, lw=1.1,
                     label=f"{label}: along interface")
        axR.plot(t, s["driverAcrossSupportL2"]
                 / np.maximum(s["driverAlongInterfaceL2"], 1e-30),
                 "-", color=c, lw=1.1, label=f"{label}: across/along")

        print(f"[split] {label}: {len(t)} steps to t={t[-1]:.4g}; "
              f"across {s['driverAcrossSupportL2'][0]:.3e} -> "
              f"{s['driverAcrossSupportL2'][-1]:.3e}, "
              f"along {s['driverAlongInterfaceL2'][0]:.3e} -> "
              f"{s['driverAlongInterfaceL2'][-1]:.3e}, "
              f"normal disagreement {s['driverNormalDisagreementMeanDeg'][0]:.2f}"
              f" -> {s['driverNormalDisagreementMeanDeg'][-1]:.2f} deg")

    # normal-disagreement history on the right axis of the bottom panel
    axN = axR.twinx()
    for i, spec in enumerate(argv[1:]):
        label, case = spec.rsplit("=", 1)
        s = _series(case)
        if not len(s["TIME"]):
            continue
        axN.semilogy(s["TIME"], np.maximum(s["driverNormalDisagreementMeanDeg"],
                                           1e-3),
                     ":", color=COLORS[i % len(COLORS)], lw=1.0)
    axN.set_ylabel("fit vs plane normal [deg] (dotted)", fontsize=8,
                   color="#555555")

    axU.set_ylabel(r"max$_c\,|U|$ [m/s]")
    axU.set_title("Where the delivered face curvature varies during the run",
                  fontsize=11)
    axD.set_ylabel(r"relative variation rate of $\kappa_f$ [1/m]")
    axR.set_ylabel("across / along")
    axR.set_xlabel(r"$t$ [s]")
    axR.axhline(1.0, color="#7f7f7f", lw=0.8, ls="-", alpha=0.7)
    for ax in (axU, axD, axR):
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.align_ylabels((axU, axD, axR))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[split] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
