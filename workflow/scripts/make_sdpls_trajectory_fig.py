#!/usr/bin/env python3
"""Band |grad psi| over the WHOLE trajectory, per SDPLS arm, with the strain.

The convergence tables report two instants, t = T/2 and t = T. That is enough to
fit an order and not enough to see a mechanism: it cannot show WHEN an arm loses
the signed-distance property, nor that `beta` is chasing a target which itself
moves with the flow. The negative-results deck has been making per-step claims
("band min |grad psi| -> 0 by t ~ 1, unchanged under R") with no figure behind
them; this is that figure.

Reads the per-time-step `gradPsiError.csv` that gradPsiErrorCSV writes in every
case directory -- which since phase 3 also carries the narrow-band statistics of
the SDPLS strain field a = n.grad(U).n (NARROW_{MIN,MEAN,MAX}_R). So the same
plot can show the measured band gradient AND the target beta - a it is chasing,
on one time axis.

Two panels:
  left   band min / mean / max |grad psi| against t, one colour per arm. The
         mean is solid, the min/max envelope is shaded. 1.0 is signed distance.
  right  the strain envelope min/max a, and the beta arm's implied target band
         beta - a. Where that band drops BELOW ZERO the fixed point is
         unreachable for a magnitude, and |grad psi| is driven toward zero
         instead -- the flattening mechanism, localised in time.

Cases are given as <label>=<case_dir>, which the curated CSV's `caseDir` column
now makes trivial to look up:

    python3 make_sdpls_trajectory_fig.py --out FIG.png --beta 1.0 \\
        noSource=studies/benchVortexEulerT2/2Dvortex_00024 \\
        R=studies/benchVortexEulerT2/2Dvortex_00026 \\
        beta=studies/benchVortexEulerT2/2Dvortex_00028
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

THEME = "sdpls-level-set"
COLS = ("TIME", "NARROW_MIN_MAG_GRAD_PSI", "NARROW_MEAN_MAG_GRAD_PSI",
        "NARROW_MAX_MAG_GRAD_PSI", "NARROW_MIN_R", "NARROW_MEAN_R",
        "NARROW_MAX_R")


def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def series(case_dir):
    """{column: np.array} from a case's per-step gradPsiError.csv."""
    path = os.path.join(case_dir, "gradPsiError.csv")
    if not os.path.isfile(path):
        return None
    out = {c: [] for c in COLS}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            for c in COLS:
                out[c].append(_f(r.get(c)))
    if not out["TIME"]:
        return None
    return {c: np.array([np.nan if v is None else v for v in out[c]], dtype=float)
            for c in out}


def label_of(case_dir, fallback):
    """N from the case sidecar, so the caption cannot drift from the data."""
    try:
        tok = json.load(open(os.path.join(case_dir, "case_params.json")))["tokens"]
        return f"{fallback} (N={tok.get('N_CELLS', '?')})"
    except (OSError, KeyError, ValueError):
        return fallback


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", metavar="LABEL=CASE_DIR")
    ap.add_argument("--out", default=None, help="output PNG (default: theme figures)")
    ap.add_argument("--beta", type=float, default=1.0,
                    help="beta of the sdplsBeta arm, for the target band")
    args = ap.parse_args(argv)

    arms = []
    for spec in args.cases:
        if "=" not in spec:
            print(f"[traj] expected LABEL=CASE_DIR, got {spec!r}")
            return 1
        lab, d = spec.split("=", 1)
        s = series(d)
        if s is None:
            print(f"[traj] no usable gradPsiError.csv in {d}; skip {lab}")
            continue
        arms.append((lab, d, s))
    if not arms:
        print("[traj] no usable cases; nothing written")
        return 1

    fig, (axg, axs) = plt.subplots(1, 2, figsize=(13, 4.6))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, (lab, d, s) in enumerate(arms):
        c = colours[i % len(colours)]
        t = s["TIME"]
        axg.plot(t, s["NARROW_MEAN_MAG_GRAD_PSI"], color=c, lw=1.6,
                 label=label_of(d, lab))
        axg.fill_between(t, s["NARROW_MIN_MAG_GRAD_PSI"],
                         s["NARROW_MAX_MAG_GRAD_PSI"], color=c, alpha=0.16,
                         linewidth=0)
    axg.axhline(1.0, color="k", lw=0.9, ls="--")
    axg.set_xlabel("t")
    axg.set_ylabel(r"band $|\nabla\psi|$  (mean, with min/max envelope)")
    axg.set_title("Signed-distance property over the trajectory", fontsize=11)
    # LOG scale, or the arms are not comparable: sdplsBeta's band maximum
    # reaches ~155 just after the flow reverses while R stays under 1.3, so on a
    # linear axis every well-behaved arm collapses onto the bottom edge.
    axg.set_yscale("log")
    axg.legend(fontsize=8)

    # The strain is a property of the FLOW: identical for every arm to plotting
    # accuracy, so show it once, from the first arm that has it.
    for lab, d, s in arms:
        if np.all(np.isnan(s["NARROW_MAX_R"])):
            continue
        t = s["TIME"]
        axs.fill_between(t, s["NARROW_MIN_R"], s["NARROW_MAX_R"],
                         color="0.6", alpha=0.35, linewidth=0,
                         label=r"strain $a$ envelope (min..max over the band)")
        axs.plot(t, s["NARROW_MEAN_R"], color="0.25", lw=1.2,
                 label=r"band mean $a$")
        # The target beta - a. Where the UPPER strain exceeds beta, the lower
        # edge of this band is negative and unreachable for a magnitude.
        lo = args.beta - s["NARROW_MAX_R"]
        hi = args.beta - s["NARROW_MIN_R"]
        axs.fill_between(t, lo, hi, color="tab:red", alpha=0.13, linewidth=0,
                         label=rf"target $\beta-a$ band ($\beta$={args.beta:g})")
        unreachable = np.nansum(lo < 0.0)
        frac = unreachable / max(1, np.sum(~np.isnan(lo)))
        axs.plot(t, np.where(lo < 0.0, lo, np.nan), color="tab:red", lw=2.0,
                 label=rf"$\beta-a<0$: unreachable ({100*frac:.0f}% of steps)")
        break
    axs.axhline(0.0, color="k", lw=0.9)
    axs.axhline(args.beta, color="tab:red", lw=0.8, ls=":")
    axs.set_xlabel("t")
    axs.set_ylabel(r"strain $a=\hat n\cdot\nabla\mathbf{u}\cdot\hat n$")
    axs.set_title(r"What $\beta$ is chasing", fontsize=11)
    axs.legend(fontsize=7.5, loc="best")

    fig.tight_layout()
    out = args.out or os.path.join(paths.figs_dir(THEME), "sdpls_trajectory.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[traj] wrote {out}  ({len(arms)} arm(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
