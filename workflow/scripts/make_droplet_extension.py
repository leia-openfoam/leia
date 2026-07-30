#!/usr/bin/env python3
"""Stationary-droplet N=128: curvature normal extension ON vs OFF.

Reads the two committed per-step time series under
``docs/.../sl-level-set-article/data/mechanism/``:

    ext_N128.csv     curvatureNewtonIters=3  -- kappa evaluated at the interface
                     foot point, extended constant along the normal (the default)
    noext_N128.csv   curvatureNewtonIters=0  -- kappa evaluated at the CELL CENTRE,
                     no normal extension

and renders ``data/figures/sl_droplet_extension.png``: max|u|(t) for both. Both use
the strict pressure solve (rel. residual 1e-12), so the only difference is the
curvature extension. The extended run diverges on this fine mesh; the un-extended
(cell-centre) run stays bounded -- the foot-point projection on the reinitialisation-
free (drifting) level set is the destabilising element.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "..", "..", "docs", "semi-lagrangian-level-set",
                   "sl-level-set-article", "data")
MECH = os.path.join(ART, "mechanism")
FIGS = os.path.join(ART, "figures")


def _read(name):
    t, u = [], []
    with open(os.path.join(MECH, name)) as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["TIME"])); u.append(float(row["maxMagU"]))
    return t, u


def main():
    te, ue = _read("ext_N128.csv")
    tn, un = _read("noext_N128.csv")
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.semilogy(te, ue, "-", color="#D55E00", lw=1.8,
                label=r"with normal extension (foot point): diverges $t\!\approx\!0.06$")
    ax.semilogy(tn, un, "-", color="#0072B2", lw=1.8,
                label=r"without extension (cell centre): diverges $t\!\approx\!0.11$")
    ax.set_xlabel(r"time $t$ [s]")
    ax.set_ylabel(r"$\max|\mathbf{u}|$ [m/s]")
    ax.set_title(r"$N=128$: normal extension delays, does not cure, the divergence")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(frameon=False, fontsize=9, loc="best")
    fig.tight_layout()
    out = os.path.join(FIGS, "sl_droplet_extension.png")
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"[extension] wrote {out}")


if __name__ == "__main__":
    main()
