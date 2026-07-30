#!/usr/bin/env python3
"""Stationary-droplet divergence-mechanism figure.

Reads the three compact per-step time series produced by
``workflow/scripts/droplet_mechanism_probe.sh`` (committed under
``docs/.../sl-level-set-article/data/mechanism/``):

    stable_N64.csv      N=64, standard transport  -> currents DECAY
    diverging_N128.csv  N=128, standard transport -> currents GROW (diverge)
    frozen_N128.csv     N=128, interface FROZEN    -> currents DECAY

and renders ``data/figures/sl_droplet_mechanism.png``:

  (a) max|u|(t) for the three runs (semilog): N=64 and the frozen N=128 both
      relax; only the standard N=128 diverges. Freezing psi -- holding it at the
      clean signed distance -- removes the divergence at the SAME resolution and
      step, so the cause is the reinitialisation-free transport, not the flow
      solver or the curvature formula.
  (b) the N=128 feedback: min|grad psi| in the band (right axis) collapses from 1
      while the band curvature error (left axis) runs away, BOTH before max|u|
      grows -- the signed-distance loss leads the current growth.

Pure-matplotlib, no seaborn; colours are colour-blind friendly.
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
    path = os.path.join(MECH, name)
    t, u, kL2, gm = [], [], [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["TIME"]))
            u.append(float(row["maxMagU"]))
            kL2.append(float(row["kErrL2Band"]))
            # min|grad psi| is set to GREAT (~1e15) once the band degenerates
            # after blow-up; clip so the axis stays meaningful.
            g = float(row["minGradPsiBand"])
            gm.append(g if g < 5.0 else float("nan"))
    return t, u, kL2, gm


def main():
    t64, u64, _, _ = _read("stable_N64.csv")
    t128, u128, k128, g128 = _read("diverging_N128.csv")
    tf, uf, _, _ = _read("frozen_N128.csv")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    # ---- (a) parasitic current vs time -----------------------------------
    axA.semilogy(t64, u64, "-", color="#0072B2", lw=1.8,
                 label=r"$N=64$ (standard)")
    axA.semilogy(t128, u128, "-", color="#D55E00", lw=1.8,
                 label=r"$N=128$ (standard)")
    axA.semilogy(tf, uf, "--", color="#009E73", lw=1.8,
                 label=r"$N=128$, interface frozen")
    axA.set_xlabel(r"time $t$ [s]")
    axA.set_ylabel(r"$\max|\mathbf{u}|$ [m/s]")
    axA.set_title(r"(a) parasitic current")
    axA.grid(True, which="both", ls=":", alpha=0.4)
    axA.legend(fontsize=8, loc="best")

    # ---- (b) the N=128 feedback: |grad psi| collapse leads ---------------
    axB.semilogy(t128, k128, "-", color="#D55E00", lw=1.8,
                 label=r"curvature error $\|\kappa-1/R\|_{2}$ (band)")
    axB.set_xlabel(r"time $t$ [s]")
    axB.set_ylabel(r"band curvature error [1/m]", color="#D55E00")
    axB.tick_params(axis="y", labelcolor="#D55E00")
    axB.set_title(r"(b) $N=128$: signed-distance loss drives curvature error")
    axB.grid(True, which="both", ls=":", alpha=0.4)

    axG = axB.twinx()
    axG.plot(t128, g128, "-", color="#0072B2", lw=1.8,
             label=r"$\min|\nabla\psi|$ (band)")
    axG.set_ylabel(r"$\min|\nabla\psi|$ (band)", color="#0072B2")
    axG.tick_params(axis="y", labelcolor="#0072B2")
    axG.set_ylim(0.0, 1.05)

    # Restrict panel (b) to the pre-blow-up window: once the interface is
    # destroyed (t > ~0.055) the band degenerates and both diagnostics become
    # meaningless noise. The clean window already shows the monotone
    # signed-distance loss LEADING the curvature-error runaway.
    axB.set_xlim(0.0, 0.055)
    axG.set_xlim(0.0, 0.055)

    lines = axB.get_lines() + axG.get_lines()
    axB.legend(lines, [ln.get_label() for ln in lines], fontsize=8, loc="upper center")

    fig.tight_layout()
    out = os.path.join(FIGS, "sl_droplet_mechanism.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"[mechanism] wrote {out}")


if __name__ == "__main__":
    main()
