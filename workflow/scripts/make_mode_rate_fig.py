#!/usr/bin/env python3
"""Figure: per-mode capillary growth rate vs psi-profile drift (N=128).

Reads docs/method-comparison/.../tables/mode_rate_vs_drift.csv (the verbatim
output of workflow/scripts/mode_rate_vs_drift.py + interface_mode_trajectory.py).

Reads the aggregate CSV from mode_rate_vs_drift.py. A point is VALID when the
trajectory fit is tight (rel residual <= 0.15) and the measured frequency is
physical (omega within 25% of the capillary dispersion relation); invalid
points are shown hollow and never connected -- they are reported, not used.
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURF, INK, INK2, MUTED = "#fcfcfb", "#22201d", "#55514b", "#8a857d"
COLORS = {2: "#2a78d6", 4: "#eb6834", 6: "#1baf7a", 8: "#eda100"}

src, out = sys.argv[1], sys.argv[2]
rows = []
for r in csv.DictReader(open(src)):
    if not r.get("mode") or r.get("N") != "128":
        continue
    d = {}
    for k, v in r.items():
        if k in ("coupling_top", "note", "snapshot_t"):
            continue
        try:
            d[k] = float(v)
        except (TypeError, ValueError):
            pass
    rows.append(d)

fig, ax = plt.subplots(figsize=(8.6, 5.6), facecolor=SURF)
ax.set_facecolor(SURF)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.8)
ax.tick_params(colors=INK2, labelsize=9)

ax.axhline(0.0, color=MUTED, lw=1.1)
modes = sorted({int(r["mode"]) for r in rows})
for m in modes:
    pts = sorted([r for r in rows if int(r["mode"]) == m],
                 key=lambda r: -r["minGradPsiBand"])
    drift = [1.0 - r["minGradPsiBand"] for r in pts]
    rm = [r["r_m"] for r in pts]
    ok = [(r["fit_rel_residual"] <= 0.15 and
           0.75 <= r["omega_ratio"] <= 1.25) for r in pts]
    c = COLORS.get(m, INK2)
    dv = [d for d, o in zip(drift, ok) if o]
    rv = [v for v, o in zip(rm, ok) if o]
    ax.plot(dv, rv, "-", color=c, lw=1.8, marker="o", ms=6,
            label=f"m = {m}")
    di = [d for d, o in zip(drift, ok) if not o]
    ri = [v for v, o in zip(rm, ok) if not o]
    if di:
        ax.plot(di, ri, "o", mfc="none", mec=c, ms=6, mew=1.4)

ax.set_xlabel(r"profile drift   $1-\min|\nabla\psi|_{\rm band}$",
              color=INK2, fontsize=10.5)
ax.set_ylabel(r"mode growth rate  $r_m$   [1/s]", color=INK2, fontsize=10.5)
ax.grid(color=MUTED, alpha=0.2, lw=0.7)
ax.set_axisbelow(True)
leg = ax.legend(frameon=False, fontsize=10, loc="upper left",
                labelcolor=INK2)
ax.set_title("Capillary-mode growth rate vs $\\psi$-profile drift — "
             "N=128, U and p reset, production delivery",
             color=INK, fontsize=11.5, loc="left", pad=10)
fig.text(0.01, 0.012,
         "hollow = fit invalid (residual > 0.15 or omega off the capillary "
         "dispersion by > 25%) — reported, not used",
         color=MUTED, fontsize=8.5)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(out, dpi=170, facecolor=SURF)
print("wrote", out)
