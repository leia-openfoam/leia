#!/usr/bin/env python3
"""The resolution ladders of the transport studies in one figure (fig_ladders.png).

Reads the exported order tables of the leia workflow (data/tables/*_orders_Co0.2.tex, never
edited by hand) and plots, against the cell size h = 1/N, the face-fraction error and the normal
error at the end of the run for the plane-clipped loop (ring normal, exact-normal control,
plicRDF) and for the dual advection with the ring-tensor paraboloid positioned by the volume
fraction (exact normal, ring quadric's normal, ring's mean normal), N = 32 to 256, Co = 0.2.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = Path(__file__).resolve().parent
TAB = HERE / "data" / "tables"

STUDIES = {
    "plane": TAB / "davof_advection_sphereTranslation3D_hex_orders_Co0.2.tex",
    "parab": TAB / "davof_advection_sphereTranslation3D_hexParaboloid_orders_Co0.2.tex",
    "ringfit": TAB / "davof_advection_sphereTranslation3D_hexParaboloidRingFit_orders_Co0.2.tex",
    "n256": TAB / "davof_advection_sphereTranslation3D_hexParaboloidN256_orders_Co0.2.tex",
}

METRIC = {"alphaf": r"$L_1(\alpha_f)$, end", "n": r"$L_1(\mathbf{n})$, end"}


def read_rows(path):
    """Return {(row, metric): [values...]} from an exported orders table."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\\texttt\{([^}]*)\}\s*&\s*(.*?)\s*&\s*(.*)\\\\", line)
        if not m:
            continue
        row, metric, rest = m.group(1), m.group(2).strip(), m.group(3)
        cells = [c.strip() for c in rest.split("&")]
        vals = []
        for c in cells:
            try:
                vals.append(float(c))
            except ValueError:
                break
        out[(row, metric)] = vals
    return out


T = {k: read_rows(v) for k, v in STUDIES.items()}
N3 = [32, 64, 128]


def series(study, row, metric, with256=None):
    vals = list(T[study][(row, METRIC[metric])][:3])
    ns = list(N3)
    if with256 is not None:
        vals.append(T["n256"][(with256, METRIC[metric])][0])
        ns.append(256)
    return [1.0 / n for n in ns], vals


SERIES = [
    # label, study, row, row at N=256, style
    ("plane-clipped loop, ring normal", "plane", "davof-ringArea", None, dict(color="#d08a5a", marker="o", ls="--")),
    ("plane-clipped loop, exact normal (control)", "plane", "davof-exact", None, dict(color="#5a7fa0", marker="s", ls="--")),
    ("plicRDF (reconstruction only)", "plane", "plicRDF", None, dict(color="#777777", marker="^", ls=":")),
    ("paraboloid by volume, ring mean normal", "parab", "davof-paraboloidRing-ringArea", "davof-paraboloidRing-ringArea", dict(color="#e0864a", marker="D", ls="-")),
    ("paraboloid by volume, ring quadric normal", "ringfit", "davof-paraboloidRing-ringFit", "davof-paraboloidRing-ringFit", dict(color="#7a3fb0", marker="p", ls="-")),
    ("paraboloid by volume, exact normal (control)", "parab", "davof-paraboloidRing-exact", "davof-paraboloidRing-exact", dict(color="#1b7d7a", marker="x", ls="-")),
]

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
for ax, metric, ylabel in [(axes[0], "alphaf", r"$L_1(\alpha_f)$ at the end of the run"),
                           (axes[1], "n", r"$L_1(\mathbf{n})$ at the end of the run")]:
    for label, study, row, row256, style in SERIES:
        key = (row, METRIC[metric])
        if key not in T[study] or not T[study][key]:
            continue
        h, v = series(study, row, metric, row256)
        ax.loglog(h, v, label=label, lw=1.6, ms=6, **style)
    h0, h1 = 1 / 256, 1 / 32
    for p, ls, lab in [(1, ":", r"$h^1$"), (2, "-.", r"$h^2$"), (3, "--", r"$h^3$")]:
        ref = 2e-5 if metric == "alphaf" else 5e-5
        ax.loglog([h0, h1], [ref * (h0 / h1) ** p, ref], color="0.6", ls=ls, lw=1, label=lab)
    ax.set_xlabel(r"$h = 1/N$")
    ax.set_ylabel(ylabel)
    ax.set_xticks([1 / 256, 1 / 128, 1 / 64, 1 / 32])
    ax.set_xticklabels(["1/256", "1/128", "1/64", "1/32"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(1 / 256 / 1.15, 1 / 32 * 1.15)
    ax.grid(True, which="major", alpha=0.3)
axes[0].set_title("face fractions", fontsize=10)
axes[1].set_title("recovered normal", fontsize=10)
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.14, 1, 1))
out = HERE / "fig_ladders.png"
fig.savefig(out, dpi=220)
print("wrote", out)
