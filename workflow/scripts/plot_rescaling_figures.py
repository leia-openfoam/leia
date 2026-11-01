#!/usr/bin/env python3
"""Figures of the rescaling studies for the gradient-control deck (STATUS 11.28).

    plot_rescaling_figures.py --fields <metrics>/staticSource_fields.csv --out <figures dir>

Reads the per-write field metrics of static_source_metrics.py (staticSourceRescale*2D) and writes
  rescale_once_convergence.png  dist3_2 after one event against N, three maps, with slope guides
  rescale_every_history.png     dist3_2 against tau, one event per step, start d, N = 256
Dark deck surface; categorical slots 1 to 3 of the dataviz reference palette (dark steps)."""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE, TEXT, MUTED, GRID = "#191919", "#f0efec", "#c3c2b7", "#3a3a38"
SLOT = ["#3987e5", "#d95926", "#199e70"]


def style(ax, fig):
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=12)
    ax.grid(True, which="major", color=GRID, linewidth=0.8)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.title.set_color(TEXT)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fields", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = list(csv.DictReader(open(a.fields)))
    val = {(r["arm"], r["init"], int(r["N"]), float(r["tau"])): float(r["dist3_2"]) for r in rows}
    Ns = [64, 128, 256, 512]

    # One event: the order of each map
    once = [("centred, second order (0b), start blend 0.1", "noSource-none+rescale-second-once-centred-SdplsLSQ", "blend01"),
            ("centred, first order (0a), start $r^2-R^2$", "noSource-none+rescale-first-once-centred-SdplsLSQ", "alg"),
            ("one-sided, first order (0a), start $d$", "noSource-none+rescale-first-once-godunov-SdplsLSQ", "sdf")]
    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=150)
    style(ax, fig)
    for (lab, arm, init), c in zip(once, SLOT):
        ys = [val.get((arm, init, N, 1.0), float("nan")) for N in Ns]
        ax.loglog(Ns, ys, "-o", color=c, lw=2, ms=8, label=lab)
    for p, y0 in ((1, 3e-2), (2, 2e-3)):
        ax.loglog([64, 512], [y0, y0 * 8.0 ** -p], "--", color=MUTED, lw=1)
        ax.text(560, y0 * 8.0 ** -p, f"order {p}", color=MUTED, fontsize=11, va="center")
    ax.set_xticks(Ns)
    ax.set_xticklabels([str(N) for N in Ns])
    ax.set_xlim(55, 700)
    ax.set_xlabel("N (cells per side)", fontsize=13)
    ax.set_ylabel(r"$L_2$ of $(\psi-d)/h$ in $|d|\leq 3h$ after one event", fontsize=13)
    leg = ax.legend(frameon=False, fontsize=11, loc="lower left")
    for t in leg.get_texts():
        t.set_color(TEXT)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "rescale_once_convergence.png"), facecolor=SURFACE)
    plt.close(fig)

    # Every step: the history from d
    every = [("centred, second order", "noSource-none+rescale-second-every-centred-SdplsLSQ"),
             ("one-sided, first order", "noSource-none+rescale-first-every-godunov-SdplsLSQ"),
             ("one-sided, second order", "noSource-none+rescale-second-every-godunov-SdplsLSQ")]
    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=150)
    style(ax, fig)
    for (lab, arm), c in zip(every, SLOT):
        pts = sorted((t, v) for (ar, i, N, t), v in val.items() if ar == arm and i == "sdf" and N == 256)
        ax.semilogy([t for t, _ in pts], [v for _, v in pts], "-o", color=c, lw=2, ms=8, label=lab)
    ax.axhline(1.0, color=MUTED, lw=1, ls="--")
    ax.text(10.1, 1.0, "1 h", color=MUTED, fontsize=11, va="center")
    ax.set_xlabel(r"pseudo-time $\tau$ (one event per step, $\Delta\tau = h/2$)", fontsize=13)
    ax.set_ylabel(r"$L_2$ of $(\psi-d)/h$ in $|d|\leq 3h$", fontsize=13)
    ax.set_title("start: the signed distance d, N = 256", fontsize=13)
    leg = ax.legend(frameon=False, fontsize=11, loc="center right")
    for t in leg.get_texts():
        t.set_color(TEXT)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "rescale_every_history.png"), facecolor=SURFACE)
    plt.close(fig)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
