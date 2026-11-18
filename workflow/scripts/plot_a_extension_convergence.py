#!/usr/bin/env python3
"""plot_a_extension_convergence.py: step A round 2, the static extension test (STATUS 11.39).

    plot_a_extension_convergence.py --study <dir with aStatic2Dr2> --out <png>

Three panels (one per input field psi with the same zero set: the distance, the scaled distance with
|grad psi| from 0.8 to 1.25, and r^2 - R^2 with |grad psi| about 0.3): the L1 error of the computed
extension velocity against u(Y), Y the Newton point of the exact psi, over the four cell layers that the
curvature reads (layers 0 to 3), against N. Two series: the Newton-point extension (projectionBlend,
distance first) and the nearest-point extension (closestPoint, newton step, layer band). A slope-2 guide.
L1 only (never L_inf). Stdlib + matplotlib.
"""
import argparse, csv, glob, json, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"    # validated categorical slots 1 and 2 (light surface)
INK, MUTED, GRID, SURF = "#0b0b0b", "#5f5e5a", "#e4e2dd", "#fcfcfb"


def load(study):
    d = {}
    for c in glob.glob(os.path.join(study, "2Dvortex_*")):
        t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        f = os.path.join(c, "leiaTestVelocityExtensionLayers.csv")
        if not (os.path.isfile(f) and os.path.getsize(f)):
            continue
        m = t["VELOCITY_EXTENSION"]
        if m == "projectionBlend" and str(t["PB_DISTANCE"]) != "first":
            continue
        r = list(csv.DictReader(open(f)))[0]
        d.setdefault((m, t["VE_TEST_INPUT"]), {})[int(t["N_CELLS"])] = float(r["E_EXTY_L1_L0TO3"])
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = load(a.study)
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True, facecolor=SURF)
    titles = {"distance": "psi = signed distance (|grad psi| = 1)",
              "scaledDistance": "psi = scaled distance (|grad psi| 0.8 to 1.25)",
              "quadratic": "psi = r^2 - R^2 (|grad psi| about 0.3)"}
    for ax, inp in zip(axs, ("distance", "scaledDistance", "quadratic")):
        ax.set_facecolor(SURF)
        for m, col, lab in (("projectionBlend", BLUE, "Newton point"), ("closestPoint", ORANGE, "nearest point")):
            v = sorted(d.get((m, inp), {}).items())
            if not v:
                continue
            ax.plot([n for n, _ in v], [e for _, e in v], color=col, lw=2, marker="o", ms=5)
            ax.annotate(lab, (v[-1][0], v[-1][1]), xytext=(6, 0), textcoords="offset points",
                        color=INK, fontsize=8, va="center")
        # slope-2 guide
        x0, y0 = 128, 2.0e-2
        ax.plot([x0, 4 * x0], [y0, y0 / 16], color=MUTED, lw=1, ls="--")
        ax.annotate("order 2", (2 * x0, y0 / 4), xytext=(4, 4), textcoords="offset points", color=MUTED, fontsize=8)
        ax.set_xscale("log", base=2); ax.set_yscale("log")
        ax.set_xticks([64, 128, 256, 512]); ax.set_xticklabels(["64", "128", "256", "512"])
        ax.set_xlim(50, 1100)
        ax.set_title(titles[inp], fontsize=9, color=INK, loc="left")
        ax.grid(True, which="major", color=GRID, lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_xlabel("cells per unit length N", color=MUTED, fontsize=8)
    axs[0].set_ylabel("L1 error of u_e against u(Y), layers 0 to 3", color=MUTED, fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=SURF)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
