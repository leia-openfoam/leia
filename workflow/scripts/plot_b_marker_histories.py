#!/usr/bin/env python3
"""plot_b_marker_histories.py: the q error on the markers against time, step B (STATUS 11.40).

    plot_b_marker_histories.py --h <clone of round 1> --h2 <clone of the controls> --out <png>

Four panels, log scale, one sequential blue ramp per panel (ordered by N or by the band width):
(a) round 1, projectionBlend w 4/8 cells, E_Q_L1; (b) the raw transport against its own exact q, E_QRAW_L1;
(c) the band width at N 256, E_Q_L1; (d) the blend at a fixed physical width (1/16 to 1/8), E_Q_L1.
L1 means over 720 markers (never L_inf). Pure strain a = 1, T = 0.5.
"""
import argparse, csv, glob, json, math, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RAMP = ["#86b6ef", "#2a78d6", "#104281"]  # validated ordinal ramp (light surface)
INK, MUTED, GRID, SURF = "#0b0b0b", "#5f5e5a", "#e4e2dd", "#fcfcfb"


def load(study_dir, **match):
    for c in sorted(glob.glob(os.path.join(study_dir, "markerTest2D_*"))):
        t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
        if all(str(t.get(k)) == str(v) for k, v in match.items()):
            f = os.path.join(c, "interfaceMarkers.csv")
            r = list(csv.DictReader(open(f)))
            return r
    return []


def panel(ax, series, col, title):
    ends = []
    for i, (lab, r) in enumerate(series):
        if not r:
            continue
        t = [float(x["TIME"]) for x in r]
        e = [float(x[col]) for x in r]
        ax.plot(t, e, color=RAMP[i], lw=2)
        ends.append((e[-1], t[-1], lab))
    # direct labels at the line ends, at least 11 points apart (positions in log10 of the error;
    # the axis spans log10(2e-1/3e-5) = 3.8 decades over about 150 points)
    ends.sort()
    pts_per_dec = 150 / 3.8
    pos = [math.log10(e) * pts_per_dec for e, _, _ in ends]
    for k in range(1, len(pos)):
        pos[k] = max(pos[k], pos[k - 1] + 11)
    for (e, t, lab), y in zip(ends, pos):
        ax.annotate(lab, (t, e), xytext=(4, y - math.log10(e) * pts_per_dec), textcoords="offset points",
                    color=INK, fontsize=8, va="center")
    ax.set_yscale("log")
    ax.set_xlim(0, 0.64)
    ax.set_ylim(3e-5, 2e-1)
    ax.set_title(title, fontsize=9, color=INK, loc="left")
    ax.grid(True, which="major", color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel("t", color=MUTED, fontsize=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", required=True); ap.add_argument("--h2", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S = lambda root, st: os.path.join(root, st)
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.4), sharey=True, facecolor=SURF)
    for ax in axs.flat:
        ax.set_facecolor(SURF)
    panel(axs[0, 0], [(f"N {N}", load(S(a.h, "bMarker2Dstrain"), VELOCITY_EXTENSION="projectionBlend", N_CELLS=N))
                      for N in (128, 256, 512)], "E_Q_L1",
          "(a) projectionBlend, band 4/8 cells: error grows under refinement")
    panel(axs[0, 1], [(f"N {N}", load(S(a.h2, "bRaw2Dstrain"), N_CELLS=N)) for N in (128, 256, 512)], "E_QRAW_L1",
          "(b) no extension, against its own exact q: order 2")
    panel(axs[1, 0], [(f"{w}/{2 * w}", load(S(a.h2, f"bBand2DstrainPB{w}"))) for w in (2, 4, 8)], "E_Q_L1",
          "(c) N 256, blend width in cells: wider band, later onset")
    panel(axs[1, 1], [(f"N {N}", load(S(a.h2, f"bEdge2DstrainN{N}"))) for N in (128, 256, 512)], "E_Q_L1",
          "(d) blend at a fixed physical width: order 1 only")
    axs[0, 0].set_ylabel("L1 mean of |q_h - q| on the markers", color=MUTED, fontsize=8)
    axs[1, 0].set_ylabel("L1 mean of |q_h - q| on the markers", color=MUTED, fontsize=8)
    fig.suptitle("Step B marker test, pure strain a = 1: the q error on 720 material markers", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=SURF)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
