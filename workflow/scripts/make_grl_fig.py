#!/usr/bin/env python3
"""GRL advection studies: per-redistancer convergence figures + order table.

Scans the geometrically-redistanced-levelset theme's data/tables/ for the
curated <study>_errors.csv files that the Snakefile report rule copies there,
and for every study with leiaRedistancedLevelSetFoam rows plots the shape /
volume / band-gradient errors over h, one line per redistancer (largest T,
each cfl). Idempotent and skip-missing: rerunning after any study only adds
or refreshes that study's outputs (same convention as make_sl_3d_fig.py).

    python3 make_grl_fig.py
"""
import csv
import glob
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

THEME = "geometrically-redistanced-levelset"

METRICS = {
    "shapeError":        r"$E_{geom}(\alpha)$ (shape)",
    "volumeError":       r"$E_{vol}(\alpha)$ rel. (volume)",
    "gradientErrorBand": r"band $L_2(||\nabla\psi|-1|)$",
}


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load(err_csv):
    with open(err_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    rows = [r for r in rows if r.get("solver") == "leiaRedistancedLevelSetFoam"
            and _num(r.get("h")) is not None]
    return rows


def studies():
    tdir = paths.tables_dir(THEME, make=False)
    for err_csv in sorted(glob.glob(os.path.join(tdir, "*_errors.csv"))):
        rows = load(err_csv)
        if rows:
            study = os.path.basename(err_csv)[:-len("_errors.csv")]
            yield study, rows


def make_figures(study, rows, figs_dir):
    ts = sorted({_num(r.get("T")) for r in rows if _num(r.get("T"))})
    if not ts:
        return
    tmax = ts[-1]
    cfls = sorted({r.get("cfl", "") for r in rows})
    fig, axes = plt.subplots(1, len(METRICS), figsize=(13.5, 4.2))
    for ax, (col, label) in zip(axes, METRICS.items()):
        drew = False
        for redist in sorted({r.get("redistancer", "") for r in rows}):
            for cfl in cfls:
                pts = sorted(
                    ((_num(r["h"]), _num(r.get(col)))
                     for r in rows
                     if r.get("redistancer") == redist
                     and r.get("cfl", "") == cfl
                     and _num(r.get("T")) == tmax
                     and _num(r.get(col)) is not None),
                )
                if not pts:
                    continue
                hs, es = zip(*pts)
                style = "o-" if cfl == cfls[0] else "s--"
                ax.loglog(hs, es, style,
                          label=f"{redist or 'n/a'} (CFL {cfl})")
                drew = True
        if drew:
            ax.set_xlabel(r"$h$")
            ax.set_title(f"{label},  T = {tmax:g}", fontsize=10)
            ax.grid(True, which="both", alpha=0.25)
            ax.legend(fontsize=6)
    fig.suptitle(f"{study}: leiaRedistancedLevelSetFoam convergence", fontsize=10)
    fig.tight_layout()
    out = os.path.join(figs_dir, f"{study}_convergence.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[make_grl_fig] wrote {out}")


def make_order_table(all_rows, tables_dir):
    lines = [
        r"\begin{tabular}{l l l r r r}",
        r"\hline",
        r"study & redistancer & cfl & $h$ & shape error & order \\",
        r"\hline",
    ]
    for study, rows in all_rows:
        ts = sorted({_num(r.get("T")) for r in rows if _num(r.get("T"))})
        if not ts:
            continue
        tmax = ts[-1]
        combos = sorted({(r.get("redistancer", ""), r.get("cfl", ""))
                         for r in rows})
        for redist, cfl in combos:
            pts = sorted(
                ((_num(r["h"]), _num(r.get("shapeError")))
                 for r in rows
                 if r.get("redistancer") == redist and r.get("cfl") == cfl
                 and _num(r.get("T")) == tmax
                 and _num(r.get("shapeError")) is not None),
                reverse=True,
            )
            prev = None
            for h, err in pts:
                if prev and err > 0 and prev[1] > 0:
                    order = math.log(prev[1]/err)/math.log(prev[0]/h)
                    otxt = f"{order:.2f}"
                else:
                    otxt = "--"
                lines.append(f"{study} & {redist} & {cfl} & {h:.5g} & "
                             f"{err:.2e} & {otxt} \\\\")
                prev = (h, err)
            lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    out = os.path.join(tables_dir, "grl_convergence_orders.tex")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[make_grl_fig] wrote {out}")


def main():
    figs = paths.figs_dir(THEME)
    tables = paths.tables_dir(THEME)
    collected = list(studies())
    if not collected:
        print("[make_grl_fig] no *_errors.csv with leiaRedistancedLevelSetFoam "
              "rows in the theme tables dir yet — nothing to do")
        return
    for study, rows in collected:
        make_figures(study, rows, figs)
    make_order_table(collected, tables)


if __name__ == "__main__":
    main()
