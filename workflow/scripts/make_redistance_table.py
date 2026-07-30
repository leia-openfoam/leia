#!/usr/bin/env python3
"""Static redistancing gate (circle): convergence figure + order table.

Reads <study>_database.csv of a leiaTestRedistance study (e.g.
config/redistanceCircle2D.yaml), plots the post-event band error norms per
redistancer over h (log-log, with O(h)/O(h^2) guides), and writes a LaTeX
order table. Outputs land in the geometrically-redistanced-levelset theme's
single data source (figures/ + tables/), consumed by BOTH the deck and the
article.

    python3 make_redistance_table.py <STUDY_DIR>
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

THEME = "geometrically-redistanced-levelset"

ERR_COLS = {
    "E_LINF_BAND_PSI": r"band $L_\infty(\psi - \psi_{exact})$",
    "E_MEAN_BAND_GRAD_PSI": r"band mean $||\nabla\psi| - 1|$",
    "E_VOL_ALPHA": r"$E_{vol}$: $|\Delta\sum\alpha V|$ (one event)",
    "E_GEOM_ALPHA": r"$E_{geom}$: $\sum V|\alpha-\alpha_0|$ (one event)",
}


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load(study_dir):
    study = os.path.basename(os.path.normpath(study_dir))
    db = os.path.join(study_dir, f"{study}_database.csv")
    with open(db, newline="") as fh:
        rows = list(csv.DictReader(fh))
    data = {}   # (surface, redistancer) -> list of (h, {col: err})
    for r in rows:
        redist = r.get("REDISTANCER", "")
        surface = r.get("SURFACE", "")   # empty for single-surface studies
        h = _num(r.get("leiaTestRedistance.H"))
        if not redist or h is None:
            continue
        errs = {c: _num(r.get(f"leiaTestRedistance.{c}")) for c in ERR_COLS}
        pre = _num(r.get("leiaTestRedistance.E_LINF_BAND_PSI_PRE"))
        data.setdefault((surface, redist), []).append((h, errs, pre))
    for key in data:
        data[key].sort(key=lambda t: t[0])
    return study, data


def _fig_name(study, surface):
    if study == "redistanceCircle2D":          # legacy names (article refs)
        return "static_redistance_convergence.png"
    sfx = f"_{surface}" if surface else ""
    return f"{study}_convergence{sfx}.png"


def make_figure(study, data, figs_dir):
    surfaces = sorted({s for s, _ in data})
    for surface in surfaces:
        cols = [c for c in ERR_COLS
                if any(e.get(c) is not None
                       for (s, _), rows in data.items() if s == surface
                       for _, e, _ in rows)]
        fig, axes = plt.subplots(1, len(cols), figsize=(4.6*len(cols), 4.0))
        if len(cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, cols):
            hs_all = set()
            for (s, redist), rows in sorted(data.items()):
                if s != surface:
                    continue
                pts = [(h, e[col]) for h, e, _ in rows
                       if e.get(col) is not None and e[col] > 0]
                if not pts:
                    continue
                hs, es = zip(*pts)
                hs_all.update(hs)
                ax.loglog(hs, es, "o-", label=redist)
            if hs_all:
                hs = sorted(hs_all)
                ref = data.get((surface, "planeFootWave"), [])
                e0 = next((e[col] for h, e, _ in ref
                           if h == hs[-1] and e.get(col)), None)
                if e0:
                    ax.loglog(hs, [e0*(h/hs[-1]) for h in hs], "k:",
                              lw=0.8, label=r"$O(h)$")
                    ax.loglog(hs, [e0*(h/hs[-1])**2 for h in hs], "k--",
                              lw=0.8, label=r"$O(h^2)$")
            ax.set_xlabel(r"$h$")
            ax.set_title(ERR_COLS[col], fontsize=10)
            ax.grid(True, which="both", alpha=0.25)
            ax.legend(fontsize=7)
        sfx = f" [{surface}]" if surface else ""
        fig.suptitle(f"static gate {study}{sfx}: one event, deviation vs input",
                     fontsize=10)
        fig.tight_layout()
        out = os.path.join(figs_dir, _fig_name(study, surface))
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"[make_redistance_table] wrote {out}")


def make_table(study, data, tables_dir):
    cols = ["E_LINF_BAND_PSI", "E_VOL_ALPHA", "E_GEOM_ALPHA"]
    lines = [
        r"\begin{tabular}{l l r r r r r}",
        r"\hline",
        r"surface & redistancer & $h$ & band $L_\infty(\psi)$ & order & "
        r"$E_{vol}$ & $E_{geom}$ \\",
        r"\hline",
    ]
    for (surface, redist), rows in sorted(data.items()):
        prev = None
        for h, e, pre in rows:
            err = e.get(cols[0])
            if err is None:
                continue
            if prev and prev[1] > 0 and err > 0:
                order = math.log(prev[1]/err)/math.log(prev[0]/h)
                otxt = f"{order:.2f}"
            else:
                otxt = "--"
            def fmt(c):
                v = e.get(c)
                return f"{v:.2e}" if v is not None else "--"
            lines.append(
                f"{surface or '--'} & {redist} & {h:.5g} & {err:.2e} & {otxt}"
                f" & {fmt('E_VOL_ALPHA')} & {fmt('E_GEOM_ALPHA')} \\\\")
            prev = (h, err)
        lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    name = ("static_redistance_orders.tex" if study == "redistanceCircle2D"
            else f"{study}_orders.tex")
    out = os.path.join(tables_dir, name)
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[make_redistance_table] wrote {out}")


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <STUDY_DIR>")
    study, data = load(sys.argv[1])
    if not data:
        sys.exit("[make_redistance_table] no leiaTestRedistance rows found")
    make_figure(study, data, paths.figs_dir(THEME))
    make_table(study, data, paths.tables_dir(THEME))


if __name__ == "__main__":
    main()
