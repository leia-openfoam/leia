#!/usr/bin/env python3
"""Method-decision benchmark report (reversed 2D vortex, all method lines).

Reads every benchVortex*_errors.csv from the sdpls-level-set theme tables dir
(the report rule copies each study's curated table there) and emits, into the
same theme's data/:
  - benchVortex_convergence_T<T>.png : per-method convergence panels --
    shape/volume/band-gradient at t = T and volume/band-gradient at t = T/2
    (maximal stretching: no shape reference there);
  - benchVortex_cost_accuracy.png    : wall clock vs shape error at t = T;
  - benchVortex_decision.tex         : the decision table (finest N per
    method: errors + total wall clock seconds).

    python3 make_method_comparison.py
"""
import csv
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import method_label  # shared method label: precise form + LaTeX escaping
import paths

THEME = "method-comparison"

PANELS = [
    ("shapeError",            r"$E_{geom}(\alpha)$ at $t=T$",              "shape_T"),
    ("volumeError",           r"$E_{vol}$ rel. at $t=T$",                  "volume_T"),
    ("gradientErrorBand",     r"band $L_2(||\nabla\psi|-1|)$ at $t=T$",    "grad_T"),
    ("volumeErrorHalf",       r"$E_{vol}$ rel. at $t=T/2$",               "volume_Thalf"),
    ("gradientErrorBandHalf", r"band $L_2(||\nabla\psi|-1|)$ at $t=T/2$", "grad_Thalf"),
]


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load():
    rows = []
    tdir = paths.tables_dir(THEME, make=False)
    for f in sorted(glob.glob(os.path.join(tdir, "benchVortex*_errors.csv"))):
        rows += [r for r in csv.DictReader(open(f))
                 if r.get("method") and _num(r.get("h")) is not None]
    return rows


# Tuning knobs (+lim, +qr) and the retired redistancing add-on / perturbed-mesh
# robustness study get their OWN slides + figures; they must not clutter the
# primary method comparison. Flux-form (+flux) is NOT a knob -- it is a distinct
# conservative advection scheme (its own RTS slScheme), so it IS a core method
# line and appears in the decision table + convergence panels.
_VARIANT_TAGS = ("+lim", "+qr", "PDEfrozen", "/pert")

# The method label now carries the DISCRETIZATION it ran with (+div:.../+dc:N),
# because two arms that differ only in div(phi,psi) are different methods. Those
# components are constant across a comparison (every Eulerian study runs one
# discretization), so for matching against a known method line we compare the
# BASE label with them stripped -- otherwise every hardcoded literal below would
# silently stop matching and its figure would vanish without an error.
_DISC_RE = re.compile(r"\+(?:div|dc):[^+]*")


def base_label(method):
    return _DISC_RE.sub("", method)


def matches(method, base):
    return base_label(method) == base


def is_core(method):
    return not any(tag in method for tag in _VARIANT_TAGS)


def core_rows(rows):
    """Core method lines only, de-duplicated by (method, T, h) so a shared
    baseline label appearing in several studies is plotted once."""
    seen, out = set(), []
    for r in rows:
        if not is_core(r["method"]):
            continue
        key = (r["method"], r["T"], r["h"])
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def _plot_panel(ax, sel, methods, col, label, guides=True):
    """One convergence panel: every method's error(col) vs h, log-log."""
    hs_all = set()
    for m in methods:
        pts = sorted((_num(r["h"]), _num(r.get(col)))
                     for r in sel if r["method"] == m
                     and _num(r.get(col)) and _num(r.get(col)) > 0)
        if pts:
            hs, es = zip(*pts)
            hs_all.update(hs)
            ax.loglog(hs, es, "o-", label=m)
    if guides and hs_all:
        hs = sorted(hs_all)
        # anchor the O(h)/O(h^2) guides to the finest SL point of this panel
        anchor = None
        for m in methods:
            if not m.startswith("SL"):
                continue
            for r in sel:
                if r["method"] == m and _num(r["h"]) == hs[0] and _num(r.get(col)):
                    anchor = _num(r[col])
        if anchor:
            ax.loglog(hs, [anchor*(h/hs[0]) for h in hs], "k:", lw=0.8,
                      label=r"$O(h)$")
            ax.loglog(hs, [anchor*(h/hs[0])**2 for h in hs], "k--", lw=0.8,
                      label=r"$O(h^2)$")
    ax.set_xlabel(r"$h$")
    ax.set_title(label, fontsize=9)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=6)


def convergence_figs(rows, figs_dir):
    for T in sorted({r["T"] for r in rows}):
        sel = [r for r in rows if r["T"] == T]
        methods = sorted({r["method"] for r in sel})
        fig, axes = plt.subplots(1, len(PANELS), figsize=(4.6*len(PANELS), 4.4))
        for ax, (col, label, _key) in zip(axes, PANELS):
            _plot_panel(ax, sel, methods, col, label)
        fig.suptitle(f"reversed vortex, T = {T}: method comparison", fontsize=11)
        fig.tight_layout()
        out = os.path.join(figs_dir, f"benchVortex_convergence_T{T}.png")
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"[method_comparison] wrote {out}")


def individual_panel_figs(rows, figs_dir):
    """One large standalone figure per (T, metric) for per-panel discussion."""
    for T in sorted({r["T"] for r in rows}):
        sel = [r for r in rows if r["T"] == T]
        methods = sorted({r["method"] for r in sel})
        for col, label, key in PANELS:
            fig, ax = plt.subplots(figsize=(7.0, 5.6))
            _plot_panel(ax, sel, methods, col, label)
            ax.set_title(f"reversed vortex, T = {T}:  {label}", fontsize=12)
            ax.legend(fontsize=8)
            fig.tight_layout()
            out = os.path.join(figs_dir, f"benchVortex_T{T}_{key}.png")
            fig.savefig(out, dpi=200)
            plt.close(fig)
            print(f"[method_comparison] wrote {out}")


def cost_accuracy_fig(rows, figs_dir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, T in zip(axes, sorted({r["T"] for r in rows})):
        sel = [r for r in rows if r["T"] == T]
        for m in sorted({r["method"] for r in sel}):
            pts = [(_num(r["totalClockTime"]), _num(r["shapeError"]),
                    round(1/_num(r["h"])))
                   for r in sel if r["method"] == m
                   and _num(r.get("totalClockTime")) and _num(r.get("shapeError"))]
            if not pts:
                continue
            cs, es, ns = zip(*sorted(pts))
            ax.loglog(cs, es, "o-", label=m)
            for c, e, n in pts:
                ax.annotate(str(n), (c, e), fontsize=6,
                            textcoords="offset points", xytext=(3, 3))
        ax.set_xlabel("total wall clock [s]")
        ax.set_ylabel(r"$E_{geom}(\alpha)$ at $t=T$")
        ax.set_title(f"T = {T}", fontsize=10)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=6)
    fig.suptitle("cost vs accuracy (labels: N)", fontsize=11)
    fig.tight_layout()
    out = os.path.join(figs_dir, "benchVortex_cost_accuracy.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[method_comparison] wrote {out}")


def decision_table(rows, tables_dir):
    lines = [r"\begin{tabular}{l l r r r r r r}", r"\hline",
             r"method & $T$ & $N$ & $E_{geom}$ & $E_{vol}$ & "
             r"band grad & $E_{vol}^{T/2}$ & clock [s] \\", r"\hline"]
    for T in sorted({r["T"] for r in rows}):
        sel = [r for r in rows if r["T"] == T]
        for m in sorted({r["method"] for r in sel}):
            mr = [r for r in sel if r["method"] == m]
            finest = min(mr, key=lambda r: _num(r["h"]))
            def g(c):
                v = _num(finest.get(c))
                return f"{v:.2e}" if v is not None else "--"
            lines.append(
                f"{method_label.latex_escape(m)} & {T} & {round(1/_num(finest['h']))} & "
                f"{g('shapeError')} & {g('volumeError')} & "
                f"{g('gradientErrorBand')} & {g('volumeErrorHalf')} & "
                f"{g('totalClockTime')} \\\\")
        lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    out = os.path.join(tables_dir, "benchVortex_decision.tex")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[method_comparison] wrote {out}")


def _series(rows, method, col, T=None):
    pts = sorted((_num(r["h"]), _num(r.get(col)))
                 for r in rows if r["method"] == method
                 and (T is None or r["T"] == T)
                 and _num(r.get(col)) and _num(r.get(col)) > 0)
    return zip(*pts) if pts else ((), ())


def improvements_fig(rows, figs_dir):
    """SL-quadratic improvements, measured: baseline vs Venkatakrishnan limiter
    vs flux-form. Left = shape error (order), right = volume error, both T = 8."""
    base = "SL:uncachedQuadraticWeightedLeastSquares"
    variants = [
        (base,           "baseline (point-value)", "o-"),
        (base + "+lim:venk", "+ Venkatakrishnan",  "s--"),
        (base + "+flux",     "+ flux-form",         "^--"),
    ]
    if not any(any(matches(r["method"], m) for r in rows) for m, _, _ in variants[1:]):
        print("[method_comparison] improvements_fig: no SL variant rows matched "
              f"{[m for m, _, _ in variants[1:]]} -- skipping")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (col, lab) in zip(
        axes,
        [("shapeError", r"$E_{geom}(\alpha)$ at $t=T$"),
         ("volumeError", r"$E_{vol}$ rel. at $t=T$")],
    ):
        for m, lg, st in variants:
            hs, es = _series(rows, m, col, T="8")
            if hs:
                ax.loglog(hs, es, st, label=lg)
        ax.set_xlabel(r"$h$"); ax.set_title(lab, fontsize=10)
        ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=7)
    fig.suptitle("SL-quadratic improvements, reversed vortex T = 8 "
                 "(the hardest case; flux-form is competitive at T = 2)", fontsize=11)
    fig.tight_layout()
    out = os.path.join(figs_dir, "benchVortex_slimprovements.png")
    fig.savefig(out, dpi=200); plt.close(fig)
    print(f"[method_comparison] wrote {out}")


def frozen_fig(rows, figs_dir):
    """Frozen-band bulk-only redistancing vs no redistancing on the Eulerian
    line, T = 8: volume error vs N (the redistancer injures the transport)."""
    pairs = [("euler", "no redistancing", "o-"),
             ("euler+RD:PDEfrozen", "frozen-band redistancing", "s--")]
    if not any(matches(r["method"], "euler+RD:PDEfrozen") for r in rows):
        print("[method_comparison] frozen_fig: no euler+RD:PDEfrozen rows "
              "-- skipping")
        return
    # Restrict to the N range the frozen study actually covers.
    Nfrozen = {round(1/_num(r["h"])) for r in rows
               if matches(r["method"], "euler+RD:PDEfrozen") and r["T"] == "8"}
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for m, lg, st in pairs:
        pts = sorted((round(1/_num(r["h"])), _num(r.get("volumeError")))
                     for r in rows if matches(r["method"], m) and r["T"] == "8"
                     and round(1/_num(r["h"])) in Nfrozen
                     and _num(r.get("volumeError")) is not None)
        if pts:
            ns, es = zip(*pts)
            ax.semilogy(ns, es, st, label=lg)
    ax.set_xlabel(r"$N$ (cells per side)")
    ax.set_ylabel(r"$E_{vol}$ rel. at $t=T$")
    ax.set_title("Frozen-band redistancing INJURES the Eulerian line "
                 "(T = 8):\nrebuilds distance to advection's spurious bulk "
                 "sign-changes", fontsize=10)
    ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(figs_dir, "benchVortex_frozen.png")
    fig.savefig(out, dpi=200); plt.close(fig)
    print(f"[method_comparison] wrote {out}")


def _order(rows, method, T, col):
    """Least-squares convergence order (slope of log e vs log 1/h) for a
    method/horizon/metric; None if < 2 usable points."""
    import math
    pts = [(_num(r["h"]), _num(r.get(col))) for r in rows
           if r["method"] == method and r["T"] == T
           and _num(r.get(col)) and _num(r.get(col)) > 0]
    pts = sorted(set(pts))
    if len(pts) < 2:
        return None
    xs = [math.log(1.0/h) for h, _ in pts]
    ys = [math.log(e) for _, e in pts]
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    d = n*sxx - sx*sx
    return (n*sxy - sx*sy)/d if abs(d) > 1e-30 else None


# Human-readable method labels for the detailed table (long RTS names -> short).
def _short(m):
    m = m.replace("SL:uncachedQuadraticWeightedLeastSquares", "SL quad (uncached)")
    m = m.replace("SL:quadraticWeightedLeastSquares", "SL quad (cached)")
    m = m.replace("+flux", " + flux-form")
    m = m.replace("+lim:venk", " + venk")
    m = m.replace("+qr", " + QR")
    m = m.replace("/pert", " [pert mesh]")
    m = m.replace("euler+VE:closestPoint", "euler + VE:closestPoint")
    m = m.replace("euler+SDPLS:R", "euler + SDPLS:R")
    m = m.replace("euler+SDPLS:beta", "euler + SDPLS:beta")
    m = m.replace("euler+RD:PDEfrozen", "euler + frozen-redist")
    # Discretization components: keep them, but compress. They are constant
    # across a comparison, so the long form only costs column width; the
    # unabbreviated string is the `method` column of the curated CSV.
    m = m.replace("+div:linearUpwind/grad(psi)", " [lU/grad(psi)]")
    m = re.sub(r"\+div:([^+\s]+)", r" [\1]", m)
    m = re.sub(r"\+dc:(\d+)", r" dc\1", m)
    return m


def _detailed_records(rows):
    """One record per (method, T) at its finest N: metrics + convergence orders,
    ordered core methods first then variants. ALL methods included."""
    recs = []
    for T in sorted({r["T"] for r in rows}):
        sel = [r for r in rows if r["T"] == T]
        methods = sorted({r["method"] for r in sel},
                         key=lambda m: (not is_core(m), m))
        for m in methods:
            mr = [r for r in sel if r["method"] == m]
            finest = min(mr, key=lambda r: _num(r["h"]))
            recs.append({
                "method": m, "T": T,
                "N": round(1/_num(finest["h"])),
                "E_geom": _num(finest.get("shapeError")),
                "p_geom": _order(rows, m, T, "shapeError"),
                "E_vol": _num(finest.get("volumeError")),
                "p_vol": _order(rows, m, T, "volumeError"),
                "band": _num(finest.get("gradientErrorBand")),
                "clock": _num(finest.get("totalClockTime")),
            })
    return recs


def detailed_table(rows, tables_dir):
    """Comprehensive LaTeX table: EVERY method x horizon, finest N, all metrics
    + convergence orders (article \\input)."""
    recs = _detailed_records(rows)
    def f(v, e=True):
        if v is None:
            return "--"
        return (f"{v:.2e}" if e else f"{v:.2f}")
    lines = [r"\begin{tabular}{l c c c c c c c}", r"\hline",
             r"method & $T$ & $N$ & $E_{geom}$ & $p_{geom}$ & $E_{vol}$ & "
             r"$p_{vol}$ & clock [s] \\", r"\hline"]
    for r in recs:
        lines.append(
            f"{method_label.latex_escape(_short(r['method']))} & {r['T']} & {r['N']} & "
            f"{f(r['E_geom'])} & {f(r['p_geom'], e=False)} & "
            f"{f(r['E_vol'])} & {f(r['p_vol'], e=False)} & "
            f"{f(r['clock'], e=False)} \\\\")
    lines += [r"\hline", r"\end{tabular}"]
    out = os.path.join(tables_dir, "benchVortex_detailed.tex")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[method_comparison] wrote {out}")


def detailed_fig(rows, figs_dir):
    """Same comprehensive table rendered as a PNG for the deck (all methods)."""
    recs = _detailed_records(rows)
    def f(v, e=True):
        if v is None:
            return "--"
        return (f"{v:.1e}" if e else f"{v:.2f}")
    header = ["method", "T", "N", "E_geom", "p", "E_vol", "p", "band|∇ψ|", "clock s"]
    cells, colours = [], []
    for r in recs:
        core = is_core(r["method"])
        cells.append([_short(r["method"]), r["T"], r["N"], f(r["E_geom"]),
                      f(r["p_geom"], e=False), f(r["E_vol"]),
                      f(r["p_vol"], e=False), f(r["band"]), f(r["clock"], e=False)])
        colours.append("#eef6ff" if core else "#fff6e8")
    fig, ax = plt.subplots(figsize=(12, 0.34*len(cells) + 0.8))
    ax.axis("off")
    tab = ax.table(cellText=cells, colLabels=header, loc="center", cellLoc="center")
    tab.auto_set_font_size(False); tab.set_fontsize(8)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#333333"); cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor(colours[r-1])
            if c == 0:
                cell.get_text().set_ha("left")
    ax.set_title("Detailed method comparison (finest N per method; "
                 "blue = core method line, tan = variant)", fontsize=11, pad=12)
    fig.tight_layout()
    out = os.path.join(figs_dir, "benchVortex_detailed.png")
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"[method_comparison] wrote {out}")


def main():
    rows = load()
    if not rows:
        print("[method_comparison] no benchVortex*_errors.csv rows yet")
        return
    figs = paths.figs_dir(THEME)
    tables = paths.tables_dir(THEME)
    # Primary comparison: the core method lines only (variants get own figures).
    core = core_rows(rows)
    convergence_figs(core, figs)
    individual_panel_figs(core, figs)
    cost_accuracy_fig(core, figs)
    decision_table(core, tables)
    # Dedicated slides for the measured improvement / redistancing results.
    improvements_fig(rows, figs)
    frozen_fig(rows, figs)
    # Comprehensive all-methods detail (table + PNG) with convergence orders.
    detailed_table(rows, tables)
    detailed_fig(rows, figs)


if __name__ == "__main__":
    main()
