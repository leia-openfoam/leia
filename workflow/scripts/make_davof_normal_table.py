#!/usr/bin/env python3
"""DAVOF static gate: order table + log-log figures of the Gauss-identity normal.

Reads, for every case of a leiaTestDavofNormal study,

    <case>/case_params.json               tokens.N_CELLS, DOMAIN_LENGTH, DAVOF_ALPHA_SOURCE
    <case>/leiaTestDavofNormalModels.csv  tidy: one row per MODEL (davof, plicRDF, ...)
    <case>/leiaTestDavofNormal.csv        the DAVOF diagnostics (consistency, closure)

and writes into the theme's data source (paths.tables_dir/figs_dir):

    tables/davof_normal_<study>.csv           the long table, every (alphaSource, model, N)
    tables/davof_normal_orders_<study>.csv    least-squares + pairwise orders per series
    tables/davof_normal_orders_<study>.tex    booktabs body of the same (\\input-able)
    figures/davof_normal_<metric>_<study>.png log-log error vs h, one per metric

Series = (alphaSource, model). Metrics = E_L1_N, E_L2_N, E_LINF_N (normal error
|n - n_exact| over the interface cells: mean, rms, max), E_L1_N_AW, E_L2_N_AW
(area-weighted) and E_AREA_REL (|sum |m_c| - 4 pi R^2| / 4 pi R^2). h = L/N.

--tpf-results <csv> merges the rows of the TwoPhaseFlow benchmark
run/benchmark/reconstruction/sphereNormal3D (results_sphereNormal3D.csv) as
model "tpf:<scheme>/<setAlphaMethod>", alphaSource "tpf", with its LNormalDiff1/2/Inf
and E_area_rel columns mapped onto E_L1_N / E_L2_N / E_LINF_N / E_AREA_REL.

Usage (from the repo root):
    python3 workflow/scripts/make_davof_normal_table.py studies/davof/sphereNormal3D \\
        [--theme davof] [--tpf-results /path/results_sphereNormal3D.csv]
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import paths

CSV_MODELS = "leiaTestDavofNormalModels.csv"
CSV_WIDE = "leiaTestDavofNormal.csv"
METRICS = ["E_L1_N", "E_L2_N", "E_LINF_N", "E_L1_N_AW", "E_L2_N_AW", "E_AREA_REL"]
FIG_METRICS = ["E_L2_N", "E_L1_N", "E_LINF_N", "E_AREA_REL"]
DIAG = ["MAX_CONSISTENCY", "SUM_M_REL", "MAX_ALPHA_DIFF_DA", "N_WISP"]

# Fixed hue per model (categorical slots in the validated order of the dataviz
# reference palette); the line style carries the alpha source. Never cycled.
COLORS = {
    "davof":     "#2a78d6",   # blue
    "plicRDF":   "#eb6834",   # orange
    "gradAlpha": "#1baf7a",   # aqua
    "isoAlpha":  "#eda100",   # yellow
}
TPF_COLORS = {"plicRDF": "#e87ba4", "gradAlpha": "#008300", "isoAlpha": "#4a3aa7"}
MARKERS = {"davof": "o", "plicRDF": "s", "gradAlpha": "D", "isoAlpha": "^"}
STYLES = {"exactSphere": "-", "linearInterpolant": "--", "planePhaseIndicator": "-.",
          "tpf": ":"}
SOURCE_LABEL = {"exactSphere": "exact sphere",
                "linearInterpolant": "linear interpolant",
                "planePhaseIndicator": "plane indicator",
                "tpf": "TwoPhaseFlow"}
METRIC_LABEL = {
    "E_L1_N": r"$L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L2_N": r"$L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_LINF_N": r"$L_\infty$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L1_N_AW": r"area-weighted $L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_L2_N_AW": r"area-weighted $L_2$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$",
    "E_AREA_REL": r"$|\sum_c |\mathbf{m}_c| - 4\pi R^2| / 4\pi R^2$",
}


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _fit(h, e):
    """Least-squares order p of e ~ h^p over the points with h, e > 0; (p, R2, n)."""
    m = [(a, b) for a, b in zip(h, e) if a and b and a > 0 and b > 0]
    if len(m) < 2:
        return None, None, len(m)
    lh, le = np.log([a for a, _ in m]), np.log([b for _, b in m])
    p = np.polyfit(lh, le, 1)
    r = le - np.polyval(p, lh)
    ss, st = float((r**2).sum()), float(((le - le.mean())**2).sum())
    return float(p[0]), (1.0 - ss/st if st > 0 else 1.0), len(m)


def _pairwise(h, e):
    out = []
    for k in range(len(h) - 1):
        a, b = e[k], e[k + 1]
        if a and b and a > 0 and b > 0 and h[k] > 0 and h[k + 1] > 0 and h[k] != h[k + 1]:
            out.append(np.log(a/b)/np.log(h[k]/h[k + 1]))
        else:
            out.append(None)
    return out


def _read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def collect(study_dir):
    """{(alphaSource, model): [record, ...]} with record = {N, h, metrics..., diag...}."""
    series = {}
    for meta in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        cdir = os.path.dirname(meta)
        mpath = os.path.join(cdir, CSV_MODELS)
        if not os.path.isfile(mpath) or os.path.getsize(mpath) == 0:
            continue
        with open(meta) as fh:
            tok = json.load(fh).get("tokens", {})
        n, L = _f(tok.get("N_CELLS")), _f(tok.get("DOMAIN_LENGTH"))
        if not n or not L:
            continue
        wide = {}
        wpath = os.path.join(cdir, CSV_WIDE)
        if os.path.isfile(wpath) and os.path.getsize(wpath) > 0:
            rows = _read_rows(wpath)
            if rows:
                wide = rows[-1]
        for row in _read_rows(mpath):
            model, src = row.get("MODEL"), row.get("ALPHA_SOURCE") or tok.get("DAVOF_ALPHA_SOURCE")
            if not model or not src:
                continue
            rec = {"N": int(n), "h": L/n, "R_OVER_H": _f(row.get("R_OVER_H")),
                   "N_INTERFACE": row.get("N_INTERFACE"), "N_WISP": row.get("N_WISP"),
                   "CPU_SECONDS": _f(row.get("CPU_SECONDS")), "case": os.path.basename(cdir)}
            for m in METRICS:
                rec[m] = _f(row.get(m))
            for d in DIAG:
                rec[d] = (wide.get(d) if model == "davof" else "")
            series.setdefault((src, model), []).append(rec)
    return series


def merge_tpf(series, path):
    """TwoPhaseFlow results_sphereNormal3D.csv rows -> series tpf:<scheme>/<setAlpha>."""
    for row in _read_rows(path):
        scheme, setA = row.get("scheme"), row.get("setAlphaMethod")
        n, h = _f(row.get("N")), _f(row.get("h"))
        if not scheme or not n or not h:
            continue
        rec = {"N": int(n), "h": h, "R_OVER_H": _f(row.get("R_over_h")),
               "N_INTERFACE": row.get("nInterface"), "N_WISP": "",
               "CPU_SECONDS": _f(row.get("k1")), "case": row.get("case", "")}
        rec.update({m: None for m in METRICS})
        rec["E_L1_N"] = _f(row.get("LNormalDiff1"))
        rec["E_L2_N"] = _f(row.get("LNormalDiff2"))
        rec["E_LINF_N"] = _f(row.get("LNormalDiffInf"))
        rec["E_AREA_REL"] = _f(row.get("E_area_rel"))
        for d in DIAG:
            rec[d] = ""
        series.setdefault(("tpf", f"tpf:{scheme}/{setA}"), []).append(rec)
    return series


def _style(src, model):
    base = model.split(":")[-1].split("/")[0] if model.startswith("tpf:") else model
    color = TPF_COLORS.get(base, "#e34948") if model.startswith("tpf:") else COLORS.get(base, "#e34948")
    return color, MARKERS.get(base, "x"), STYLES.get(src, "-")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="davof")
    ap.add_argument("--tpf-results", default=None,
                    help="results_sphereNormal3D.csv of the TwoPhaseFlow benchmark")
    a = ap.parse_args(argv)

    series = collect(a.study_dir)
    if a.tpf_results and os.path.isfile(a.tpf_results):
        merge_tpf(series, a.tpf_results)
    if not series:
        print(f"[davof] no completed cases under {a.study_dir}")
        return 1
    for recs in series.values():
        recs.sort(key=lambda r: -r["h"])          # coarse -> fine

    study = os.path.basename(os.path.normpath(a.study_dir))
    tables, figs = paths.tables_dir(a.theme), paths.figs_dir(a.theme)

    # --- long table ---------------------------------------------------------
    long_rows = []
    for (src, model), recs in sorted(series.items()):
        for r in recs:
            row = {"alphaSource": src, "model": model, "N": r["N"], "h": r["h"],
                   "R_over_h": r["R_OVER_H"], "N_interface": r["N_INTERFACE"],
                   "N_wisp": r["N_WISP"]}
            row.update({m: r[m] for m in METRICS})
            row.update({d: r[d] for d in DIAG})
            row["cpu_s"] = r["CPU_SECONDS"]
            row["case"] = r["case"]
            long_rows.append(row)
    lpath = os.path.join(tables, f"davof_normal_{study}.csv")
    with open(lpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(long_rows[0]))
        w.writeheader(); w.writerows(long_rows)
    print(f"[davof] wrote {lpath} ({len(long_rows)} rows)")

    # --- orders ---------------------------------------------------------------
    order_rows = []
    for (src, model), recs in sorted(series.items()):
        h = [r["h"] for r in recs]
        row = {"alphaSource": src, "model": model, "N_coarsest": recs[0]["N"],
               "N_finest": recs[-1]["N"], "n_levels": len(recs)}
        for m in METRICS:
            e = [r[m] for r in recs]
            p, r2, n = _fit(h, e)
            row[f"{m}_finest"] = e[-1]
            row[f"p_{m}"] = p
            row[f"R2_{m}"] = r2
            row[f"pairwise_{m}"] = " ".join(
                "--" if q is None else f"{q:.2f}" for q in _pairwise(h, e))
        order_rows.append(row)
    opath = os.path.join(tables, f"davof_normal_orders_{study}.csv")
    with open(opath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(order_rows[0]))
        w.writeheader(); w.writerows(order_rows)
    print(f"[davof] wrote {opath} ({len(order_rows)} rows)")

    def _num(x, fmt="{:.2f}"):
        return "--" if x is None else fmt.format(x)

    tpath = os.path.join(tables, f"davof_normal_orders_{study}.tex")
    with open(tpath, "w") as fh:
        fh.write("% Auto-generated by workflow/scripts/make_davof_normal_table.py -- do not edit.\n")
        fh.write("\\begin{tabular}{llrcccccc}\n\\toprule\n")
        fh.write("$\\alpha$ source & model & $N$ & $L_2(\\mathbf{n})$ & $p(L_2)$ & $R^2$ & "
                 "$p(L_1)$ & $p(L_\\infty)$ & $p(A)$ \\\\\n\\midrule\n")
        for r in order_rows:
            fh.write(f"{SOURCE_LABEL.get(r['alphaSource'], r['alphaSource'])} & "
                     f"{r['model'].replace('_', chr(92) + '_')} & {r['N_finest']} & "
                     f"{_num(r['E_L2_N_finest'], '{:.3e}')} & {_num(r['p_E_L2_N'])} & "
                     f"{_num(r['R2_E_L2_N'], '{:.3f}')} & {_num(r['p_E_L1_N'])} & "
                     f"{_num(r['p_E_LINF_N'])} & {_num(r['p_E_AREA_REL'])} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")
    print(f"[davof] wrote {tpath}")

    # --- figures: one log-log panel per metric --------------------------------
    hs = sorted({r["h"] for recs in series.values() for r in recs})
    for metric in FIG_METRICS:
        fig, ax = plt.subplots(figsize=(8.4, 5.0))
        ymin, ymax = np.inf, -np.inf
        for (src, model), recs in sorted(series.items()):
            h = [r["h"] for r in recs]
            e = [r[metric] for r in recs]
            pts = [(a, b) for a, b in zip(h, e) if a and b and b > 0]
            if not pts:
                continue
            p, _, _ = _fit(h, e)
            color, marker, ls = _style(src, model)
            lab = f"{model}, {SOURCE_LABEL.get(src, src)}"
            if p is not None:
                lab += rf"  ($\propto h^{{{p:.2f}}}$)"
            ax.loglog([q[0] for q in pts], [q[1] for q in pts], ls, color=color,
                      marker=marker, ms=5, lw=1.6, label=lab)
            ymin = min(ymin, min(q[1] for q in pts))
            ymax = max(ymax, max(q[1] for q in pts))
        if np.isfinite(ymin) and len(hs) >= 2:
            h0, h1 = hs[0], hs[-1]
            y1 = 0.6*ymax
            ax.loglog([h0, h1], [y1*h0/h1, y1], color="#7f7f7f", ls="-", lw=0.8, alpha=0.7)
            ax.annotate(r"$h^1$", xy=(h0*1.05, 1.15*y1*h0/h1), color="#7f7f7f", fontsize=9)
            y2 = 0.3*ymax
            ax.loglog([h0, h1], [y2*(h0/h1)**2, y2], color="#7f7f7f", ls="-", lw=0.8, alpha=0.7)
            ax.annotate(r"$h^2$", xy=(h0*1.05, 1.15*y2*(h0/h1)**2), color="#7f7f7f", fontsize=9)
        ax.set_xlabel(r"$h$ [m]")
        ax.set_ylabel(METRIC_LABEL.get(metric, metric))
        ax.set_title("DAVOF Gauss-identity normal on a sphere: " + metric.replace("_", " "),
                     fontsize=11)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()
        fpath = os.path.join(figs, f"davof_normal_{metric.lower()}_{study}.png")
        fig.savefig(fpath, dpi=200); plt.close(fig)
        print(f"[davof] wrote {fpath}")

    for r in order_rows:
        print(f"[davof] {r['alphaSource']:20s} {r['model']:28s} N={r['N_coarsest']}..{r['N_finest']}"
              f"  p(L1)={_num(r['p_E_L1_N'])}  p(L2)={_num(r['p_E_L2_N'])}"
              f"  p(Linf)={_num(r['p_E_LINF_N'])}  p(A)={_num(r['p_E_AREA_REL'])}"
              f"  L2@finest={_num(r['E_L2_N_finest'], '{:.3e}')}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
