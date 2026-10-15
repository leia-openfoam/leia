#!/usr/bin/env python3
"""Tables and figures of the TwoPhaseFlow advection studies of the DAVOF face
update (solver advectorVoF, cases/davof/{plane,sphere}Translation3D,
config/davof/*Translation3D_*.yaml; 2026-10-13).

Reads, for every case of the study,

    <case>/case_params.json   tokens TPF_ROW, TPF_VELOCITY_ARM, N_CELLS,
                              TARGET_ADVECTIVE_CO, DOMAIN_LENGTH, END_TIME
                              (h = DOMAIN_LENGTH/N_CELLS), the mesh class
    <case>/advectorVoF.csv    one row per step (tpf_advection_metrics.py)

and writes into the davof theme's data source (paths.tables_dir / figs_dir):

    tables/davof_advection_<study>.csv
        the long table, one row per case: the end values (the last non-blank
        row of every column), the run maxima, min alpha over the run, the
        interface cells at the first and the last step, the sliver and
        fit-fallback and ring-fallback maxima, the cpu time and the cpu
        time per step (cpu / nSteps), the status
    tables/davof_advection_history_<study>.csv
        the davof rows: every step's face-fraction errors and counts
    tables/davof_advection_growth_<study>.csv
        the davof rows: L1(alphaf) of the first 20 steps, the per-step
        factors g_k = L1_k/L1_{k-1}, their maximum and geometric mean
    tables/davof_advection_<study>_orders_Co<co>.{csv,tex}   per Courant number
        rows x N: L1(alphaf) at the end, max over the run, L1(n), L1(x),
        Eshape, Ebound at the end, min alpha, status; the least-squares and
        the pairwise orders over N per row; a diverged (aborted) rung prints
        "div." and is excluded from the fit
    tables/davof_advection_<study>_co_N<n>.{csv,tex}
        the Courant-number dependence at one N (--co-table-N; default the
        middle rung)
    tables/davof_advection_<study>_gate.{csv,tex}          plane studies only
        per row, N, Co, arm: max L1(alphaf), max Linf(alphaf), max L1(n) over
        the run and PASS/FAIL at 1e-12 for the davof rows ("-" for the
        comparators, which carry no face fractions)
    figures/davof_advection_<metric>_<study>_Co<co>.png
        error against h, log-log, one hue per row (tpf_rows.COLORS), the
        velocity arm by line style; metrics L1_alphaf_end, L1_n_end,
        L1_x_end, Eshape_end
    figures/davof_advection_history_<study>_N<n>_Co<co>.png
        L1(alphaf) against the step, log y, the davof rows

The .tex files are table BODIES (booktabs tabular only, the header comment
first); \\input them with their full name including the extension (the Co
name carries a dot). The number of a comparator's face-fraction error is -1
in the CSV and "--" in the tables.

Usage (from the repo root):
    python3 workflow/scripts/make_davof_advection_table.py studies/davof/sphereTranslation3D_hex \\
        [--theme davof] [--co-table-N 64]
"""
import argparse
import csv
import glob
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402
import tpf_rows  # noqa: E402
from make_davof_normal_table import _fit, _pairwise, _p, _num, _write_csv, COLORS  # noqa: E402,F401

CSV = "advectorVoF.csv"
NHIST = 20
GATE_TOL = 1e-12

# end-value / run-maximum columns of the long table: name -> (CSV column, reduction)
END_COLS = {
    "L1_alphaf_end": "L1_alphaf", "L2_alphaf_end": "L2_alphaf", "Linf_alphaf_end": "Linf_alphaf",
    "nInterfaceCells_end": "nInterfaceCells", "consistency_end": "consistency",
    "L1_n_end": "LNormalDiff1", "L2_n_end": "LNormalDiff2", "Linf_n_end": "LNormalDiffInf",
    "L1_x_end": "LCentreDiff1", "L2_x_end": "LCentreDiff2", "Linf_x_end": "LCentreDiffInf",
    "areaSum_end": "areaSum", "areaExact_end": "areaExact",
    "Eshape_end": "Eshape", "Emass_end": "Emass", "Ebound_end": "Ebound",
    "sumAlphaV_end": "sumAlphaV", "devPercent_end": "devPercent",
    "oneMinusMaxAlpha_end": "oneMinusMaxAlpha", "minAlpha_end": "minAlpha",
}
MAX_COLS = {
    "L1_alphaf_max": "L1_alphaf", "L2_alphaf_max": "L2_alphaf", "Linf_alphaf_max": "Linf_alphaf",
    "L1_n_max": "LNormalDiff1", "Linf_n_max": "LNormalDiffInf",
    "Eshape_max": "Eshape", "Ebound_max": "Ebound", "Emass_max": "Emass",
    "nInterfaceCells_max": "nInterfaceCells", "nSliver_max": "nSliver",
    "nReplaced_max": "nReplaced", "nFitFallback_max": "nFitFallback",
    "nRingFallback_max": "nRingFallback",
}
# the order tables: (long-table column, tex label, has an order)
ORDER_METRICS = [
    ("L1_alphaf_end", r"$L_1(\alpha_f)$, end", True),
    ("L1_alphaf_max", r"$\max_t L_1(\alpha_f)$", True),
    ("L1_n_end", r"$L_1(\mathbf{n})$, end", True),
    ("L1_x_end", r"$L_1(\mathbf{x})$, end", True),
    ("Eshape_end", r"$E_{\mathrm{shape}}$, end", True),
    ("Ebound_end", r"$E_{\mathrm{bound}}$, end", False),
    ("minAlpha_run", r"$\min_t \alpha$", False),
    ("status", "status", False),
]
FIG_METRICS = {
    "L1_alphaf_end": r"$L_1$ of $|\alpha_f - \alpha_{f,\mathrm{exact}}|$ at the end time",
    "L1_n_end": r"$L_1$ of $|\mathbf{n}_c - \mathbf{n}_{\mathrm{exact}}|$ at the end time",
    "L1_x_end": r"$L_1$ of the centre error at the end time",
    "Eshape_end": r"$E_{\mathrm{shape}} = \sum_c |\alpha_c - \alpha_{\mathrm{exact}}| V_c$ at the end time",
}
ARM_STYLE = {"normal": "-", "oblique": "--", "axial": "-", "": "-"}
STATUS_SHORT = {"diverged": "div.", "crashed": "crash", "incomplete": "inc.", "completed": "ok"}


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _last_nonblank(rows, col):
    for r in reversed(rows):
        v = r.get(col, "")
        if v not in ("", None):
            return v
    return ""


def _col_floats(rows, col):
    out = []
    for r in rows:
        v = _f(r.get(col))
        if v is not None:
            out.append(v)
    return out


def _row_rank(row):
    names = list(tpf_rows.ROWS)
    return names.index(row) if row in names else len(names)


def _fmt_co(co):
    return f"{float(co):g}"


def collect(study_dir):
    """One record per case; the per-step rows kept for the histories."""
    recs = []
    for params in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        case = os.path.dirname(params)
        meta = json.load(open(params))
        tok = meta.get("tokens", {})
        n = _f(tok.get("N_CELLS"))
        length = _f(tok.get("DOMAIN_LENGTH"))
        if not n or not length:
            continue
        rec = {
            "case": os.path.basename(case), "mesh": meta.get("mesh", ""),
            "row": tok.get("TPF_ROW", ""), "arm": tok.get("TPF_VELOCITY_ARM", ""),
            "N": int(round(n)), "Co": _fmt_co(tok.get("TARGET_ADVECTIVE_CO", "nan")),
            "h": length / n, "L": length, "T": _f(tok.get("END_TIME")),
            "dt": _f(tok.get("MAX_DELTA_T")), "scheme": tok.get("TPF_RECON_SCHEME", ""),
            "normalModel": tok.get("TPF_NORMAL_MODEL", ""),
            "gitCommit": meta.get("gitCommit", ""),
            "tpfCommit": "",   # filled from the case CSV below (the fork's commit)
        }
        path = os.path.join(case, CSV)
        rows = []
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            with open(path, newline="") as fh:
                rows = list(csv.DictReader(fh))
        rec["_steps"] = rows
        rec["nSteps"] = len(rows)
        rec["lastTime"] = _last_nonblank(rows, "TIME")
        rec["status"] = _last_nonblank(rows, "status") or ("absent" if not rows else "incomplete")
        rec["abortStep"] = _last_nonblank(rows, "abortStep")
        rec["abortReason"] = _last_nonblank(rows, "abortReason")
        # provenance of the solver: the TwoPhaseFlow fork's commit and version stamp the
        # converter recorded (cases before 2026-10-15 carry none: the fork commit per study
        # is then the one STATUS 10.17 names)
        rec["tpfCommit"] = _last_nonblank(rows, "tpfCommit") or "unrecorded"
        rec["tpfVersion"] = _last_nonblank(rows, "tpfVersion") or "unrecorded"
        for name, col in END_COLS.items():
            rec[name] = _last_nonblank(rows, col)
        for name, col in MAX_COLS.items():
            vals = _col_floats(rows, col)
            rec[name] = max(vals) if vals else ""
        vals = _col_floats(rows, "minAlpha")
        rec["minAlpha_run"] = min(vals) if vals else ""
        vals = _col_floats(rows, "minAlphaBound")
        rec["minAlphaBound_run"] = min(vals) if vals else ""
        cells = _col_floats(rows, "nInterfaceAdv") or _col_floats(rows, "nInterfaceCells")
        rec["nInterfaceAdv_first"] = int(cells[0]) if cells else ""
        rec["nInterfaceAdv_last"] = int(cells[-1]) if cells else ""
        rec["cpu"] = _last_nonblank(rows, "executionTime")
        # the cpu time per step (the timing studies): the solver's cumulative
        # ExecutionTime at the last step over the number of steps
        cpu = _f(rec["cpu"])
        rec["cpu_per_step"] = cpu / rec["nSteps"] if (cpu is not None and rec["nSteps"]) else ""
        recs.append(rec)
    recs.sort(key=lambda r: (_row_rank(r["row"]), r["row"], r["arm"], float(r["Co"]), r["N"]))
    return recs


def long_rows(recs):
    keys = [k for k in recs[0] if not k.startswith("_")]
    return [{k: r.get(k, "") for k in keys} for r in recs]


def history_rows(recs):
    out = []
    for r in recs:
        if r["scheme"] != "davof":
            continue
        for s in r["_steps"]:
            out.append({"row": r["row"], "arm": r["arm"], "N": r["N"], "Co": r["Co"],
                        "step": s.get("step", ""), "TIME": s.get("TIME", ""),
                        "L1_alphaf": s.get("L1_alphaf", ""), "L2_alphaf": s.get("L2_alphaf", ""),
                        "Linf_alphaf": s.get("Linf_alphaf", ""),
                        "nInterfaceCells": s.get("nInterfaceCells", ""),
                        "nInterfaceAdv": s.get("nInterfaceAdv", ""),
                        "nSliver": s.get("nSliver", ""), "nReplaced": s.get("nReplaced", ""),
                        "nFitFallback": s.get("nFitFallback", ""),
                        "nRingFallback": s.get("nRingFallback", ""),
                        "minAlpha": s.get("minAlpha", ""), "Eshape": s.get("Eshape", ""),
                        "status": s.get("status", "")})
    return out


def growth_rows(recs):
    out = []
    for r in recs:
        if r["scheme"] != "davof":
            continue
        l1 = [_f(s.get("L1_alphaf")) for s in r["_steps"][:NHIST]]
        row = {"row": r["row"], "arm": r["arm"], "N": r["N"], "Co": r["Co"], "status": r["status"]}
        for k in range(1, NHIST + 1):
            v = l1[k - 1] if k <= len(l1) else None
            row[f"L1_{k}"] = "" if v is None else v
        gs = []
        for k in range(2, NHIST + 1):
            g = None
            if k <= len(l1) and l1[k - 2] and l1[k - 1] and l1[k - 2] > 0 and l1[k - 1] > 0:
                g = l1[k - 1] / l1[k - 2]
            row[f"g_{k}"] = "" if g is None else g
            if g:
                gs.append(g)
        row["maxGrowth"] = max(gs) if gs else ""
        row["geoMeanGrowth"] = math.exp(sum(math.log(g) for g in gs) / len(gs)) if gs else ""
        out.append(row)
    return out


def _series_key(r):
    return (r["row"], r["arm"])


def _label(row, arm, with_arm):
    return f"{row}, {arm}" if with_arm and arm else row


def _tex_escape(s):
    return str(s).replace("_", r"\_")


def _value_cell(rec, metric):
    """The printed cell of a metric at one rung and its numeric value for the fit
    (None where nothing is fitted: a comparator's -1, a diverged rung)."""
    if metric == "status":
        return STATUS_SHORT.get(rec["status"], rec["status"]), None
    v = _f(rec.get(metric))
    if rec["status"] in ("diverged", "crashed", "incomplete", "absent"):
        return STATUS_SHORT.get(rec["status"], rec["status"]), None
    if v is None:
        return "--", None
    if metric.startswith("L1_alphaf") and v < 0:
        return "--", None
    if metric == "minAlpha_run":
        return f"{v:.1e}", None
    return f"{v:.2e}", v


def write_orders_tables(recs, study, tables, with_arm):
    """Per Courant number: rows x N with the orders over N."""
    cos = sorted({r["Co"] for r in recs}, key=float)
    ns = sorted({r["N"] for r in recs})
    written = []
    for co in cos:
        sub = [r for r in recs if r["Co"] == co]
        series = {}
        for r in sub:
            series.setdefault(_series_key(r), {})[r["N"]] = r
        csv_rows, tex_lines = [], []
        head = ["row"] + (["arm"] if with_arm else []) + ["metric"] + [f"N={n}" for n in ns] + ["p_ls", "p_pairwise"]
        tex_lines.append("% Auto-generated by workflow/scripts/make_davof_advection_table.py -- do not edit.")
        tex_lines.append(r"\begin{tabular}{l" + ("l" if with_arm else "") + "l" + "r" * len(ns) + "rl}")
        tex_lines.append(r"\toprule")
        tex_lines.append(" & ".join(["row"] + (["arm"] if with_arm else []) + ["metric"]
                                    + [f"$N={n}$" for n in ns] + ["$p$", "pairwise"]) + r" \\")
        tex_lines.append(r"\midrule")
        first_block = True
        for key in sorted(series, key=lambda k: (_row_rank(k[0]), k[0], k[1])):
            row, arm = key
            per_n = series[key]
            if not first_block:
                tex_lines.append(r"\midrule")
            first_block = False
            for metric, label, has_order in ORDER_METRICS:
                cells, hs, es = [], [], []
                for n in ns:
                    rec = per_n.get(n)
                    if rec is None:
                        cells.append("")
                        continue
                    cell, v = _value_cell(rec, metric)
                    cells.append(cell)
                    if has_order and v is not None and v > 0:
                        hs.append(rec["h"]); es.append(v)
                p = pw = None
                if has_order and len(hs) >= 2:
                    p = _fit(hs, es)[0]
                    pw = _pairwise(hs, es)
                pw_str = " / ".join(_p(x) for x in pw) if pw else ""
                crow = {"row": row, "metric": metric}
                if with_arm:
                    crow["arm"] = arm
                for n, c in zip(ns, cells):
                    crow[f"N={n}"] = c
                crow["p_ls"] = _p(p) if has_order else ""
                crow["p_pairwise"] = pw_str
                csv_rows.append({k: crow.get(k, "") for k in head})
                tex_lines.append(" & ".join(
                    [r"\texttt{" + _tex_escape(row) + "}"] + ([arm] if with_arm else []) + [label]
                    + cells + [(_p(p) if has_order else ""), pw_str]) + r" \\")
        tex_lines.append(r"\bottomrule")
        tex_lines.append(r"\end{tabular}")
        base = os.path.join(tables, f"davof_advection_{study}_orders_Co{co}")
        if csv_rows:
            _write_csv(base + ".csv", csv_rows)
            with open(base + ".tex", "w") as fh:
                fh.write("\n".join(tex_lines) + "\n")
            written.append(base)
    return written


def write_co_table(recs, study, tables, with_arm, n_pick):
    ns = sorted({r["N"] for r in recs})
    if not ns:
        return None
    if n_pick is None or n_pick not in ns:
        n_pick = ns[len(ns) // 2]
    sub = [r for r in recs if r["N"] == n_pick]
    cos = sorted({r["Co"] for r in sub}, key=float)
    series = {}
    for r in sub:
        series.setdefault(_series_key(r), {})[r["Co"]] = r
    head = ["row"] + (["arm"] if with_arm else []) + ["metric"] + [f"Co={co}" for co in cos]
    csv_rows, tex_lines = [], []
    tex_lines.append("% Auto-generated by workflow/scripts/make_davof_advection_table.py -- do not edit.")
    tex_lines.append(r"\begin{tabular}{l" + ("l" if with_arm else "") + "l" + "r" * len(cos) + "}")
    tex_lines.append(r"\toprule")
    tex_lines.append(" & ".join(["row"] + (["arm"] if with_arm else []) + ["metric"]
                                + [f"Co $= {co}$" for co in cos]) + r" \\")
    tex_lines.append(r"\midrule")
    first_block = True
    for key in sorted(series, key=lambda k: (_row_rank(k[0]), k[0], k[1])):
        row, arm = key
        if not first_block:
            tex_lines.append(r"\midrule")
        first_block = False
        for metric, label, _ in ORDER_METRICS:
            cells = []
            for co in cos:
                rec = series[key].get(co)
                cells.append(_value_cell(rec, metric)[0] if rec else "")
            crow = {"row": row, "metric": metric}
            if with_arm:
                crow["arm"] = arm
            for co, c in zip(cos, cells):
                crow[f"Co={co}"] = c
            csv_rows.append({k: crow.get(k, "") for k in head})
            tex_lines.append(" & ".join(
                [r"\texttt{" + _tex_escape(row) + "}"] + ([arm] if with_arm else []) + [label] + cells) + r" \\")
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    base = os.path.join(tables, f"davof_advection_{study}_co_N{n_pick}")
    if csv_rows:
        _write_csv(base + ".csv", csv_rows)
        with open(base + ".tex", "w") as fh:
            fh.write("\n".join(tex_lines) + "\n")
    return base


def write_gate_table(recs, study, tables):
    """Plane studies: the round-off gate per (row, N, Co, arm)."""
    csv_rows, tex_lines = [], []
    tex_lines.append("% Auto-generated by workflow/scripts/make_davof_advection_table.py -- do not edit.")
    tex_lines.append(r"\begin{tabular}{llrrrrrl}")
    tex_lines.append(r"\toprule")
    tex_lines.append(r"row & arm & $N$ & Co & $\max_t L_1(\alpha_f)$ & $\max_t L_\infty(\alpha_f)$ & "
                     r"$\max_t L_1(\mathbf{n})$ & gate \\")
    tex_lines.append(r"\midrule")
    last_row = None
    for r in recs:
        davof = r["scheme"] == "davof"
        l1 = _f(r.get("L1_alphaf_max")); li = _f(r.get("Linf_alphaf_max")); ln = _f(r.get("L1_n_max"))
        if not davof:
            gate = "-"
        elif r["status"] != "completed" or l1 is None or ln is None:
            gate = "FAIL (" + STATUS_SHORT.get(r["status"], r["status"]) + ")"
        else:
            gate = "PASS" if (l1 < GATE_TOL and ln < GATE_TOL) else "FAIL"
        csv_rows.append({"row": r["row"], "arm": r["arm"], "N": r["N"], "Co": r["Co"],
                         "L1_alphaf_max": "" if l1 is None else l1,
                         "Linf_alphaf_max": "" if li is None else li,
                         "L1_n_max": "" if ln is None else ln,
                         "status": r["status"], "gate": gate})
        if last_row is not None and r["row"] != last_row:
            tex_lines.append(r"\midrule")
        last_row = r["row"]
        tex_lines.append(" & ".join([
            r"\texttt{" + _tex_escape(r["row"]) + "}", r["arm"], str(r["N"]), r["Co"],
            "--" if (l1 is None or l1 < 0) else f"{l1:.1e}",
            "--" if (li is None or li < 0) else f"{li:.1e}",
            "--" if ln is None else f"{ln:.1e}", gate]) + r" \\")
    tex_lines.append(r"\bottomrule")
    tex_lines.append(r"\end{tabular}")
    base = os.path.join(tables, f"davof_advection_{study}_gate")
    if csv_rows:
        _write_csv(base + ".csv", csv_rows)
        with open(base + ".tex", "w") as fh:
            fh.write("\n".join(tex_lines) + "\n")
    return base


def plot_convergence(recs, study, figs, with_arm):
    cos = sorted({r["Co"] for r in recs}, key=float)
    out = []
    for co in cos:
        sub = [r for r in recs if r["Co"] == co]
        for metric, ylabel in FIG_METRICS.items():
            series = {}
            for r in sub:
                v = _f(r.get(metric))
                if v is None or v <= 0 or r["status"] != "completed":
                    continue
                series.setdefault(_series_key(r), []).append((r["h"], v))
            if not series:
                continue
            fig, ax = plt.subplots(figsize=(6.4, 4.6))
            for key in sorted(series, key=lambda k: (_row_rank(k[0]), k[0], k[1])):
                row, arm = key
                pts = sorted(series[key])
                ax.loglog([p[0] for p in pts], [p[1] for p in pts],
                          ARM_STYLE.get(arm, "-"), marker=tpf_rows.MARKERS.get(row, "o"),
                          color=tpf_rows.COLORS.get(row, "#777777"),
                          label=_label(tpf_rows.LABELS.get(row, row), arm, with_arm))
            hs = sorted({r["h"] for r in sub})
            if len(hs) >= 2:
                h0, h1 = hs[0], hs[-1]
                ymin = min(min(v for _, v in s) for s in series.values())
                for p, ls in ((1, ":"), (2, "-.")):
                    ax.loglog([h0, h1], [ymin, ymin * (h1 / h0) ** p], ls, color="#999999",
                              label=f"$h^{p}$")
            ax.set_xlabel("$h$")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{study}, Co = {co}")
            ax.grid(True, which="both", alpha=0.3)
            ax.legend(fontsize=7)
            fig.tight_layout()
            path = os.path.join(figs, f"davof_advection_{metric}_{study}_Co{co}.png")
            fig.savefig(path, dpi=150)
            plt.close(fig)
            out.append(path)
    return out


def plot_histories(recs, study, figs, with_arm):
    out = []
    combos = sorted({(r["N"], r["Co"]) for r in recs if r["scheme"] == "davof"},
                    key=lambda t: (t[0], float(t[1])))
    for n, co in combos:
        sub = [r for r in recs if r["scheme"] == "davof" and r["N"] == n and r["Co"] == co]
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        drawn = False
        for r in sub:
            xs, ys = [], []
            for s in r["_steps"]:
                v = _f(s.get("L1_alphaf")); k = _f(s.get("step"))
                if v is not None and v > 0 and k is not None:
                    xs.append(k); ys.append(v)
            if not xs:
                continue
            drawn = True
            lab = _label(tpf_rows.LABELS.get(r["row"], r["row"]), r["arm"], with_arm)
            if r["status"] == "diverged":
                lab += f" (abort at step {r['abortStep']})"
            ax.semilogy(xs, ys, ARM_STYLE.get(r["arm"], "-"),
                        color=tpf_rows.COLORS.get(r["row"], "#777777"), label=lab)
        if not drawn:
            plt.close(fig)
            continue
        ax.set_xlabel("step")
        ax.set_ylabel(r"$L_1$ of $|\alpha_f - \alpha_{f,\mathrm{exact}}|$")
        ax.set_title(f"{study}, N = {n}, Co = {co}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
        fig.tight_layout()
        path = os.path.join(figs, f"davof_advection_history_{study}_N{n}_Co{co}.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        out.append(path)
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="davof")
    ap.add_argument("--co-table-N", type=int, default=None,
                    help="the N of the Courant-number table (default: the middle rung)")
    a = ap.parse_args(argv)
    study_dir = os.path.abspath(a.study_dir.rstrip("/"))
    study = os.path.basename(study_dir)
    recs = collect(study_dir)
    if not recs:
        print(f"[davof-advection] no case with case_params.json under {study_dir}", file=sys.stderr)
        return 1
    tables, figs = paths.tables_dir(a.theme), paths.figs_dir(a.theme)
    with_arm = len({r["arm"] for r in recs}) > 1
    is_plane = any(r["case"].startswith("plane") for r in recs)
    _write_csv(os.path.join(tables, f"davof_advection_{study}.csv"), long_rows(recs))
    n_cases = sum(1 for r in recs if r["_steps"])
    print(f"[davof-advection] {study}: {len(recs)} cases ({n_cases} with {CSV}), "
          f"rows {sorted({r['row'] for r in recs}, key=_row_rank)}, N {sorted({r['N'] for r in recs})}, "
          f"Co {sorted({r['Co'] for r in recs}, key=float)}"
          + (f", arms {sorted({r['arm'] for r in recs})}" if with_arm else ""))
    hist = history_rows(recs)
    if hist:
        _write_csv(os.path.join(tables, f"davof_advection_history_{study}.csv"), hist)
    grow = growth_rows(recs)
    if grow:
        _write_csv(os.path.join(tables, f"davof_advection_growth_{study}.csv"), grow)
    for base in write_orders_tables(recs, study, tables, with_arm):
        print(f"[davof-advection] wrote {os.path.basename(base)}.{{csv,tex}}")
    base = write_co_table(recs, study, tables, with_arm, a.co_table_N)
    if base:
        print(f"[davof-advection] wrote {os.path.basename(base)}.{{csv,tex}}")
    if is_plane:
        base = write_gate_table(recs, study, tables)
        print(f"[davof-advection] wrote {os.path.basename(base)}.{{csv,tex}}")
    for path in plot_convergence(recs, study, figs, with_arm) + plot_histories(recs, study, figs, with_arm):
        print(f"[davof-advection] wrote {os.path.relpath(path, paths.REPO)}")
    # the console summary: the end values per case
    print(f"[davof-advection] {'row':22s} {'arm':8s} {'N':>4s} {'Co':>5s} {'L1(af)end':>10s} {'L1(af)max':>10s} "
          f"{'L1(n)':>10s} {'L1(x)':>10s} {'Eshape':>10s} {'minAlpha':>10s} {'cells':>11s} {'steps':>5s} status")
    for r in recs:
        def s(k, fmt="{:.2e}"):
            v = _f(r.get(k))
            return "--" if v is None else fmt.format(v)
        print(f"[davof-advection] {r['row']:22s} {r['arm']:8s} {r['N']:4d} {r['Co']:>5s} {s('L1_alphaf_end'):>10s} "
              f"{s('L1_alphaf_max'):>10s} {s('L1_n_end'):>10s} {s('L1_x_end'):>10s} {s('Eshape_end'):>10s} "
              f"{s('minAlpha_run', '{:.1e}'):>10s} {str(r['nInterfaceAdv_first']) + '->' + str(r['nInterfaceAdv_last']):>11s} "
              f"{r['nSteps']:5d} {r['status']}" + (f" (step {r['abortStep']})" if r["abortStep"] else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
