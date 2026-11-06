#!/usr/bin/env python3
"""The comparison tables of the swept-face advection with plicRDF + isoAdvector at OpenFOAM
Courant numbers 1 and 0.2 (sec:swept-cmp): one table body per flow from data/swept/cmp_table.csv
(the per-case summary of the timed comparison, copied from the cluster). Per Courant number and
mesh the three rows (swept-face advection with the ring normal, with the exact normal, plicRDF
with isoAdvector) with the error norms at the end of the run, the bounds, the total wall-clock
time of the case at the stated rank count (median of three repetitions) and its ratio to
plicRDF's; after the three meshes a row with the mean of the two pairwise orders over
N = 32, 64, 128 per method. Writes table_cmp_linex.tex, table_cmp_linediag.tex,
table_cmp_rotation.tex. Run on Lichtenberg, not on the laptop:  python3 make_table_cmp.py
"""
import csv
import math

SRC = "data/swept/cmp_table.csv"
FLOWS = {"LineX": "table_cmp_linex.tex", "LineDiag": "table_cmp_linediag.tex", "Rot": "table_cmp_rotation.tex"}
METHODS = [("ring", "swept-face, ring normal"), ("exact", "swept-face, exact normal"), ("rdf", "plicRDF + isoAdvector")]
NS = [32, 64, 128]


def f2(x):
    return "--" if x is None or math.isnan(x) else "%.2e" % x


def f1(x):
    return "--" if x is None or math.isnan(x) else "%.1e" % x


def order(vals):
    """Mean of the pairwise orders log2(E(N)/E(2N)) over the three meshes."""
    if any(v is None or math.isnan(v) or v <= 0 for v in vals):
        return "--"
    p = [math.log2(vals[i] / vals[i + 1]) for i in range(len(vals) - 1)]
    return "%.2f" % (sum(p) / len(p))


rows = {}
with open(SRC) as fh:
    for r in csv.DictReader(fh):
        co = 1 if float(r["ofCo"]) > 0.5 else 0.2
        key = (r["flow"], co, int(r["N"]), r["method"])
        rows[key] = {k: (float(v) if v not in ("", "nan") else float("nan")) for k, v in r.items()
                     if k in ("L1x", "Eshape", "L1n", "L1af", "minAlpha", "excess", "wall_med", "wall_ratio_to_rdf", "np")}

for flow, out in FLOWS.items():
    lines = ["%% written by make_table_cmp.py from %s -- do not edit" % SRC,
             r"\begin{tabular}{llrrrrrrrr}", r"\toprule",
             r"Co & $N$ (ranks) & method & $L_1(\mathbf{x})$ & $E_{\mathrm{shape}}$ & $L_1(\mathbf{n})$ & $L_1(\alpha_f)$ & $\min_t\alpha$ / excess & wall s & ratio \\",
             r"\midrule"]
    for co in (1, 0.2):
        for n in NS:
            for m, name in METHODS:
                r = rows[(flow, co, n, m)]
                first = ("%g" % co) if (n == NS[0] and m == METHODS[0][0]) else ""
                nn = ("%d (%d)" % (n, int(r["np"]))) if m == METHODS[0][0] else ""
                lines.append(" & ".join([first, nn, name, f2(r["L1x"]), f2(r["Eshape"]), f2(r["L1n"]), f2(r["L1af"]),
                                         "%s / %s" % (f1(r["minAlpha"]), f1(r["excess"])),
                                         "%.1f" % r["wall_med"], "%.1f" % r["wall_ratio_to_rdf"]]) + r" \\")
        for m, name in METHODS:
            vals = [rows[(flow, co, n, m)] for n in NS]
            lines.append(" & ".join(["", "order" if m == METHODS[0][0] else "", name,
                                     order([v["L1x"] for v in vals]), order([v["Eshape"] for v in vals]),
                                     order([v["L1n"] for v in vals]), order([v["L1af"] for v in vals]), "", "", ""]) + r" \\")
        if co == 1:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("wrote", out)
