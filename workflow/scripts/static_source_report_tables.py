#!/usr/bin/env python3
"""LaTeX tables of the static source report (docs/gradient-controlled-level-set/gcls-static-source-report)
from the curated tables of static_source_metrics.py (analyse + summary).

    static_source_report_tables.py <tables dir> <out dir>

<tables dir> holds staticSource_{representation,survival,growth,fields,verdicts}.csv (the curated copies
under docs/.../gcls-level-set-article/data/tables). Writes into <out dir>:
  ss_representation.tex  tau = 0, no source: volume error, profile volume difference, orders
  ss_growth.tex          growth rate of the band gradient error against G(0), per law and region
  ss_survival.tex        first write with dist3_2 > 1, per law, start and N
  ss_shift.tex           zero-crossing shift of d in the band at tau = 80, per law and N, with the order
  ss_verdicts.tex        the count of every verdict per pre-registered line
Norms are L1 and L2 only."""
import csv
import collections
import math
import os
import sys

import numpy as np

LAWS = ["linearQ", "linearZ", "saturatedLinearZ", "twoThirdsZReg", "cubicQ", "cubicZ", "softWall", "beta"]
INITS = ["sdf", "blendH2", "blendH1", "blend01", "algNorm", "alg", "modM2", "modM4"]
INIT_TEX = {"sdf": "$d$", "blendH2": "$\\varepsilon=(h/R)^2$", "blendH1": "$\\varepsilon=h/R$",
            "blend01": "$\\varepsilon=0.1$", "algNorm": "$(r^2-R^2)/(2R)$", "alg": "$r^2-R^2$",
            "modM2": "$m=2$", "modM4": "$m=4$"}
NS = (64, 128, 256, 512)


def G0(law):
    return {"linearQ": 1.0, "beta": 1.0, "linearZ": 0.5, "saturatedLinearZ": 1 / 3,
            "twoThirdsZReg": 0.5 / (0.25 + 0.0004) ** (1 / 6), "cubicQ": 1.0, "cubicZ": 0.125,
            "softWall": 1.0}[law]


def arm(law, moll):
    return f"beta-{moll}" if law == "beta" else f"gc-{law}-none-{moll}"


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def order(hs, vs):
    pts = [(h, abs(v)) for h, v in zip(hs, vs) if math.isfinite(v) and v != 0]
    if len(pts) < 3:
        return float("nan")
    return float(np.polyfit(np.log([p[0] for p in pts]), np.log([p[1] for p in pts]), 1)[0])


def rd(d, n):
    p = os.path.join(d, f"staticSource_{n}.csv")
    return list(csv.DictReader(open(p))) if os.path.isfile(p) else []


def e(x, fmt="{:.2e}"):
    return "--" if not math.isfinite(x) else fmt.format(x)


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    R, S, Gr, F, V = (rd(src, n) for n in ("representation", "survival", "growth", "fields", "verdicts"))

    # representation
    sdf = {int(r["N"]): r for r in R if r["init"] == "sdf"}
    lines = ["\\begin{tabular}{lrrrrrr}", "\\toprule",
             "$\\psi_0$ & volEx $N{=}64$ & volEx $N{=}512$ & order & $\\Delta V$ $N{=}64$ & $\\Delta V$ $N{=}512$ & order \\\\",
             "\\midrule"]
    for init in INITS:
        rs = sorted((r for r in R if r["init"] == init), key=lambda r: int(r["N"]))
        if not rs:
            continue
        hs = [1 / int(r["N"]) for r in rs]
        v = [f(r["volEx"]) for r in rs]
        dv = [(f(r["VDA"]) - f(sdf[int(r["N"])]["VDA"])) / f(r["Vex"]) for r in rs]
        lines.append(f"{INIT_TEX[init]} & {e(v[0])} & {e(v[-1])} & {e(order(hs, v), '{:.2f}')} & "
                     f"{e(dv[0]) if init != 'sdf' else '0'} & {e(dv[-1]) if init != 'sdf' else '0'} & "
                     f"{e(order(hs, dv), '{:.2f}') if init != 'sdf' else '--'} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(out, "ss_representation.tex"), "w").write("\n".join(lines) + "\n")

    # growth
    lines = ["\\begin{tabular}{lrll}", "\\toprule", "law & $G(0)$ & rate, band & rate, whole domain \\\\", "\\midrule"]
    for law in LAWS:
        cells = []
        for moll in ("band", "none"):
            rates = [f(g["rate"]) for g in Gr if g["arm"] == arm(law, moll) and g["rate"] not in ("", None)]
            rates = [x for x in rates if math.isfinite(x)]
            cells.append("no growth phase" if not rates else f"{min(rates):.3f} to {max(rates):.3f} ({len(rates)})")
        lines.append(f"\\texttt{{{law}}} & {G0(law):.3f} & {cells[0]} & {cells[1]} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(out, "ss_growth.tex"), "w").write("\n".join(lines) + "\n")

    # survival
    for moll, name in (("band", "band"), ("none", "none")):
        lines = ["\\begin{tabular}{l" + "c" * len(INITS) + "}", "\\toprule",
                 "law & " + " & ".join(INIT_TEX[i] for i in INITS) + " \\\\", "\\midrule"]
        for law in LAWS:
            cells = []
            for init in INITS:
                s = []
                for N in NS:
                    r = next((r for r in S if r["arm"] == arm(law, moll) and r["init"] == init and int(r["N"]) == N), None)
                    s.append("" if r is None else ("-" if r["tauFail"] == "none" else f"{f(r['tauFail']):.0f}"))
                cells.append("/".join(s))
            lines.append(f"\\texttt{{{law}}} & " + " & ".join(cells) + " \\\\")
        lines += ["\\bottomrule", "\\end{tabular}"]
        open(os.path.join(out, f"ss_survival_{name}.tex"), "w").write("\n".join(lines) + "\n")

    # shift of d in the band at tau = 80
    by = {(r["arm"], r["init"], int(r["N"]), f(r["tau"])): r for r in F}
    lines = ["\\begin{tabular}{lrrrrr}", "\\toprule",
             "law & $N{=}64$ & $N{=}128$ & $N{=}256$ & $N{=}512$ & order (length) \\\\", "\\midrule"]
    for law in LAWS:
        vals = [f(by[(arm(law, 'band'), 'sdf', N, 80.0)]["shift2"]) if (arm(law, 'band'), 'sdf', N, 80.0) in by
                else float("nan") for N in NS]
        o = order([1 / N for N in NS], [v / N for v, N in zip(vals, NS)])
        lines.append(f"\\texttt{{{law}}} & " + " & ".join(e(v, "{:.1e}") for v in vals) + f" & {e(o, '{:.2f}')} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(out, "ss_shift.tex"), "w").write("\n".join(lines) + "\n")

    # verdicts
    c = collections.defaultdict(collections.Counter)
    for v in V:
        c[v["item"]][v["verdict"]] += 1
    lines = ["\\begin{tabular}{ll}", "\\toprule", "line & verdicts \\\\", "\\midrule"]
    for k in sorted(c, key=lambda s: [int(x) if x.isdigit() else x for x in s.replace("P", "").replace("S", "0.").split(".")]):
        lines.append(f"{k} & " + ", ".join(f"{n} {w}" for w, n in sorted(c[k].items())) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(out, "ss_verdicts.tex"), "w").write("\n".join(lines) + "\n")
    print("report tables written to", out)


if __name__ == "__main__":
    main()
