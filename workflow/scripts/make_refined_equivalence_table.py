#!/usr/bin/env python3
"""Refined vs uniform at matched N_CELLS: equivalence rows, cost ratio, fitted orders.

For every (refined study, uniform study) pair, arms are joined on N_CELLS -- on the refined
mesh that is the FINE count, so matched N_CELLS means matched h, matched capillary dt and,
asserted here, an equal number of steps and equal read-out instants. Read at t = T/2 and
t = T: meanMagUPrime, l2MagUPrime (L1 and L2 of the parasitic velocity -- never the
maximum, which does not converge), phaseVolumeRelError, zeroSetRadialL2 (shape),
pLaplace relative to the exact 2 sigma/R, kErrL2Band and A2hL2Band; each with the
relative difference refined vs uniform. Core-hours from the solver log's last ClockTime
times np, and their ratio. Then, per ladder, least-squares orders p with the correlation
R of log(metric at T) against log(h), refined and uniform side by side.

Pairs whose two arms are not BOTH complete (last TIME within 0.1 % of END_TIME) are
reported as incomplete and excluded; a number nobody has seen land is not reported.

Writes refined_equivalence.csv and refined_convergence_orders.tex under the article's
data/tables and prints both.

Usage:  python3 workflow/scripts/make_refined_equivalence_table.py [--root .]
            [--pairs refinedStudy:uniformStudy[,...]] [--controls study[,...]]
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import sys

SOLVER_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
SOLVER_LOG = "log.leiaSemiLagrangianLevelSetTwoPhaseFoam"
DEFAULT_PAIRS = ("stationaryDroplet3DrefinedWB:stationaryDroplet3DuniformWB",
                 "stationaryDroplet3Drefined:stationaryDroplet3Duniform",
                 "stationaryDroplet3Drefined:stationaryDroplet3Duniform120",
                 "polyDroplet3Drefined_r13p8:polyDroplet3D_r13p8")
# Controls are compared with the refined ladder arm of the same N_CELLS.
DEFAULT_CONTROLS = ("stationaryDroplet3DrefinedL2:stationaryDroplet3Drefined",
                    "stationaryDroplet3DrefinedBall:stationaryDroplet3Drefined")
METRICS = ("meanMagUPrime", "l2MagUPrime", "phaseVolumeRelError", "zeroSetRadialL2",
           "pLaplaceRelError", "kErrL2Band", "A2hL2Band")
# Studies that are rungs of the SAME ladder for the order fit (the N = 120 uniform twin runs
# as its own study because it needs np 64 and ~90 core-hours; it is the fourth uniform rung).
LADDER_ALIASES = {"stationaryDroplet3Duniform120": "stationaryDroplet3Duniform"}


def _f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def load_arms(root, study):
    arms = {}
    for d in sorted(glob.glob(os.path.join(root, "studies", study, "*_[0-9]*"))):
        cp = os.path.join(d, "case_params.json")
        f = os.path.join(d, SOLVER_CSV)
        if not (os.path.isfile(cp) and os.path.isfile(f)):
            continue
        meta = json.load(open(cp))
        tok = meta["tokens"]
        rows = list(csv.DictReader(open(f)))
        if not rows:
            continue
        key = (int(_f(tok.get("N_CELLS"))), tok.get("REFINE_SOURCE", "interface"),
               int(_f(tok.get("REFINE_LEVELS"), 0)))
        arms[key] = {"dir": d, "meta": meta, "tok": tok, "rows": rows}
    return arms


def complete(arm):
    t_end = _f(arm["tok"].get("END_TIME"))
    t_last = _f(arm["rows"][-1].get("TIME"))
    return t_end and t_last and abs(t_last - t_end) <= 1e-3 * t_end


def at_time(rows, t):
    return min(rows, key=lambda r: abs(_f(r.get("TIME"), float("inf")) - t))


def metric(row, name, tok):
    if name == "pLaplaceRelError":
        sigma, R = _f(tok.get("SIGMA"), 0.07274), _f(tok.get("DROPLET_RADIUS"), 1e-3)
        exact = 2.0 * sigma / R                       # 3D sphere: kappa = 2/R
        p = _f(row.get("pLaplace"))
        return abs(p - exact) / exact if p is not None else None
    return _f(row.get(name))


def core_hours(arm):
    path = os.path.join(arm["dir"], SOLVER_LOG)
    if not os.path.isfile(path):
        return None
    last = None
    with open(path, errors="replace") as fh:
        for line in fh:
            if "ClockTime = " in line:
                last = line
    if not last:
        return None
    m = re.search(r"ClockTime = ([\d.]+) s", last)
    return float(m.group(1)) * int(arm["meta"].get("np", 1)) / 3600.0 if m else None


def fit(hs, es):
    pts = [(math.log(h), math.log(e)) for h, e in zip(hs, es) if h and e and e > 0]
    if len(pts) < 2:
        return None, None
    n = len(pts)
    mx = sum(x for x, _ in pts) / n
    my = sum(y for _, y in pts) / n
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    syy = sum((y - my) ** 2 for _, y in pts)
    if sxx == 0 or syy == 0:
        return None, None
    return sxy / sxx, sxy / math.sqrt(sxx * syy)


def compare(a, b, label_a, label_b, pair, ncell, source="", levels=""):
    """Rows for one matched pair (a = refined/control, b = reference)."""
    if len(a["rows"]) != len(b["rows"]):
        raise RuntimeError(f"{pair} N={ncell}: unequal step counts {len(a['rows'])} vs "
                           f"{len(b['rows'])} -- endpoint metrics are not comparable")
    t_end = _f(a["tok"]["END_TIME"])
    out = []
    ch_a, ch_b = core_hours(a), core_hours(b)
    # maximum over time of the two parasitic-velocity norms (the well-balanced gate's
    # read-out: the peak of a round-off-level signal, not its endpoint)
    for m in ("meanMagUPrime", "l2MagUPrime"):
        va = max(_f(r.get(m), 0.0) for r in a["rows"]); vb = max(_f(r.get(m), 0.0) for r in b["rows"])
        out.append({"pair": pair, "N_CELLS": ncell, "source": source, "levels": levels,
                    "t": "max_t", "metric": m, label_a: va, label_b: vb,
                    "relDiff": ((va - vb) / abs(vb)) if vb else "", "steps": len(a["rows"]),
                    "coreHours_" + label_a: ch_a if ch_a is not None else "",
                    "coreHours_" + label_b: ch_b if ch_b is not None else "",
                    "coreHoursRatio": (ch_b / ch_a) if (ch_a and ch_b) else ""})
    for t in (0.5 * t_end, t_end):
        ra, rb = at_time(a["rows"], t), at_time(b["rows"], t)
        if abs(_f(ra["TIME"]) - _f(rb["TIME"])) > 1e-9 * max(1.0, t):
            raise RuntimeError(f"{pair} N={ncell}: read-out instants differ at t={t}")
        for m in METRICS:
            va, vb = metric(ra, m, a["tok"]), metric(rb, m, b["tok"])
            out.append({"pair": pair, "N_CELLS": ncell, "source": source, "levels": levels,
                        "t": _f(ra["TIME"]), "metric": m,
                        label_a: va, label_b: vb,
                        "relDiff": ((va - vb) / abs(vb)) if (va is not None and vb) else "",
                        "steps": len(a["rows"]),
                        "coreHours_" + label_a: ch_a if ch_a is not None else "",
                        "coreHours_" + label_b: ch_b if ch_b is not None else "",
                        "coreHoursRatio": (ch_b / ch_a) if (ch_a and ch_b) else ""})
    return out


def orders(arms, tok_h):
    """Fitted p (R) at t = T per metric over the arms of one ladder."""
    res = {}
    keys = sorted(k for k in arms if complete(arms[k]))
    if len(keys) < 2:
        return res
    hs = [tok_h(arms[k]) for k in keys]
    for m in METRICS:
        es = [metric(arms[k]["rows"][-1], m, arms[k]["tok"]) for k in keys]
        p, r = fit(hs, es)
        res[m] = (p, r, len(keys))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--pairs", default=",".join(DEFAULT_PAIRS))
    ap.add_argument("--controls", default=",".join(DEFAULT_CONTROLS))
    ap.add_argument("--outdir", default="docs/semi-lagrangian-level-set/sl-level-set-article/data/tables")
    a = ap.parse_args()

    rows, order_rows, notes = [], [], []
    ladders = {}
    for spec, kind in [(s, "pair") for s in a.pairs.split(",") if s] + \
                      [(s, "control") for s in a.controls.split(",") if s]:
        left, right = spec.split(":")
        A, B = load_arms(a.root, left), load_arms(a.root, right)
        if kind == "pair":
            ladders.setdefault(LADDER_ALIASES.get(left, left), {}).update(
                {k: v for k, v in A.items() if k[1] == "interface"})
            ladders.setdefault(LADDER_ALIASES.get(right, right), {}).update(B)
        for ka, arm_a in sorted(A.items()):
            n = ka[0]
            # a control is matched with the interface, one-level ladder arm at the same N
            kb = next((k for k in B if k[0] == n and (kind == "pair" or (k[1] == "interface" and k[2] == 1))), None)
            if kb is None:
                notes.append(f"{spec}: no partner for N={n} ({ka[1]}, L={ka[2]})")
                continue
            arm_b = B[kb]
            if not (complete(arm_a) and complete(arm_b)):
                notes.append(f"{spec} N={n}: INCOMPLETE -- refined last t={arm_a['rows'][-1]['TIME']}, "
                             f"reference last t={arm_b['rows'][-1]['TIME']}; excluded")
                continue
            la = "refined" if kind == "pair" else f"control_{ka[1]}_L{ka[2]}"
            lb = "uniform" if kind == "pair" else "refined"
            rows.extend(compare(arm_a, arm_b, la, lb, spec, n, ka[1], ka[2]))

    # orders per ladder (only complete arms; needs >= 2 of them)
    for study, arms in sorted(ladders.items()):
        res = orders(arms, lambda arm: _f(arm["tok"]["DOMAIN_LENGTH"]) / _f(arm["tok"]["N_CELLS"]))
        for m, (p, r, npts) in res.items():
            order_rows.append({"ladder": study, "metric": m, "p": p, "R": r, "nPoints": npts})

    os.makedirs(os.path.join(a.root, a.outdir), exist_ok=True)
    if rows:
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(os.path.join(a.root, a.outdir, "refined_equivalence.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        for r in rows:
            vals = [f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys if k in r]
            print("  ".join(vals))
    # LaTeX: one row per production-force rung (refined vs uniform at matched N_CELLS),
    # relative differences at t = T for the metric vector plus the peak velocity norms and
    # the core-hour ratio. The well-balanced pair and the controls are not in this table.
    prod = [r for r in rows if r["pair"].startswith("stationaryDroplet3Drefined:stationaryDroplet3Duniform")
            and r.get("source") == "interface"]
    if prod:
        cols = [("max_t", "meanMagUPrime", "$\\max_t\\overline{|\\mathbf u|}$"),
                ("max_t", "l2MagUPrime", "$\\max_t\\lVert\\mathbf u\\rVert_2$"),
                ("T", "meanMagUPrime", "$\\overline{|\\mathbf u|}(T)$"),
                ("T", "l2MagUPrime", "$\\lVert\\mathbf u\\rVert_2(T)$"),
                ("T", "phaseVolumeRelError", "volume"), ("T", "zeroSetRadialL2", "shape"),
                ("T", "pLaplaceRelError", "$\\Delta p$"), ("T", "kErrL2Band", "$\\kappa$ err."),
                ("T", "A2hL2Band", "$A_{2h}$")]
        ns = sorted({r["N_CELLS"] for r in prod})
        with open(os.path.join(a.root, a.outdir, "refined_equivalence.tex"), "w") as fh:
            fh.write("% regenerated by workflow/scripts/make_refined_equivalence_table.py -- do not edit\n")
            fh.write("\\begin{tabular}{r" + "c" * (len(cols) + 2) + "}\n\\hline\n")
            fh.write("$N_{\\mathrm{fine}}$ & steps & " + " & ".join(c[2] for c in cols) + " & core-h ratio \\\\\n\\hline\n")
            for n in ns:
                rr = [r for r in prod if r["N_CELLS"] == n]
                # each rung's own final instant (they differ in the last digits between rungs)
                t_end = max(r["t"] for r in rr if isinstance(r["t"], float))
                def pick(when, m):
                    for r in rr:
                        if r["metric"] == m and ((when == "max_t" and r["t"] == "max_t") or
                                                 (when == "T" and isinstance(r["t"], float) and abs(r["t"] - t_end) < 1e-12)):
                            return r["relDiff"]
                    return None
                cells = []
                for when, m, _ in cols:
                    v = pick(when, m)
                    cells.append(f"{100 * v:+.2f}\\,\\%" if isinstance(v, float) else "--")
                steps = rr[0]["steps"]; ratio = rr[0].get("coreHoursRatio", "")
                fh.write(f"{n} & {steps} & " + " & ".join(cells) +
                         (f" & {ratio:.1f}" if isinstance(ratio, float) else " & --") + " \\\\\n")
            fh.write("\\hline\n\\end{tabular}\n")
    if order_rows:
        by_metric = {}
        for r in order_rows:
            by_metric.setdefault(r["metric"], {})[r["ladder"]] = r
        ladder_names = sorted({r["ladder"] for r in order_rows})
        with open(os.path.join(a.root, a.outdir, "refined_convergence_orders.tex"), "w") as fh:
            fh.write("% regenerated by workflow/scripts/make_refined_equivalence_table.py -- do not edit\n")
            fh.write("\\begin{tabular}{l" + "c" * len(ladder_names) + "}\n\\hline\n")
            fh.write("metric & " + " & ".join(n.replace("stationaryDroplet3D", "").replace("_", "\\_")
                                             for n in ladder_names) + " \\\\\n\\hline\n")
            for m, d in by_metric.items():
                cells = []
                for ln in ladder_names:
                    r = d.get(ln)
                    cells.append(f"{r['p']:.2f} ({r['R']:.3f}, n={r['nPoints']})" if r and r["p"] is not None else "--")
                fh.write(m.replace("_", "\\_") + " & " + " & ".join(cells) + " \\\\\n")
            fh.write("\\hline\n\\end{tabular}\n")
        for r in order_rows:
            print(f"order  {r['ladder']:34s} {r['metric']:20s} p={r['p']:.3f} R={r['R']:.4f} n={r['nPoints']}"
                  if r["p"] is not None else f"order  {r['ladder']} {r['metric']}: not fittable")
    for n in notes:
        print("[equivalence] " + n)
    if not rows and not order_rows:
        raise SystemExit("[equivalence] nothing complete to compare yet")


if __name__ == "__main__":
    main()
