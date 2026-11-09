#!/usr/bin/env python3
"""The tables of the iterative face-fraction corrector (leia davof-iter, 2026-11-09;
docs/davof/REPORT-swept.md "The iterative face-fraction corrector", predictions P1-P5).

Per case (a study directory's <case>_NNNNN with case_params.json, advectorVoF.csv and
the solver log): the tokens (N, Co, arm, corrector, iterations, scheme, np), the end
values of the converter's CSV (L1(alphaf), L1(n) = LNormalDiff1, L1(x) = LCentreDiff1,
E_shape, min alpha, max alpha, the mass deviation, the interface cells, status) and from
the log the corrector's per-step lines

    davof corrector (mean): iterations K of N, interface cells (alpha) I, faces F,
    cells without a plane C; max |dalphaf| per iteration (a b c ...), L1 |dalphaf|
    (a b c ...); identity |sum_k m_k - sum_b alphaf S_f| / sum_k |m_k| = X

(iterations used: median and max over the steps; the change at iteration 1: the median
and the max over the steps of max |dalphaf| and L1 |dalphaf|; the contraction factors
change(i)/change(i+1) of the max change per iteration: the median and the worst over the
steps; the iterations to 1e-6 and 1e-10 of the max change: the median over the steps,
">k" when not reached within the k evaluated; the identity defect: the max over the
steps), the wall clock (the last ClockTime) and the steps. Per group (study, Co, arm,
corrector) the least-squares orders over N of L1(alphaf), E_shape, L1(n) and of the
first correction's max and L1 change.

    python3 workflow/scripts/iter_tables.py --out DIR STUDYDIR [STUDYDIR ...]

Writes DIR/iter_cases.csv, DIR/iter_orders.csv and DIR/iter_tables.md.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import statistics

CORR_RE = re.compile(
    r"davof corrector \((\w+)\): iterations (\d+) of (\d+), interface cells \(alpha\) (\d+), faces (\d+),"
    r" cells without a plane (\d+); max \|dalphaf\| per iteration \(([^)]*)\), L1 \|dalphaf\| \(([^)]*)\);"
    r" identity [^=]*= (\S+)")
IC_RE = re.compile(r"^davof: (\d+) interface cells")
CLOCK_RE = re.compile(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s")
TIME_RE = re.compile(r"^Time = (\S+)")


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def order(ns, es):
    pts = [(1.0 / n, e) for n, e in zip(ns, es) if e is not None and e > 0]
    if len(pts) < 2:
        return None
    lh = [math.log(p[0]) for p in pts]
    le = [math.log(p[1]) for p in pts]
    mh, me = sum(lh) / len(lh), sum(le) / len(le)
    sxx = sum((a - mh) ** 2 for a in lh)
    return sum((a - mh) * (b - me) for a, b in zip(lh, le)) / sxx if sxx > 0 else None


def med(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def parse_log(logf):
    """The corrector lines per step, the interface-cell count of the first reconstruction
    per step, the wall clock, the step count."""
    steps = []          # one dict per step with a corrector line
    icounts = []        # the first 'davof: N interface cells' per step
    clock = None
    nsteps = 0
    cur_ic = None
    with open(logf, errors="replace") as fh:
        for line in fh:
            if TIME_RE.match(line):
                nsteps += 1
                cur_ic = None
                continue
            m = IC_RE.match(line)
            if m and cur_ic is None and nsteps > 0:
                cur_ic = int(m.group(1))
                icounts.append(cur_ic)
                continue
            m = CORR_RE.search(line)
            if m:
                mx = [fnum(x) for x in m.group(7).split()]
                l1 = [fnum(x) for x in m.group(8).split()]
                steps.append(dict(used=int(m.group(2)), nmax=int(m.group(3)), icells=int(m.group(4)),
                                  faces=int(m.group(5)), noplane=int(m.group(6)), mx=mx, l1=l1,
                                  ident=fnum(m.group(9))))
                continue
            m = CLOCK_RE.search(line)
            if m:
                clock = fnum(m.group(2))
    return steps, icounts, clock, nsteps


def iters_to(mx, tol):
    """The first evaluation index (1-based) whose max change is below tol; None when
    not reached within the evaluated ones."""
    for i, v in enumerate(mx):
        if v is not None and v <= tol:
            return i   # i replacements made before the evaluation below tol
    return None


def measure(c):
    tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) else []
    last = rows[-1] if rows else {}
    col = lambda k: [fnum(r.get(k)) for r in rows if fnum(r.get(k)) is not None]
    logf = next((os.path.join(c, n) for n in ("log.davofAdvectorFoam", "log.advectorVoF")
                 if os.path.isfile(os.path.join(c, n))), None)
    steps, icounts, clock, nsteps = parse_log(logf) if logf else ([], [], None, 0)
    corr = str(tok.get("TPF_FACE_FRACTION_CORRECTOR", "none"))
    niter = int(corr[4:]) if corr.startswith("mean") and len(corr) > 4 else (int(tok.get("TPF_CORRECTOR_ITERATIONS", 1)) if corr != "none" else 0)
    out = dict(
        study=os.path.basename(os.path.dirname(c)), case=os.path.basename(c),
        scheme=tok.get("DAVOF_ADVECTION_SCHEME", "davofSwept"), row=tok.get("TPF_ROW"),
        arm=tok.get("TPF_VELOCITY_ARM", "-"), N=int(float(tok["N_CELLS"])),
        h=float(tok.get("DOMAIN_LENGTH", 1)) / float(tok["N_CELLS"]),
        Co=tok.get("TARGET_ADVECTIVE_CO"), corrector=corr, iterations=niter, np=tok.get("NP", "-"),
        grazingCut=tok.get("TPF_GRAZING_CUT", "-"),
        steps=nsteps, status=last.get("status", "-"), tEnd=fnum(last.get("TIME")),
        L1a_end=fnum(last.get("L1_alphaf")), L1n_end=fnum(last.get("LNormalDiff1")),
        L1x_end=fnum(last.get("LCentreDiff1")), Eshape_end=fnum(last.get("Eshape")),
        L1a_max=max(col("L1_alphaf")) if col("L1_alphaf") else None,
        L1n_max=max(col("LNormalDiff1")) if col("LNormalDiff1") else None,
        minAlpha=min(col("minAlphaBound")) if col("minAlphaBound") else None,
        maxAlphaExcess=max(col("maxAlphaBound")) if col("maxAlphaBound") else None,
        massDev=fnum(last.get("devPercent")), nInterface_end=fnum(last.get("nInterfaceAdv")),
        clock_s=clock,
    )
    # the corrector's statistics
    if steps:
        used = [s["used"] for s in steps]
        out.update(
            it_median=med(used), it_max=max(used),
            d1max_median=med([s["mx"][0] for s in steps if s["mx"]]),
            d1max_max=max(s["mx"][0] for s in steps if s["mx"]),
            d1L1_median=med([s["l1"][0] for s in steps if s["l1"]]),
            d1L1_max=max(s["l1"][0] for s in steps if s["l1"]),
            ident_max=max(s["ident"] for s in steps if s["ident"] is not None),
            noplane_max=max(s["noplane"] for s in steps),
            icells_median=med([s["icells"] for s in steps]),
            faces_median=med([s["faces"] for s in steps]),
        )
        # the contraction factors per iteration: mx[i]/mx[i+1], i = 0..3
        for i in range(4):
            fac = []
            for s in steps:
                mx = s["mx"]
                if len(mx) > i + 1 and mx[i] is not None and mx[i + 1] is not None and mx[i + 1] > 0 and mx[i] > 1e-14:
                    fac.append(mx[i] / mx[i + 1])
            out["fac%d_median" % (i + 1)] = med(fac)
            out["fac%d_worst" % (i + 1)] = min(fac) if fac else None
            out["fac%d_n" % (i + 1)] = len(fac)
        for tol, name in ((1e-6, "to1em6"), (1e-10, "to1em10")):
            vals = [iters_to(s["mx"], tol) for s in steps]
            reached = [v for v in vals if v is not None]
            out[name + "_median"] = med(reached)
            out[name + "_notReached"] = len(vals) - len(reached)
            out[name + "_max"] = max(reached) if reached else None
        # the alpha-based interface-cell count against the reconstruction's count
        out["icells_eq_log"] = all(s["icells"] == i for s, i in zip(steps, icounts[-len(steps):])) if icounts else None
    else:
        for k in ("it_median", "it_max", "d1max_median", "d1max_max", "d1L1_median", "d1L1_max", "ident_max",
                  "noplane_max", "icells_median", "faces_median"):
            out[k] = None
        for i in range(4):
            out["fac%d_median" % (i + 1)] = out["fac%d_worst" % (i + 1)] = None
            out["fac%d_n" % (i + 1)] = 0
        for name in ("to1em6", "to1em10"):
            out[name + "_median"] = out[name + "_max"] = None
            out[name + "_notReached"] = None
        out["icells_eq_log"] = None
    return out


def g(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        if x == 0:
            return "0"
        if abs(x) >= 1000 or abs(x) < 1e-3:
            return "%.*e" % (nd - 1, x)
        return "%.*g" % (nd + 1, x)
    return str(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cases = []
    for st in a.studies:
        for c in sorted(glob.glob(os.path.join(st, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            if os.path.isfile(os.path.join(c, "case_params.json")):
                cases.append(measure(c))
    cases.sort(key=lambda r: (r["study"], r["scheme"], r["arm"], str(r["Co"]), r["corrector"], r["N"]))
    keys = list(cases[0].keys()) if cases else []
    with open(os.path.join(a.out, "iter_cases.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in cases:
            w.writerow(r)

    # the orders over N per (case family, scheme, arm, Co, corrector)
    fam = lambda r: r["study"].rsplit("_N", 1)[0]
    groups = {}
    for r in cases:
        groups.setdefault((fam(r), r["scheme"], r["arm"], str(r["Co"]), r["corrector"]), []).append(r)
    orders = []
    for k, rs in sorted(groups.items()):
        rs = sorted(rs, key=lambda r: r["N"])
        ns = [r["N"] for r in rs]
        if len(set(ns)) < 2:
            continue
        orders.append(dict(family=k[0], scheme=k[1], arm=k[2], Co=k[3], corrector=k[4], Ns=" ".join(map(str, ns)),
                           order_L1a=order(ns, [r["L1a_end"] for r in rs]),
                           order_Eshape=order(ns, [r["Eshape_end"] for r in rs]),
                           order_L1n=order(ns, [r["L1n_end"] for r in rs]),
                           order_L1x=order(ns, [r["L1x_end"] for r in rs]),
                           order_d1max=order(ns, [r["d1max_median"] for r in rs]),
                           order_d1L1=order(ns, [r["d1L1_median"] for r in rs])))
    with open(os.path.join(a.out, "iter_orders.csv"), "w", newline="") as fh:
        if orders:
            w = csv.DictWriter(fh, fieldnames=list(orders[0].keys()))
            w.writeheader()
            for r in orders:
                w.writerow(r)

    # markdown
    md = []
    md.append("Notation: V_k cell volume, A_k cell surface area, |Sigma_k| interface area, m_k = -sum_f alpha_f S_f, "
              "alpha_k, alpha_f, S_f, x_k, x_f, p_k. d1 = the change of alpha_f at iteration 1 (the median over the steps; "
              "max over the faces / L1 over the set faces); it = iterations used per step (median / max); fac_i = the "
              "contraction factor max|dalpha_f|(i) / max|dalpha_f|(i+1) (median / worst over the steps); to 1e-6 / 1e-10 = "
              "the replacements before the max change fell below it (median; the steps not reaching it); ident = the "
              "identity defect |sum_k m_k - sum_b alpha_f S_f| / sum_k |m_k| (max over the steps); wall = ClockTime at "
              "the end (shared nodes unless stated).\n")
    md.append("### The cases\n")
    md.append("| study | scheme | arm | Co | N | h | np | corrector | it med/max | d1 max | d1 L1 | fac1 med/worst | fac2 | fac3 | "
              "to 1e-6 | to 1e-10 | L1(alphaf) end | L1(alphaf) max | E_shape | L1(x) | L1(n) | min alpha | max alpha - 1 | "
              "mass dev % | interface cells | no plane | ident | steps | wall s | status |")
    md.append("|" + "---|" * 31)
    for r in cases:
        md.append("| %s | %s | %s | %s | %d | %s | %s | %s | %s/%s | %s | %s | %s/%s | %s/%s | %s/%s | %s (%s) | %s (%s) | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d | %s | %s |" % (
            r["study"], r["scheme"], r["arm"], r["Co"], r["N"], g(r["h"]), r["np"], r["corrector"],
            g(r["it_median"]), g(r["it_max"]), g(r["d1max_median"]), g(r["d1L1_median"]),
            g(r["fac1_median"]), g(r["fac1_worst"]), g(r["fac2_median"]), g(r["fac2_worst"]), g(r["fac3_median"]), g(r["fac3_worst"]),
            g(r["to1em6_median"]), g(r["to1em6_notReached"]), g(r["to1em10_median"]), g(r["to1em10_notReached"]),
            g(r["L1a_end"]), g(r["L1a_max"]), g(r["Eshape_end"]), g(r["L1x_end"]), g(r["L1n_end"]),
            g(r["minAlpha"]), g(r["maxAlphaExcess"]), g(r["massDev"]), g(r["nInterface_end"]), g(r["noplane_max"]),
            g(r["ident_max"]), r["steps"], g(r["clock_s"]), r["status"]))
    md.append("\n### The orders over N\n")
    md.append("| family | scheme | arm | Co | corrector | N | order L1(alphaf) | order E_shape | order L1(n) | order L1(x) | order d1 max | order d1 L1 |")
    md.append("|" + "---|" * 12)
    for r in orders:
        md.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["family"], r["scheme"], r["arm"], r["Co"], r["corrector"], r["Ns"], g(r["order_L1a"]), g(r["order_Eshape"]),
            g(r["order_L1n"]), g(r["order_L1x"]), g(r["order_d1max"]), g(r["order_d1L1"])))
    with open(os.path.join(a.out, "iter_tables.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
