#!/usr/bin/env python3
"""The tables of the normal-order study (leia davof-normal-order, 2026-11-12;
docs/davof/REPORT-swept.md "Normal order versus normal bias", predictions N-1 .. N-5),
and the 1 - n.n_e column for any finished DAVOF advection study.

Per case (a study directory's <case>_NNNNN with case_params.json and advectorVoF.csv),
the end row of the converter's CSV:
  the angle theta = 2 asin(|n - n_e|/2): L1 and Linf (LNormalAngle1/Inf);
  1 - n.n_e: L1 and Linf (LNormalOneMinusDot1/Inf). For a case run before 2026-11-12
    the two are DERIVED exactly from the stored |n - n_e| columns: 1 - n.n_e =
    |n - n_e|^2/2 for unit vectors, so L1(1 - n.n_e) = L2(|n - n_e|)^2/2 and
    Linf(1 - n.n_e) = Linf(|n - n_e|)^2/2; Linf(theta) = 2 asin(Linf|n - n_e|/2); L1(theta)
    is bounded, L1|n - n_e| <= L1(theta) <= L1|n - n_e| * 2 asin(m/2)/m with m the Linf
    (the column `angleSource` says direct or derived);
  |n - n_e| (LNormalDiff1/2/Inf, the existing L1(n) column);
  the plicRDF benchmark's L_inf^n (Scheufler and Roenby 2019, Sec. 4.1, Eq. (21)) =
    Linf(1 - n.n_e);
  the position error L1(d) (LCentreDiff1);
  the per-step flux error (fluxErrL1: the max and the median over the steps; fluxErrInf the
    max) where the run logged it;
  E_shape, the mass deviation %, min alpha and the excess over 1 over the run, the
  interface cells, the steps, the ClockTime at the end, the status.
Per group (study family, row, Co, np) the least-squares and the pairwise orders over N of
every error, and the four-order table of the predictions with the verdicts.

    python3 workflow/scripts/nor_tables.py --out DIR STUDYDIR [STUDYDIR ...]

Writes DIR/nor_cases.csv, DIR/nor_orders.csv, DIR/nor_tables.md.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import statistics


def f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def ls_order(ns, es):
    pts = [(1.0 / n, e) for n, e in zip(ns, es) if e is not None and e > 0]
    if len(pts) < 2:
        return None
    lh = [math.log(p[0]) for p in pts]
    le = [math.log(p[1]) for p in pts]
    mh, me = sum(lh) / len(lh), sum(le) / len(le)
    sxx = sum((a - mh) ** 2 for a in lh)
    return sum((a - mh) * (b - me) for a, b in zip(lh, le)) / sxx if sxx > 0 else None


def pair_orders(ns, es):
    out = []
    for i in range(len(ns) - 1):
        a, b = es[i], es[i + 1]
        if a and b and a > 0 and b > 0 and ns[i + 1] != ns[i]:
            out.append(math.log(a / b) / math.log(ns[i + 1] / ns[i]))
        else:
            out.append(None)
    return out


def family(study):
    """The study name without its resolution: standing_shear3D_N32_variants ->
    standing_shear3D_variants, port_ladder_sphereN128_leia -> port_ladder_sphere_leia,
    nor_trans_N32_np1 -> nor_trans_np1."""
    return re.sub(r"_?N\d+(?=_|$)", "", study)


def measure(c):
    tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
    rows = list(csv.DictReader(open(os.path.join(c, "advectorVoF.csv"))))
    if not rows:
        return None
    last = rows[-1]
    # the end row with the reconstruction error (the last row that has it)
    rec = next((r for r in reversed(rows) if f(r.get("LNormalDiff1")) is not None), {})
    col = lambda k: [f(r.get(k)) for r in rows if f(r.get(k)) is not None]
    d1, d2, dinf = f(rec.get("LNormalDiff1")), f(rec.get("LNormalDiff2")), f(rec.get("LNormalDiffInf"))
    a1, ainf = f(rec.get("LNormalAngle1")), f(rec.get("LNormalAngleInf"))
    o1, oinf = f(rec.get("LNormalOneMinusDot1")), f(rec.get("LNormalOneMinusDotInf"))
    src = "direct"
    a1_hi = None
    if o1 is None and d2 is not None:
        src = "derived"
        o1 = 0.5 * d2 * d2
        oinf = 0.5 * dinf * dinf if dinf is not None else None
        ainf = 2 * math.asin(min(0.5 * dinf, 1.0)) if dinf is not None else None
        if d1 is not None and dinf:
            a1 = d1
            a1_hi = d1 * (2 * math.asin(min(0.5 * dinf, 1.0)) / dinf)
    elif o1 is None:
        src = "missing"
    # the identity check of the direct columns: L1(1 - n.n_e) = L2^2/2
    ident = None
    if src == "direct" and d2 and o1:
        ident = abs(o1 - 0.5 * d2 * d2) / o1
    fe = col("fluxErrL1")
    clock = col("clockTime")
    fw = col("fluxErrWall")
    out = dict(
        study=os.path.basename(os.path.dirname(c)), case=os.path.basename(c),
        family=family(os.path.basename(os.path.dirname(c))),
        row=tok.get("TPF_ROW"), arm=tok.get("TPF_VELOCITY_ARM", "-"),
        corrector=tok.get("TPF_FACE_FRACTION_CORRECTOR", "none"),
        enright=tok.get("TPF_ENRIGHT_SETTING", "-"),
        Co=tok.get("TARGET_ADVECTIVE_CO"), N=int(float(tok["N_CELLS"])), np=tok.get("NP", "1"),
        rep=tok.get("TPF_REPETITION", "1"),
        steps=len(rows), status=last.get("status"),
        angle1=a1, angle1_hi=a1_hi, angleInf=ainf, omd1=o1, omdInf=oinf, angleSource=src,
        srLinf=oinf, nd1=d1, nd2=d2, ndInf=dinf, identity=ident,
        pos1=f(rec.get("LCentreDiff1")), posInf=f(rec.get("LCentreDiffInf")),
        flux1Max=max(fe) if fe else None, flux1Med=statistics.median(fe) if fe else None,
        fluxInfMax=max(col("fluxErrInf")) if col("fluxErrInf") else None,
        fluxWall=sum(fw) if fw else None,
        Eshape=f(next((r.get("Eshape") for r in reversed(rows) if f(r.get("Eshape")) is not None), None)),
        L1af=f(next((r.get("L1_alphaf") for r in reversed(rows) if f(r.get("L1_alphaf")) is not None), None)),
        massDev=f(last.get("devPercent")),
        minAlpha=min(col("minAlphaBound")) if col("minAlphaBound") else (min(col("minAlpha")) if col("minAlpha") else None),
        maxAlphaExcess=max(col("maxAlphaBound")) if col("maxAlphaBound") else None,
        nInterface=f(rec.get("nInterface")),
        clock=clock[-1] if clock else None,
    )
    return out


QUANT = [("angle1", "L1(theta)"), ("angleInf", "Linf(theta)"), ("omd1", "L1(1-n.ne)"),
         ("omdInf", "Linf(1-n.ne) = SR L_inf^n"), ("nd1", "L1|n-ne|"), ("nd2", "L2|n-ne|"),
         ("ndInf", "Linf|n-ne|"), ("pos1", "L1(d)"), ("flux1Max", "flux L1 max/step"),
         ("Eshape", "E_shape"), ("L1af", "L1(alphaf)")]


def g(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    if x == 0:
        return "0"
    if abs(x) >= 1000 or abs(x) < 1e-2:
        return "%.*e" % (nd - 1, x)
    return "%.*g" % (nd, x)


def go(x):
    return "-" if x is None else "%.2f" % x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cases = []
    for st in a.studies:
        for c in sorted(glob.glob(os.path.join(st, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            if os.path.isfile(os.path.join(c, "advectorVoF.csv")) and os.path.isfile(os.path.join(c, "case_params.json")):
                m = measure(c)
                if m:
                    cases.append(m)
    keys = list(cases[0].keys()) if cases else []
    cases.sort(key=lambda r: (r["family"], r["row"], r["arm"], r["enright"], r["corrector"], str(r["Co"]), r["np"], r["N"]))
    with open(os.path.join(a.out, "nor_cases.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in cases:
            w.writerow(r)

    groups = {}
    for r in cases:
        groups.setdefault((r["family"], r["row"], r["arm"], r["enright"], r["corrector"], str(r["Co"]), r["np"]), []).append(r)
    orders = []
    for k, rs in sorted(groups.items()):
        # one case per N: the repetition 1 (the repetitions of the standing matrix
        # are identical in their errors; their wall clock is in the case table)
        byN = {}
        for r in sorted(rs, key=lambda r: (r["N"], str(r["rep"]))):
            byN.setdefault(r["N"], r)
        rs = [byN[n] for n in sorted(byN)]
        ns = [r["N"] for r in rs]
        if len(set(ns)) < 2:
            continue
        o = dict(family=k[0], row=k[1], arm=k[2], enright=k[3], corrector=k[4], Co=k[5], np=k[6],
                 Ns=" ".join(map(str, ns)), statuses=" ".join(r["status"] or "-" for r in rs))
        for q, _ in QUANT:
            es = [r[q] for r in rs]
            o["ls_" + q] = ls_order(ns, es)
            o["pw_" + q] = " ".join(go(x) for x in pair_orders(ns, es))
        # the verdicts of N-1 .. N-4 (p = the L1 order of the angle)
        p = o["ls_angle1"]
        if p is not None:
            pos, flx, sh, omd = o["ls_pos1"], o["ls_flux1Max"], o["ls_Eshape"], o["ls_omd1"]
            o["N1"] = "-" if pos is None else ("held" if abs(pos - (p + 1)) <= 0.3 else "failed")
            o["N2"] = "-" if flx is None else ("held" if abs(flx - (p + 1)) <= 0.3 else "failed")
            if sh is None:
                o["N3"] = "-"
            else:
                d = sh - p
                cls = "near p+1" if d >= 0.6 else ("near p" if d <= 0.4 else "undecided")
                want = "near p" if "youngs" in (o["row"] or "") else "near p+1"
                o["N3"] = "%s (q-p %.2f): %s" % (cls, d, "held" if cls == want else ("undecided" if cls == "undecided" else "failed"))
            o["N4"] = "-" if omd is None else ("held" if abs(omd - 2 * p) <= 0.3 else "failed")
        orders.append(o)
    with open(os.path.join(a.out, "nor_orders.csv"), "w", newline="") as fh:
        if orders:
            fields = []
            for o in orders:
                for kk in o:
                    if kk not in fields:
                        fields.append(kk)
            w = csv.DictWriter(fh, fieldnames=fields, restval="")
            w.writeheader()
            for r in orders:
                w.writerow(r)

    md = ["Notation: V_k cell volume, A_k cell surface area, |Sigma_k| interface area, m_k = -sum_f alpha_f S_f, "
          "alpha_k, alpha_f, S_f, x_k, x_f, p_k; n_k the unit normal of the positioned plane, n_e the exact unit "
          "normal at its centroid x_k, theta_k = 2 asin(|n_k - n_e|/2); 1 - n.n_e = |n - n_e|^2/2. SR L_inf^n = the "
          "plicRDF benchmark's normal error (Scheufler and Roenby 2019, Sec. 4.1, Eq. (21)) = Linf(1 - n.n_e). "
          "Derived: the case predates the direct columns; 1 - n.n_e from the stored L2 and Linf of |n - n_e| "
          "(exact), L1(theta) between the two printed bounds. Wall = ClockTime at the end (shared nodes).\n",
          "### The cases\n",
          "| family | row | arm | Co | np | N | L1(theta) | Linf(theta) | L1(1-n.ne) | Linf(1-n.ne) = SR | L1|n-ne| | L2|n-ne| | L1(d) | flux L1 max (med) | E_shape | L1(alphaf) | mass dev % | min alpha | max alpha - 1 | steps | wall s | status | source |",
          "|" + "---|" * 23]
    for r in cases:
        a1 = g(r["angle1"]) if r["angle1_hi"] is None else "%s..%s" % (g(r["angle1"]), g(r["angle1_hi"]))
        fl = "-" if r["flux1Max"] is None else "%s (%s)" % (g(r["flux1Max"]), g(r["flux1Med"]))
        md.append("| %s | %s | %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d | %s | %s | %s |" % (
            r["family"], r["row"], r["arm"], r["Co"], r["np"], r["N"], a1, g(r["angleInf"]), g(r["omd1"]), g(r["omdInf"]),
            g(r["nd1"]), g(r["nd2"]), g(r["pos1"]), fl, g(r["Eshape"]), g(r["L1af"]), g(r["massDev"]), g(r["minAlpha"]),
            g(r["maxAlphaExcess"]), r["steps"], g(r["clock"]), r["status"], r["angleSource"]))
    md.append("\n### The orders over N (least squares; pairwise in parentheses)\n")
    md.append("| family | row | Co | np | N | " + " | ".join(lbl for _, lbl in QUANT) + " |")
    md.append("|" + "---|" * (5 + len(QUANT)))
    for o in orders:
        md.append("| %s | %s | %s | %s | %s | " % (o["family"], o["row"], o["Co"], o["np"], o["Ns"]) +
                  " | ".join("%s (%s)" % (go(o["ls_" + q]), o["pw_" + q]) for q, _ in QUANT) + " |")
    md.append("\n### The four orders and the verdicts (p = the L1 order of the angle)\n")
    md.append("| family | row | Co | np | p | position (N-1) | flux (N-2) | E_shape (N-3) | 1-n.ne (N-4) |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    for o in orders:
        if "N1" not in o:
            continue
        md.append("| %s | %s | %s | %s | %s | %s %s | %s %s | %s %s | %s %s |" % (
            o["family"], o["row"], o["Co"], o["np"], go(o["ls_angle1"]), go(o["ls_pos1"]), o["N1"],
            go(o["ls_flux1Max"]), o["N2"], go(o["ls_Eshape"]), o["N3"], go(o["ls_omd1"]), o["N4"]))
    with open(os.path.join(a.out, "nor_tables.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
