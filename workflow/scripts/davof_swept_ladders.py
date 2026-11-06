#!/usr/bin/env python3
"""The swept-face advection in the preprint's framing (2026-11-02; predictions
in the headers of config/davof/planeTranslation3D_hexSwept.yaml (G-1),
sphereRotation3D_hexRotFinal.yaml (R-0..R-4) and
sphereTranslation3D_hexSweptLadder.yaml (T-0..T-4)).

Per case, from the converter's advectorVoF.csv (one row per step): the end
values of L1(alphaf), L1(n) (LNormalDiff1), L1(x) (LCentreDiff1) and Eshape
(volumeFractionError; the fixed one when the .dat header names exactState
exactSphere), the run maxima of L1(alphaf) and L1(n) (the plane gate,
GATE_TOL 1e-12), min alpha and the largest excess over 1 after the bounding,
the mass deviation at the end, the rank fallbacks; per (study group, row,
arm, Co) the least-squares orders over N.

    python3 workflow/scripts/davof_swept_ladders.py --group NAME STUDY [STUDY ...] [--group ...] [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os


def f(x):
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


def measure(c):
    tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
    rows = list(csv.DictReader(open(os.path.join(c, "advectorVoF.csv"))))
    last = rows[-1]
    col = lambda k: [f(r.get(k)) for r in rows if f(r.get(k)) is not None]
    dat = os.path.join(c, "postProcessing", "volumeFractionError", "0", "volumeFractionError.dat")
    fixed = os.path.isfile(dat) and "exactState exactSphere" in open(dat).readline()
    logf = next((os.path.join(c, n) for n in ("log.advectorVoF", "log.davofAdvectorFoam") if os.path.isfile(os.path.join(c, n))), os.path.join(c, "log.advectorVoF"))  # either route (2026-11-06)
    ex = [l for l in open(logf) if l.startswith("ExecutionTime")]
    L1a, L1n = col("L1_alphaf"), col("LNormalDiff1")
    return dict(
        study=os.path.basename(os.path.dirname(c)), case=os.path.basename(c), row=tok.get("TPF_ROW"),
        arm=tok.get("TPF_VELOCITY_ARM"), N=int(float(tok["N_CELLS"])), Co=tok.get("TARGET_ADVECTIVE_CO"),
        np=tok.get("NP", "-"), steps=len(rows), status=last.get("status"), tEnd=f(last.get("TIME")),
        L1a_end=f(last.get("L1_alphaf")), L1n_end=f(last.get("LNormalDiff1")), L1x_end=f(last.get("LCentreDiff1")),
        Eshape_end=f(last.get("Eshape")), Eshape_fixed=fixed,
        L1a_max=max(L1a) if L1a else None, L1n_max=max(L1n) if L1n else None,
        minAlpha=min(col("minAlphaBound")) if col("minAlphaBound") else None,
        maxExcess=max(col("maxAlphaBound")) if col("maxAlphaBound") else None,
        devPercent_end=f(last.get("devPercent")),
        rankFallbacks=int(sum(col("nSweptRankFallback"))) if col("nSweptRankFallback") else 0,
        cpu_s=f(ex[-1].split()[2]) if ex else None,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", nargs="+", action="append", required=True)
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    out = []
    for g in a.group:
        name, studies = g[0], g[1:]
        recs = []
        for st in studies:
            for c in sorted(glob.glob(os.path.join(st, "*_0*"))):
                if "." in os.path.basename(c) or not os.path.isfile(os.path.join(c, "advectorVoF.csv")):
                    continue
                r = measure(c)
                r["group"] = name
                recs.append(r)
        keyf = lambda r: (r["row"], r["arm"], r["Co"])
        for key in sorted({keyf(r) for r in recs}, key=str):
            sel = sorted([r for r in recs if keyf(r) == key], key=lambda r: r["N"])
            ok = [r for r in sel if r["status"] == "completed"]
            ns = [r["N"] for r in ok]
            for m in ("L1a_end", "L1n_end", "L1x_end", "Eshape_end"):
                p = order(ns, [r[m] for r in ok]) if len(ok) >= 2 else None
                for r in sel:
                    r["p_" + m] = p
        out += recs
    keys = ["group", "row", "arm", "Co", "N", "status", "steps", "tEnd", "L1a_end", "p_L1a_end", "L1n_end", "p_L1n_end",
            "L1x_end", "p_L1x_end", "Eshape_end", "p_Eshape_end", "Eshape_fixed", "L1a_max", "L1n_max", "minAlpha",
            "maxExcess", "devPercent_end", "rankFallbacks", "cpu_s", "study", "case"]
    fmt = lambda v: ("%.4e" % v) if isinstance(v, float) else str(v)
    print("\t".join(keys))
    for r in sorted(out, key=lambda r: (r["group"], str(r["row"]), str(r["arm"]), str(r["Co"]), r["N"])):
        print("\t".join(fmt(r.get(k)) for k in keys))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
