#!/usr/bin/env python3
"""The straight translations of sphereRotation3D against the confound of
the per-axis displacement per step (config/davof/sphereRotation3D_hexFineLine*,
2026-11-02).

Per case: the direction, dt, the displacement per step along each mesh axis
in cells (|U_k| dt N, the largest over k), the steps, the LATE SLOPE (the
least-squares slope of L1(alphaf) against t over [0.25 T, T], T = END_TIME;
the column rate of davof_rotation_history.py), the mean of L1(alphaf) over
t in [0.1, 0.45], L1(alphaf), L1(x) and Eshape at the end, min(alpha), the
largest excess over 1 and the mass deviation at the end (percent), and the
log-position of the late slope between the reference x-line and diagonal of
the same N,

    x = ln(s / s_X) / ln(s_D / s_X),

with s_X, s_D the late slopes of sphereRotation3D_hexFineLineX and
_hexFineLineDiag (--ref-x, --ref-d): x >= 0.75 D-like (accumulates as the
diagonal), x <= 0.25 X-like (levels off as the x-line), between AMBIGUOUS.

    python3 workflow/scripts/davof_line_confound.py --ref-x studies/davof/sphereRotation3D_hexFineLineX --ref-d studies/davof/sphereRotation3D_hexFineLineDiag studies/davof/<study> [...] [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os

import numpy as np


def f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def cases(study):
    for c in sorted(glob.glob(os.path.join(study, "*_0*"))):
        if "." in os.path.basename(c) or not os.path.isfile(os.path.join(c, "case_params.json")):
            continue
        p = os.path.join(c, "advectorVoF.csv")
        if not (os.path.isfile(p) and os.path.getsize(p)):
            continue
        yield c


def measure(case):
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    rows = list(csv.DictReader(open(os.path.join(case, "advectorVoF.csv"))))
    # a step whose time prints differently in two postProcessing files (an
    # irrational fixed dt) is split by the converter into two rows, each
    # with the other's columns empty: every column is read from its
    # non-empty entries
    col = lambda k: [f(r.get(k)) for r in rows if f(r.get(k)) is not None]
    lastOf = lambda k: (col(k) or [None])[-1]
    tl = [(f(r["TIME"]), f(r["L1_alphaf"])) for r in rows if f(r.get("L1_alphaf")) is not None]
    t = np.array([a for a, _ in tl])
    L = np.array([b for _, b in tl])
    T = f(tok["END_TIME"])
    N = int(float(tok["N_CELLS"]))
    dt = f(tok["MAX_DELTA_T"])
    U = [f(tok.get(k)) or 0.0 for k in ("TPF_UX", "TPF_UY", "TPF_UZ")]
    m = (t >= 0.25 * T) & (t <= T * (1 + 1e-9))
    slope = float(np.polyfit(t[m], L[m], 1)[0])
    w = (t >= 0.1 - 1e-12) & (t <= 0.45 + 1e-12)
    return dict(
        study=os.path.basename(os.path.dirname(case)), case=os.path.basename(case), N=N,
        arm=tok.get("TPF_VELOCITY_ARM"), Co=f(tok.get("TARGET_ADVECTIVE_CO")), dt=dt,
        axisDisp=max(abs(u) for u in U) * dt * N, speedDisp=math.sqrt(sum(u * u for u in U)) * dt * N,
        steps=len(t), tEnd=float(t[-1]), lateSlope=slope, meanL1=float(L[w].mean()),
        L1_end=float(L[-1]), L1_max=float(L.max()), L1x_end=lastOf("LCentreDiff1"),
        Es_end=lastOf("Eshape"),
        minAlpha=min(col("minAlphaBound")),
        maxExcess=max(col("maxAlphaBound")),
        devPercent_end=lastOf("devPercent"),
        fallbacks=int(sum(col("nSweptRankFallback"))),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-x", required=True)
    ap.add_argument("--ref-d", required=True)
    ap.add_argument("--csv", default="")
    ap.add_argument("studies", nargs="+")
    a = ap.parse_args()
    sX = {r["N"]: r["lateSlope"] for r in map(measure, cases(a.ref_x))}
    sD = {r["N"]: r["lateSlope"] for r in map(measure, cases(a.ref_d))}
    out = []
    for st in [a.ref_x, a.ref_d] + a.studies:
        for c in cases(st):
            r = measure(c)
            if r["N"] in sX and r["N"] in sD and r["lateSlope"] > 0:
                x = math.log(r["lateSlope"] / sX[r["N"]]) / math.log(sD[r["N"]] / sX[r["N"]])
                r["x"] = x
                r["class"] = "D-like" if x >= 0.75 else ("X-like" if x <= 0.25 else "AMBIGUOUS")
            out.append(r)
    keys = ["study", "N", "arm", "Co", "dt", "axisDisp", "speedDisp", "steps", "tEnd", "lateSlope", "x",
            "class", "meanL1", "L1_end", "L1_max", "L1x_end", "Es_end", "minAlpha", "maxExcess",
            "devPercent_end", "fallbacks"]
    fmt = lambda v: ("%.4g" % v) if isinstance(v, float) else str(v)
    print("\t".join(keys))
    for r in out:
        print("\t".join(fmt(r.get(k)) for k in keys))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
