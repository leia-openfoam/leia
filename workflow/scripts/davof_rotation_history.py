#!/usr/bin/env python3
"""The growth of the face-fraction error over the rotation turn
(config/davof/sphereRotation3D_hexHist*.yaml, 2026-11-01): per case the
step history of advectorVoF.csv (tpf_advection_metrics.py), and

  - L1(alphaf), L1(x), Eshape at 25 / 50 / 75 / 100 % of the run's end time
    (the turn; for a translation reference its own end time),
  - the growth law: the increment D(t) = L1(t) - L1(t_1) fitted as b t^p by
    least squares in log-log over t in [0.1 T, T] and over the second half
    [0.5 T, T] (p ~ 1: a fixed remainder accumulates; p ~ 0.5: random;
    p ~ 0: saturation),
  - the mean increment per unit time over [0.25 T, T] (the slope of the
    least-squares line of L1 against t) and its h-scaling between the two
    resolutions of a study (log2 of the ratio N = 32 over N = 64),
  - the bounds and the mass: min(alpha) after the bounding, the largest
    excess over 1, the deviation of sum(alpha V) at the end,
  - the abort step where the run diverged.

Usage (on the cluster, from the leia clone):
    python3 workflow/scripts/davof_rotation_history.py studies/davof/<study> [...] [--csv out.csv]
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


def load(case):
    t = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    path = os.path.join(case, "advectorVoF.csv")
    rows = list(csv.DictReader(open(path))) if os.path.isfile(path) and os.path.getsize(path) else []
    return t, rows


def series(rows, col):
    tt, vv = [], []
    for r in rows:
        a, b = f(r.get("TIME")), f(r.get(col))
        if a is not None and b is not None and b >= 0:
            tt.append(a); vv.append(b)
    return np.array(tt), np.array(vv)


def at_fraction(t, v, T, frac):
    if len(t) == 0 or t[-1] < frac * T * (1 - 1e-9):
        return None
    return float(v[np.argmin(np.abs(t - frac * T))])


def growth(t, v, T, lo):
    """b, p of D(t) = v - v[0] ~ b t^p over t in [lo T, T] (log-log least squares)."""
    if len(t) < 5:
        return None, None
    D = v - v[0]
    m = (t >= lo * T) & (t <= T * (1 + 1e-9)) & (D > 0)
    if m.sum() < 5:
        return None, None
    p, lnb = np.polyfit(np.log(t[m]), np.log(D[m]), 1)
    return math.exp(lnb), p


def rate(t, v, T):
    m = (t >= 0.25 * T) & (t <= T * (1 + 1e-9))
    if m.sum() < 5:
        return None
    return float(np.polyfit(t[m], v[m], 1)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    out = []
    for st in a.studies:
        for case in sorted(glob.glob(os.path.join(st, "*_0*"))):
            if not os.path.isfile(os.path.join(case, "case_params.json")) or "." in os.path.basename(case):
                continue
            tok, rows = load(case)
            if not rows:
                continue
            T = f(tok.get("END_TIME"))
            t, L1 = series(rows, "L1_alphaf")
            tx, Lx = series(rows, "LCentreDiff1")
            te, Es = series(rows, "Eshape")
            mins = [f(r.get("minAlphaBound")) for r in rows if f(r.get("minAlphaBound")) is not None]
            maxs = [f(r.get("maxAlphaBound")) for r in rows if f(r.get("maxAlphaBound")) is not None]
            dev = f(rows[-1].get("devPercent"))
            status = rows[-1].get("status", "")
            b1, p1 = growth(t, L1, T, 0.1)
            b2, p2 = growth(t, L1, T, 0.5)
            rec = dict(study=os.path.basename(st), case=os.path.basename(case), row=tok.get("TPF_ROW"),
                       flux=tok.get("TPF_FLUX_INTEGRATION"), plane=tok.get("TPF_PLANE_SOURCE", ""),
                       flow=tok.get("TPF_SETFLOW_MODE", "rotation"),
                       N=int(float(tok["N_CELLS"])), Co=tok.get("TARGET_ADVECTIVE_CO"), T=T,
                       steps=len(rows), status=status, abortStep=rows[-1].get("abortStep", ""),
                       L1_first=float(L1[0]) if len(L1) else None,
                       **{f"L1_{int(100*q)}": at_fraction(t, L1, T, q) for q in (0.25, 0.5, 0.75, 1.0)},
                       **{f"L1x_{int(100*q)}": at_fraction(tx, Lx, T, q) for q in (0.25, 0.5, 0.75, 1.0)},
                       **{f"Eshape_{int(100*q)}": at_fraction(te, Es, T, q) for q in (0.25, 0.5, 0.75, 1.0)},
                       p_full=p1, b_full=b1, p_half=p2, b_half=b2, rate=rate(t, L1, T),
                       minAlpha=min(mins) if mins else None, maxExcess=max(maxs) if maxs else None,
                       devPercent_end=dev)
            out.append(rec)
    # the h-scaling of the rate within a study (and row/flux/Co)
    for r in out:
        r["rate_order"] = None
        if r["N"] == 32:
            mate = [s for s in out if s["study"] == r["study"] and s["row"] == r["row"] and s["flux"] == r["flux"]
                    and s["Co"] == r["Co"] and s["plane"] == r["plane"] and s["flow"] == r["flow"]
                    and s["N"] == 64]
            if mate and r["rate"] and mate[0]["rate"] and r["rate"] > 0 and mate[0]["rate"] > 0:
                r["rate_order"] = math.log2(r["rate"] / mate[0]["rate"])
    fmt = lambda v: ("%.3e" % v) if isinstance(v, float) else str(v)
    keys = ["study", "N", "Co", "flow", "flux", "plane", "steps", "status", "L1_first", "L1_25", "L1_50", "L1_75", "L1_100",
            "p_full", "p_half", "rate", "rate_order", "L1x_100", "Eshape_100", "minAlpha", "maxExcess", "devPercent_end"]
    print("\t".join(keys))
    for r in out:
        print("\t".join(fmt(r.get(k)) for k in keys))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
