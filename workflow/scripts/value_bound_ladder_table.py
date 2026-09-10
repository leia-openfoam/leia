#!/usr/bin/env python3
"""Resolution ladder for the semi-Lagrangian value bound: error AND observed order.

WHY THE ORDER COLUMN IS THE POINT. An error that is lower at one resolution is not a fix.
MEASURED 2026-09-10: the distance-cone bound lowered every interface error at N = 64 and its
eikonal error then moved 7.119e-03 -> 7.066e-03 across a 2x refinement -- an order of 0.01.
That is a FLOOR, not a converging error: the unbounded run converges at order 1.10 and
reaches the same value at about N = 256, after which the bound is the worse choice. A table
of errors alone would have shown a win at both rungs and hidden the crossover.

Arms are keyed by the (valueBound, lipschitzMode) pair READ FROM case_params.json, never by
the arm index: an index changes when an axis is added and the rows then silently swap.

Reports L1/L2-type metrics only. NEVER L_inf, which does not converge here.
Reports the geometric error, the volume error and the diagnostics TOGETHER -- a single
headline metric has already misread a candidate in this repository.

Usage:
  python3 workflow/scripts/value_bound_ladder_table.py \
      <study-at-N1> <study-at-N2> [<study-at-N3> ...] [--case popinetTranslating2D]
      [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os

METRICS = ["phaseVolumeRelError", "zeroSetRadialL2", "centroidError",
           "l2MagUPrime", "meanMagUPrime", "gradPsiL2ErrorBand"]


def arms(study, case):
    out = {}
    for d in sorted(glob.glob(os.path.join("studies", study, f"{case}_*/"))):
        pj = os.path.join(d, "case_params.json")
        if not os.path.exists(pj):
            continue
        t = json.load(open(pj))["tokens"]
        csvs = glob.glob(os.path.join(d, "leiaSemiLagrangian*TwoPhaseFoam.csv"))
        if not csvs:
            continue
        rows = list(csv.DictReader(open(csvs[0])))
        if not rows:
            continue
        key = (t.get("SL_VALUE_BOUND"), t.get("SL_CONE_L_MODE"))
        out[key] = (rows[-1], len(rows), int(t["N_CELLS"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--case", default="popinetTranslating2D")
    ap.add_argument("--csv")
    a = ap.parse_args()

    rungs = {}
    for s in a.studies:
        d = arms(s, a.case)
        if not d:
            print(f"  WARNING: no arms found in studies/{s}")
            continue
        N = next(iter(d.values()))[2]
        rungs[N] = d
    Ns = sorted(rungs)
    if len(Ns) < 2:
        raise SystemExit("a ladder needs at least two rungs")

    # An unequal step count makes the endpoints incomparable at a rung.
    for N in Ns:
        counts = {k: v[1] for k, v in rungs[N].items()}
        if len(set(counts.values())) != 1:
            print(f"  WARNING: N = {N} arms have UNEQUAL step counts {counts}"
                  f" -- endpoint metrics are not comparable there")

    keys = [k for k in rungs[Ns[0]] if all(k in rungs[N] for N in Ns)]
    out = []
    for m in METRICS:
        print(f"\n== {m} ==")
        print(f"{'arm':30s}" + "".join(f"{'N=' + str(N):>14s}" for N in Ns)
              + f"{'order':>9s}{'vs none':>11s}")
        base = {N: float(rungs[N][("none", "unity")][0][m]) for N in Ns} \
            if ("none", "unity") in keys else {}
        for k in keys:
            vals = [float(rungs[N][k][0][m]) for N in Ns]
            # order from the two FINEST rungs, which is where it matters
            o = (math.log(vals[-2] / vals[-1]) / math.log(Ns[-1] / Ns[-2])
                 if vals[-1] > 0 else float("nan"))
            rel = ((vals[-1] - base[Ns[-1]]) / abs(base[Ns[-1]]) * 100.0) if base else 0.0
            print(f"{k[0] + '/' + k[1]:30s}" + "".join(f"{v:14.4e}" for v in vals)
                  + f"{o:9.2f}{rel:+10.1f}%")
            out.append({"metric": m, "valueBound": k[0], "lipschitzMode": k[1],
                        **{f"N{N}": float(rungs[N][k][0][m]) for N in Ns},
                        "order": o, "relToNonePercent": rel})

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
