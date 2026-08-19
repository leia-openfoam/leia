#!/usr/bin/env python3
"""Curate the coupled stationary-droplet ladders (2D and 3D) into one CSV.

Reads the per-arm droplet-metrics CSV and case_params.json of every study dir
given, takes the LAST row of each arm, and emits one row per arm with the four
metrics that must always be reported together -- max|U|, volume error, shape
error and min|grad psi| -- plus the WP8.1 delivered-curvature-variation split at
t = 0 (the cellCentreInverse property under refinement) and at the end.

Reporting all four is not stylistic. In this campaign a single headline metric has
already pointed the wrong way twice: the domain-size control's rejected 4R box had
a 26x BETTER volume error than the reference while its max|U| was 15x worse, and
here the 2D N=512 arm improves in volume and shape while its max|U| order falls
from 2.02 to 1.00.

Usage: make_wide_ladder_table.py --out <csv> STUDY_DIR [STUDY_DIR ...]
"""
import argparse
import csv
import json
import math
import os
import sys

METRICS = (("maxMagU_end", "maxMagU"),
           ("volumeRelError_end", "phaseVolumeRelError"),
           ("shapeL2_end", "zeroSetRadialL2"),
           ("minGradPsiBand_end", "minGradPsiBand"),
           ("pLaplace_end", "pLaplace"),
           ("driverAcrossSupportL2_end", "driverAcrossSupportL2"),
           ("driverAlongInterfaceL2_end", "driverAlongInterfaceL2"),
           ("A2hL2Band_end", "A2hL2Band"))
T0 = (("driverAcrossSupportL2_t0", "driverAcrossSupportL2"),
      ("driverAlongInterfaceL2_t0", "driverAlongInterfaceL2"),
      ("maxMagU_t0", "maxMagU"),
      ("A2hL2Band_t0", "A2hL2Band"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    for sdir in args.studies:
        study = os.path.basename(os.path.normpath(sdir))
        for name in sorted(os.listdir(sdir)):
            cdir = os.path.join(sdir, name)
            mp = os.path.join(cdir, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
            pp = os.path.join(cdir, "case_params.json")
            if not (os.path.isfile(mp) and os.path.isfile(pp)):
                continue
            tok = json.load(open(pp))["tokens"]
            hist = list(csv.DictReader(open(mp)))
            # A DIVERGED ARM'S LAST LINE IS TRUNCATED MID-WRITE: the solver dies
            # on a floating-point exception part-way through the row, so the final
            # record has None for the columns after the break. Scoring an arm on a
            # half-written row would silently mix a real value with a missing one,
            # so keep only fully-parseable rows and remember how many were dropped.
            need = [c for _, c in METRICS] + ["TIME"]
            full = [r for r in hist
                    if all(r.get(c) not in (None, "") for c in need)]
            truncated = len(hist) - len(full)
            hist = full
            if not hist:
                continue
            # Arms predating the DOMAIN_LENGTH token all ran on the 0.01 m
            # reference box -- that is what the token was defaulted to when it
            # was introduced, so the fallback is the historical value and not a
            # guess.
            L = float(tok.get("DOMAIN_LENGTH", 0.01))
            N = int(float(tok["N_CELLS"]))
            dim = 3 if "3D" in tok.get("SOLVER_CASE", study) or "3D" in name else 2
            h = L/N
            rec = dict(study=study, case=name, dim=dim, domainLength=L, LoverR=L/1e-3,
                       N=N, cells=N**dim, h=h, RoverH=1e-3/h,
                       dt=float(tok.get("MAX_DELTA_T", "nan")),
                       curvatureExtension=tok.get("CURVATURE_EXTENSION", ""),
                       psiFilter=tok.get("PSI_FILTER", ""),
                       psiFilterTheta=tok.get("PSI_FILTER_THETA", ""),
                       # The K switch is part of the SCHEME, so it must key the
                       # grouping: without it the K-on and K-off 3D ladders merge
                       # into one series with duplicated h and the order fit
                       # divides by a zero log-ratio.
                       gaussianCurvature=tok.get("CURVATURE_INVERSE_GAUSSIAN", "yes"),
                       np=tok.get("np", ""), nSteps=len(hist),
                       tEnd=float(hist[-1]["TIME"]),
                       reachedEndTime=int(float(hist[-1]["TIME"]) >= 0.999*float(tok["END_TIME"])),
                       truncatedRows=truncated)
            for out, col in METRICS:
                rec[out] = abs(float(hist[-1][col])) if col in hist[-1] else ""
            for out, col in T0:
                rec[out] = abs(float(hist[0][col])) if col in hist[0] else ""
            rows.append(rec)

    rows.sort(key=lambda r: (r["dim"], r["h"]*-1))
    cols = list(rows[0].keys())
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"[wide-ladder] wrote {args.out}: {len(rows)} arms")

    # orders per dimension, printed for the record (not written -- h ratios differ)
    # Group by (dim, delivery, filter, theta): an order fit across different
    # filter strengths would be meaningless, and theta = 0.05 and 0.2 arms of the
    # same ladder sit side by side in these studies.
    key = lambda r: (r["dim"], r["curvatureExtension"], r["psiFilter"],
                     r["psiFilterTheta"], r["gaussianCurvature"])
    for k in sorted({key(r) for r in rows}):
        sub = [r for r in rows if key(r) == k]
        if len(sub) < 2:
            continue
        dim, ext, filt, th, gk = k
        print(f"\n  {dim}D  {ext} + {filt} theta={th}  K={gk}  orders:")
        for a, b in zip(sub, sub[1:]):
            lr = math.log(a["h"]/b["h"])
            f = lambda k: (math.log(a[k]/b[k])/lr
                           if a[k] and b[k] else float("nan"))
            print(f"    N {a['N']:>4} -> {b['N']:<4} (R/h {a['RoverH']:.1f} -> {b['RoverH']:.1f})"
                  f"  max|U| {f('maxMagU_end'):+6.2f}  volume {f('volumeRelError_end'):+6.2f}"
                  f"  shape {f('shapeL2_end'):+6.2f}"
                  f"  driver(t=0) {f('driverAcrossSupportL2_t0'):+6.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
