#!/usr/bin/env python3
"""Global (least-squares) convergence orders from a curated ladder CSV.

WHY GLOBAL AND NOT PAIRWISE. A pairwise order log(e_i/e_{i+1})/log(h_i/h_{i+1})
uses two points and inherits the full noise of both, so it is unstable and can
exceed the scheme's formal order outright -- this ladder produced +2.43 and +7.59
for a method that is at most second order. Those are not orders; they are the
signature of a coarse point that is not yet in the asymptotic regime, where the
error is dominated by a term that dies faster than the asymptotic one. A
least-squares fit of log(e) = log(C) + p log(h) over the whole ladder uses every
point once, reports a residual that says whether a single power law even
describes the data, and cannot be inflated by one bad pair.

Both are printed, together with R^2 and the fit excluding the coarsest point --
the honest way to show that the exponent depends on which points you trust.

Note on reading orders: p is defined against the CELL SIZE h, and halving h
multiplies the error by 2^-p. A factor-2 error drop per h-halving is p = 1; a
factor-4 drop is p = 2. Dimension does not enter this. It matters only if one
indexes by cell COUNT instead: cells ~ h^-d, so a second-order method drops the
error 4x per 4x cells in 2D and 4x per 8x cells in 3D.

Usage: convergence_order_fit.py <curated.csv> [--metric COL ...] [--drop-coarsest]
"""
import argparse
import csv
import math
import sys

DEFAULT_METRICS = ("maxMagU_end", "volumeRelError_end", "shapeL2_end",
                   "driverAcrossSupportL2_t0")


def fit(hs, es):
    """Least squares p in log(e) = log(C) + p log(h); returns (p, R^2, n)."""
    pts = [(math.log(h), math.log(e)) for h, e in zip(hs, es) if e > 0]
    n = len(pts)
    if n < 2:
        return float("nan"), float("nan"), n
    mx = sum(x for x, _ in pts)/n
    my = sum(y for _, y in pts)/n
    sxx = sum((x-mx)**2 for x, _ in pts)
    sxy = sum((x-mx)*(y-my) for x, y in pts)
    if sxx == 0:
        return float("nan"), float("nan"), n
    p = sxy/sxx
    c = my - p*mx
    ss_res = sum((y - (c + p*x))**2 for x, y in pts)
    ss_tot = sum((y - my)**2 for _, y in pts)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else float("nan")
    return p, r2, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--metric", action="append", default=None)
    args = ap.parse_args()
    metrics = args.metric or list(DEFAULT_METRICS)

    rows = list(csv.DictReader(open(args.csv)))
    groups = {}
    for r in rows:
        key = (int(r["dim"]), r["curvatureExtension"], r["psiFilter"],
               r["psiFilterTheta"], r["domainLength"])
        groups.setdefault(key, []).append(r)

    for key, sub in sorted(groups.items()):
        sub.sort(key=lambda r: -float(r["h"]))
        if len(sub) < 3:
            continue
        dim, ext, filt, th, L = key
        print(f"\n{dim}D  {ext} + {filt} theta={th}  L={L}   "
              f"N = {', '.join(r['N'] for r in sub)}  "
              f"(R/h = {', '.join(f'{float(r[chr(82)+chr(111)+chr(118)+chr(101)+chr(114)+chr(72)]):.1f}' for r in sub)})")
        hs = [float(r["h"]) for r in sub]
        print(f"  {'metric':>26} {'global p':>9} {'R^2':>7} {'n':>3} "
              f"{'p w/o coarsest':>15} {'R^2':>7}   pairwise")
        for m in metrics:
            try:
                es = [abs(float(r[m])) for r in sub]
            except (KeyError, ValueError):
                continue
            p, r2, n = fit(hs, es)
            p2, r22, n2 = fit(hs[1:], es[1:])
            pw = []
            for i in range(len(hs)-1):
                if es[i] > 0 and es[i+1] > 0:
                    pw.append(math.log(es[i]/es[i+1])/math.log(hs[i]/hs[i+1]))
            print(f"  {m:>26} {p:>+9.2f} {r2:>7.4f} {n:>3} {p2:>+15.2f} {r22:>7.4f}   "
                  + " ".join(f"{v:+5.2f}" for v in pw))
    return 0


if __name__ == "__main__":
    sys.exit(main())
