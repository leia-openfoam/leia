#!/usr/bin/env python3
"""Compare two solver metrics CSVs step by step (serial vs parallel, before vs after).

For every common numeric column the maximum over rows of |a - b| / max(|a|, |b|, floor)
is reported; PASS when every column's maximum is below the tolerance. Row counts must
agree (unequal step counts are not comparable). Columns that are identically zero in both
files pass trivially; `floor` keeps a round-off-level column (e.g. a parasitic velocity of
1e-10 m/s) from being judged on its relative noise when both values are below it.

Usage:  compare_metrics_csv.py A.csv B.csv [--tol 1e-10] [--floor 1e-14] [--skip TIME,...]
Exit 0 on PASS, 2 on FAIL, 1 on a structural mismatch.
"""
import argparse
import csv
import sys


def load(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"[compare] {path}: empty")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--floor", type=float, default=1e-14,
                    help="absolute floor below which differences are not judged relatively")
    ap.add_argument("--skip", default="", help="comma-separated columns to ignore")
    ap.add_argument("--top", type=int, default=8, help="how many worst columns to print")
    a = ap.parse_args()
    A, B = load(a.a), load(a.b)
    if len(A) != len(B):
        print(f"[compare] FAIL: {len(A)} rows vs {len(B)} rows -- unequal step counts")
        return 1
    skip = {c for c in a.skip.split(",") if c}
    cols = [c for c in A[0] if c in B[0] and c not in skip]
    worst = {}
    for c in cols:
        m, at = 0.0, -1
        for i, (ra, rb) in enumerate(zip(A, B)):
            try:
                va, vb = float(ra[c]), float(rb[c])
            except (TypeError, ValueError):
                continue
            den = max(abs(va), abs(vb), a.floor)
            d = abs(va - vb) / den
            if d > m:
                m, at = d, i
        worst[c] = (m, at)
    ranked = sorted(worst.items(), key=lambda kv: -kv[1][0])
    bad = [c for c, (m, _) in ranked if m > a.tol]
    print(f"[compare] {len(A)} rows, {len(cols)} columns, tol {a.tol:g}, floor {a.floor:g}")
    for c, (m, at) in ranked[:a.top]:
        print(f"  {c:28s} max rel diff {m:.3e}" + (f"  (row {at}, t={A[at].get('TIME')})" if at >= 0 else ""))
    if bad:
        print(f"[compare] FAIL: {len(bad)} column(s) above tolerance: " + ", ".join(bad[:12]))
        return 2
    print("[compare] PASS: every column within tolerance at every step")
    return 0


if __name__ == "__main__":
    sys.exit(main())
