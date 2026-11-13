#!/usr/bin/env python3
"""Per-column differences of two metric CSVs with equal row counts.

    column_differences.py A.csv B.csv [--skip COLS] [--only-diff]

For every numeric column outside --skip the script prints the first row and TIME where the
two files differ, the largest relative difference |a - b| / max(|a|, |b|) with its row and
TIME, and the two values of the last row. It is a measurement (for example serial against
4 ranks); it has no pass or fail. Use compare_metrics_csv.py for a verdict.
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--skip", default="ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME")
    ap.add_argument("--only-diff", action="store_true", help="omit the columns that are identical")
    a = ap.parse_args()
    A = list(csv.DictReader(open(a.a, newline="")))
    B = list(csv.DictReader(open(a.b, newline="")))
    if len(A) != len(B):
        print("row counts differ: %d against %d" % (len(A), len(B)))
        return 1
    skip = {c for c in a.skip.split(",") if c}
    print("| column | first row | first TIME | largest rel. diff. | row | TIME | last value A | last value B |")
    print("|---|---|---|---|---|---|---|---|")
    n_diff = 0
    for c in A[0]:
        if c in skip or c not in B[0]:
            continue
        first, worst, at = None, 0.0, None
        for i, (ra, rb) in enumerate(zip(A, B)):
            if ra[c] == rb[c]:
                continue
            try:
                va, vb = float(ra[c]), float(rb[c])
                rel = abs(va - vb) / max(abs(va), abs(vb), 1e-300)
            except ValueError:
                rel = float("inf")
            if first is None:
                first = i
            if rel > worst:
                worst, at = rel, i
        if first is None:
            if not a.only_diff:
                print("| %s | identical | | 0 | | | %s | %s |" % (c, A[-1][c], B[-1][c]))
            continue
        n_diff += 1
        print("| %s | %d | %s | %.3e | %d | %s | %s | %s |" % (
            c, first, A[first].get("TIME", "?"), worst, at, A[at].get("TIME", "?"), A[-1][c], B[-1][c]))
    print("%d rows; %d column(s) differ" % (len(A), n_diff))
    return 0


if __name__ == "__main__":
    sys.exit(main())
