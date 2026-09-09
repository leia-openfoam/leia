#!/usr/bin/env python3
"""Amplification census of the semi-Lagrangian fit: Lambda by cell class AND by fit order.

Lambda_c = |1 - sum_j g_j| + sum_j |g_j| is the Lebesgue constant of the fit at the departure
displacement d, with g_j = b(d)^T M^-1 w_j^2 b(d_j). Lambda = 1 means the update is a convex
combination of the stencil values and cannot create a new extremum; Lambda > 1 is the
necessary condition for the checkerboard growth that destroys the far field of a translating
polyhedral case. It depends on the stencil GEOMETRY alone, so it costs one mesh pass.

sl_fit_pivot_census.py already tabulates the Cholesky PIVOT by cell class. This script asks
the two questions the pivot alone cannot answer:

  1. WHICH MESHERS PRODUCE AMPLIFYING CELLS? Lambda side by side for blockMesh, cartesianMesh,
     pMesh and snappyHexMesh on the same box at the same nominal cell size and the same probe
     displacement. cartesianMesh is the sharpest arm: cfMesh's HEXAHEDRAL mesher on the same
     meshDict and surface file as pMesh, so it separates hex-versus-polyhedral from
     mesher-versus-mesher.

  2. IS THERE HEADROOM IN RANK REDUCTION? Lambda cross-tabulated against slFitOrder
     (2 = quadratic, 1 = linear, 0 = the cell value only). Rank reduction shrinks the
     pseudo-inverse in exactly the directions that make ||g||_1 large, so it moves Lambda
     toward 1, and the rank-1 limit is the weighted mean: every g_j >= 0, sum g = 1,
     Lambda = 1 exactly. The question is whether the high-Lambda cells are ALREADY demoted to
     a linear fit by quadraticPivotTol. If they are, a truncated-SVD rank choice has little
     left to remove and the lever is weak. If they are still quadratic, it has a lot.

Cell classes follow sl_fit_pivot_census.py so the two tables are read together:

    band          |psi| <  6 h_band   (h_band = median cell size of the 0 < alpha < 1 cells)
    near          |psi| < 12 h_band
    far-interior  everything else with cell size >= 0.7 h_band
    far-small     everything else smaller: cfMesh boundary-layer / grading-transition cells

Report L1 and L2 style summaries -- median, p99, p99.9, max -- never a bare maximum alone.

Recipe (one mesh pass per mesher, serial, minutes):

    workflow/scripts/fit_amplification_probe.sh <poly-study-case> cartesianMesh /tmp/probe_cart
    python3 workflow/scripts/sl_fit_amplification_census.py /tmp/probe_cart --label cartesianMesh \\
        --csv docs/method-comparison/method-comparison-article/data/tables/sl_fit_amplification.csv

Usage:  python3 workflow/scripts/sl_fit_amplification_census.py <probe-dir> [--label L]
            [--csv out.csv] [--tol 0.3] [--top 8]
"""
import argparse
import csv
import glob
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr  # noqa: E402

ORDER_NAME = {2: "quadratic", 1: "linear", 0: "cell value"}


def pct(sorted_vals, q):
    """Percentile of an already sorted list, linear interpolation."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    i = q * (len(sorted_vals) - 1)
    lo = int(i)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (i - lo) * (sorted_vals[hi] - sorted_vals[lo])


def summarise(vals):
    v = sorted(vals)
    return {
        "n": len(v),
        "median": pct(v, 0.5),
        "p99": pct(v, 0.99),
        "p999": pct(v, 0.999),
        "max": v[-1] if v else float("nan"),
    }


def load(case):
    cwd = os.getcwd()
    os.chdir(case)
    try:
        n = lr.n_cells(".")
        def field(name, required=True):
            hits = sorted(glob.glob(f"*/{name}"))
            if not hits:
                if required:
                    raise SystemExit(
                        f"{case}: no {name} -- run fit_amplification_probe.sh first")
                return None
            return lr.read_scalar_field(hits[0], n)
        amp = field("slFitAmplification")
        order = field("slFitOrder", required=False)
        pivot = field("slFitPivot", required=False)
        if not os.path.exists("0/V"):
            raise SystemExit(f"{case}: no 0/V -- postProcess -func writeCellVolumes")
        V = lr.read_scalar_field("0/V", n)
        psi = lr.read_scalar_field("0/psi", n)
        alpha_f = [f for f in glob.glob("0/alpha*") if not f.endswith(".template")]
        alpha = lr.read_scalar_field(alpha_f[0], n)
    finally:
        os.chdir(cwd)
    return n, amp, order, pivot, V, psi, alpha


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--label", default=None)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--tol", type=float, default=0.3)
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()
    label = args.label or os.path.basename(os.path.abspath(args.case))

    n, amp, order, pivot, V, psi, alpha = load(args.case)
    size = [v ** (1.0 / 3.0) for v in V]
    band = [i for i in range(n) if 1e-9 < alpha[i] < 1 - 1e-9]
    if not band:
        raise SystemExit(f"{args.case}: no interface cells -- is 0/alpha* initialised?")
    h = st.median(size[i] for i in band)

    # A singular normal matrix is written as -1 by the diagnostic; it is not a Lambda.
    singular = [i for i in range(n) if amp[i] < 0]
    ok = [i for i in range(n) if amp[i] >= 0]

    def cls(i):
        if abs(psi[i]) < 6 * h:
            return "band"
        if abs(psi[i]) < 12 * h:
            return "near"
        return "far-interior" if size[i] >= 0.7 * h else "far-small"

    print(f"{label}: {n} cells, h_band {h:.4e} m, "
          f"{len(singular)} cells with a singular normal matrix")
    if order:
        counts = {}
        for i in range(n):
            counts[int(round(order[i]))] = counts.get(int(round(order[i])), 0) + 1
        parts = [f"{ORDER_NAME.get(k, k)} {v} ({100*v/n:.2f} %)"
                 for k, v in sorted(counts.items(), reverse=True)]
        print("  fit order: " + ", ".join(parts))
    print()

    # ---- 1. Lambda by cell class ---------------------------------------
    print(f"  {'class':<14}{'n':>9}{'median':>10}{'p99':>10}{'p99.9':>10}{'max':>10}"
          f"{'> 1.10':>9}")
    by_class = {}
    for i in ok:
        by_class.setdefault(cls(i), []).append(amp[i])
    for k in ("band", "near", "far-interior", "far-small"):
        if k not in by_class:
            continue
        s = summarise(by_class[k])
        over = sum(1 for v in by_class[k] if v > 1.10)
        print(f"  {k:<14}{s['n']:>9}{s['median']:>10.4f}{s['p99']:>10.4f}"
              f"{s['p999']:>10.4f}{s['max']:>10.4f}{over:>9}")

    # ---- 2. Lambda by FIT ORDER: is there rank-reduction headroom? -----
    row_head = None
    if order:
        print()
        print(f"  {'fit order':<14}{'n':>9}{'median':>10}{'p99':>10}{'p99.9':>10}"
              f"{'max':>10}{'> 1.10':>9}")
        by_order = {}
        for i in ok:
            by_order.setdefault(int(round(order[i])), []).append(amp[i])
        for k in sorted(by_order, reverse=True):
            s = summarise(by_order[k])
            over = sum(1 for v in by_order[k] if v > 1.10)
            print(f"  {ORDER_NAME.get(k, k):<14}{s['n']:>9}{s['median']:>10.4f}"
                  f"{s['p99']:>10.4f}{s['p999']:>10.4f}{s['max']:>10.4f}{over:>9}")

        # The headroom question, stated as a number: of the cells that carry the
        # worst Lambda, how many are ALREADY demoted below a quadratic fit?
        worst = sorted(ok, key=lambda i: -amp[i])[:max(1, len(ok) // 1000)]
        dem = sum(1 for i in worst if int(round(order[i])) < 2)
        print()
        print(f"  RANK-REDUCTION HEADROOM: of the {len(worst)} cells with the highest Lambda"
              f" (top 0.1 %), {dem} ({100*dem/len(worst):.1f} %) are ALREADY demoted below a"
              f" quadratic fit.")
        print(f"    Low share  -> the amplifiers still carry a full quadratic, so choosing the"
              f" rank per cell has room to lower Lambda.")
        print(f"    High share -> quadraticPivotTol has already demoted them and Lambda is"
              f" still > 1, so rank reduction is a weak lever here.")
        row_head = 100.0 * dem / len(worst)

    # ---- the worst cells, with the properties that a rule could use -----
    print()
    print(f"  the {args.top} worst cells:")
    print(f"    {'cell':>9}{'Lambda':>9}{'size/h':>9}{'|psi|/h':>10}"
          f"{'order':>7}{'pivot':>10}{'class':>14}")
    for i in sorted(ok, key=lambda j: -amp[j])[:args.top]:
        o = int(round(order[i])) if order else -1
        p = pivot[i] if pivot else float("nan")
        print(f"    {i:>9}{amp[i]:>9.4f}{size[i]/h:>9.3f}{abs(psi[i])/h:>10.2f}"
              f"{o:>7}{p:>10.2e}{cls(i):>14}")

    # ---- Lambda against cell size: the recorded monotone trend ---------
    print()
    print(f"  Lambda by cell size (the excess falls monotonically with size):")
    print(f"    {'size/h':<14}{'n':>9}{'median':>10}{'max':>10}")
    edges = [(0, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.1), (1.1, 1e9)]
    for lo, hi in edges:
        vals = [amp[i] for i in ok if lo <= size[i] / h < hi]
        if not vals:
            continue
        s = summarise(vals)
        tag = f"[{lo}, {hi})" if hi < 1e9 else f">= {lo}"
        print(f"    {tag:<14}{s['n']:>9}{s['median']:>10.4f}{s['max']:>10.4f}")

    if args.csv:
        alls = summarise([amp[i] for i in ok])
        small = summarise(by_class.get("far-small", [])) if "far-small" in by_class else None
        row = {
            "mesh": label, "nCells": n, "hBand": f"{h:.4e}",
            "lambdaMedian": f"{alls['median']:.4f}",
            "lambdaP999": f"{alls['p999']:.4f}",
            "lambdaMax": f"{alls['max']:.4f}",
            "nOver110": sum(1 for i in ok if amp[i] > 1.10),
            "fracOver110": f"{sum(1 for i in ok if amp[i] > 1.10)/n:.6f}",
            "smallN": small["n"] if small else 0,
            "smallMax": f"{small['max']:.4f}" if small else "",
            "nSingular": len(singular),
            "worstAlreadyDemotedPct": f"{row_head:.1f}" if row_head is not None else "",
        }
        rows = []
        if os.path.exists(args.csv):
            with open(args.csv) as fh:
                rows = [r for r in csv.DictReader(fh) if r.get("mesh") != label]
        rows.append(row)
        os.makedirs(os.path.dirname(os.path.abspath(args.csv)), exist_ok=True)
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row))
            w.writeheader()
            w.writerows(rows)
        print(f"\n  wrote {args.csv}")


if __name__ == "__main__":
    main()
