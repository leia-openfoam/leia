#!/usr/bin/env python3
"""Pivot census of the semi-Lagrangian quadratic fit, by cell class.

The reconstruction decides once per mesh, from stencil POSITIONS alone, whether a cell's
quadratic fit is admissible: the smallest Cholesky pivot of the Jacobi-scaled weighted
normal matrix must exceed `quadraticPivotTol` (token SL_QUAD_PIVOT_TOL, default 0.3).
This script reads the per-cell pivot the solver writes with SL_WRITE_FIT_ORDER true
(field slFitPivot, -1 where no full quadratic stencil exists) and tabulates it by class:

    band          |psi| <  6 h_band   (h_band = median cell size of the 0 < alpha < 1 cells)
    near          |psi| < 12 h_band
    far-interior  everything else with cell size >= 0.7 h_band
    far-small     everything else smaller: cfMesh boundary-layer / grading-transition cells

Recipe for the one-step run (serial, in a COPY of a rendered case, never in the study):

    cp -r <case>/{0,0.org,constant,system} <copy>; cd <copy>
    foamDictionary -entry levelSet/semiLagrangian/writeFitOrder -set true system/fvSolution
    foamDictionary -entry endTime -set <deltaT> system/controlDict
    leiaSemiLagrangianLevelSetTwoPhaseFoam > log.census
    python3 workflow/scripts/sl_fit_pivot_census.py <copy> --label "<mesh>" [--csv out.csv]

The `--csv` row per mesh feeds docs/semi-lagrangian-level-set/sl-level-set-article/
data/tables/sl_fit_pivot_census.csv. The threshold is read off the printed histogram:
it must sit in the gap between the degenerate tail and the healthy population, and no
band cell may fall below it (that would cost second order where the interface is).
"""
import argparse, csv, glob, os, statistics as st, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr

EDGES = [0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]


def census(case, label, tol):
    cwd = os.getcwd(); os.chdir(case)
    try:
        n = lr.n_cells(".")
        piv_f = sorted(glob.glob("*/slFitPivot"))
        if not piv_f:
            raise SystemExit(f"{case}: no slFitPivot field -- run one step with writeFitOrder true")
        piv = lr.read_scalar_field(piv_f[0], n)
        if not os.path.exists("0/V"):
            subprocess.run("postProcess -func writeCellVolumes -time 0 > log.writeCellVolumes 2>&1", shell=True)
        V = lr.read_scalar_field("0/V", n)
        psi = lr.read_scalar_field("0/psi", n)
        alpha_f = [f for f in glob.glob("0/alpha*") if not f.endswith(".template")]
        alpha = lr.read_scalar_field(alpha_f[0], n)
    finally:
        os.chdir(cwd)
    size = [v ** (1 / 3) for v in V]
    band = [i for i in range(n) if 1e-9 < alpha[i] < 1 - 1e-9]
    h = st.median(size[i] for i in band)

    def cls(i):
        if abs(psi[i]) < 6 * h: return "band"
        if abs(psi[i]) < 12 * h: return "near"
        return "far-interior" if size[i] >= 0.7 * h else "far-small"

    groups, nofull = {}, {}
    for i in range(n):
        k = cls(i)
        if piv[i] >= 0: groups.setdefault(k, []).append(piv[i])
        else: nofull[k] = nofull.get(k, 0) + 1
    print(f"{label}: {n} cells, h_band {h:.4e}, tolerance {tol}")
    print(f"{'class':13s} {'n':>7s} {'min':>8s} " + " ".join(f"{'<' + str(e):>6s}" for e in EDGES[1:]))
    for k in ("band", "near", "far-interior", "far-small"):
        vals = groups.get(k, [])
        if not vals: continue
        hist = [0] * (len(EDGES) - 1)
        for v in vals:
            for b in range(len(EDGES) - 1):
                if EDGES[b] <= v < EDGES[b + 1]: hist[b] += 1; break
        print(f"{k:13s} {len(vals):7d} {min(vals):8.2e} " + " ".join(f"{c:6d}" for c in hist))
    if nofull: print("cells without a full quadratic stencil (size fallback):", nofull)
    demoted = sum(1 for v in piv if 0 <= v < tol)
    row = {"mesh": label, "nCells": n, "hBand": f"{h:.4e}",
           "minBand": f"{min(groups['band']):.3f}", "minNear": f"{min(groups['near']):.3f}",
           "minInterior": f"{min(groups['far-interior']):.3f}",
           "nSmall": len(groups.get("far-small", [])),
           "minSmall": f"{min(groups['far-small']):.1e}" if groups.get("far-small") else "",
           "nDemoted": demoted, "fracDemoted": f"{demoted / n:.4f}",
           "nGap03to04": sum(1 for v in piv if 0.3 <= v < 0.4), "tol": tol}
    print("demoted at tol:", demoted, f"({100 * demoted / n:.2f} %), cells in [0.3, 0.4):", row["nGap03to04"])
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case"); ap.add_argument("--label", default=None)
    ap.add_argument("--tol", type=float, default=0.3); ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    row = census(a.case, a.label or os.path.basename(os.path.abspath(a.case)), a.tol)
    if a.csv:
        new = not os.path.exists(a.csv)
        rows = [] if new else list(csv.DictReader(open(a.csv)))
        rows = [r for r in rows if r["mesh"] != row["mesh"]] + [row]
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(row)); w.writeheader(); w.writerows(rows)
