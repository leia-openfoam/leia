#!/usr/bin/env python3
"""Curate the Popinet (2009) Section 6.2.2 translating-droplet reproduction.

Popinet introduced this test case in JCP 228 (2009) 5838-5866 for exactly the coupling this
work studies -- a droplet advected across a fixed grid, where advection error perturbs the
curvature and hence, through surface tension, the velocity. His configuration, taken from
the paper: D = 0.4 in the unit square at 64x64 effective (R/h = 12.8), We = 0.4,
La = sigma D/(rho nu^2) over {120, 1200, 12000, infinity}, and CONSTANT density AND
viscosity -- both ratios are 1.

REPORTING CONVENTION IS HIS: the maximum OVER TIME of the spatial norm, with the velocity
relative to U. That is not the endpoint value used elsewhere in this repository, and the
two differ substantially, so the column names say "max_t" explicitly.

L_inf is included here ONLY because it is the quantity Popinet reports and the comparison
would otherwise be impossible. It is not used for any verdict: it converges at half order
against L2's near-first order, which is his finding as well as ours.

Usage:  python3 workflow/scripts/make_popinet_table.py [--root .] [--out FILE]
"""
import argparse, csv, glob, math, os, subprocess, sys

U, TU = 1.0, 0.4          # reference velocity and timescale T_U = D/U


def load(root, study):
    out = []
    for d in sorted(glob.glob(os.path.join(root, "studies", study, "popinetTranslating2D_*"))):
        f = os.path.join(d, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
        if not os.path.isfile(f):
            continue
        rows = list(csv.DictReader(open(f)))
        if not rows:
            continue
        nu = None
        tp = os.path.join(d, "constant", "transportProperties")
        for line in open(tp):
            if line.strip().startswith("nu"):
                nu = float(line.split()[1].rstrip(";"))
                break
        out.append((nu, rows))
    return out


def maxcol(rows, col):
    vals = []
    for r in rows:
        try:
            vals.append(float(r[col]))
        except (TypeError, ValueError):
            pass
    return max(vals) if vals else float("nan")


def fit(hs, es):
    x = [math.log(h) for h in hs]
    y = [math.log(e) for e in es]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / sxx, sxy / math.sqrt(sxx * syy)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=".")
    p.add_argument("--out", default="docs/semi-lagrangian-level-set/sl-level-set-article/"
                                    "data/tables/popinetTranslating.csv")
    a = p.parse_args()
    try:
        commit = subprocess.check_output(["git", "-C", a.root, "rev-parse", "--short", "HEAD"],
                                         text=True).strip()
    except Exception:
        commit = "unknown"

    recs = []
    for N in (64, 128, 256):
        got = load(a.root, f"popinet2D_La12000_N{N}")
        if not got:
            sys.stderr.write(f"warning: popinet2D_La12000_N{N} missing\n")
            continue
        nu, rows = got[0]
        recs.append({"series": "convergence", "N": N, "RoverH": 0.2 * N, "La": 12000,
                     "nu": nu, "steps": len(rows),
                     "maxt_L2": maxcol(rows, "l2MagUPrime") / U,
                     "maxt_Linf": maxcol(rows, "maxMagUPrime") / U,
                     "maxt_L1": maxcol(rows, "meanMagUPrime") / U,
                     "maxt_shape": maxcol(rows, "zeroSetRadialL2"),
                     "gitCommit": commit})
    for nu, rows in sorted(load(a.root, "popinet2D_LaSweep_N64"),
                           key=lambda t: -(t[0] or 0.0)):
        recs.append({"series": "LaSweep", "N": 64, "RoverH": 12.8,
                     "La": (0.4 / nu ** 2 if nu > 0 else float("inf")), "nu": nu,
                     "steps": len(rows),
                     "maxt_L2": maxcol(rows, "l2MagUPrime") / U,
                     "maxt_Linf": maxcol(rows, "maxMagUPrime") / U,
                     "maxt_L1": maxcol(rows, "meanMagUPrime") / U,
                     "maxt_shape": maxcol(rows, "zeroSetRadialL2"),
                     "gitCommit": commit})
    if not recs:
        sys.exit("no Popinet studies found")

    path = os.path.join(a.root, a.out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader()
        w.writerows(recs)
    print(f"wrote {path}  ({len(recs)} rows)")

    conv = [r for r in recs if r["series"] == "convergence"]
    if len(conv) >= 3:
        H = [1.0 / r["N"] for r in conv]
        print("\nconvergence orders (least squares over N = "
              + "/".join(str(r["N"]) for r in conv) + "):")
        for key, label, claim in (("maxt_L2", "L2  (his RMS)", "close to first order"),
                                  ("maxt_Linf", "Linf (his Max)", "less than first order"),
                                  ("maxt_shape", "shape", "roughly first order")):
            pp, rr = fit(H, [r[key] for r in conv])
            print(f"  {label:<15} p = {pp:5.2f}  R = {rr:6.3f}   Popinet: {claim}")


if __name__ == "__main__":
    main()
