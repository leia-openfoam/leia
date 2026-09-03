#!/usr/bin/env python3
"""Curate the repaired translating-droplet matrix (config/translatingRepaired{,EqualRho}2D).

Both studies are 8 arms: MASS_FLUX (rhoLENT, geometricFaceDensity) crossed with
MOMENTUM_DIV_SCHEME (upwind, limitedLinearV 1, vanLeerV, linearUpwind gradU), at density
ratio 838.824 and at ratio 1 with rho1+rho2 held at 999.39.

Every row is read at the SAME step index, because endpoint estimators are not comparable
across unequal horizons (CLAUDE.md, research loop step 5). That index is NOT the largest
one every arm reached: at ratio 838.8 the arms terminate between 8427 and 9987 steps, so
the common maximum (8427) is set by the FIRST ARM TO DIE and samples every other
ratio-838.8 arm inside its own blow-up. Read there, max|U-U0|/U0 ranges from 4.5 to 89 --
i.e. up to 8900 percent of the translation speed -- which is a divergence transient, not a
parasitic-current level, and reporting it as one was a curation error made on 2026-09-02.

--read-step therefore defaults to 5000 (t = 0.0375 s), where every arm in both studies is
still healthy; the script refuses to run if any arm is shorter than that. Velocity errors
are reported RELATIVE to the translation speed, since a dimensionless 0.2 reads as "twenty
percent of the speed the droplet is carried at" and 1e-2 m/s does not. The per-arm
termination step and final time are carried alongside, so the divergence outcome is still
visible without contaminating the comparison.

Usage:  python3 workflow/scripts/make_translating_matrix_table.py [--root .] [--out FILE]
                [--read-step 5000]
"""
import argparse, csv, glob, os, subprocess, sys

SCHEMES = ["upwind", "limitedLinearV 1", "vanLeerV", "linearUpwind gradU"]
MASSFLUX = ["rhoLENT", "geometricFaceDensity"]
STUDIES = [("translatingRepaired2D", 838.824), ("translatingRepairedEqualRho2D", 1.0)]
U0, X0 = 0.05, 0.0025002
COLS = ["maxMagUPrime", "meanMagUPrime", "phaseVolumeRelError", "zeroSetRadialL2",
        "gradPsiL2ErrorBand", "rhoMassResidualRelL1", "centroidX"]


def num(row, col):
    try:
        return float(row.get(col))
    except (TypeError, ValueError):
        return float("nan")


def load(root, study):
    arms = {}
    for d in sorted(glob.glob(os.path.join(root, "studies", study, "translatingDroplet2D_*"))):
        f = os.path.join(d, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
        if not os.path.isfile(f):
            continue
        rows = list(csv.DictReader(open(f)))
        if rows:
            arms[int(d[-5:])] = rows
    return arms


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=".")
    p.add_argument("--out", default="docs/method-comparison/method-comparison-article/"
                                   "data/tables/translatingRepairedMatrix.csv")
    p.add_argument("--read-step", type=int, default=5000,
                   help="1-based step index every arm is read at; must be inside the "
                        "healthy window of EVERY arm, not merely reached by all of them")
    a = p.parse_args()

    loaded = [(s, r, load(a.root, s)) for s, r in STUDIES]
    missing = [s for s, _, arms in loaded if len(arms) != 8]
    if missing:
        sys.exit(f"incomplete studies: {missing}")
    shortest = min(len(rows) for _, _, arms in loaded for rows in arms.values())
    common = a.read_step
    if shortest < common:
        sys.exit(f"--read-step {common} exceeds the shortest arm ({shortest} rows); an arm "
                 f"diverged before the read point, so no comparison is valid there")

    try:
        commit = subprocess.check_output(["git", "-C", a.root, "rev-parse", "--short", "HEAD"],
                                         text=True).strip()
    except Exception:
        commit = "unknown"

    out = []
    for study, ratio, arms in loaded:
        for k in sorted(arms):
            rows = arms[k]
            r, t = rows[common - 1], float(rows[common - 1]["TIME"])
            rec = {"study": study, "densityRatio": ratio, "massFlux": MASSFLUX[k // 4],
                   "momentumDivScheme": SCHEMES[k % 4], "arm": f"translatingDroplet2D_{k:05d}",
                   "stepsReached": len(rows), "tEnd": float(rows[-1]["TIME"]),
                   "completed": int(len(rows) >= 13334), "readStep": common, "readTime": t,
                   "travelFraction": (num(r, "centroidX") - X0) / (U0 * t) if t else float("nan"),
                   # velocity errors are reported RELATIVE to the translation speed: a
                   # dimensionless 0.2 reads as "twenty percent of the speed the droplet
                   # is carried at", which 1e-2 m/s does not.
                   "maxMagUPrimeRel": num(r, "maxMagUPrime") / U0,
                   "meanMagUPrimeRel": num(r, "meanMagUPrime") / U0,
                   "gitCommit": commit}
            rec.update({c: num(r, c) for c in COLS})
            out.append(rec)

    os.makedirs(os.path.dirname(os.path.join(a.root, a.out)), exist_ok=True)
    path = os.path.join(a.root, a.out)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"wrote {path}  ({len(out)} rows, read at step {common}; "
          f"shortest arm {shortest} rows)")


if __name__ == "__main__":
    main()
