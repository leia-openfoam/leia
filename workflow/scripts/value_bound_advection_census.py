#!/usr/bin/env python3
"""Tabulate the value-bound arms of an ADVECTION case, at the LAST COMMON step.

WHY THE LAST COMMON STEP. An endpoint estimator is not comparable across unequal horizons,
and the whole point of this table is that some arms DIVERGE while others complete. Reading
each arm at its own last row would compare a run that died at step 198 against runs that
reached 967 -- and the comparison would be meaningless in the direction that matters. Every
row is therefore read at the last step for which EVERY arm has data, and the per-arm row
count is printed so a short arm is visible rather than hidden.

WHY THE MESH DIGEST IS A COLUMN. cfMesh is not bit-reproducible, so arms whose meshes differ
are comparing the mesher. The digest column must be IDENTICAL down a case's block; the script
says NOT COMPARABLE when it is not.

NEVER L_inf: it does not converge for these cases and a verdict built on it had to be
retracted. This reports the geometric error, the volume error and the boundedness error
together, which is the rule for an interface method -- a single headline metric has already
misread a candidate in this repository.

Usage:
  python3 workflow/scripts/value_bound_advection_census.py <arm-dir> [<arm-dir> ...]
      [--label L] [--csv out.csv]

Each <arm-dir> is an output directory of workflow/scripts/advect_bound_arm.sh.
"""
import argparse
import csv
import hashlib
import os
import sys

METRICS = ["E_GEOM_ALPHA_REL", "E_VOL_ALPHA_REL", "E_BOUND_ALPHA"]


def digest(arm):
    p = os.path.join(arm, "constant", "polyMesh", "points")
    if not os.path.exists(p):
        return "NO-MESH"
    with open(p, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()[:12]


def rows(arm, solver="leiaSemiLagrangeLevelSetFoam"):
    p = os.path.join(arm, f"{solver}.csv")
    if not os.path.exists(p):
        return []
    out = []
    for r in csv.DictReader(open(p)):
        try:
            out.append({k: (float(v) if v not in (None, "") else None)
                        for k, v in r.items() if k is not None})
        except (ValueError, TypeError):
            break          # a truncated final row from a killed run: keep what is whole
    return [r for r in out if r.get(METRICS[0]) is not None]


def state(arm):
    """Outcome from the log, via the committed classifier -- never from a return code."""
    log = os.path.join(arm, "log.solve")
    if not os.path.exists(log):
        return "NO-LOG"
    import subprocess
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        r = subprocess.run([os.path.join(here, "foam_log_state.sh"), log],
                           capture_output=True, text=True, timeout=60)
        return r.stdout.strip().split()[0] if r.stdout.strip() else "UNKNOWN"
    except Exception:
        return "UNKNOWN"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+")
    ap.add_argument("--label", default="")
    ap.add_argument("--csv")
    a = ap.parse_args()

    data = {arm: rows(arm) for arm in a.arms}
    have = {k: v for k, v in data.items() if v}
    if not have:
        sys.exit("no arm produced a metrics CSV")
    common = min(len(v) for v in have.values())

    digs = {arm: digest(arm) for arm in a.arms}
    uniq = set(digs.values())
    comparable = len(uniq) == 1 and "NO-MESH" not in uniq

    print(f"\n=== value-bound advection census {a.label} ===")
    print(f"  arms                {len(a.arms)}")
    print(f"  last COMMON step    {common}")
    print(f"  mesh digest         {'IDENTICAL ' + uniq.pop() if comparable else 'DIFFERS -- NOT COMPARABLE: ' + str(sorted(uniq))}")
    if not comparable:
        print("  A metric difference between arms on DIFFERENT meshes measures the mesher.")
    hdr = f"{'arm':38s}{'state':>11s}{'rows':>6s}" + "".join(f"{m[:14]:>15s}" for m in METRICS)
    print(hdr)
    print("-" * len(hdr))
    out = []
    for arm in a.arms:
        st = state(arm)
        v = data[arm]
        name = os.path.basename(arm.rstrip("/"))
        if not v:
            print(f"{name:38s}{st:>11s}{0:6d}" + "".join(f"{'--':>15s}" for _ in METRICS))
            continue
        r = v[common - 1]
        print(f"{name:38s}{st:>11s}{len(v):6d}"
              + "".join(f"{r[m]:15.4e}" for m in METRICS))
        out.append({"arm": name, "state": st, "rows": len(v),
                    "commonStep": common, "meshDigest": digs[arm],
                    **{m: r[m] for m in METRICS}})

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
