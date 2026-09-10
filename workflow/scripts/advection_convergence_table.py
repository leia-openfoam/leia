#!/usr/bin/env python3
"""Mesh convergence table for an ADVECTION study: error AND observed order per arm.

WHY THIS SCRIPT EXISTS. CLAUDE.md: "Touch advection, run a mesh convergence study --
before stating anything." A single resolution gives an error, not a result, and the ORDER
is the only column that separates a real gain from a floor. MEASURED 2026-09-10: a bound
lowered every error at one rung and its eikonal error then moved 7.119e-03 -> 7.066e-03
across a 2x refinement, an order of 0.01.

TWO THINGS IT REFUSES TO DO.

1. It does NOT plot against MAX_CELL_SIZE. That token is a PIN, not a mesh size: cfMesh
   treats it as a request and the mesh it returns has its own spacing. The effective
   spacing comes from the cell count, h_eff = (V_domain/nCells)^(1/3), read from the mesh
   itself, and the order is computed against that. A ladder plotted against a request is
   not a ladder.

2. It does NOT report L_inf. L_INF_E_PSI does not converge for these cases and a verdict
   built on it had to be retracted. The geometric error, the volume error and the
   boundedness error are reported TOGETHER, which is the rule for an interface method.

An arm that DIVERGED is reported as such and excluded from its order, because an endpoint
estimator cannot be compared across unequal horizons. The state comes from the committed
log classifier, never from a return code.

Usage:
  python3 workflow/scripts/advection_convergence_table.py <study> [<study> ...]
      [--root .] [--case CASE] [--solver leiaSemiLagrangeLevelSetFoam] [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess

METRICS = ["E_GEOM_ALPHA_REL", "E_VOL_ALPHA_REL", "E_BOUND_ALPHA"]
HERE = os.path.dirname(os.path.abspath(__file__))


def ncells(case):
    """Cell count from the mesh header note, so h_eff comes from the MESH."""
    for f in ("owner", "faces"):
        p = os.path.join(case, "constant", "polyMesh", f)
        if not os.path.exists(p):
            continue
        with open(p, "rb") as fh:
            head = fh.read(4096).decode("utf-8", "replace")
        m = re.search(r"nCells:\s*(\d+)", head)
        if m:
            return int(m.group(1))
    return None


def state(case, solver):
    log = os.path.join(case, f"log.{solver}")
    if not os.path.exists(log):
        return "NO-LOG"
    try:
        r = subprocess.run([os.path.join(HERE, "foam_log_state.sh"), log],
                           capture_output=True, text=True, timeout=60)
        return r.stdout.split()[0] if r.stdout.strip() else "UNKNOWN"
    except Exception:
        return "UNKNOWN"


def rows(case, solver):
    p = os.path.join(case, f"{solver}.csv")
    if not os.path.exists(p):
        return []
    out = []
    for r in csv.DictReader(open(p)):
        try:
            out.append({k: (float(v) if v not in (None, "") else None)
                        for k, v in r.items() if k is not None})
        except (ValueError, TypeError):
            break
    return [r for r in out if r.get(METRICS[0]) is not None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--root", default=".")
    ap.add_argument("--case", default=None)
    ap.add_argument("--solver", default="leiaSemiLagrangeLevelSetFoam")
    ap.add_argument("--csv")
    a = ap.parse_args()

    arms = []
    for s in a.studies:
        pat = os.path.join(a.root, "studies", s, (a.case or "*") + "_*/")
        for d in sorted(glob.glob(pat)):
            pj = os.path.join(d, "case_params.json")
            if not os.path.exists(pj):
                continue
            t = json.load(open(pj))["tokens"]
            nc = ncells(d)
            arms.append(dict(
                dir=d, study=s,
                bound=t.get("SL_VALUE_BOUND", "?"),
                lmode=t.get("SL_CONE_L_MODE", "?"),
                pin=t.get("N_CELLS") or t.get("MAX_CELL_SIZE") or "?",
                nCells=nc,
                heff=(1.0 / nc ** (1.0 / 3.0)) if nc else None,   # unit-box domain
                state=state(d, a.solver), rows=rows(d, a.solver)))

    if not arms:
        raise SystemExit("no arms found")

    out = []
    bounds = sorted({x["bound"] for x in arms})
    for m in METRICS:
        print(f"\n=== {m} ===")
        print(f"{'bound':16s}{'pin':>10s}{'nCells':>10s}{'h_eff':>11s}"
              f"{'value':>13s}{'order':>8s}{'state':>11s}")
        for b in bounds:
            sel = sorted([x for x in arms if x["bound"] == b],
                         key=lambda x: (x["nCells"] or 0))
            prev = None
            for x in sel:
                v = x["rows"][-1][m] if x["rows"] else None
                o = float("nan")
                if (prev and v and prev[1] and x["heff"] and prev[0]
                        and x["state"] == "COMPLETED" and prev[2] == "COMPLETED"
                        and v > 0 and prev[1] > 0 and x["heff"] != prev[0]):
                    o = math.log(prev[1] / v) / math.log(prev[0] / x["heff"])
                print(f"{b:16s}{str(x['pin']):>10s}{str(x['nCells'] or '?'):>10s}"
                      f"{(f'{x[chr(104)+chr(101)+chr(102)+chr(102)]:.4e}' if x['heff'] else '?'):>11s}"
                      f"{(f'{v:.4e}' if v is not None else '--'):>13s}"
                      f"{(f'{o:.2f}' if o == o else '--'):>8s}{x['state']:>11s}")
                out.append(dict(metric=m, bound=b, pin=x["pin"], nCells=x["nCells"],
                                hEff=x["heff"], value=v, order=(o if o == o else ""),
                                state=x["state"], study=x["study"]))
                if v and x["heff"]:
                    prev = (x["heff"], v, x["state"])
        print("  order is computed between consecutive COMPLETED rungs, against h_eff")

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader(); w.writerows(out)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
