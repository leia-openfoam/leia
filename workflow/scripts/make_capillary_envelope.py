#!/usr/bin/env python3
"""Curate the OPERATING ENVELOPE: transport dissipation vs capillary survival.

Reads the per-case output of config/capillaryEnvelope.yaml (the capillary arm)
and writes one curated CSV + a readable table. Optionally also reads
config/capillaryEnvelopeKinematic.yaml (the transport-cost arm) so both halves
of the trade appear in one place.

WHAT THE ENVELOPE IS FOR. The campaign established by elimination that the
parasitic-current instability is not systematic error in any operator: the
pressure-velocity coupling is exact (imposed constant curvature leaves an m=2
perturbation frozen to 1.0000 over 3200 steps), so the instability is entirely
curvature-driven, and curvature being a second derivative of psi amplifies psi
roughness by 1/h^2 -- which is why refinement hurts. Reinitialization-free means
nothing removes that roughness; the only damping anywhere is the SL
interpolation error, and a quadratic reconstruction has almost none, which is
precisely why the transport is accurate. Virtue and vice are the same property.
So the deliverable is not a winning scheme but a MAPPED trade-off: how much
transport accuracy must be spent to buy how much capillary survival.

COLUMNS, and why each is here:
  t_blow            last time reached; with solverFailed it is the blow-up time,
                    otherwise a survival lower bound (possibly wall-clock capped
                    -- flagged separately, never silently conflated)
  reachedEndTime    whether the horizon was reached, so a censored run is never
                    read as a blow-up
  maxMagU_final     the parasitic current at the end of the run
  volumeRelError,   BOTH are reported at a MATCHED time across arms, because the
  shapeL2           campaign has a standing case where they disagreed (a delivery
                    improved volume order while failing the gradient) and a
                    single-metric view misread it as stability
  minGradPsiBand    the profile drift, the quantity the delivered-curvature error
                    is slaved to (measured slope 1.00, r^2 = 0.98)
  A2hL2Band         the grid-scale corrugation the filter is designed to damp
  driverAcrossSupportL2  the non-absorbable capillary driver (WP8.1)

Usage:
  make_capillary_envelope.py <capillary_study_dir> <out_dir>
                             [--kinematic <kinematic_study_dir>]
"""

import argparse
import csv
import glob
import json
import math
import os
import sys

SOLVER_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"

# tokens that identify an arm
AXES = ("N_CELLS", "SL_RECONSTRUCTION", "PSI_FILTER_THETA")

SHORT = {
    "uncachedQuadraticWeightedLeastSquares": "quadratic",
    "linearTaylor": "linearTaylor",
    "linearWeightedLeastSquares": "linearWLSQ",
}


def _rows(path):
    """Per-step rows, skipping any line torn by a mid-write blow-up."""
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if any(v is None or v == "" for v in r.values()):
                continue
            try:
                out.append({k: float(v) for k, v in r.items()})
            except ValueError:
                continue
    return out


def _tokens(case_dir):
    p = os.path.join(case_dir, "case_params.json")
    if not os.path.exists(p):
        return {}
    try:
        return json.load(open(p)).get("tokens", {})
    except Exception:
        return {}


def collect(study):
    cases = []
    for d in sorted(glob.glob(os.path.join(study, "*_0*"))):
        if not os.path.isdir(d):
            continue
        tk = _tokens(d)
        rows = _rows(os.path.join(d, SOLVER_CSV))
        if not tk or not rows:
            continue
        failed = os.path.exists(os.path.join(d, ".leia_solver_failed"))
        last = rows[-1]
        end_time = None
        try:
            end_time = float(tk.get("END_TIME"))
        except (TypeError, ValueError):
            pass
        reached = (end_time is not None
                   and last["TIME"] >= 0.999 * end_time)
        cases.append({
            "case_dir": os.path.basename(d),
            "N": int(tk.get("N_CELLS", 0)),
            "transport": SHORT.get(tk.get("SL_RECONSTRUCTION", ""),
                                   tk.get("SL_RECONSTRUCTION", "")),
            "theta": float(tk.get("PSI_FILTER_THETA", 0.0)),
            "t_last": last["TIME"],
            "diverged": failed,
            "reachedEndTime": reached,
            "rows": rows,
        })
    return cases


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("outdir")
    ap.add_argument("--kinematic", default=None)
    a = ap.parse_args(argv)

    cases = collect(a.study)
    if not cases:
        print(f"[envelope] no usable cases under {a.study}", file=sys.stderr)
        return 1

    # A matched instant for the error comparison: the earliest final time over
    # ALL arms, so volume and shape are read at the same physical instant
    # everywhere rather than at each arm's own end (which would compare a
    # survivor against a blow-up).
    t_match = min(c["t_last"] for c in cases)

    def at(rows, t):
        return min(rows, key=lambda r: abs(r["TIME"] - t))

    out = []
    for c in cases:
        m = at(c["rows"], t_match)
        f = c["rows"][-1]
        out.append({
            "N": c["N"], "transport": c["transport"], "theta": c["theta"],
            "t_last": f"{c['t_last']:.6g}",
            "diverged": int(c["diverged"]),
            "reachedEndTime": int(c["reachedEndTime"]),
            "t_matched": f"{t_match:.6g}",
            "maxMagU_final": f"{f['maxMagU']:.6g}",
            "maxMagU_matched": f"{m['maxMagU']:.6g}",
            "volumeRelError_matched": f"{m['phaseVolumeRelError']:.6g}",
            "shapeL2_matched": f"{m['zeroSetRadialL2']:.6g}",
            "minGradPsiBand_matched": f"{m['minGradPsiBand']:.6g}",
            "A2hL2Band_matched": f"{m.get('A2hL2Band', float('nan')):.6g}",
            "driverAcrossSupportL2_matched":
                f"{m.get('driverAcrossSupportL2', float('nan')):.6g}",
            "case_dir": c["case_dir"],
        })
    out.sort(key=lambda r: (r["transport"], r["theta"], r["N"]))

    os.makedirs(a.outdir, exist_ok=True)
    p = os.path.join(a.outdir, "capillary_envelope.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"[envelope] wrote {p} ({len(out)} arms)")

    # ---- readable table -----------------------------------------------------
    print(f"\nCAPILLARY SURVIVAL vs TRANSPORT and FILTER "
          f"(errors at the matched instant t = {t_match:.5g} s)\n")
    hdr = (f"{'transport':<13} {'theta':>6} {'N':>5} {'t_last':>10} "
           f"{'fate':>10} {'maxU_end':>11} {'vol err':>10} {'shape L2':>10} "
           f"{'minGradPsi':>11}")
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        fate = ("DIVERGED" if r["diverged"]
                else ("reached" if r["reachedEndTime"] else "censored"))
        print(f"{r['transport']:<13} {r['theta']:>6} {r['N']:>5} "
              f"{float(r['t_last']):>10.5g} {fate:>10} "
              f"{float(r['maxMagU_final']):>11.4g} "
              f"{float(r['volumeRelError_matched']):>10.4g} "
              f"{float(r['shapeL2_matched']):>10.4g} "
              f"{float(r['minGradPsiBand_matched']):>11.5f}")

    # ---- the h-trend, which is the whole question ---------------------------
    print("\nh-TREND per (transport, theta): does refinement still blow first?")
    print("  (positive slope = finer survives LONGER = the trend has inverted)")
    keys = sorted({(r["transport"], r["theta"]) for r in out})
    for tr, th in keys:
        pts = sorted([(r["N"], float(r["t_last"]), r["diverged"])
                      for r in out if r["transport"] == tr and r["theta"] == th])
        if len(pts) < 2:
            continue
        cens = any(not d for _, _, d in pts)
        ts = "  ".join(f"N={n}:{t:.4g}{'' if d else '*'}" for n, t, d in pts)
        p0, p1 = pts[0], pts[-1]
        slope = (math.log(p1[1] / p0[1]) / math.log(p1[0] / p0[0])
                 if p0[1] > 0 and p1[1] > 0 else float("nan"))
        print(f"  {tr:<13} theta={th:<5} {ts}   slope={slope:+.3f}"
              + ("   (* = not a blow-up: survived or censored)" if cens else ""))

    if a.kinematic:
        print("\nTRANSPORT COST (pure advection, prescribed velocity) -- "
              "an UPPER BOUND on the coupled cost, since coupled interface")
        print("motion is capped at 1st order in time by Euler momentum on a "
              "dt ~ h^1.5 step:")
        kp = os.path.join(a.kinematic, os.path.basename(a.kinematic.rstrip("/"))
                          + "_errors.csv")
        if os.path.exists(kp):
            for r in csv.DictReader(open(kp)):
                rec = r.get("reconstruction", "")
                if rec:
                    print(f"  {SHORT.get(rec, rec):<13} h={r.get('h','?'):>10} "
                          f"shape={r.get('shapeError','?'):>12} "
                          f"vol={r.get('volumeError','?'):>12}")
        else:
            print(f"  (no {kp} yet)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
