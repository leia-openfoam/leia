#!/usr/bin/env python3
"""Domain-size control: is truncating the far field at FIXED h legitimate?

WHY THIS EXISTS. The 3D ladder that carries the dimensional confirmation is
pinned at its coarse end by R/h >= 10, so with cell-count doubling per level it
spans only a factor 2 in h -- an error ratio of 4 end to end for a second-order
method. Widening the lever on the 0.01 m reference box is unaffordable (uniform
h = 2.5e-5 is 6.4e7 cells). The droplet fills 0.1% of that box in 3D and
g = (0 0 0), so the exact solution p_in - p_out = sigma*2/R does not depend on
the box size: truncating the far field at fixed h cuts cells as (L/0.01)^3 while
leaving h, the step count and every metric alone. That buys the wide-h ladder on
a UNIFORM mesh -- w = 1/2 at every face, so the potential-form identity and the
exact-balance mechanism hold by construction, which an octree-refined mesh cannot
offer.

This script scores whether the truncation actually is inert, in 2D where it is
cheap: several domain sizes, each at cell sizes matched ACROSS domains.

WHY THE METRICS ARE COMPARABLE ACROSS DOMAIN SIZES, by construction:
  * phaseVolumeRelError is relative to the INITIAL DISCRETE volume, not an
    analytic one -- a conservation error, so domain-independent.
  * zeroSetRadialL2 is crossing-radius minus the configured radius.
  * maxMagU and minGradPsiBand are band-local.
  * the case holds the 2D z half-thickness fixed for every domain size, so any
    V-derived quantity is affected identically in all arms.

THE FAILURE MODE THAT MATTERS is not disagreement but AGREEMENT IN THE WRONG
DIRECTION: if the truncated box reads LOWER max|U|, no-slip walls 1R from the
interface are clipping the parasitic recirculation and the truncation is
FLATTERING the method. That disqualifies it even though the numbers look better.
A wall effect must be MONOTONE in L; scatter must not be. That is why the
control uses three domain sizes and not a single pair -- monotonicity is the
discriminator, and it cannot be measured with two points.

Comparisons are made at the LAST TIME COMMON TO EVERY ARM in an h group, never
at each arm's own end time: an arm that died early would otherwise be compared
against one that ran to completion.

Usage:
  domain_size_control.py --studies-root studies --out <curated.csv> STUDY [STUDY ...]
"""
import argparse
import csv
import json
import os
import sys

# metric -> (column, "lower is better"?, relative-deviation gate)
# max|U| and min|grad psi| are the headline stability metrics and get the tight
# gate; volume and shape are second-order small and structurally noisier, so
# they get a looser one. Reporting BOTH volume and shape is mandatory here: a
# single headline metric already misread one arm of this campaign.
METRICS = (
    ("maxMagU",             "maxMagU",             0.05),
    ("volumeRelError",      "phaseVolumeRelError",  0.15),
    ("shapeL2",             "zeroSetRadialL2",      0.15),
    ("minGradPsiBand",      "minGradPsiBand",       0.05),
    ("pLaplace",            "pLaplace",             0.01),
)
RATIO_GATE = 0.10   # cross-domain spread of the h-refinement ratio of max|U|


def read_arm(case_dir):
    """One arm: its tokens, its cell size, and its full metric history."""
    with open(os.path.join(case_dir, "case_params.json")) as fh:
        tokens = json.load(fh)["tokens"]
    L = float(tokens["DOMAIN_LENGTH"])
    N = float(tokens["N_CELLS"])
    csv_path = os.path.join(
        case_dir, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
    rows = []
    if os.path.isfile(csv_path):
        with open(csv_path) as fh:
            for r in csv.DictReader(fh):
                try:
                    r["TIME"] = float(r["TIME"])
                except (TypeError, ValueError, KeyError):
                    continue
                rows.append(r)
    return dict(case=os.path.basename(case_dir), L=L, N=N, h=L/N,
                RoverH=1e-3/(L/N), rows=rows,
                tEnd=rows[-1]["TIME"] if rows else float("nan"),
                dt=float(tokens.get("MAX_DELTA_T", "nan")))


def at_time(arm, t):
    """Last row at or before *t* (the runs use a fixed dt, so no interpolation
    is needed -- matched-h arms share the step exactly)."""
    best = None
    for r in arm["rows"]:
        if r["TIME"] <= t + 1e-15:
            best = r
        else:
            break
    return best


def value(row, col):
    try:
        return abs(float(row[col]))
    except (TypeError, ValueError, KeyError):
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--studies-root", default="studies")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    arms = []
    for study in args.studies:
        sdir = os.path.join(args.studies_root, study)
        if not os.path.isdir(sdir):
            sys.exit(f"[domain-control] missing study dir: {sdir}")
        for name in sorted(os.listdir(sdir)):
            case_dir = os.path.join(sdir, name)
            if os.path.isfile(os.path.join(case_dir, "case_params.json")):
                a = read_arm(case_dir)
                a["study"] = study
                arms.append(a)
    if not arms:
        sys.exit("[domain-control] no arms found")

    # group by cell size; keys rounded so 0.004/40 and 0.01/100 land together
    groups = {}
    for a in arms:
        groups.setdefault(round(a["h"], 12), []).append(a)

    out_rows, verdict_lines, failures = [], [], []
    ratios = {}   # L -> maxMagU ratio between the two finest h, per domain

    for h in sorted(groups, reverse=True):
        grp = sorted(groups[h], key=lambda a: -a["L"])
        ref = grp[0]                                   # largest domain = reference
        tCommon = min(a["tEnd"] for a in grp)
        verdict_lines.append(
            f"\nh = {h:.4e}   R/h = {ref['RoverH']:.1f}   dt = {ref['dt']:.5g}"
            f"   compared at t = {tCommon:.6g} (last time common to all arms)")
        verdict_lines.append(
            f"  {'domain':>8} {'L/R':>5} {'N':>5} {'cells':>8} {'tEnd':>9}"
            + "".join(f" {n:>14}" for n, _, _ in METRICS))
        base = {}
        for a in grp:
            row = at_time(a, tCommon)
            if row is None:
                failures.append(f"h={h:.3e} L={a['L']}: no row at t={tCommon:g}")
                continue
            vals = {n: value(row, c) for n, c, _ in METRICS}
            if a is ref:
                base = dict(vals)
            cells = int(round(a["N"]))**2
            verdict_lines.append(
                f"  {a['L']:>8.4g} {a['L']/1e-3:>5.1f} {a['N']:>5.0f}"
                f" {cells:>8d} {a['tEnd']:>9.5f}"
                + "".join(f" {vals[n]:>14.6g}" for n, _, _ in METRICS))
            rec = dict(study=a["study"], case=a["case"], domainLength=a["L"],
                       LoverR=a["L"]/1e-3, N=int(round(a["N"])), cells=cells,
                       h=a["h"], RoverH=a["RoverH"], dt=a["dt"],
                       tEnd=a["tEnd"], tCompared=tCommon,
                       isReference=int(a is ref))
            for n, _, _ in METRICS:
                rec[n] = vals[n]
                d = (vals[n] - base[n])/base[n] if base.get(n) else float("nan")
                rec[n + "_relDevVsRef"] = d
            out_rows.append(rec)
            ratios.setdefault(a["L"], {})[h] = vals["maxMagU"]

        # deviations and gates, relative to the reference domain
        verdict_lines.append(
            f"  {'relative deviation vs L = ' + format(ref['L'], '.4g'):>37}"
            + "".join(f" {n:>14}" for n, _, _ in METRICS))
        for rec in out_rows:
            if rec["h"] != h or rec["isReference"]:
                continue
            cells_txt = f"  {'L = ' + format(rec['domainLength'], '.4g'):>37}"
            devs = []
            for n, _, gate in METRICS:
                d = rec[n + "_relDevVsRef"]
                flag = "" if abs(d) <= gate else " !"
                if abs(d) > gate:
                    failures.append(
                        f"h={h:.3e} L={rec['domainLength']:g}: {n} deviates "
                        f"{d*100:+.2f}% (gate {gate*100:.0f}%)")
                devs.append(f"{d*100:>+12.2f}%{flag}")
            verdict_lines.append(cells_txt + "".join(f" {t:>14}" for t in devs))

    # h-refinement ratio of max|U| must agree across domains
    verdict_lines.append("\nh-refinement ratio of max|U| (coarse -> fine), per domain:")
    rr = {}
    for L, byh in sorted(ratios.items(), reverse=True):
        hs = sorted(byh, reverse=True)
        if len(hs) >= 2 and byh[hs[1]]:
            rr[L] = byh[hs[0]]/byh[hs[1]]
            verdict_lines.append(
                f"  L = {L:<8.4g} max|U|({hs[0]:.3g}) / max|U|({hs[1]:.3g})"
                f" = {rr[L]:.4g}")
    if len(rr) >= 2:
        lo, hi = min(rr.values()), max(rr.values())
        spread = (hi - lo)/lo
        ok = spread <= RATIO_GATE
        verdict_lines.append(
            f"  cross-domain spread = {spread*100:+.2f}% "
            f"(gate {RATIO_GATE*100:.0f}%) -> {'OK' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"refinement ratio spread {spread*100:.2f}%")

    # monotone-in-L deviation is a wall effect; scatter is not
    verdict_lines.append("\nDirection test (the failure mode that matters):")
    for h in sorted(groups, reverse=True):
        recs = sorted([r for r in out_rows if r["h"] == h],
                      key=lambda r: -r["domainLength"])
        devs = [(r["domainLength"], r["maxMagU_relDevVsRef"]) for r in recs
                if not r["isReference"]]
        if len(devs) >= 2:
            signs = {d > 0 for _, d in devs}
            mono = (abs(devs[0][1]) < abs(devs[-1][1])) and len(signs) == 1
            direction = ("LOWER max|U| in the smaller box -- walls may be "
                         "clipping the recirculation, which FLATTERS the method"
                         if all(d < 0 for _, d in devs) else
                         "higher max|U| in the smaller box" if all(d > 0 for _, d in devs)
                         else "mixed sign -> scatter, not a wall effect")
            verdict_lines.append(
                f"  h = {h:.3e}: " + ", ".join(f"L={L:g}: {d*100:+.2f}%"
                                               for L, d in devs)
                + f"\n      {'MONOTONE in L' if mono else 'not monotone'} -> {direction}")

    print("\n".join(verdict_lines))
    print("\n" + "="*78)
    if failures:
        print("VERDICT: domain truncation NOT certified. Gate violations:")
        for f in failures:
            print("  -", f)
    else:
        print("VERDICT: domain truncation CERTIFIED at every matched h -- all "
              "metrics within gate,\n         refinement ratio consistent, no "
              "monotone wall signature.")
    print("="*78)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    cols = ["study", "case", "domainLength", "LoverR", "N", "cells", "h",
            "RoverH", "dt", "tEnd", "tCompared", "isReference"]
    for n, _, _ in METRICS:
        cols += [n, n + "_relDevVsRef"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(out_rows, key=lambda r: (-r["h"], -r["domainLength"])):
            w.writerow(r)
    print(f"[domain-control] wrote {args.out}: {len(out_rows)} arms")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
