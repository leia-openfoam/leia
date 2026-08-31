#!/usr/bin/env python3
"""Curate the projectedFlux-vs-cellCentred refinement ladders.

Studies: config/projFluxStationary2D (N = 32..256) and
         config/projFluxOscillating2D (N = 32..128, moving interface).

Both were run with each rung carrying its OWN matched cellCentred baseline on
one commit, because the t_blow baseline was retracted 2026-08-31 and no
historical number is a valid reference.

Reported per arm, never max|U| alone (standing rule): the parasitic level, the
volume and shape errors, the level-set gradient health and the grid-scale
corrugation A2h.

TWO NORMALISATION TRAPS this script handles explicitly, both of which produced
a wrong reading before they were understood:

  * min|grad psi| IS NOT COMPARABLE TO 1 IN EVERY CASE. The stationary case
    initialises psi as a signed distance (|grad psi| ~ 1), but the oscillating
    case initialises it as an ALGEBRAIC ELLIPSE, (x/a)^2 + (y/b)^2 - 1, whose
    gradient is ~2R/a^2 ~ 1650. A raw value of 1035 there is healthy, not
    broken. The metric is therefore reported as a RATIO to the arm's own t = 0
    value.
  * A DIVERGED arm's last CSV row is garbage, and a live run's last row is
    half-written. Rows are filtered to those with every scored field finite,
    and the arm's status and reached time are reported alongside so a short run
    is never silently scored as if it were a full one.

Usage: make_trace_refinement_table.py [--studies-root DIR] [-o OUT_CSV]
"""
import argparse, csv, json, math, os, sys

METRICS = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"
SCORED = ["TIME", "maxMagU", "phaseVolumeRelError", "zeroSetRadialL2",
          "minGradPsiBand", "A2hMaxBand"]
STUDIES = ["projFluxStationary2D", "projFluxOscillating2D"]
T_END = 0.1


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def read_arm(d):
    mp, cp = os.path.join(d, METRICS), os.path.join(d, "case_params.json")
    if not (os.path.exists(mp) and os.path.exists(cp)):
        return None
    raw = list(csv.DictReader(open(mp)))
    rows = [r for r in raw
            if all(r.get(k) not in (None, "") and math.isfinite(fnum(r[k]))
                   for k in SCORED)]
    if len(rows) < 10:
        return None
    meta = json.load(open(cp))
    tok = meta.get("tokens", meta)
    first, last = rows[0], rows[-1]
    g0, g1 = fnum(first["minGradPsiBand"]), fnum(last["minGradPsiBand"])
    a0, a1 = fnum(first["A2hMaxBand"]), fnum(last["A2hMaxBand"])
    t_last = fnum(last["TIME"])
    return dict(
        arm=os.path.basename(d), commit=meta.get("gitCommit", "?"),
        N=int(fnum(tok.get("N_CELLS", "nan"))),
        trace=tok.get("SL_TRACE_VELOCITY", "?"),
        dt=fnum(tok.get("MAX_DELTA_T", "nan")),
        steps=len(rows) - 1, t_last=t_last,
        reached=t_last / T_END,
        complete=t_last >= 0.995 * T_END,
        maxU=fnum(last["maxMagU"]),
        volErr=fnum(last["phaseVolumeRelError"]),
        shapeL2=fnum(last["zeroSetRadialL2"]),
        gradPsi0=g0, gradPsiEnd=g1,
        gradPsiRatio=(g1 / g0) if g0 else float("nan"),
        A2h0=a0, A2hEnd=a1,
        A2hGrowth=(a1 / a0) if a0 else float("nan"),
    )


def order(fine, coarse):
    """Observed order between two rungs whose h differs by 2."""
    if not (fine > 0 and coarse > 0):
        return float("nan")
    return math.log(coarse / fine) / math.log(2.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies-root", default="studies")
    ap.add_argument("-o", "--out", default=os.path.join(
        "docs/method-comparison/method-comparison-article/data/tables",
        "trace_refinement.csv"))
    a = ap.parse_args()

    all_arms = []
    for study in STUDIES:
        sd = os.path.join(a.studies_root, study)
        if not os.path.isdir(sd):
            print("skip (absent): %s" % sd, file=sys.stderr)
            continue
        for name in sorted(os.listdir(sd)):
            d = os.path.join(sd, name)
            if os.path.isdir(d):
                r = read_arm(d)
                if r:
                    r["study"] = study
                    all_arms.append(r)
    if not all_arms:
        sys.exit("no scoreable arms under %s" % a.studies_root)

    all_arms.sort(key=lambda x: (x["study"], x["N"], x["trace"]))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    cols = ["study", "arm", "N", "trace", "dt", "steps", "t_last", "reached",
            "complete", "maxU", "volErr", "shapeL2", "gradPsi0", "gradPsiEnd",
            "gradPsiRatio", "A2h0", "A2hEnd", "A2hGrowth", "commit"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for x in all_arms:
            w.writerow(x)

    for study in STUDIES:
        arms = [x for x in all_arms if x["study"] == study]
        if not arms:
            continue
        print("\n" + "=" * 100)
        print(study)
        print("=" * 100)
        print("%-5s %-15s %8s %9s %8s %11s %11s %11s %9s %8s" % (
            "N", "trace", "steps", "reached", "status", "max|U|", "volErr",
            "shapeL2", "gPsi/g0", "A2h x"))
        for x in arms:
            print("%-5d %-15s %8d %8.1f%% %8s %11.4e %11.4e %11.4e %9.3f %8.1f"
                  % (x["N"], x["trace"], x["steps"], 100 * x["reached"],
                     "full" if x["complete"] else "SHORT",
                     x["maxU"], x["volErr"], x["shapeL2"],
                     x["gradPsiRatio"], x["A2hGrowth"]))

        # Convergence is only meaningful between arms that reached the horizon.
        for tr in sorted({x["trace"] for x in arms}):
            s = sorted([x for x in arms if x["trace"] == tr], key=lambda x: x["N"])
            full = [x for x in s if x["complete"]]
            if len(full) < 2:
                print("  %-15s only %d full-horizon rung(s): no order fitted"
                      % (tr, len(full)))
                continue
            print("  %-15s observed order between full-horizon rungs:" % tr)
            for i in range(len(full) - 1):
                c, f = full[i], full[i + 1]
                print("    N %3d->%3d   max|U| %6.2f   volume %6.2f   shape %6.2f"
                      % (c["N"], f["N"],
                         order(f["maxU"], c["maxU"]),
                         order(f["volErr"], c["volErr"]),
                         order(f["shapeL2"], c["shapeL2"])))
            short = [x for x in s if not x["complete"]]
            if short:
                print("    EXCLUDED (did not reach the horizon): "
                      + ", ".join("N=%d at %.1f%%" % (x["N"], 100 * x["reached"])
                                  for x in short))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
