#!/usr/bin/env python3
"""The per-step amplification g = r*dt: the parasitic instability's order parameter.

WHY THIS AND NOT A GROWTH RATE, AND NOT t_blow.

The capillary law ties the step to the cell size, dt = CAPILLARY_DT_COEFF/nRef^1.5,
so a refinement ladder changes the mesh AND the step together. Measured on the 3D
L = 6R ladder, the fitted growth rate of max|U| over t = 0.05..0.1 rises
40.47 -> 52.32 -> 90.55 1/s at R/h = 10.0/12.7/15.8 -- a factor 2.2 that reads as
"finer is less stable". But multiplied by each arm's own step:

    r*dt = 4.395e-04, 3.986e-04, 4.936e-04     constant to 24%

The amplification PER STEP does not move. r ~ 1/dt, and the e-folds accumulated by a
fixed physical time (4.05, 5.23, 9.05) grow only because a finer mesh takes more steps
to get there. So the rate is not a property of the discretisation; the per-step gain
is. An h-independent per-step gain also cannot be a spatial truncation error -- it is
the signature of something that does not refine: an operator splitting, a
once-per-step lag, a fixed relative error, or a solver tolerance.

t_blow was abandoned earlier for a related reason: the e-fold count at blow-up varies
5.4-13.3 across the matrix, so it measures the horizon, not the mechanism.

WHAT IT IS FOR. Two things:

  1. Scoring. g is dt-normalised and dimension-agnostic, so it compares 2D against
     3D, stationary against oscillating, and any delivery against any other -- which
     a growth rate cannot.
  2. A GATE. --gate G exits non-zero if g exceeds G, so a candidate can be rejected
     from a ~200-step run instead of a 12-hour ladder. The specification comes from
     working backwards: at 3D N_L = 120, t = 0.1 is 26042 steps, so keeping the
     accumulated amplification O(1) needs g <~ 4e-5, against the measured 4.4e-4 --
     a 10x reduction. A lever that improves the delivered force by two orders and
     leaves g alone is worthless for stability, which is precisely what the
     cell-centre inverse did (+2.03 on the delivered non-gradient content, no
     stability gain).

TWO ESTIMATORS, both reported, because max|U| is not monotone. The 2D N=64 arm
excursions to 3.9e-2 at t = 0.05 and recovers to 2.2e-3, which makes a single
exponential fit meaningless there (its R^2 is 0.73 and the fitted rate is negative).
  * gFit   -- least squares on ln(y) against t over the late window, with R^2 so a
              bad fit is visible rather than silent.
  * gMed   -- median of the per-step increment of ln(y), over the same window.
              Insensitive to excursions; this is the one to trust when R^2 is low.
  * gAvg   -- ln(y_end/y_0) divided by the TOTAL number of steps. No window, no fit,
              no model. This is the estimator to quote, because it is the one that
              predicts the end state exactly and without assumption:
                  y_end = y_0 * exp(gAvg * nSteps)
              which is an identity, so all the physics sits in how y_0 and gAvg
              depend on h. Verified against the ladder: 3D N_L = 95 starts at
              max|U| = 4.091e-05 with 7.45 e-folds -> 7.06e-02 against 7.04e-02
              measured; 2D N = 512 starts at 4.672e-07 with 4.69 e-folds ->
              5.08e-05 against 5.10e-05 measured.
A large gap between gFit and gMed means the trajectory is not exponential, and
neither should be quoted without the other. gAvg is immune to that.

SCORED ON TWO SIGNALS. max|U| is the headline, but the quantity the plan identifies
as causal is the delivered non-gradient content alpha_f*snGrad(kappa_c) -- the only
part of the capillary force the pressure solve cannot absorb (dropping it takes
max|U| from 3.6e-5 to 1.3e-9). Where the solver writes it (driverAcrossSupportL2),
g is reported for it too.

Usage:
  per_step_gain.py STUDY_DIR [STUDY_DIR ...] [--window 0.5] [--gate 4e-5]
                             [--signal maxMagU] [--out curated.csv]
"""
import argparse
import csv
import glob
import json
import math
import os
import sys
import statistics

METRIC_CSVS = ("leiaSemiLagrangianLevelSetTwoPhaseFoam.csv", "leiaLevelSetFoam.csv",
               "interFoam.csv", "interFlow.csv")
SIGNALS = ("maxMagU", "driverAcrossSupportL2")


def find_metrics(case_dir):
    for n in METRIC_CSVS:
        p = os.path.join(case_dir, n)
        if os.path.isfile(p) and os.path.getsize(p) > 200:
            return p
    return None


def read_rows(path, need):
    """Rows with every needed column parseable. A DIVERGED arm's final row is
    truncated mid-write, so a half-written record must never be scored."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if any(r.get(c) in (None, "") for c in need):
                continue
            try:
                out.append({c: float(r[c]) for c in need})
            except (TypeError, ValueError):
                continue
    return out


def fit_rate(ts, ys):
    """Least squares r in ln(y) = c + r t. Returns (r, R^2, n)."""
    pts = [(t, math.log(y)) for t, y in zip(ts, ys) if y > 0]
    n = len(pts)
    if n < 3:
        return float("nan"), float("nan"), n
    mx = sum(t for t, _ in pts)/n
    my = sum(v for _, v in pts)/n
    sxx = sum((t-mx)**2 for t, _ in pts)
    if sxx <= 0:
        return float("nan"), float("nan"), n
    r = sum((t-mx)*(v-my) for t, v in pts)/sxx
    c = my - r*mx
    ssr = sum((v-(c+r*t))**2 for t, v in pts)
    sst = sum((v-my)**2 for _, v in pts)
    return r, (1-ssr/sst if sst > 0 else float("nan")), n


def median_step_gain(ts, ys):
    """Median per-step increment of ln(y). Robust to excursions."""
    inc = [math.log(ys[i+1]/ys[i])
           for i in range(len(ys)-1) if ys[i] > 0 and ys[i+1] > 0]
    return statistics.median(inc) if inc else float("nan")


def score(case_dir, window, signal):
    mp = find_metrics(case_dir)
    pp = os.path.join(case_dir, "case_params.json")
    if mp is None or not os.path.isfile(pp):
        return None
    tok = json.load(open(pp))["tokens"]
    rows = read_rows(mp, ["TIME", signal])
    if len(rows) < 6:
        return None
    dt = float(tok.get("MAX_DELTA_T", "nan"))
    tEnd = rows[-1]["TIME"]
    t0 = tEnd*(1.0 - window)
    late = [r for r in rows if r["TIME"] >= t0]
    if len(late) < 6:
        late = rows
    ts = [r["TIME"] for r in late]
    ys = [abs(r[signal]) for r in late]
    r, r2, n = fit_rate(ts, ys)
    gmed = median_step_gain(ts, ys)
    y0, yE = abs(rows[0][signal]), abs(rows[-1][signal])
    L = float(tok.get("DOMAIN_LENGTH", 0.01))
    N = float(tok.get("N_CELLS", "nan"))
    return dict(
        case=os.path.basename(case_dir), signal=signal,
        N=int(N) if N == N else -1, L=L, h=(L/N if N == N and N else float("nan")),
        RoverH=(1e-3/(L/N) if N == N and N else float("nan")),
        dt=dt, nSteps=len(rows), tEnd=tEnd,
        rFit=r, R2=r2, nWindow=n,
        gFit=r*dt if r == r else float("nan"),
        gMed=gmed,
        gAvg=((math.log(yE/y0)/len(rows)) if y0 > 0 and yE > 0 and len(rows) else float("nan")),
        eFolds=(math.log(yE/y0) if y0 > 0 and yE > 0 else float("nan")),
        yStart=y0, yEnd=yE,
        delivery=tok.get("CURVATURE_EXTENSION", ""),
        filt=tok.get("PSI_FILTER", ""), theta=tok.get("PSI_FILTER_THETA", ""),
        K=tok.get("CURVATURE_INVERSE_GAUSSIAN", ""),
        outerCorr=tok.get("PSI_OUTER_CORRECTORS", ""),
        stfModel=tok.get("STF_MODEL", ""),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--window", type=float, default=0.5,
                    help="trailing fraction of the run to fit (default 0.5)")
    ap.add_argument("--gate", type=float, default=None,
                    help="exit non-zero if any arm's gMed exceeds this")
    ap.add_argument("--signal", default=None, help="score only this column")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    signals = [args.signal] if args.signal else list(SIGNALS)

    recs = []
    for sd in args.studies:
        for cand in sorted(glob.glob(os.path.join(sd, "*"))):
            if not os.path.isdir(cand):
                continue
            for sig in signals:
                r = score(cand, args.window, sig)
                if r:
                    r["study"] = os.path.basename(os.path.normpath(sd))
                    recs.append(r)
    if not recs:
        sys.exit("[gain] no scorable arms found")

    for sig in signals:
        sub = [r for r in recs if r["signal"] == sig]
        if not sub:
            continue
        print(f"\n=== per-step gain on {sig}   (window = last {args.window:.0%} of each run)")
        print(f"  {'study':<30} {'N':>4} {'R/h':>5} {'dt':>10} {'steps':>7} "
              f"{'tEnd':>8} {'gFit':>10} {'R2':>6} {'gMed':>10} {'gAvg':>10} {'e-folds':>8}  cfg")
        for r in sorted(sub, key=lambda r: (r["study"], r["N"])):
            cfg = " ".join(filter(None, [
                r["delivery"], (f"{r['filt']}@{r['theta']}" if r["filt"] not in ("", "none") else ""),
                (f"K={r['K']}" if r["K"] not in ("", "yes") else ""),
                (f"outer={r['outerCorr']}" if r["outerCorr"] not in ("", "no") else ""),
                r["stfModel"]]))
            print(f"  {r['study']:<30} {r['N']:>4} {r['RoverH']:>5.1f} {r['dt']:>10.4g} "
                  f"{r['nSteps']:>7} {r['tEnd']:>8.5f} {r['gFit']:>+10.3e} {r['R2']:>6.3f} "
                  f"{r['gMed']:>+10.3e} {r['gAvg']:>+10.3e} {r['eFolds']:>+8.2f}  {cfg}")

    if args.out:
        cols = list(recs[0].keys())
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(recs)
        print(f"\n[gain] wrote {args.out}: {len(recs)} rows")

    if args.gate is not None:
        # Gate on gAvg: window-free, fit-free, and the estimator that predicts the
        # end state, so it cannot be gamed by a lucky fit window.
        bad = [r for r in recs if r["signal"] == "maxMagU" and r["gAvg"] > args.gate]
        print(f"\n[gain] GATE g <= {args.gate:.3g} on maxMagU: "
              f"{len(recs)//max(len(signals),1)} arms, {len(bad)} over the gate")
        for r in bad:
            print(f"    OVER  {r['study']}/{r['case']}  gAvg={r['gAvg']:+.3e}")
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
