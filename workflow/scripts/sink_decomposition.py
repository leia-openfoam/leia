#!/usr/bin/env python3
"""Split the net per-step gain into a measured sink and an inferred amplifier.

The stationary-droplet net gain is

    A(h) = -sink(h) + amplifier(h)

and it changes sign under refinement in both dimensions, which makes every local
logarithmic slope of A an inflated estimate of the amplifier's exponent and makes
a 3-parameter additive fit to 4 points unidentifiable (it degenerates outright in
2D). `sink_decay_restart.sh` measures sink(h) directly by restarting each arm
from its written field -- so the leading eigenmode is already established -- with
sigma = 0, which removes the interface-to-momentum feedback and leaves only
viscous diffusion plus the pressure projection.

This script reads the parent study and the two restart studies and reports

  * the CONTROL GATE: the sigma-unchanged restart's gain against the parent's own
    gain over the SAME time window. A restart loses the multi-level time history
    (OpenFOAM's `backward` restarts on a single level, so the first step after a
    restart is Euler), so this gate is what licenses the sigma = 0 numbers. If it
    fails, nothing else here means anything.
  * sink(h), with both an endpoint and a least-squares estimator, and the fitted
    exponent against the two mechanistic predictions h^+1.5 (geometry-set mode,
    k ~ 1/R) and h^-0.5 (mesh-set mode, k ~ 1/h).
  * amplifier(h) = A(h) + sink(h) on the parent's own window, and its exponent --
    the quantity a fix has to change.
  * volume AND shape error over the restart window, never a single metric. Note
    that phaseVolumeRelError is referenced to the RESTART state (the solver takes
    the startup volume as its reference), so `zeroSetRadialL2`, which is
    referenced to the analytic shape, is the trustworthy geometric column here.

Usage:
  sink_decomposition.py --parent studies/filterOffAmplifier3D \
                        --sink studies/sinkDecay3D --control studies/sinkControl3D \
                        [--radius 1e-3] [--out <curated.csv>]
"""

import argparse
import csv
import json
import math
import os


def arms(study):
    """Load every arm of a study as (tokens, rows), keyed by cell count."""
    out = {}
    if not study or not os.path.isdir(study):
        return out
    for name in sorted(os.listdir(study)):
        case = os.path.join(study, name)
        csvf = os.path.join(case, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
        parf = os.path.join(case, "case_params.json")
        if not (os.path.isfile(csvf) and os.path.isfile(parf)):
            continue
        tok = json.load(open(parf))["tokens"]
        need = ("TIME", "maxMagU", "phaseVolumeRelError", "zeroSetRadialL2")
        rows = [r for r in csv.DictReader(open(csvf))
                if all(r.get(k) not in (None, "") for k in need)]
        if len(rows) < 20:
            continue
        out.setdefault(int(float(tok["N_CELLS"])), []).append((name, tok, rows))
    return out


def gains(rows, tlo=None, thi=None):
    """Endpoint and least-squares per-step gain of max|U| over a time window."""
    sel = [r for r in rows
           if (tlo is None or float(r["TIME"]) >= tlo - 1e-14)
           and (thi is None or float(r["TIME"]) <= thi + 1e-14)]
    if len(sel) < 10:
        return None
    u = [abs(float(r["maxMagU"])) for r in sel]
    if min(u) <= 0:
        return None
    n = len(u)
    y = [math.log(v) for v in u]
    # least squares slope of ln u against step index
    xm = (n - 1) / 2.0
    ym = sum(y) / n
    sxy = sum((i - xm) * (y[i] - ym) for i in range(n))
    sxx = sum((i - xm) ** 2 for i in range(n))
    gfit = sxy / sxx if sxx > 0 else float("nan")
    return dict(
        n=n, t0=float(sel[0]["TIME"]), t1=float(sel[-1]["TIME"]),
        u0=u[0], u1=u[-1],
        gAvg=(y[-1] - y[0]) / (n - 1), gFit=gfit,
        vol=abs(float(sel[-1]["phaseVolumeRelError"])),
        shape=abs(float(sel[-1]["zeroSetRadialL2"])),
    )


def order(pairs):
    """Fitted exponent of a quantity against h, pairwise and least squares."""
    out = []
    for i in range(1, len(pairs)):
        (h0, v0), (h1, v1) = pairs[i - 1], pairs[i]
        if v0 > 0 and v1 > 0 and h0 != h1:
            out.append(math.log(v0 / v1) / math.log(h0 / h1))
    ls = float("nan")
    pos = [(h, v) for h, v in pairs if v > 0]
    if len(pos) >= 2:
        lx = [math.log(h) for h, _ in pos]
        ly = [math.log(v) for _, v in pos]
        n = len(pos)
        xm, ym = sum(lx) / n, sum(ly) / n
        sxx = sum((x - xm) ** 2 for x in lx)
        if sxx > 0:
            ls = sum((lx[i] - xm) * (ly[i] - ym) for i in range(n)) / sxx
    return out, ls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", required=True)
    ap.add_argument("--sink", required=True)
    ap.add_argument("--control", default=None)
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    P, S, C = arms(a.parent), arms(a.sink), arms(a.control)
    if not S:
        raise SystemExit("no sink arms found in %s" % a.sink)

    print("=" * 100)
    print("CONTROL GATE -- sigma-unchanged restart vs the parent over the same window")
    if not C:
        print("  no control study given: the sigma = 0 numbers below are UNGATED.")
    for N in sorted(C):
        for name, tok, rows in C[N]:
            g = gains(rows)
            if not g:
                continue
            par = None
            for pname, ptok, prows in P.get(N, []):
                par = gains(prows, tlo=g["t0"], thi=g["t1"])
                if par:
                    break
            h = float(tok.get("DOMAIN_LENGTH", 0.01)) / N
            if par:
                d = (g["gFit"] - par["gFit"]) / abs(par["gFit"]) if par["gFit"] else float("nan")
                flag = "PASS" if abs(d) < 0.15 else "**FAIL**"
                print("  R/h %5.1f  parent gFit %+.4e (%d rows)  restart gFit %+.4e (%d rows)"
                      "  delta %+7.1f%%  %s"
                      % (a.radius / h, par["gFit"], par["n"], g["gFit"], g["n"],
                         100.0 * d, flag))
            else:
                print("  R/h %5.1f  restart gFit %+.4e over t=[%.6g, %.6g] -- no overlapping "
                      "parent rows, gate cannot be evaluated" % (a.radius / h, g["gFit"], g["t0"], g["t1"]))

    print()
    print("=" * 100)
    print("SINK, measured (sigma = 0 restart) -- and the amplifier it implies")
    hdr = ("  %5s %7s %13s %13s %13s %13s %11s %11s"
           % ("R/h", "steps", "A parent", "sink gAvg", "sink gFit",
              "amplifier=A+sink", "volume", "shape"))
    print(hdr)
    table, curated = [], []
    for N in sorted(S):
        for name, tok, rows in S[N]:
            g = gains(rows)
            if not g:
                continue
            h = float(tok.get("DOMAIN_LENGTH", 0.01)) / N
            # parent net gain over ITS OWN full record (that is the A we published)
            A = None
            for pname, ptok, prows in P.get(N, []):
                pg = gains(prows)
                if pg:
                    A = pg["gFit"]
                    break
            sink = -g["gFit"]          # decay: gFit < 0, sink is its magnitude
            amp = (A + sink) if A is not None else float("nan")
            print("  %5.1f %7d %+13.4e %+13.4e %+13.4e %+13.4e %11.3e %11.3e"
                  % (a.radius / h, g["n"], A if A is not None else float("nan"),
                     -g["gAvg"], sink, amp, g["vol"], g["shape"]))
            table.append((h, A, sink, amp))
            curated.append(dict(
                cells=N, h=h, RoverH=a.radius / h, restartSteps=g["n"],
                tStart=g["t0"], tEnd=g["t1"], u0=g["u0"], uEnd=g["u1"],
                netGainParent=A, sinkGAvg=-g["gAvg"], sinkGFit=sink,
                amplifier=amp, volumeRelErrorVsRestart=g["vol"],
                zeroSetRadialL2=g["shape"]))

    print()
    print("  sink(h) exponent -- predictions: h^+1.5 if the damped mode is set by the")
    print("  geometry (k ~ 1/R, sink ~ nu k^2 dt ~ dt ~ h^1.5), h^-0.5 if it is set by")
    print("  the mesh (k ~ 1/h). A positive exponent means damping WEAKENS on refinement.")
    pw, ls = order([(h, s) for h, _, s, _ in table if s is not None and s > 0])
    print("    pairwise %s   least squares h^%+.2f"
          % (", ".join("h^%+.2f" % v for v in pw) or "n/a", ls))
    print("  amplifier(h) exponent -- THE quantity a fix has to change:")
    pw, ls = order([(h, m) for h, _, _, m in table if m == m and m > 0])
    print("    pairwise %s   least squares h^%+.2f"
          % (", ".join("h^%+.2f" % v for v in pw) or "n/a", ls))
    print("  net A(h) local slopes, for comparison (inflated near the sign change):")
    pw, ls = order([(h, A) for h, A, _, _ in table if A is not None and A > 0])
    print("    pairwise %s   least squares h^%+.2f"
          % (", ".join("h^%+.2f" % v for v in pw) or "n/a", ls))

    if a.out and curated:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(curated[0].keys()))
            w.writeheader()
            w.writerows(curated)
        print("\n  wrote %s" % a.out)


if __name__ == "__main__":
    main()
