#!/usr/bin/env python3
"""The summary of the fork's tilt diagnostic of the curved dual advection on
the translated plane (TwoPhaseFlow branch davof-tiltdiag,
isoAdvectionTiltDiagnostics.C, key tiltDiagnostics; study
config/davof/planeTranslation3D_hexParaboloidTiltDiag.yaml, 2026-10-26).

Reads, in the case directory:

    postProcessing/tiltDiagnostics/cells_<step>.csv
        one line per interface cell (0 < alpha < 1) at the end of the step:
        step, celli, alpha, A_over_h2 (the exact plane's wetted patch area in
        the cell over h^2), a_over_h2 (|a_k|/h^2, a_k = -sum_f alpha_f S_f),
        nFaces, nTraced (faces whose pre-image carries a trace of the exact
        plane at t^n), nForeign (traced faces set by a cell other than k),
        nNoSurface (traced faces no surface set), sumAbsW, maxW (the weights
        w_f = ell_f |m_f - x*_u|/A_k), tauMeas (|tau_meas|), tauPred
        (|tau_pred|), tauDiff (|tau_meas - tau_pred|), cosine, maxNeighbourTilt
        (radians), maxNeighbourOffset (over h), selfTilt (the cell's own
        surface at t^n, -1 without one), selfOffset, exactCut (1 when the
        exact plane cuts the cell at t^{n+1}).
    postProcessing/tiltDiagnostics/ownSetsTrace.csv
        step, celli, ownSetsTrace: 1 when a face set by the cell's surface in
        that step carries a trace. The cell line of step s (the fractions at
        n+1) is joined with the flag of step s+1 (the next step's face
        update, where the surface built from those fractions sets faces).
    advectorVoF.csv (tpf_advection_metrics.py): step, L1_alphaf.

Writes, under --out (the case directory by default):

    tiltdiag_per_step.csv      analysis 1: per step the interface cells, the
                               fraction with sumAbsW > 1, max and p99 of
                               sumAbsW, median and p90 of |tau_meas|, the
                               log-log Pearson correlation and the median ratio
                               of |tau_meas| to |tau_pred| over the cells with
                               |tau_pred| > 1e-13, the fraction with cosine
                               above 0.9 (and below -0.9: the flipped sign),
                               max |tau_meas|, L1(alphaf)
    tiltdiag_amplification.csv analysis 2: per step and pooled, the median of
                               |tau_meas|/maxNeighbourTilt per bin of sumAbsW
                               (< 0.5, 0.5-1, 1-2, 2-5, > 5)
    tiltdiag_growth.csv        analysis 3: per step max |tau_meas| over all
                               cells and over the cells with sumAbsW > 1,
                               whether the maximum sits in that set, the growth
                               factor per step of max |tau_meas| and of
                               L1(alphaf), and the range of sumAbsW of the
                               cells with sumAbsW > 1
    tiltdiag_deadends.csv      analysis 4: per step and pooled, among the
                               cells with sumAbsW > 1 the fraction whose own
                               surface sets a traced face at the next step

Standard library only.

    python3 tilt_diagnostics_summary.py --case . --out .
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

BINS = [("lt0.5", 0.0, 0.5), ("0.5-1", 0.5, 1.0), ("1-2", 1.0, 2.0),
        ("2-5", 2.0, 5.0), ("gt5", 5.0, float("inf"))]


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def percentile(values, q):
    """The q-th percentile (0..100) by linear interpolation; nan when empty."""
    v = sorted(x for x in values if not math.isnan(x))
    if not v:
        return float("nan")
    if len(v) == 1:
        return v[0]
    pos = (len(v) - 1)*q/100.0
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo])*(pos - lo)


def median(values):
    return percentile(values, 50)


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x)/n
    my = sum(y)/n
    sxx = sum((a - mx)**2 for a in x)
    syy = sum((b - my)**2 for b in y)
    sxy = sum((a - mx)*(b - my) for a, b in zip(x, y))
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy/math.sqrt(sxx*syy)


def read_cells(tdir):
    """{step: [row dict with floats]} from cells_<step>.csv."""
    steps = {}
    for path in glob.glob(os.path.join(tdir, "cells_*.csv")):
        m = re.search(r"cells_(\d+)\.csv$", path)
        if not m:
            continue
        step = int(m.group(1))
        rows = []
        with open(path) as fh:
            for r in csv.DictReader(fh):
                rows.append({k: (int(v) if k in ("step", "celli", "nFaces", "nTraced",
                                                 "nForeign", "nNoSurface", "exactCut")
                                 else fnum(v)) for k, v in r.items()})
        steps[step] = rows
    return steps


def read_own_sets(tdir):
    """{(step, celli): flag}"""
    out = {}
    path = os.path.join(tdir, "ownSetsTrace.csv")
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                out[(int(r["step"]), int(r["celli"]))] = int(r["ownSetsTrace"])
            except (KeyError, ValueError):
                continue
    return out


def read_l1(case):
    """{step: L1_alphaf} from advectorVoF.csv (the converter's output); the
    areaFractionError .dat as the fallback, its rows numbered from 1."""
    out = {}
    path = os.path.join(case, "advectorVoF.csv")
    if os.path.isfile(path):
        with open(path) as fh:
            for r in csv.DictReader(fh):
                try:
                    out[int(r["step"])] = fnum(r.get("L1_alphaf", ""))
                except (KeyError, ValueError):
                    continue
        if out:
            return out
    fo = os.path.join(case, "postProcessing", "areaFractionError")
    tdirs = [d for d in glob.glob(os.path.join(fo, "*")) if os.path.isdir(d)]
    if not tdirs:
        return out
    latest = max(tdirs, key=lambda d: fnum(os.path.basename(d)))
    dat = os.path.join(latest, "areaFractionError.dat")
    if not os.path.isfile(dat):
        return out
    names = None
    i = 0
    with open(dat) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                parts = s[1:].split()
                if parts and parts[0] == "Time":
                    names = parts
                continue
            if names is None:
                continue
            parts = s.split()
            d = dict(zip(names, parts))
            i += 1
            out[i] = fnum(d.get("L1", ""))
    return out


def bin_of(w):
    for name, lo, hi in BINS:
        if lo <= w < hi:
            return name
    return BINS[-1][0]


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--case", default=".", help="case directory")
    ap.add_argument("--out", default=None, help="output directory (default: the case)")
    a = ap.parse_args(argv)
    case = a.case
    out = a.out or case
    tdir = os.path.join(case, "postProcessing", "tiltDiagnostics")
    steps = read_cells(tdir)
    if not steps:
        print(f"[tilt_diagnostics_summary] no cells_<step>.csv under {tdir}: nothing to do",
              file=sys.stderr)
        return 0
    own = read_own_sets(tdir)
    l1 = read_l1(case)
    os.makedirs(out, exist_ok=True)
    order = sorted(steps)

    # ---- analysis 1: per step ------------------------------------------------
    per_step = []
    prev_max = None
    prev_l1 = None
    growth_rows = []
    amp_rows = []
    dead_rows = []
    pooled_amp = {name: [] for name, _, _ in BINS}
    pooled_dead = [0, 0]
    for s in order:
        rows = steps[s]
        n = len(rows)
        sumw = [r["sumAbsW"] for r in rows]
        tm = [r["tauMeas"] for r in rows]
        high = [r for r in rows if r["sumAbsW"] > 1]
        sel = [r for r in rows if r["tauPred"] > 1e-13 and r["tauMeas"] > 0]
        lx = [math.log10(r["tauPred"]) for r in sel]
        ly = [math.log10(r["tauMeas"]) for r in sel]
        corr = pearson(lx, ly)
        ratio = median([r["tauMeas"]/r["tauPred"] for r in sel])
        fcos = (sum(1 for r in sel if r["cosine"] > 0.9)/len(sel)) if sel else float("nan")
        fcosneg = (sum(1 for r in sel if r["cosine"] < -0.9)/len(sel)) if sel else float("nan")
        max_tm = max(tm) if tm else float("nan")
        per_step.append({
            "step": s, "nCells": n,
            "fracSumAbsWgt1": (len(high)/n) if n else float("nan"),
            "maxSumAbsW": max(sumw) if sumw else float("nan"),
            "p99SumAbsW": percentile(sumw, 99),
            "medTauMeas": median(tm), "p90TauMeas": percentile(tm, 90),
            "maxTauMeas": max_tm,
            "nPred": len(sel), "pearsonLogLog": corr, "medianRatioMeasOverPred": ratio,
            "fracCosGt0.9": fcos, "fracCosLtMinus0.9": fcosneg,
            "medCosine": median([r["cosine"] for r in sel]),
            "L1_alphaf": l1.get(s, float("nan")),
        })

        # ---- analysis 2: the amplification against sumAbsW ------------------
        bins = {name: [] for name, _, _ in BINS}
        for r in rows:
            if r["maxNeighbourTilt"] > 0:
                amp = r["tauMeas"]/r["maxNeighbourTilt"]
                bins[bin_of(r["sumAbsW"])].append(amp)
                pooled_amp[bin_of(r["sumAbsW"])].append(amp)
        for name, _, _ in BINS:
            amp_rows.append({"step": s, "bin": name, "n": len(bins[name]),
                             "medianAmplification": median(bins[name]),
                             "p90Amplification": percentile(bins[name], 90)})

        # ---- analysis 3: the growth -----------------------------------------
        imax = max(range(n), key=lambda i: tm[i]) if n else -1
        max_high = max((r["tauMeas"] for r in high), default=float("nan"))
        g = (max_tm/prev_max) if (prev_max is not None and prev_max > 0) else float("nan")
        gl1 = (l1.get(s, float("nan"))/prev_l1) if (prev_l1 is not None and prev_l1 > 0) else float("nan")
        wlo = min((r["sumAbsW"] for r in high), default=float("nan"))
        whi = max((r["sumAbsW"] for r in high), default=float("nan"))
        growth_rows.append({
            "step": s, "maxTauMeasAll": max_tm, "maxTauMeasHighW": max_high,
            "argmaxSumAbsW": rows[imax]["sumAbsW"] if imax >= 0 else float("nan"),
            "argmaxCell": rows[imax]["celli"] if imax >= 0 else -1,
            "maxInHighW": int(imax >= 0 and rows[imax]["sumAbsW"] > 1),
            "growthMaxTau": g, "nHighW": len(high),
            "minSumAbsWHigh": wlo, "maxSumAbsWHigh": whi,
            "growthInRange": int(not math.isnan(g) and not math.isnan(wlo) and wlo <= g <= whi),
            "L1_alphaf": l1.get(s, float("nan")), "growthL1": gl1,
        })
        prev_max = max_tm
        prev_l1 = l1.get(s, None)

        # ---- analysis 4: the dead ends --------------------------------------
        nh = 0
        nset = 0
        nknown = 0
        for r in high:
            nh += 1
            flag = own.get((s + 1, r["celli"]))
            if flag is None:
                continue
            nknown += 1
            nset += flag
        dead_rows.append({"step": s, "nHighW": nh, "nJoined": nknown,
                          "fracOwnSetsTraceNextStep": (nset/nknown) if nknown else float("nan")})
        pooled_dead[0] += nknown
        pooled_dead[1] += nset

    for name, _, _ in BINS:
        amp_rows.append({"step": "all", "bin": name, "n": len(pooled_amp[name]),
                         "medianAmplification": median(pooled_amp[name]),
                         "p90Amplification": percentile(pooled_amp[name], 90)})
    dead_rows.append({"step": "all", "nHighW": sum(r["nHighW"] for r in dead_rows),
                      "nJoined": pooled_dead[0],
                      "fracOwnSetsTraceNextStep": (pooled_dead[1]/pooled_dead[0]) if pooled_dead[0] else float("nan")})

    def write(name, rows):
        if not rows:
            return
        path = os.path.join(out, name)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in rows:
                w.writerow({k: (f"{v:.6e}" if isinstance(v, float) else v) for k, v in r.items()})
        print(f"[tilt_diagnostics_summary] wrote {path} ({len(rows)} rows)")

    write("tiltdiag_per_step.csv", per_step)
    write("tiltdiag_amplification.csv", amp_rows)
    write("tiltdiag_growth.csv", growth_rows)
    write("tiltdiag_deadends.csv", dead_rows)

    # the one-screen summary for the log
    print("[tilt_diagnostics_summary] steps", len(order),
          "cells/step", per_step[0]["nCells"], "..", per_step[-1]["nCells"])
    corr = [r["pearsonLogLog"] for r in per_step[1:] if not math.isnan(r["pearsonLogLog"])]
    rat = [r["medianRatioMeasOverPred"] for r in per_step[1:] if not math.isnan(r["medianRatioMeasOverPred"])]
    fc = [r["fracCosGt0.9"] for r in per_step[1:] if not math.isnan(r["fracCosGt0.9"])]
    fcn = [r["fracCosLtMinus0.9"] for r in per_step[1:] if not math.isnan(r["fracCosLtMinus0.9"])]
    if corr:
        print(f"  P-TD1 correlation log-log: min {min(corr):.3f} median {median(corr):.3f} max {max(corr):.3f};"
              f" median ratio meas/pred: min {min(rat):.3g} median {median(rat):.3g} max {max(rat):.3g};"
              f" frac cos>0.9: min {min(fc):.3f} median {median(fc):.3f}; frac cos<-0.9 median {median(fcn):.3f}")
    pooled = [r for r in amp_rows if r["step"] == "all"]
    print("  P-TD2 pooled median amplification per bin:",
          ", ".join(f"{r['bin']}: {r['medianAmplification']:.3g} (n={r['n']})" for r in pooled))
    inset = [r["maxInHighW"] for r in growth_rows]
    inrange = [r["growthInRange"] for r in growth_rows[1:]]
    print(f"  P-TD3 max in high-W set: {sum(inset)}/{len(inset)} steps;"
          f" growth of max|tau| within the sumAbsW range of those cells: {sum(inrange)}/{len(inrange)} steps")
    d = dead_rows[-1]
    print(f"  P-TD4 dead ends: among {d['nJoined']} joined high-W cell-steps the fraction whose own surface"
          f" sets a traced face at the next step is {d['fracOwnSetsTraceNextStep']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
