#!/usr/bin/env python3
"""Curate the trace-amplifier dt sweep (config/traceAmplifierDt2D.yaml).

Order parameter, on the fixed physical window [T_REF, T_END] pre-registered in
the config header (chosen to start after the startup kick every arm shows near
t = 0.0022):

    G_x(dt)  = ln( max|U|_x(T_END) / max|U|_x(T_REF) )
    dr_x(dt) = [ G_x(dt) - G_projectedFlux(dt) ] / (T_END - T_REF)

dr subtracts the physical (viscous) damping the arms share, leaving the
NUMERICAL amplifier alone.  Two questions are then read off:

  * does dr_cellCentred scale with dt   -> is the amplifier the explicit
    capillary time coupling (r = r0 + c*dt)?
  * rho = dr_reconstructedU / dr_cellCentred -> is the mechanism SOLENOIDALITY
    (rho ~ 1, since reconstructedU is smoothed but not divergence-free) or the
    reconstruct SMOOTHING (rho ~ 0)?

Usage:  make_trace_amplifier_table.py [STUDY_DIR] [-o OUT_CSV]
"""
import argparse, csv, json, math, os, sys

T_REF, T_END = 0.02, 0.1
REF_TRACE = "projectedFlux"
METRICS = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"


def fnum(v):
    """float() that tolerates a missing/short field (nan), never raises."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def complete(rows, keys):
    """Drop truncated rows.

    A live run's last CSV line is usually half-written -- csv.DictReader then
    yields None for the fields past the truncation and float() raises. Scoring
    a partial row would also be wrong even when it parses, so require every
    field this script reads to be present and numeric.
    """
    ok = []
    for r in rows:
        if all(r.get(k) not in (None, "") and math.isfinite(fnum(r[k]))
               for k in keys):
            ok.append(r)
    return ok


def at(rows, t, key="maxMagU"):
    """Value of `key` at the sample nearest time t, with that sample's time."""
    i = min(range(len(rows)), key=lambda j: abs(fnum(rows[j]["TIME"]) - t))
    return fnum(rows[i][key]), fnum(rows[i]["TIME"])


def read_arm(d):
    mp = os.path.join(d, METRICS)
    cp = os.path.join(d, "case_params.json")
    if not (os.path.exists(mp) and os.path.exists(cp)):
        return None
    NEEDED = ["TIME", "maxMagU", "phaseVolumeRelError", "zeroSetRadialL2",
              "minGradPsiBand"]
    raw = list(csv.DictReader(open(mp)))
    rows = complete(raw, NEEDED)
    if len(rows) < 10:
        return None
    tok = json.load(open(cp))
    u_ref, t_ref = at(rows, T_REF)
    u_end, t_end = at(rows, T_END)
    # A blown-up arm is a result: report it, do not fit through it.
    diverged = not all(math.isfinite(fnum(r["maxMagU"])) for r in rows)
    t_last = fnum(rows[-1]["TIME"])
    partial = t_last < 0.995*T_END
    G = math.log(u_end / u_ref) if (u_ref > 0 and u_end > 0) else float("nan")
    return dict(
        arm=os.path.basename(d),
        trace=tok.get("SL_TRACE_VELOCITY", "?"),
        dt=float(tok.get("MAX_DELTA_T", "nan")),
        steps=len(rows) - 1,
        t_ref=t_ref, t_end=t_end,
        u_ref=u_ref, u_end=u_end, G=G,
        r=G / (t_end - t_ref) if t_end > t_ref else float("nan"),
        volErr=fnum(rows[-1]["phaseVolumeRelError"]),
        shapeL2=fnum(rows[-1]["zeroSetRadialL2"]),
        minGradPsi=fnum(rows[-1]["minGradPsiBand"]),
        diverged=diverged, partial=partial, t_last=t_last,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study_dir", nargs="?", default="studies/traceAmplifierDt2D")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args()

    arms = []
    for name in sorted(os.listdir(a.study_dir)):
        d = os.path.join(a.study_dir, name)
        if os.path.isdir(d):
            r = read_arm(d)
            if r:
                arms.append(r)
    if not arms:
        sys.exit("no completed arms in %s" % a.study_dir)

    # dr is defined per dt against the reference trace at the SAME dt.
    ref = {round(x["dt"], 12): x["G"] for x in arms if x["trace"] == REF_TRACE}
    for x in arms:
        g0 = ref.get(round(x["dt"], 12))
        x["dr"] = ((x["G"] - g0) / (T_END - T_REF)
                   if g0 is not None and math.isfinite(x["G"]) else float("nan"))

    arms.sort(key=lambda x: (x["trace"], -x["dt"]))
    out = a.out or os.path.join(
        "docs/method-comparison/method-comparison-article/data/tables",
        "trace_amplifier_dt.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cols = ["arm", "trace", "dt", "steps", "u_ref", "u_end", "G", "r", "dr",
            "volErr", "shapeL2", "minGradPsi", "diverged", "partial", "t_last"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for x in arms:
            w.writerow(x)

    print("window [%g, %g] s, reference trace %s\n" % (T_REF, T_END, REF_TRACE))
    print("%-16s %12s %8s %10s %10s %10s" %
          ("trace", "dt", "steps", "G", "r [1/s]", "dr [1/s]"))
    for x in arms:
        print("%-16s %12.6g %8d %10.4f %10.2f %10.2f%s" %
              (x["trace"], x["dt"], x["steps"], x["G"], x["r"], x["dr"],
               "  DIVERGED" if x["diverged"]
               else ("  PARTIAL t=%.4f" % x["t_last"]) if x["partial"] else ""))

    # PREDICTION 1: dr_cellCentred halves with dt.
    print("\n-- PREDICTION 1: dr proportional to dt (ratio 2.0 per rung) --")
    for trace in sorted({x["trace"] for x in arms}):
        if trace == REF_TRACE:
            continue
        s = [x for x in arms if x["trace"] == trace]
        s.sort(key=lambda x: -x["dt"])
        rs = ["%.2f" % (s[i]["dr"] / s[i + 1]["dr"])
              for i in range(len(s) - 1)
              if s[i + 1]["dr"] not in (0,) and math.isfinite(s[i + 1]["dr"])]
        print("  %-16s dr = %s" % (trace, ", ".join("%.2f" % x["dr"] for x in s)))
        print("  %-16s ratios per halving: %s  (predicted 2.00)"
              % ("", ", ".join(rs) if rs else "n/a"))

    # PREDICTION 2: rho decides the mechanism.
    print("\n-- PREDICTION 2: rho = dr_reconstructedU / dr_cellCentred --")
    cc = {round(x["dt"], 12): x["dr"] for x in arms if x["trace"] == "cellCentred"}
    ru = {round(x["dt"], 12): x["dr"] for x in arms if x["trace"] == "reconstructedU"}
    for k in sorted(set(cc) & set(ru), reverse=True):
        rho = ru[k] / cc[k] if cc[k] else float("nan")
        verdict = ("SOLENOIDALITY (hypothesis confirmed)" if rho >= 0.7 else
                   "SMOOTHING (hypothesis falsified)" if rho <= 0.3 else
                   "MIXED -- claim neither")
        print("  dt=%-12.6g rho=%6.3f  -> %s" % (k, rho, verdict))
    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
