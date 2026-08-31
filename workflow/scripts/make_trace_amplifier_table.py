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
    # case_params.json nests the materialized tokens under "tokens"; fall back
    # to the top level so the script also reads a flat params file.
    meta = json.load(open(cp))
    tok = meta.get("tokens", meta)
    u_ref, t_ref = at(rows, T_REF)
    u_end, t_end = at(rows, T_END)
    # A blown-up arm is a result: report it, do not fit through it.
    diverged = not all(math.isfinite(fnum(r["maxMagU"])) for r in rows)
    t_last = fnum(rows[-1]["TIME"])
    partial = t_last < 0.995*T_END
    G = math.log(u_end / u_ref) if (u_ref > 0 and u_end > 0) else float("nan")
    return dict(
        arm=os.path.basename(d),
        commit=meta.get("gitCommit", "?"),
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
    bad = [x for x in arms if x["trace"] == "?" or not math.isfinite(x["dt"])]
    if bad:
        sys.exit("could not read SL_TRACE_VELOCITY/MAX_DELTA_T for %d arm(s) "
                 "(e.g. %s) -- a table of '?' and nan is not a result"
                 % (len(bad), bad[0]["arm"]))

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
            "volErr", "shapeL2", "minGradPsi", "diverged", "partial", "t_last",
            "commit"]
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
    # rho is only meaningful where the amplifier it normalises by is LARGE.
    # dr_cc -> 0 as dt -> 0 by construction (all three traces converge to the
    # same physics), so at the fine end rho is a ratio of two small numbers and
    # its verdict is noise, not evidence. Score it at the coarsest step -- the
    # one where the instability actually manifests -- and mark the rest.
    ks = sorted(set(cc) & set(ru), reverse=True)
    dr_max = max((abs(cc[k]) for k in ks), default=0.0)
    for k in ks:
        rho = ru[k] / cc[k] if cc[k] else float("nan")
        weak = abs(cc[k]) < 0.25*dr_max
        if weak:
            note = "dr_cc only %.1f%% of its coarsest value -- rho is noise here" \
                   % (100*abs(cc[k])/dr_max if dr_max else float("nan"))
        else:
            note = ("SOLENOIDALITY (hypothesis confirmed)" if rho >= 0.7 else
                    "SMOOTHING (hypothesis falsified)" if rho <= 0.3 else
                    "MIXED -- claim neither")
        print("  dt=%-12.6g rho=%6.3f  dr_cc=%7.2f  -> %s" % (k, rho, cc[k], note))
    if ks:
        k0 = ks[0]
        rho0 = ru[k0] / cc[k0] if cc[k0] else float("nan")
        print("\n  SCORED AT THE COARSEST STEP dt=%g (largest signal):" % k0)
        print("    rho = %.3f -- solenoidality accounts for %.0f%% of the"
              % (rho0, 100*rho0))
        print("    cellCentred amplifier; the other %.0f%% is removed by"
              % (100*(1 - rho0)))
        print("    dropping the extension and/or by the reconstruct smoothing,")
        print("    which THIS matrix does not separate from each other.")
    # ---- LADDER: split the amplifier when all four traces are present -------
    # Every step of the ladder changes exactly one thing, so the three gaps in
    # dr attribute the cellCentred amplifier to three mechanisms:
    #   cellCentred -> reconstructedU        : the reconstruct SMOOTHING
    #   reconstructedU -> reconstructedMomentum : the off-interface EXTENSION
    #   reconstructedMomentum -> projectedFlux  : SOLENOIDALITY
    by = {}
    for x in arms:
        by.setdefault(round(x["dt"], 12), {})[x["trace"]] = x["dr"]
    need = ["cellCentred", "reconstructedU", "reconstructedMomentum"]
    rung = [k for k in sorted(by, reverse=True) if all(t in by[k] for t in need)]
    if rung:
        print("\n-- LADDER: where the cellCentred amplifier actually goes --")
        print("   (shares of dr_cellCentred; each ladder step changes ONE thing)")
        print("   %-12s %9s %11s %11s %11s" %
              ("dt", "dr_cc", "smoothing", "extension", "solenoidal"))
        for k in rung:
            cc_, ru_, rm_ = (by[k]["cellCentred"], by[k]["reconstructedU"],
                             by[k]["reconstructedMomentum"])
            if not cc_:
                continue
            sm, ex, so = (cc_ - ru_)/cc_, (ru_ - rm_)/cc_, rm_/cc_
            mono = cc_ >= ru_ >= rm_ >= 0
            print("   %-12.6g %9.2f %10.1f%% %10.1f%% %10.1f%%%s"
                  % (k, cc_, 100*sm, 100*ex, 100*so,
                     "" if mono else "   NON-MONOTONE -- ladder void at this dt"))
        k0 = rung[0]
        cc_, ru_, rm_ = (by[k0]["cellCentred"], by[k0]["reconstructedU"],
                         by[k0]["reconstructedMomentum"])
        if cc_:
            order = sorted(
                [("the reconstruct smoothing", (cc_ - ru_)/cc_),
                 ("dropping the velocity extension", (ru_ - rm_)/cc_),
                 ("solenoidality of the traced field", rm_/cc_)],
                key=lambda t: -t[1])
            print("\n   AT THE BASE STEP dt=%g, largest share first:" % k0)
            for name, share in order:
                print("     %5.1f%%  %s" % (100*share, name))

    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
