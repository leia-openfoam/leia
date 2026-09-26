#!/usr/bin/env python3
"""Summary, orders and verdict of one method-gate candidate (CLAUDE.md "Method gates").

Reads the arm studies listed in <summaryDir>/<candidate>/manifest.json (written by
render_gate_configs.py), extracts the error vector of every case at the instants the gate
defines, and writes into the same directory:

  summary.csv     one row per (arm, rung): the whole vector, rank count, wall clock, state
  orders.csv      one row per (arm, metric): pairwise and least-squares orders; for a
                  quantity without an exact value the Celik (2008) Richardson estimate
  seam.csv        the coarsest shear rung at other decompositions against np = gate np
  vsBaseline.csv  (with --baseline) candidate against baseline at every rung, and orders
  verdict.txt     (with --baseline) the pre-registered criteria, PASS or FAIL each

Instants: reversed kinematic flow -- gradient at T/2, shape at T, volume at T/2 and T
(CLAUDE.md research loop, step 5); droplets -- at T, volume also at T/2. No L_inf metric is
in the vector or gets an order. A case whose run did not reach END_TIME has blank final
metrics: a mid-run row is never published as final.

Usage:
  make_gate_summary.py --gate config/gates/methodGate2D.yaml --summary-dir <dir>/<cand>
      [--baseline <dir>/baseline] [--exact1d-only] [--docs]
Exit code: 0; with --exact1d-only, 1 if the exact1D arm fails; with --baseline, 1 if the
verdict is FAIL.
"""
import argparse
import csv
import glob
import json
import math
import os
import shutil
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import aggregate    # noqa: E402  (_launch_record, _reached, _last_step)
import paths        # noqa: E402
import richardson   # noqa: E402

CLASSIFIER = os.path.join(HERE, "foam_log_state.sh")
CMP = os.path.join(HERE, "compare_metrics_csv.py")
SKIP_CLOCK = "ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME"
# Decomposition invariance: max over rows of |a - b| divided by the column's largest |value|.
SEAM_TOL = 1e-10
# Every vector metric is an error: lower is better, the exact value is 0.
ERROR_METRICS = ["shapeError", "gradientBandError", "volumeError", "volumeErrorHalf",
                 "boundsError", "rhoClipFraction", "maxMagU", "l2MagUPrime",
                 "pressureJumpError", "curvatureError", "travelledFractionError", "qError"]
# Quantities without an exact value: Celik/Richardson, reported, not scored.
REFERENCE_FREE = ["oscPeriod", "oscDampingRate"]
KINEMATIC_CSV = {"leiaSemiLagrangeLevelSetFoam": "leiaSemiLagrangeLevelSetFoam.csv",
                 "leiaLevelSetFoam": "leiaLevelSetFoam.csv"}


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def read_rows(path):
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        return [{(k or "").strip(): num(v) for k, v in r.items()} for r in csv.DictReader(fh)]


def at(rows, t, dt=None):
    """Row nearest time t, if the run reached it; {} otherwise."""
    if not rows:
        return {}
    t_last = rows[-1].get("TIME")
    if dt is None and len(rows) >= 2:
        dt = rows[-1]["TIME"] - rows[-2]["TIME"]
    if not aggregate._reached(t_last, t, dt):
        return {}
    return min(rows, key=lambda r: abs((r.get("TIME") or 0.0) - t))


def colmax(rows, col):
    vals = [r.get(col) for r in rows if r.get(col) is not None]
    return max(vals) if vals else None


def classify(case_dir, solver):
    log = os.path.join(case_dir, f"log.{solver}")
    r = subprocess.run(["bash", CLASSIFIER, log], capture_output=True, text=True)
    parts = r.stdout.split()
    out = {"state": parts[0] if parts else "UNKNOWN"}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k] = v
    return out


def osc_period_damping(rows, col="m2CosCoefficient"):
    """Period from the zero crossings of the signed mode-2 coefficient, damping rate from
    a least-squares fit of ln|peak| between successive crossings."""
    ts = [(r["TIME"], r.get(col)) for r in rows if r.get("TIME") is not None and r.get(col) is not None]
    cross = []
    for (t0, a0), (t1, a1) in zip(ts, ts[1:]):
        if a0 == 0 or a0 * a1 < 0:
            cross.append(t0 - a0 * (t1 - t0) / (a1 - a0) if a1 != a0 else t0)
    period = (2.0 * (cross[-1] - cross[0]) / (len(cross) - 1)) if len(cross) >= 3 else None
    peaks = []
    for c0, c1 in zip(cross, cross[1:]):
        seg = [(t, abs(a)) for t, a in ts if c0 <= t <= c1]
        if seg:
            peaks.append(max(seg, key=lambda x: x[1]))
    damping = None
    pts = [(t, math.log(a)) for t, a in peaks if a > 0]
    if len(pts) >= 3:
        n = len(pts)
        mt = sum(p[0] for p in pts) / n
        ml = sum(p[1] for p in pts) / n
        stt = sum((p[0] - mt) ** 2 for p in pts)
        damping = -sum((p[0] - mt) * (p[1] - ml) for p in pts) / stt if stt else None
    return period, damping, len(cross)


def lamb_period(dims, R, sigma, rho_d, rho_a, n=2):
    if dims == 2:
        w2 = (n ** 3 - n) * sigma / ((rho_d + rho_a) * R ** 3)
    else:
        w2 = n * (n - 1) * (n + 2) * sigma / ((n * rho_d + (n + 1) * rho_a) * R ** 3)
    return 2.0 * math.pi / math.sqrt(w2)


def case_vector(case_dir, arm, gate, cand, solver):
    meta = json.load(open(os.path.join(case_dir, "case_params.json")))
    tok = meta.get("tokens", {})
    dims = 1 if arm["case"].startswith("1D") else int(gate["dims"])
    N = int(float(tok.get("N_CELLS")))
    L = float(arm.get("DOMAIN_LENGTH", tok.get("DOMAIN_LENGTH", 1.0)))
    h = L / N
    T = float(tok.get("END_TIME"))
    row = {"arm": arm["_name"], "case": arm["case"], "N": N, "nCells": N ** dims, "h": h,
           "R_over_h": (float(arm["R"]) / h) if arm.get("R") else "",
           "caseDir": os.path.relpath(case_dir, REPO), "gitCommit": meta.get("gitCommit", "")}
    st = classify(case_dir, solver)
    launch = aggregate._launch_record(case_dir)
    row.update({"state": st.get("state"), "steps": st.get("steps", ""),
                "nRanks": st.get("nprocs", launch.get("nRanks", "")),
                "wallClockSolve": launch.get("wallClockSolve", "")})
    stamp = os.path.join(case_dir, "leia.version")
    row["libStamps"] = ";".join(l.strip() for l in open(stamp) if l.strip()) if os.path.isfile(stamp) else ""

    if arm["kind"] == "kinematic":
        adv = read_rows(os.path.join(case_dir, KINEMATIC_CSV.get(solver, f"{solver}.csv")))
        gp = read_rows(os.path.join(case_dir, "gradPsiError.csv"))
    else:
        adv = read_rows(os.path.join(case_dir, "leiaLevelSetFoam.csv"))
        drop = read_rows(os.path.join(case_dir, f"{solver}.csv"))
    fin_adv, half_adv = at(adv, T), at(adv, 0.5 * T)
    row["endTimeReached"] = "1" if fin_adv else "0"
    row["solverClockTime"] = fin_adv.get("ELAPSED_CLOCK_TIME", "") if fin_adv else ""
    row["boundsError"] = colmax(adv, "E_BOUND_ALPHA")
    try:
        wc, nr, stp = float(row["wallClockSolve"]), float(row["nRanks"]), float(row["steps"])
        row["coreSeconds"] = wc * nr
        row["secondsPerStep"] = wc / stp if stp else ""
    except (TypeError, ValueError):
        row["coreSeconds"], row["secondsPerStep"] = "", ""

    if arm["kind"] == "kinematic":
        fin_gp, half_gp = at(gp, T), at(gp, 0.5 * T)
        row["staticGradientFloor"] = gp[0].get("E_NARROW_L2_GRAD_PSI") if gp else None
        if arm["_name"] == "exact1D":
            # u = alpha x has div u = alpha, and q = exp(-alpha t) is the exact solution: the
            # volume change and |q - 1| are physics here, not errors. Only the closed form is.
            row["staticGradientFloor"] = None
            q = fin_gp.get("NARROW_MEAN_MAG_GRAD_PSI")
            row["qBandMean"] = q
            plain = not any(cand["tokens"].get(k, v) != v for k, v in
                            (("VELOCITY_EXTENSION", "none"), ("SL_SOURCE", "none"),
                             ("SDPLS_SOURCE", "noSource")))
            alpha = float(tok.get("STRAIN_RATE", 1.0))
            q_exact = math.exp(-alpha * T) if plain else None
            row["qExact"] = q_exact if q_exact is not None else ""
            row["qError"] = abs(q - q_exact) / q_exact if (q is not None and q_exact) else None
        else:
            row["volumeError"] = abs(fin_adv["E_VOL_ALPHA_REL"]) if fin_adv.get("E_VOL_ALPHA_REL") is not None else None
            row["volumeErrorHalf"] = abs(half_adv["E_VOL_ALPHA_REL"]) if half_adv.get("E_VOL_ALPHA_REL") is not None else None
            row["gradientBandError"] = (half_gp if arm.get("reversed") else fin_gp).get("E_NARROW_L2_GRAD_PSI")
            row["shapeError"] = fin_adv.get("E_GEOM_ALPHA_REL")
        return row

    # two-phase droplet arms
    fin, half = at(drop, T), at(drop, 0.5 * T)
    R = float(arm["R"])
    row["staticGradientFloor"] = drop[0].get("gradPsiL2ErrorBand") if drop else None
    row["gradientBandError"] = fin.get("gradPsiL2ErrorBand")
    row["volumeError"] = abs(fin["phaseVolumeRelError"]) if fin.get("phaseVolumeRelError") is not None else None
    row["volumeErrorHalf"] = abs(half["phaseVolumeRelError"]) if half.get("phaseVolumeRelError") is not None else None
    row["rhoClipFraction"] = colmax(drop, "rhoClipFraction")
    row["l2MagUPrime"] = fin.get("l2MagUPrime")
    flow = arm.get("flow")
    if flow in ("stationary", "translating"):
        rad = fin.get("zeroSetRadialL2")
        row["shapeError"] = rad / R if rad is not None else None
        sigma = float(tok.get("SIGMA", 0.07274))
        jump = sigma * (int(gate["dims"]) - 1) / R
        pl = fin.get("pLaplace")
        row["pressureJumpError"] = abs(pl - jump) / jump if pl is not None else None
        kerr = fin.get("kErrL2Band")          # absolute [1/m], against the exact (dims-1)/R
        row["curvatureError"] = kerr * R / (int(gate["dims"]) - 1) if kerr is not None else None
        row["maxMagU"] = colmax(drop, "maxMagU" if flow == "stationary" else "maxMagUPrime")
    if flow == "translating" and drop and fin:
        x0, x1 = drop[0].get("centroidX"), fin.get("centroidX")
        U = float(arm["U"])
        if x0 is not None and x1 is not None:
            row["travelledFractionError"] = abs((x1 - x0) / (U * (fin["TIME"] - drop[0]["TIME"])) - 1.0)
    if flow == "oscillating":
        period, damping, ncross = osc_period_damping(drop if fin else [])
        row["oscPeriod"], row["oscDampingRate"], row["oscZeroCrossings"] = period, damping, ncross
        tl = lamb_period(int(gate["dims"]), float(arm["R"]), float(tok.get("SIGMA", 0.07274)),
                         float(tok.get("RHO1", 998.2)), float(tok.get("RHO2", 1.19)))
        row["lambPeriod"] = tl
        row["periodVsLamb"] = (period / tl - 1.0) if period else None
    return row


def arm_cases(studies_dir, study):
    return sorted(d for d in glob.glob(os.path.join(studies_dir, study, "*_0*")) if os.path.isdir(d))


def write_csv(path, rows, first=()):
    cols = list(first)
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def orders_of(summary):
    out = []
    for arm in sorted({r["arm"] for r in summary}):
        rows = sorted([r for r in summary if r["arm"] == arm], key=lambda r: -r["h"])
        h = [r["h"] for r in rows]
        for m in ERROR_METRICS + REFERENCE_FREE + ["staticGradientFloor"]:
            vals = [num(r.get(m)) for r in rows]
            if all(v is None for v in vals):
                continue
            rec = {"arm": arm, "metric": m, "N": " ".join(str(r["N"]) for r in rows),
                   "values": " ".join("" if v is None else f"{v:.6g}" for v in vals)}
            if m in REFERENCE_FREE:
                c = richardson.celik(h, vals) if all(v is not None for v in vals) and len(vals) == 3 else {}
                rec.update({"celikOrder": c.get("p"), "extrapolated": c.get("f_ext"),
                            "gciFine": c.get("gci_fine"), "asymptoticRatio": c.get("asymptotic_ratio"),
                            "convergence": c.get("conv_type", "insufficient")})
            else:
                po = richardson.pairwise_orders(h, [v if v is not None else None for v in vals])
                rec["pairwiseOrders"] = " ".join("" if p is None else f"{p:.3f}" for _, _, p in po)
                rec["lsqOrder"] = richardson.lsq_order(h, vals)
            out.append(rec)
    return out


def seam_rows(man, studies_dir):
    out = []
    ref_arm = next((a for a in man["arms"] if a["arm"] == "shear"), None)
    if not ref_arm:
        return out
    ref_case = (arm_cases(studies_dir, ref_arm["study"]) or [None])[0]
    for a in man["arms"]:
        if not a.get("seamOf") or not ref_case:
            continue
        case = (arm_cases(studies_dir, a["study"]) or [None])[0]
        if not case:
            out.append({"decomposition": a["arm"], "np": a["np"], "verdict": "MISSING"})
            continue
        for name in ("gradPsiError.csv", "leiaSemiLagrangeLevelSetFoam.csv", "leiaLevelSetFoam.csv"):
            pa, pb = os.path.join(ref_case, name), os.path.join(case, name)
            if not (os.path.isfile(pa) and os.path.isfile(pb)):
                continue
            A, B = read_rows(pa), read_rows(pb)
            if len(A) != len(B):
                out.append({"decomposition": a["arm"], "np": a["np"], "csv": name, "verdict": "FAIL",
                            "detail": f"{len(A)} rows vs {len(B)} rows"})
                continue
            worst, wcol = 0.0, ""
            for col in A[0]:
                if col in ("ELAPSED_CPU_TIME", "ELAPSED_CLOCK_TIME") or col == "":
                    continue
                va = [r.get(col) for r in A]; vb = [r.get(col) for r in B]
                if any(x is None for x in va + vb):
                    continue
                scale = max(max(abs(x) for x in va), max(abs(x) for x in vb), 1e-300)
                d = max(abs(x - y) for x, y in zip(va, vb)) / scale
                if d > worst:
                    worst, wcol = d, col
            out.append({"decomposition": a["arm"], "np": a["np"], "csv": name,
                        "verdict": "PASS" if worst <= SEAM_TOL else "FAIL",
                        "maxColumnScaledDiff": worst, "worstColumn": wcol})
    return out


def identical_to(man, base_man, studies_dir):
    """True if every per-case CSV of the candidate equals the baseline's at tolerance 0."""
    pairs = 0
    for a in man["arms"]:
        b = next((x for x in base_man["arms"] if x["arm"] == a["arm"]), None)
        if not b:
            continue
        for ca, cb in zip(arm_cases(studies_dir, a["study"]), arm_cases(studies_dir, b["study"])):
            for p in glob.glob(os.path.join(ca, "*.csv")):
                q = os.path.join(cb, os.path.basename(p))
                if not os.path.isfile(q):
                    return False
                r = subprocess.run([sys.executable, CMP, p, q, "--tol", "0", "--skip", SKIP_CLOCK],
                                   capture_output=True, text=True)
                if r.returncode != 0:
                    return False
                pairs += 1
    return pairs > 0


def summarize(gate, man, studies_dir, cand):
    rows = []
    for a in man["arms"]:
        if a.get("seamOf"):
            continue
        arm = dict(gate["arms"][a["arm"]]); arm["_name"] = a["arm"]
        solver = gate["solvers"][man["line"]][arm["kind"]]
        for d in arm_cases(studies_dir, a["study"]):
            rows.append(case_vector(d, arm, gate, cand, solver))
    return rows


def verdict(gate, cand, summ, base_summ, orders, base_orders, no_effect, seam):
    tol = float(gate["verdict"]["regressionTolerance"])
    otol = float(gate["verdict"]["orderTolerance"])
    lines, ok, cmp_rows = [], True, []
    key = lambda r: (r["arm"], int(r["N"]))
    base = {key(r): r for r in base_summ}
    for r in summ:
        b = base.get(key(r))
        if not b:
            continue
        if b["state"] == "COMPLETED" and r["state"] != "COMPLETED":
            lines.append(f"FAIL completion: {r['arm']} N={r['N']} is {r['state']} where the baseline COMPLETED")
            ok = False
        for m in ERROR_METRICS:
            cv, bv = num(r.get(m)), num(b.get(m))
            if cv is None or bv is None:
                continue
            ratio = (cv / bv) if bv else (0.0 if cv == 0 else float("inf"))
            cmp_rows.append({"arm": r["arm"], "N": r["N"], "metric": m, "candidate": cv,
                             "baseline": bv, "ratio": ratio})
    finest = {}
    for r in summ:
        if r["arm"] not in finest or r["h"] < finest[r["arm"]]["h"]:
            finest[r["arm"]] = r
    for arm, r in finest.items():
        b = base.get(key(r))
        if not b:
            continue
        for m in ERROR_METRICS:
            cv, bv = num(r.get(m)), num(b.get(m))
            if cv is None or bv is None:
                continue
            bad = (cv > (1 + tol) * bv) if bv > 0 else (cv > 1e-12)
            if bad:
                lines.append(f"FAIL regression: {arm} {m} at N={r['N']}: {cv:.4g} vs baseline {bv:.4g}")
                ok = False
    bo = {(o["arm"], o["metric"]): o for o in base_orders}
    for o in orders:
        b = bo.get((o["arm"], o["metric"]))
        if not b or o.get("lsqOrder") is None or b.get("lsqOrder") is None:
            continue
        if o["lsqOrder"] < b["lsqOrder"] - otol:
            lines.append(f"FAIL order: {o['arm']} {o['metric']} order {o['lsqOrder']:.2f} vs baseline {b['lsqOrder']:.2f}")
            ok = False
    if no_effect and cand["candidate"] not in ("baseline", "baselineEulerian"):
        lines.append("FAIL no effect: every CSV equals the baseline's; the method tokens were not consumed")
        ok = False
    tgt = cand.get("target") or {}
    if tgt:
        for arm in tgt.get("arms", list(finest)):
            r = finest.get(arm)
            b = base.get(key(r)) if r else None
            cv, bv = (num(r.get(tgt["metric"])) if r else None), (num(b.get(tgt["metric"])) if b else None)
            if cv is None or bv is None or bv == 0:
                lines.append(f"FAIL target: {arm} {tgt['metric']} not available")
                ok = False
            elif cv / bv > float(tgt["ratioToBaselineAtFinest"]):
                lines.append(f"FAIL target: {arm} {tgt['metric']} ratio {cv / bv:.3f} > {tgt['ratioToBaselineAtFinest']}")
                ok = False
            else:
                lines.append(f"PASS target: {arm} {tgt['metric']} ratio {cv / bv:.3f}")
    for s in seam:
        if s.get("verdict") != "PASS":
            lines.append(f"FAIL seam: {s.get('decomposition')} {s.get('csv', '')} "
                         f"{s.get('worstColumn', '')} {s.get('maxColumnScaledDiff', s.get('detail', ''))}")
            ok = False
    if not cand.get("preRegistered", True):
        lines.append("NOTE: ad-hoc candidate (SET=), NOT pre-registered; exploratory only")
    lines.insert(0, f"VERDICT {'PASS' if ok else 'FAIL'}: candidate {cand['candidate']} against baseline")
    return ok, lines, cmp_rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gate", required=True)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--baseline", default=None, help="summary dir of the matching baseline")
    ap.add_argument("--studies-dir", default=os.path.join(REPO, "studies"))
    ap.add_argument("--exact1d-only", action="store_true")
    ap.add_argument("--docs", action="store_true", help="copy the tables into the theme")
    a = ap.parse_args(argv)
    gate = yaml.safe_load(open(a.gate))
    man = json.load(open(os.path.join(a.summary_dir, "manifest.json")))
    cand = json.load(open(os.path.join(a.summary_dir, "candidate.json")))

    if a.exact1d_only:
        ex = next(x for x in man["arms"] if x["arm"] == "exact1D")
        arm = dict(gate["arms"]["exact1D"]); arm["_name"] = "exact1D"
        solver = gate["solvers"][man["line"]][arm["kind"]]
        rows = [case_vector(d, arm, gate, cand, solver) for d in arm_cases(a.studies_dir, ex["study"])]
        bad = [r for r in rows if r["state"] != "COMPLETED"]
        finest = min(rows, key=lambda r: r["h"]) if rows else None
        qerr = num(finest.get("qError")) if finest else None
        ok = bool(rows) and not bad and (qerr is None or qerr <= 1e-3)
        write_csv(os.path.join(a.summary_dir, "exact1D.csv"), rows, first=("arm", "N"))
        print(f"[gate] exact1D {cand['candidate']}: {len(rows)} cases, not completed {len(bad)}, "
              f"qError(finest) = {qerr}" + ("" if qerr is not None else " (oracle pending: source or extension active)")
              + f" -> {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 1

    summ = summarize(gate, man, a.studies_dir, cand)
    write_csv(os.path.join(a.summary_dir, "summary.csv"), summ, first=("arm", "N", "nCells", "h", "R_over_h", "state"))
    orders = orders_of(summ)
    write_csv(os.path.join(a.summary_dir, "orders.csv"), orders, first=("arm", "metric"))
    seam = seam_rows(man, a.studies_dir)
    write_csv(os.path.join(a.summary_dir, "seam.csv"), seam or [{"decomposition": "none"}])
    outputs = ["summary.csv", "orders.csv", "seam.csv"]
    rc = 0
    if a.baseline:
        bman = json.load(open(os.path.join(a.baseline, "manifest.json")))
        bcand = json.load(open(os.path.join(a.baseline, "candidate.json")))
        bsumm = summarize(gate, bman, a.studies_dir, bcand)
        border = orders_of(bsumm)
        no_eff = identical_to(man, bman, a.studies_dir)
        ok, lines, cmp_rows = verdict(gate, cand, summ, bsumm, orders, border, no_eff, seam)
        write_csv(os.path.join(a.summary_dir, "vsBaseline.csv"), cmp_rows or [{"arm": "none"}])
        with open(os.path.join(a.summary_dir, "verdict.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print("\n".join(lines))
        outputs += ["vsBaseline.csv", "verdict.txt"]
        rc = 0 if ok else 1
    print(f"[gate] {man['gate']} {cand['candidate']}: {len(summ)} cases -> {a.summary_dir}")
    if a.docs and not cand.get("smoke"):
        tdir = paths.tables_dir(gate["theme"])
        for f in outputs:
            if f.endswith(".csv"):
                dst = os.path.join(tdir, f"{man['gate']}_{cand['candidate']}_{f}")
                shutil.copy(os.path.join(a.summary_dir, f), dst)
                print(f"[gate] -> {os.path.relpath(dst, REPO)}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
