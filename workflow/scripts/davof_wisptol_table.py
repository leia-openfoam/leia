#!/usr/bin/env python3
"""The table of the interface-cell threshold study (Tomislav, 2026-11-08; the configs
config/davof/standing/deformation3D_N*_wisptol*.yaml): one row per run of the studies
studies/davof/standing_deformation3D_N*_wisptol*, the DAVOF swept-face advection in its
recommended setting on the 3D deformation with the reconstruction's wispTol tau varied
(TPF_RECON_WISP_TOL) and the positioning direction unit | areaNormal
(TPF_POSITIONING_DIRECTION).

    python3 workflow/scripts/davof_wisptol_table.py [--studies-dir DIR] [--out DIR] [--N 32]

Writes, under docs/davof/data/wisptol/:
  wisptol_deformation3D.csv         one row per run: the end-of-run errors, the bounds, the
                                    mass, the counts, the wall clock, the stamp
  wisptol_deformation3D_writes.csv  one row per run and write time (t = 0.75, 1.5, 2.25, 3):
                                    the errors and counts at the write
  wisptol_deformation3D_steps.csv   one row per run and step: the reconstruction's
                                    interface-cell count, the rule's counter, the wisps, the
                                    min alpha before and after the bounding, the fragment
                                    log's new seeds (0 where the log is off)
and prints the table and the verdicts W-0 .. W-7b of the configs' pre-registration.

The measures: E_shape, L1(x) = LCentreDiff1, L1(n) = LNormalDiff1, L1(alphaf) from the
converter's CSV (tpf_advection_metrics.py) at the writes; from the log per step: the
first `davof: N interface cells` line after `Time =` (the advection step's
reconstruction), its `K kept m` clause (cells positioned along their own m_k because
the ring sum gave no usable normal) and `tiny directions (...) T` clause (areaNormal),
`interface-cell rule ...: cells with |m_c| above the threshold and alpha outside
(tol, 1 - tol) N, with a plane M`, the curved-clipping line's `wisps without a neighbour
surface W`, `before the bounding: min(alpha) = a, max(alpha) = 1 + b`, `After conservative
bounding: min(alpha) = a, max(alpha) = 1 + b`, `sweptFragmentLog: new face-fraction
interface cells with alpha within 1e-6 of 0 or 1: S (fallback-adjacent a, second-layer
bounded b, wisp without a surface c, bounded face d, none e)`, `ExecutionTime`; the
solve's wall clock from .leia_launch.
"""
import argparse, csv, glob, json, os, re, statistics, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
WRITES = (0.75, 1.5, 2.25, 3.0)


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def tau_of(t):
    return float(t.get("TPF_RECON_WISP_TOL", t.get("TPF_WISP_TOL", "1e-3")))


def read_log(path):
    """Per-step records from the solver log."""
    log = open(path, errors="replace").read()
    chunks = re.split(r"\nTime = ", log)[1:]
    steps = []
    rx = {
        "nInterface": re.compile(r"davof: (\d+) interface cells \(wispTol"),
        "nKept": re.compile(r"(\d+) kept m"),
        "nTiny": re.compile(r"tiny directions \([^)]*\) (\d+)"),
        "ruleOutside": re.compile(r"alpha outside \(tol, 1 - tol\) (\d+), with a plane (\d+)"),
        "wisps": re.compile(r"wisps without a neighbour surface (\d+)"),
        "before": re.compile(r"before the bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)"),
        "after": re.compile(r"After conservative bounding: min\(alpha\) = (\S+), max\(alpha\) = 1 \+ (\S+)"),
        "seeds": re.compile(r"sweptFragmentLog: new face-fraction interface cells with alpha within 1e-6 of 0 or 1: (\d+) \(fallback-adjacent (\d+), second-layer bounded (\d+), wisp without a surface (\d+), bounded face (\d+), none (\d+)\)"),
        "exec": re.compile(r"ExecutionTime = (\S+) s"),
        "stepTime": re.compile(r"timeStepTime = (\S+) s"),
        "posCuts": re.compile(r"cuts per cell (\S+) \(max (\d+)"),
        "posFloor": re.compile(r"floor exits (\d+), maxIter exits (\d+)"),
        "rankFallback": re.compile(r"rank fallbacks (\d+)"),
    }
    for i, ch in enumerate(chunks):
        t = num(ch.split("\n", 1)[0].strip())
        rec = dict(step=i + 1, t=t)
        m = rx["nInterface"].search(ch); rec["nInterface"] = int(m.group(1)) if m else None
        m = rx["nKept"].search(ch); rec["nKept"] = int(m.group(1)) if m else 0
        m = rx["nTiny"].search(ch); rec["nTiny"] = int(m.group(1)) if m else 0
        m = rx["ruleOutside"].search(ch); rec["ruleOutside"] = int(m.group(1)) if m else None; rec["ruleOutsidePlane"] = int(m.group(2)) if m else None
        m = rx["wisps"].search(ch); rec["wisps"] = int(m.group(1)) if m else None
        m = rx["before"].search(ch); rec["minBefore"] = num(m.group(1)) if m else None; rec["excessBefore"] = num(m.group(2)) if m else None
        m = rx["after"].search(ch); rec["minAfter"] = num(m.group(1)) if m else None; rec["excessAfter"] = num(m.group(2)) if m else None
        m = rx["seeds"].search(ch)
        rec["seeds"] = int(m.group(1)) if m else 0
        for k, j in (("seedsFallback", 2), ("seedsBound2", 3), ("seedsWisp", 4), ("seedsBoundFace", 5), ("seedsNone", 6)):
            rec[k] = int(m.group(j)) if m else 0
        m = rx["exec"].search(ch); rec["exec"] = num(m.group(1)) if m else None
        m = rx["stepTime"].search(ch); rec["stepTime"] = num(m.group(1)) if m else None
        m = rx["posCuts"].search(ch); rec["cutsPerCell"] = num(m.group(1)) if m else None
        m = rx["posFloor"].search(ch); rec["floorExits"] = int(m.group(1)) if m else None; rec["maxIterExits"] = int(m.group(2)) if m else None
        m = rx["rankFallback"].search(ch); rec["rankFallbacks"] = int(m.group(1)) if m else None
        steps.append(rec)
    hasLog = "sweptFragmentLog: new" in log or "sweptFragmentLog true" in log
    return steps, ("ABORT at step" in log), hasLog


def read_case(c):
    p = os.path.join(c, "case_params.json")
    if not os.path.isfile(p):
        return None
    cp = json.load(open(p)); t = cp["tokens"]
    rec = dict(case=os.path.relpath(c, REPO), study=os.path.basename(os.path.dirname(c)),
               N=int(float(t.get("N_CELLS", 0))), co=float(t.get("TARGET_ADVECTIVE_CO", 0)),
               tau=tau_of(t), direction=t.get("TPF_POSITIONING_DIRECTION", "unit"),
               setting=t.get("TPF_ENRIGHT_SETTING", "reference"), grazingCut=t.get("TPF_GRAZING_CUT", ""),
               dt=num(t.get("MAX_DELTA_T")), T=num(t.get("END_TIME")), np=int(cp.get("np", 1)))
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) and os.path.getsize(csvp) else []
    last = rows[-1] if rows else {}
    rec.update(status=last.get("status", "noCSV"), abortStep=last.get("abortStep", ""), tLast=num(last.get("TIME")),
               Eshape=num(last.get("Eshape")), L1x=num(last.get("LCentreDiff1")), L1n=num(last.get("LNormalDiff1")),
               L1af=num(last.get("L1_alphaf")), massDevPct=num(last.get("devPercent")),
               nInterfaceFO=num(last.get("nInterfaceCells")), nInterfaceRecFO=num(last.get("nInterface")))
    devs = [abs(x) for x in (num(r.get("devPercent")) for r in rows) if x is not None]
    rec["massDevPctMax"] = max(devs) if devs else None
    # the writes
    writes = []
    for r in rows:
        tt = num(r.get("TIME"))
        if tt is None or num(r.get("Eshape")) is None:
            continue
        writes.append(dict(t=tt, Eshape=num(r.get("Eshape")), L1x=num(r.get("LCentreDiff1")), L1n=num(r.get("LNormalDiff1")),
                           L1af=num(r.get("L1_alphaf")), nInterfaceFO=num(r.get("nInterfaceCells")), nInterfaceRecFO=num(r.get("nInterface")),
                           massDevPct=num(r.get("devPercent"))))
    rec["writes"] = writes
    logp = os.path.join(c, "log.davofAdvectorFoam")
    steps, aborted, hasLog = read_log(logp) if os.path.isfile(logp) else ([], False, False)
    rec["steps"] = steps; rec["nSteps"] = len(steps); rec["fragmentLog"] = hasLog
    if aborted and rec["status"] == "noCSV":
        rec["status"] = "aborted"
    def col(k):
        return [s[k] for s in steps if s.get(k) is not None]
    rec["minBefore"] = min(col("minBefore")) if col("minBefore") else None
    rec["minAfter"] = min(col("minAfter")) if col("minAfter") else None
    rec["excessBefore"] = max(col("excessBefore")) if col("excessBefore") else None
    rec["excessAfter"] = max(col("excessAfter")) if col("excessAfter") else None
    rec["ruleOutsideMax"] = max(col("ruleOutside")) if col("ruleOutside") else None
    rec["ruleOutsideSum"] = sum(col("ruleOutside")) if col("ruleOutside") else None
    rec["wispsMax"] = max(col("wisps")) if col("wisps") else None
    rec["nInterfaceMax"] = max(col("nInterface")) if col("nInterface") else None
    rec["nInterfaceStep1"] = steps[0]["nInterface"] if steps else None
    rec["nKeptMax"] = max(col("nKept")) if col("nKept") else 0
    rec["nTinyMax"] = max(col("nTiny")) if col("nTiny") else 0
    rec["seedsTotal"] = sum(col("seeds")) if hasLog else None
    rec["seedsFallback"] = sum(col("seedsFallback")) if hasLog else None
    rec["seedsBound2"] = sum(col("seedsBound2")) if hasLog else None
    rec["seedsWisp"] = sum(col("seedsWisp")) if hasLog else None
    rec["seedsBoundFace"] = sum(col("seedsBoundFace")) if hasLog else None
    rec["seedsNone"] = sum(col("seedsNone")) if hasLog else None
    rec["firstSeedStep"] = next((s["step"] for s in steps if s["seeds"] > 0), None) if hasLog else None
    rec["rankFallbacksSum"] = sum(col("rankFallbacks")) if col("rankFallbacks") else None
    rec["cutsPerCellMean"] = statistics.mean(col("cutsPerCell")) if col("cutsPerCell") else None
    rec["floorExitsSum"] = sum(col("floorExits")) if col("floorExits") else None
    rec["maxIterExitsSum"] = sum(col("maxIterExits")) if col("maxIterExits") else None
    st = [x for x in col("stepTime") if x > 0]
    rec["stepMedian"] = statistics.median(st) if st else None
    rec["execTime"] = col("exec")[-1] if col("exec") else None
    # the per-step counts at the writes (the step whose time is the write time)
    for w in writes:
        s = next((s for s in steps if s["t"] is not None and abs(s["t"] - w["t"]) <= 1e-6 * max(1.0, w["t"])), None)
        for k in ("nInterface", "ruleOutside", "wisps", "nKept", "nTiny"):
            w[k] = s[k] if s else None
        lo = [x["seeds"] for x in steps if x["t"] is not None and w["t"] - 0.75 - 1e-9 < x["t"] <= w["t"] + 1e-9]
        w["seedsQuarter"] = sum(lo) if hasLog else None
    launch = {}
    lp = os.path.join(c, ".leia_launch")
    if os.path.isfile(lp):
        for l in open(lp):
            if "=" in l:
                k, v = l.strip().split("=", 1); launch[k] = v
    try:
        rec["wall"] = float(launch["end"]) - float(launch["start"])
    except (KeyError, ValueError):
        rec["wall"] = None
    rec["host"] = launch.get("host", "")
    sp = os.path.join(c, ".leia_stamp.json")
    s = json.load(open(sp)) if os.path.isfile(sp) else {}
    bins = s.get("trigger", {}).get("binaries", {})
    rec["stampId"] = s.get("id", "")
    rec["stampBinary"] = (next(iter(bins.values()), {}).get("sha256") or "")[:12] if bins else ""
    rec["stampLeia"] = s.get("record", {}).get("leiaAtRun") or s.get("record", {}).get("leia", "")
    return rec


def at_write(rec, t, key):
    w = next((w for w in rec["writes"] if abs(w["t"] - t) < 1e-6), None)
    return w.get(key) if w else None


def rel(a, b):
    if a is None or b is None:
        return None
    if a == b:
        return 0.0
    d = max(abs(a), abs(b))
    return abs(a - b) / d if d > 0 else 0.0


def f(x, nd=3):
    v = num(x)
    return "-" if v is None else ("%.*e" % (nd, v))


def fi(x):
    return "-" if x is None else str(int(x))


def label(r):
    return ("%g" % r["tau"]) + ("" if r["direction"] == "unit" else " " + r["direction"]) + ("" if r["setting"] == "recommendedLog" else " [" + r["setting"] + "]")


def verdicts(recs, N):
    """The pre-registered verdicts for one rung."""
    out = []
    series = sorted([r for r in recs if r["N"] == N and r["direction"] == "unit" and r["setting"] == "recommendedLog"], key=lambda r: -r["tau"])
    ref = next((r for r in series if r["tau"] == 1e-3), None)
    rep = next((r for r in recs if r["N"] == N and r["direction"] == "unit" and r["setting"] == "recommended" and r["tau"] == 1e-3), None)
    half = (ref["T"] or 3.0) / 2.0 if ref else 1.5
    def say(tag, ok, text):
        out.append("%s %s: %s" % (tag, "HELD" if ok else ("FALSIFIED" if ok is False else "undecided"), text))
    if ref is None:
        out.append("no reference run (tau 1e-3, unit, recommendedLog) at N %d" % N); return out
    # W-0 the replica against the logged run (the result columns at the writes)
    if rep:
        cols = ("Eshape", "L1x", "L1n", "L1af", "nInterfaceFO", "nInterfaceRecFO", "massDevPct")
        same = all(at_write(rep, t, k) == at_write(ref, t, k) for t in WRITES for k in cols) and rep["status"] == ref["status"] \
            and rep["minAfter"] == ref["minAfter"] and rep["nInterfaceMax"] == ref["nInterfaceMax"]
        say("W-0", same, "the replica (recommended, no log) and the logged tau 1e-3 run: result columns at the writes, min alpha and the interface-cell counts %s" % ("identical" if same else "DIFFER"))
    # W-1
    small = [r for r in series if r["tau"] <= 1e-8]
    ok = bool(small) and all(r["status"] == "completed" and r["minAfter"] is not None and r["minAfter"] >= -1e-6 for r in small)
    say("W-1", ok if small else None, "tau <= 1e-8: " + ", ".join("%g: %s, min alpha after %s" % (r["tau"], r["status"], f(r["minAfter"], 1)) for r in small))
    # W-2
    nh = [(r["tau"], at_write(r, half, "nInterface")) for r in series]
    vals = [v for _, v in nh]
    mono = all(v is not None for v in vals) and all(b >= a for a, b in zip(vals, vals[1:]))
    n3 = at_write(ref, half, "nInterface"); n10 = next((v for tt, v in nh if tt == 1e-10), None)
    rise = (n10 - n3) / n3 if (n3 and n10 is not None) else None
    say("W-2", (mono and rise is not None and rise < 0.2) if all(v is not None for v in vals) and rise is not None else None,
        "interface cells at T/2 over tau: " + ", ".join("%g: %s" % (tt, fi(v)) for tt, v in nh) + "; monotone %s; rise 1e-3 -> 1e-10 %s" % (mono, ("%.1f %%" % (100 * rise)) if rise is not None else "-"))
    # W-3
    es = [(r["tau"], r["Eshape"]) for r in series if r["status"] == "completed" and r["Eshape"]]
    ch = (max(e for _, e in es) / min(e for _, e in es) - 1) if len(es) >= 2 else None
    say("W-3", (ch < 0.05) if ch is not None else None, "E_shape at T: " + ", ".join("%g: %s" % (tt, f(e)) for tt, e in es) + "; max/min - 1 = %s" % (("%.2f %%" % (100 * ch)) if ch is not None else "-"))
    # W-4
    r0 = ref["ruleOutsideMax"]
    grow = [r for r in series if r["tau"] < 1e-3 and r["ruleOutsideMax"] is not None and r0 is not None and r["ruleOutsideMax"] > 1.1 * r0 + 1]
    say("W-4", (not grow) if r0 is not None else None, "the rule's counter, max over the steps: " + ", ".join("%g: %s (sum %s)" % (r["tau"], fi(r["ruleOutsideMax"]), fi(r["ruleOutsideSum"])) for r in series))
    # W-5
    w0 = ref["wall"]
    ws = [(r["tau"], r["wall"]) for r in series]
    bad = [tt for tt, w in ws if w is None or w0 is None or not (0.9 * w0 <= w <= 1.1 * w0)]
    say("W-5", (not bad) if w0 else None, "wall clock s: " + ", ".join("%g: %s" % (tt, ("%.0f" % w) if w else "-") for tt, w in ws) + ("; outside 10 %%: %s" % bad if bad else ""))
    # W-7
    lim = [r for r in series if r["tau"] <= 1e-10]
    if len(lim) >= 2:
        okc = all(r["status"] == "completed" and r["minAfter"] is not None and r["minAfter"] >= -1e-6 for r in lim)
        nI = [at_write(r, half, "nInterface") for r in lim]; eI = [r["Eshape"] for r in lim]
        dn = max(rel(a, b) or 0 for a in nI for b in nI) if all(v is not None for v in nI) else None
        de = max(rel(a, b) or 0 for a in eI for b in eI) if all(v is not None for v in eI) else None
        held = okc and dn is not None and de is not None and dn < 0.01 and de < 0.01
        fals = (not okc) or (dn is not None and dn > 0.05) or (de is not None and de > 0.05)
        say("W-7", True if held else (False if fals else None), "tau <= 1e-10 (%s): interface cells at T/2 within %s, E_shape within %s, all completed bounded %s" % (", ".join("%g" % r["tau"] for r in lim), ("%.2f %%" % (100 * dn)) if dn is not None else "-", ("%.2f %%" % (100 * de)) if de is not None else "-", okc))
    # W-6 and W-7b: the areaNormal runs against the unit runs at the same tau
    for tag, tt in (("W-6", 1e-3), ("W-7b", 0.0)):
        a = next((r for r in recs if r["N"] == N and r["direction"] == "areaNormal" and r["tau"] == tt), None)
        u = next((r for r in series if r["tau"] == tt), None)
        if not (a and u):
            continue
        ds = []
        for t in WRITES:
            for k in ("Eshape", "L1x", "L1n", "L1af"):
                d = rel(at_write(a, t, k), at_write(u, t, k))
                if d is not None:
                    ds.append(d)
        dn = max((rel(at_write(a, t, "nInterface"), at_write(u, t, "nInterface")) or 0) for t in WRITES)
        dmax = max(ds) if ds else None
        held = a["status"] == u["status"] and dmax is not None and dmax <= 1e-4 and dn <= 0.005 and (tag == "W-6" or a["nTinyMax"] == 0)
        fals = a["status"] != u["status"] or (dmax is not None and dmax > 1e-2) or (tag == "W-7b" and a["nTinyMax"] > 0)
        say(tag, True if held else (False if fals else None), "areaNormal vs unit at tau %g: max relative difference of E_shape, L1(x), L1(n), L1(alphaf) over the writes %s, interface cells %s, status %s / %s, tiny directions max %d" % (tt, ("%.1e" % dmax) if dmax is not None else "-", "%.2e" % dn, a["status"], u["status"], a["nTinyMax"]))
    return out


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies-dir", default=os.path.join(REPO, "studies", "davof"))
    ap.add_argument("--out", default=os.path.join(REPO, "docs", "davof", "data", "wisptol"))
    ap.add_argument("--extra-case", action="append", default=[], help="a case directory to add (the matrix's reference row)")
    a = ap.parse_args(argv)
    recs = []
    for st in sorted(glob.glob(os.path.join(a.studies_dir, "standing_deformation3D_N*_wisptol*"))):
        for c in sorted(glob.glob(os.path.join(st, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = read_case(c)
            if r:
                recs.append(r)
    for c in a.extra_case:
        r = read_case(c)
        if r:
            r["study"] = "extra:" + os.path.basename(os.path.dirname(c)); recs.append(r)
    recs.sort(key=lambda r: (r["N"], r["direction"] != "unit", r["setting"] != "recommendedLog", -r["tau"]))
    os.makedirs(a.out, exist_ok=True)
    base = os.path.join(a.out, "wisptol_deformation3D")
    cols = ["study", "N", "np", "co", "tau", "direction", "setting", "grazingCut", "status", "abortStep", "nSteps", "tLast",
            "Eshape", "L1x", "L1n", "L1af", "massDevPct", "massDevPctMax", "minBefore", "minAfter", "excessBefore", "excessAfter",
            "nInterfaceStep1", "nInterfaceMax", "nInterfaceFO", "nInterfaceRecFO", "ruleOutsideMax", "ruleOutsideSum", "wispsMax",
            "nKeptMax", "nTinyMax", "fragmentLog", "seedsTotal", "seedsFallback", "seedsBound2", "seedsWisp", "seedsBoundFace", "seedsNone",
            "firstSeedStep", "rankFallbacksSum", "cutsPerCellMean", "floorExitsSum", "maxIterExitsSum", "wall", "execTime", "stepMedian",
            "host", "stampId", "stampBinary", "stampLeia", "case"]
    with open(base + ".csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(recs)
    wcols = ["study", "N", "tau", "direction", "setting", "t", "Eshape", "L1x", "L1n", "L1af", "massDevPct", "nInterface", "nInterfaceFO",
             "nInterfaceRecFO", "ruleOutside", "wisps", "nKept", "nTiny", "seedsQuarter"]
    with open(base + "_writes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wcols, extrasaction="ignore"); w.writeheader()
        for r in recs:
            for wr in r["writes"]:
                w.writerow(dict(r, **wr))
    scols = ["study", "N", "tau", "direction", "setting", "step", "t", "nInterface", "nKept", "nTiny", "ruleOutside", "ruleOutsidePlane", "wisps",
             "minBefore", "excessBefore", "minAfter", "excessAfter", "seeds", "seedsFallback", "seedsBound2", "seedsWisp", "seedsBoundFace", "seedsNone",
             "rankFallbacks", "cutsPerCell", "floorExits", "maxIterExits", "stepTime"]
    with open(base + "_steps.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols, extrasaction="ignore"); w.writeheader()
        for r in recs:
            for s in r["steps"]:
                w.writerow(dict(r, **s))
    # the table
    print("| N (np) | tau | status | steps | E_shape(T) | L1(x) | L1(n) | L1(alphaf) | mass % | min a before | min a after | excess | n_if step 1 | n_if T/2 | n_if T | n_if max | rule max (sum) | wisps max | kept m | tiny | seeds (first) | wall s | s/step |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in recs:
        half = (r["T"] or 3.0) / 2.0
        print("| %d (%d) | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s (%s) | %s | %s | %s | %s (%s) | %s | %s |" % (
            r["N"], r["np"], label(r), r["status"] + ((" " + str(r["abortStep"])) if r["status"] == "aborted" else ""), fi(r["nSteps"]),
            f(r["Eshape"]), f(r["L1x"]), f(r["L1n"]), f(r["L1af"]), f(r["massDevPct"], 1), f(r["minBefore"], 1), f(r["minAfter"], 1), f(r["excessAfter"], 1),
            fi(r["nInterfaceStep1"]), fi(at_write(r, half, "nInterface")), fi(at_write(r, r["T"] or 3.0, "nInterface")), fi(r["nInterfaceMax"]),
            fi(r["ruleOutsideMax"]), fi(r["ruleOutsideSum"]), fi(r["wispsMax"]), fi(r["nKeptMax"]), fi(r["nTinyMax"]),
            fi(r["seedsTotal"]), fi(r["firstSeedStep"]), ("%.0f" % r["wall"]) if r["wall"] else "-", f(r["stepMedian"], 2)))
    print()
    for N in sorted({r["N"] for r in recs}):
        print("N %d:" % N)
        for v in verdicts(recs, N):
            print("  " + v)
    print("written:", base + ".csv", base + "_writes.csv", base + "_steps.csv", "(%d runs)" % len(recs))


if __name__ == "__main__":
    main(sys.argv[1:])
