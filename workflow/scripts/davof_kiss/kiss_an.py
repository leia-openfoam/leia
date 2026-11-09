#!/usr/bin/env python3
"""The tables and verdicts of the KISS variant (sweptClip onePlane; Tomislav, 2026-11-09;
docs/davof/REPORT-swept.md "The KISS variant: one plane per face").

    python3 workflow/scripts/davof_kiss/kiss_an.py plane     [--studies-dir D] [--out D]
    python3 workflow/scripts/davof_kiss/kiss_an.py identity  ...
    python3 workflow/scripts/davof_kiss/kiss_an.py rotation  ...
    python3 workflow/scripts/davof_kiss/kiss_an.py standing  ...   (after make_standing_tables.py)

Reads every case of the kiss studies under studies/davof (case_params.json tokens, the
converter's advectorVoF.csv last row, the solver log's per-step lines: the bounds after the
conservative bounding, the Courant number, the one-plane clause of the clipping line and the
rule 2 line, the abort guard; .leia_launch for the wall clock and the ranks) and writes,
under <out> (default studies/davof/kiss_tables/), one CSV per mode and prints the tables with
the verdicts of the pre-registered predictions:

  plane      kiss_plane_np1 / np8: per (setting, arm, N, Co) L1(alphaf), L1(n), E_shape, L1(x),
             min/max alpha, mass deviation, wall clock, steps, the donor classes; P1 per onePlane
             cell (L1(alphaf) <= 1e-12 and L1(n) <= 1e-12); the np1-np8 identity per cell
  identity   kiss_identity_rotation_np1 / np8: the same numbers paired, plus the per-step bounds
  rotation   kiss_rotation_N* and kiss_rotationRDF_N*: per (row, N, Co) with the orders over N;
             P2: E_shape(onePlane)/E_shape(recommended) in [0.9, 1.1], min alpha >= plicRDF's
  standing   the standing_<family>.csv of make_standing_tables.py (the kiss row beside the
             record's rows): P3 (completion, min alpha vs plicRDF, E_shape <= 1.3 x the
             recommended / signedTets row, mass), P4 (step time <= 0.7 x signedTets), P5 (the
             seeds of deformation3D N 64 at T/2 and over the run, from the kiss logs)
"""
import argparse, csv, glob, json, math, os, re, statistics, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
IDENT_TOL = 1e-12
P1_TOL = 1e-12

ONEPLANE_RE = re.compile(
    r"; one plane: faces (\d+), clipped (\d+), clipped triangles (\d+); donors: upwind (\d+), downwind (\d+), "
    r"ring (\d+), remote (\d+), beyond \(upwind taken\) (\d+), outside \(boundary value\) (\d+); "
    r"without a plane: full (\d+), empty (\d+), wisps (\d+), wisp planes (\d+); zero-flux faces (\d+), clamped (\d+)"
    r"(?:; donor content: donors scaled liquid (\d+), gas (\d+), min ratio (\S+), faces with a remote donor not bounded (\d+))?")
ONEPLANE_KEYS = ("faces", "clipped", "clipTris", "donorUp", "donorDown", "donorRing", "donorRemote", "donorBeyond",
                 "donorOutside", "full", "empty", "wisps", "wispPlanes", "zeroFlux", "clamped",
                 "boundLiq", "boundGas", "boundMin", "boundSkipped")
RULE2_RE = re.compile(r"one plane, rule 2: faces reset to 0 (\d+), to 1 (\d+); before the reset: seeds \(alpha <= \S+ with a wet face\) (\d+), "
                      r"slivers \(eps < alpha < 1e-3 with a wet face\) (\d+)")


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def f(x, nd=3):
    v = num(x)
    return "-" if v is None else "%.*e" % (nd, v)


def g(x, nd=3):
    v = num(x)
    return "-" if v is None else "%.*g" % (nd, v)


def read_case(c):
    p = os.path.join(c, "case_params.json")
    if not os.path.isfile(p):
        return None
    t = json.load(open(p))["tokens"]
    rec = dict(case=os.path.basename(c), study=os.path.basename(os.path.dirname(c)),
               row=t.get("TPF_ROW", ""), setting=t.get("TPF_ENRIGHT_SETTING", "reference"),
               arm=t.get("TPF_VELOCITY_ARM", ""), N=int(float(t.get("N_CELLS", 0))),
               co=float(t.get("TARGET_ADVECTIVE_CO", 0)), rep=int(float(t.get("TPF_REPETITION", 1))),
               dt=num(t.get("MAX_DELTA_T")), T=num(t.get("END_TIME")), grazingCut=t.get("TPF_GRAZING_CUT", ""))
    rec["label"] = rec["row"] if rec["setting"] in ("reference", None) else rec["row"] + " [" + rec["setting"] + "]"
    csvp = os.path.join(c, "advectorVoF.csv")
    rows = list(csv.DictReader(open(csvp))) if os.path.isfile(csvp) and os.path.getsize(csvp) else []
    last = rows[-1] if rows else {}
    rec.update(Eshape=num(last.get("Eshape")), L1x=num(last.get("LCentreDiff1")), L1n=num(last.get("LNormalDiff1")),
               L1af=num(last.get("L1_alphaf")) if rec["row"].startswith("davof") else None,
               massDevPct=num(last.get("devPercent")), nInterfaceEnd=num(last.get("nInterfaceCells")),
               tLast=num(last.get("TIME")), status=last.get("status", "noCSV"), abortStep=last.get("abortStep", ""))
    st = [num(r.get("timeStepTime")) for r in rows]
    st = [x for x in st if x is not None and x > 0]
    rec["stepMedian"] = statistics.median(st) if st else None
    # the log
    logp = next((os.path.join(c, n) for n in ("log.davofAdvectorFoam", "log.advectorVoF") if os.path.isfile(os.path.join(c, n))), None)
    rec.update(minAlpha=None, excess=None, ofCo=None, steps=None, bounds=[], onePlane={}, rule2={}, seedsHalf=None,
               seedsMax=None, seedsSteps=None, sliversMax=None)
    if logp:
        log = open(logp, errors="replace").read()
        mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
        exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)]
        cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
        rec.update(minAlpha=min(mins) if mins else None, excess=max(exc) if exc else None, ofCo=max(cos) if cos else None,
                   steps=len(re.findall(r"\nTime = ", log)), bounds=list(zip(mins, exc)))
        if "ABORT at step" in log:
            rec["status"] = "aborted"
            m = re.search(r"ABORT at step (\d+)", log)
            rec["abortStep"] = m.group(1) if m else rec["abortStep"]
        # the one-plane clause per step: the sums over the run, the maxima per step
        ops = ONEPLANE_RE.findall(log)
        if ops:
            agg = {}
            for k, i in zip(ONEPLANE_KEYS, range(len(ONEPLANE_KEYS))):
                vals = [num(o[i]) for o in ops if o[i] != ""]
                vals = [v for v in vals if v is not None]
                if not vals:
                    continue
                if k == "boundMin":
                    agg[k] = min(vals)
                else:
                    agg[k + "Sum"] = sum(vals); agg[k + "MaxStep"] = max(vals)
            agg["stepsWithClause"] = len(ops)
            rec["onePlane"] = agg
        r2 = RULE2_RE.findall(log)
        if r2:
            reset0 = [int(x[0]) for x in r2]; reset1 = [int(x[1]) for x in r2]
            seeds = [int(x[2]) for x in r2]; sliv = [int(x[3]) for x in r2]
            rec["rule2"] = dict(reset0Sum=sum(reset0), reset1Sum=sum(reset1), reset0Max=max(reset0), reset1Max=max(reset1))
            rec["seedsMax"] = max(seeds); rec["seedsSteps"] = sum(1 for s in seeds if s > 0)
            rec["sliversMax"] = max(sliv); rec["sliversSteps"] = sum(1 for s in sliv if s > 0)
            # the step at T/2: the rule 2 line of the time block nearest END_TIME/2
            if rec["T"]:
                blocks = log.split("\nTime = ")[1:]
                best = None
                for b in blocks:
                    tm = num(b.split("\n", 1)[0].strip())
                    m = RULE2_RE.search(b)
                    if tm is None or not m:
                        continue
                    d = abs(tm - rec["T"]/2)
                    if best is None or d < best[0]:
                        best = (d, tm, int(m.group(3)), int(m.group(4)))
                if best:
                    rec["seedsHalf"] = best[2]; rec["sliversHalf"] = best[3]; rec["tHalf"] = best[1]
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
    rec["np"] = int(launch.get("np", 0) or 0)
    rec["rc"] = launch.get("rc", "")
    sp = os.path.join(c, ".leia_stamp.json")
    s = json.load(open(sp)) if os.path.isfile(sp) else {}
    rec["stampId"] = s.get("id", "")
    rec["stampLeia"] = s.get("record", {}).get("leiaAtRun") or s.get("record", {}).get("leia", "")
    return rec


def collect(studies_dir, pattern):
    recs = []
    for st in sorted(glob.glob(os.path.join(studies_dir, pattern))):
        for c in sorted(glob.glob(os.path.join(st, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = read_case(c)
            if r:
                recs.append(r)
    return recs


def flat(r):
    d = {k: v for k, v in r.items() if k not in ("bounds", "onePlane", "rule2")}
    for k, v in r.get("onePlane", {}).items():
        d["op_" + k] = v
    for k, v in r.get("rule2", {}).items():
        d["r2_" + k] = v
    return d


def write_csv(path, recs):
    if not recs:
        print("no records for", path); return
    keys = []
    for r in recs:
        for k in r:
            if k not in keys:
                keys.append(k)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(recs)
    print("written:", path, f"({len(recs)} rows)")


def order(pairs):
    """The least-squares order of E over h: pairs (N, E); None below three points or E <= 0."""
    pts = [(math.log(1.0/N), math.log(E)) for N, E in pairs if E is not None and E > 0]
    if len(pts) < 2:
        return None
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    den = sum((x - mx)**2 for x in xs)
    return sum((x - mx)*(y - my) for x, y in zip(xs, ys))/den if den > 0 else None


def rel(a, b):
    if a is None or b is None:
        return None
    s = max(abs(a), abs(b))
    return abs(a - b)/s if s > 0 else 0.0


def identity_pairs(r1, r8, key):
    """The np1 and np8 records paired by key; the largest relative difference over the
    converter numbers and the per-step bounds."""
    out = []
    d8 = {key(r): r for r in r8}
    for a in r1:
        b = d8.get(key(a))
        if not b:
            continue
        diffs = {}
        for k in ("Eshape", "L1x", "L1n", "L1af", "massDevPct", "minAlpha", "excess"):
            diffs[k] = rel(a[k], b[k])
        bd = None
        if a["bounds"] and b["bounds"] and len(a["bounds"]) == len(b["bounds"]):
            bd = max(max(rel(x[0], y[0]) or 0, rel(x[1], y[1]) or 0) for x, y in zip(a["bounds"], b["bounds"]))
        worst = max([v for v in diffs.values() if v is not None] + ([bd] if bd is not None else [0]))
        out.append(dict(key=key(a), steps1=a["steps"], steps8=b["steps"], boundsRel=bd, worstRel=worst,
                        verdict=("identical to round-off" if worst <= IDENT_TOL else "DIFFERS"), **{"rel_" + k: v for k, v in diffs.items()}))
    return out


def table(rows, cols, title):
    print("\n" + title)
    widths = [max(len(c), max((len(str(r.get(c, ""))) for r in rows), default=0)) for c in cols]
    print("  " + " | ".join(c.ljust(w) for c, w in zip(cols, widths)))
    print("  " + "-+-".join("-"*w for w in widths))
    for r in rows:
        print("  " + " | ".join(str(r.get(c, "")).ljust(w) for c, w in zip(cols, widths)))


def mode_plane(args):
    recs = collect(args.studies_dir, "kiss_plane_np*")
    write_csv(os.path.join(args.out, "kiss_plane_cases.csv"), [flat(r) for r in recs])
    rows = []
    worstAf = worstN = 0.0; nCells = 0; fails = []
    for r in sorted(recs, key=lambda r: (r["study"], r["setting"], r["arm"], r["co"], r["N"])):
        op = r["onePlane"]
        v = ""
        if r["setting"] == "onePlane":
            ok = (r["L1af"] is not None and r["L1n"] is not None and r["L1af"] <= P1_TOL and r["L1n"] <= P1_TOL and r["status"] != "aborted")
            v = "P1 ok" if ok else "P1 FAIL"
            nCells += 1
            if r["L1af"] is not None: worstAf = max(worstAf, r["L1af"])
            if r["L1n"] is not None: worstN = max(worstN, r["L1n"])
            if not ok: fails.append(r)
        rows.append(dict(study=r["study"], setting=r["setting"], arm=r["arm"], N=r["N"], Co=r["co"], ofCo=g(r["ofCo"]),
                         L1af=f(r["L1af"]), L1n=f(r["L1n"]), Eshape=f(r["Eshape"]), L1x=f(r["L1x"]),
                         minAlpha=g(r["minAlpha"], 6), excess=g(r["excess"]), massDev=g(r["massDevPct"]),
                         steps=r["steps"], wall_s=g(r["wall"], 4), np=r["np"], status=r["status"],
                         donors="up %s dn %s ring %s rem %s bey %s out %s" % tuple(int(op.get(k + "Sum", 0)) for k in ("donorUp", "donorDown", "donorRing", "donorRemote", "donorBeyond", "donorOutside")) if op else "",
                         noPlane="full %s empty %s wisp %s" % tuple(int(op.get(k + "Sum", 0)) for k in ("full", "empty", "wisps")) if op else "",
                         bound="liq %s gas %s min %s" % (int(op.get("boundLiqSum", 0)), int(op.get("boundGasSum", 0)), g(op.get("boundMin"))) if op else "",
                         zeroFlux=int(op.get("zeroFluxSum", 0)) if op else "", clamped=int(op.get("clampedSum", 0)) if op else "",
                         reset="0:%s 1:%s" % (r["rule2"].get("reset0Sum", ""), r["rule2"].get("reset1Sum", "")) if r["rule2"] else "",
                         seedsMax=r["seedsMax"], verdict=v))
    table(rows, ["study", "setting", "arm", "N", "Co", "ofCo", "L1af", "L1n", "Eshape", "L1x", "minAlpha", "excess", "massDev", "steps", "wall_s", "np", "status", "verdict"], "P1 the plane gate")
    table(rows, ["study", "setting", "arm", "N", "Co", "donors", "noPlane", "bound", "zeroFlux", "clamped", "reset", "seedsMax"], "the one-plane counters (sums over the run)")
    write_csv(os.path.join(args.out, "kiss_plane_table.csv"), rows)
    print(f"\nP1: {nCells} onePlane cells, max L1(alphaf) {worstAf:.3e}, max L1(n) {worstN:.3e}, failing {len(fails)}: "
          + ("HELD" if nCells and not fails else "FALSIFIED" if fails else "no cells"))
    r1 = [r for r in recs if r["study"].endswith("np1")]; r8 = [r for r in recs if r["study"].endswith("np8")]
    pairs = identity_pairs(r1, r8, lambda r: (r["setting"], r["arm"], r["N"], r["co"]))
    for p in pairs:
        p["setting"], p["arm"], p["N"], p["Co"] = p.pop("key")
    prow = [dict(setting=p["setting"], arm=p["arm"], N=p["N"], Co=p["Co"], steps=f"{p['steps1']}/{p['steps8']}",
                 Eshape=g(p["rel_Eshape"]), L1x=g(p["rel_L1x"]), L1n=g(p["rel_L1n"]), L1af=g(p["rel_L1af"]), mass=g(p["rel_massDevPct"]),
                 minAlpha=g(p["rel_minAlpha"]), bounds=g(p["boundsRel"]), worst=g(p["worstRel"]), verdict=p["verdict"]) for p in pairs]
    table(prow, ["setting", "arm", "N", "Co", "steps", "Eshape", "L1x", "L1n", "L1af", "mass", "minAlpha", "bounds", "worst", "verdict"], "the plane identity np1 vs np8 (relative differences)")
    write_csv(os.path.join(args.out, "kiss_plane_identity.csv"), prow)
    if pairs:
        w = max(p["worstRel"] for p in pairs)
        print(f"identity (plane): {len(pairs)} pairs, worst relative difference {w:.3e}: " + ("HELD" if w <= IDENT_TOL else "FALSIFIED"))


def mode_identity(args):
    recs = collect(args.studies_dir, "kiss_identity_rotation_np*")
    write_csv(os.path.join(args.out, "kiss_identity_cases.csv"), [flat(r) for r in recs])
    rows = [dict(study=r["study"], N=r["N"], Co=r["co"], ofCo=g(r["ofCo"]), steps=r["steps"], Eshape=f(r["Eshape"], 12), L1x=f(r["L1x"], 12),
                 L1n=f(r["L1n"], 12), L1af=f(r["L1af"], 12), minAlpha=g(r["minAlpha"], 15), excess=g(r["excess"], 15), massDev=g(r["massDevPct"], 15),
                 wall_s=g(r["wall"], 4), np=r["np"], status=r["status"]) for r in sorted(recs, key=lambda r: (r["co"], r["study"]))]
    table(rows, ["study", "N", "Co", "ofCo", "steps", "Eshape", "L1x", "L1n", "L1af", "minAlpha", "excess", "massDev", "wall_s", "np", "status"], "the rotation identity cases")
    r1 = [r for r in recs if r["study"].endswith("np1")]; r8 = [r for r in recs if r["study"].endswith("np8")]
    pairs = identity_pairs(r1, r8, lambda r: (r["N"], r["co"]))
    prow = [dict(N=p["key"][0], Co=p["key"][1], steps=f"{p['steps1']}/{p['steps8']}", Eshape=g(p["rel_Eshape"]), L1x=g(p["rel_L1x"]), L1n=g(p["rel_L1n"]),
                 L1af=g(p["rel_L1af"]), mass=g(p["rel_massDevPct"]), minAlpha=g(p["rel_minAlpha"]), bounds=g(p["boundsRel"]), worst=g(p["worstRel"]), verdict=p["verdict"]) for p in pairs]
    table(prow, ["N", "Co", "steps", "Eshape", "L1x", "L1n", "L1af", "mass", "minAlpha", "bounds", "worst", "verdict"], "the rotation identity np1 vs np8 (relative differences)")
    write_csv(os.path.join(args.out, "kiss_identity_table.csv"), prow)
    if pairs:
        w = max(p["worstRel"] for p in pairs)
        print(f"identity (rotation): {len(pairs)} pairs, worst relative difference {w:.3e}: " + ("HELD" if w <= IDENT_TOL else "FALSIFIED"))


def mode_rotation(args):
    recs = collect(args.studies_dir, "kiss_rotation*_N*")
    write_csv(os.path.join(args.out, "kiss_rotation_cases.csv"), [flat(r) for r in recs])
    by = {}
    for r in recs:
        by[(r["label"], r["N"], r["co"])] = r
    labels = sorted({k[0] for k in by}); cos = sorted({k[2] for k in by}); Ns = sorted({k[1] for k in by})
    rows = []
    for lab in labels:
        for co in cos:
            pts = [(N, by[(lab, N, co)]) for N in Ns if (lab, N, co) in by]
            oE = order([(N, r["Eshape"]) for N, r in pts]); oA = order([(N, r["L1af"]) for N, r in pts]); oX = order([(N, r["L1x"]) for N, r in pts]); oN = order([(N, r["L1n"]) for N, r in pts])
            for N, r in pts:
                rec = by.get(("davof-planeSurf-ringArea [recommended]", N, co)); rdf = by.get(("plicRDF", N, co))
                ratio = (r["Eshape"]/rec["Eshape"]) if (rec and rec["Eshape"] and r["Eshape"] is not None) else None
                v = ""
                if lab.endswith("[onePlane]"):
                    okE = ratio is not None and 0.9 <= ratio <= 1.1
                    okA = (rdf is not None and r["minAlpha"] is not None and rdf["minAlpha"] is not None and r["minAlpha"] >= rdf["minAlpha"])
                    v = ("P2 E_shape " + ("ok" if okE else "FAIL")) + "; min alpha " + ("ok" if okA else "FAIL")
                    if r["status"] == "aborted": v += "; ABORTED"
                rows.append(dict(row=lab, N=N, Co=co, ofCo=g(r["ofCo"]), Eshape=f(r["Eshape"]), orderE=g(oE), ratioRec=g(ratio), L1af=f(r["L1af"]), orderAf=g(oA),
                                 L1x=f(r["L1x"]), orderX=g(oX), L1n=f(r["L1n"]), orderN=g(oN), minAlpha=g(r["minAlpha"], 6), excess=g(r["excess"]),
                                 massDev=g(r["massDevPct"]), steps=r["steps"], stepMedian_s=g(r["stepMedian"]), wall_s=g(r["wall"], 4), np=r["np"], status=r["status"], verdict=v))
    table(rows, ["row", "N", "Co", "ofCo", "Eshape", "orderE", "ratioRec", "L1af", "orderAf", "L1x", "L1n", "orderN", "minAlpha", "excess", "massDev", "steps", "stepMedian_s", "wall_s", "np", "status", "verdict"], "P2 the rotation (one turn)")
    write_csv(os.path.join(args.out, "kiss_rotation_table.csv"), rows)
    one = [r for r in rows if r["row"].endswith("[onePlane]")]
    print(f"\nP2: {len(one)} onePlane cells; E_shape within [0.9, 1.1] of recommended: {sum('E_shape ok' in r['verdict'] for r in one)}; "
          f"min alpha >= plicRDF: {sum('min alpha ok' in r['verdict'] for r in one)}")


def mode_standing(args):
    tab_dir = args.standing_tables or os.path.join(REPO, "docs", "davof", "data", "standing")
    out = []
    for fam in ("shear2D", "shear3D", "deformation3D"):
        p = os.path.join(tab_dir, f"standing_{fam}.csv")
        if not os.path.isfile(p):
            print("missing", p); continue
        cells = list(csv.DictReader(open(p)))
        by = {(c["row"], int(float(c["N"])), float(c["coTarget"])): c for c in cells}
        for (row, N, co), c in sorted(by.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][0])):
            if not row.endswith("[onePlane]"):
                continue
            rdf = by.get(("plicRDF", N, co)); rec = by.get(("davof-planeSurf-ringArea [recommended]", N, co)); st = by.get(("davof-planeSurf-ringArea [signedTets]", N, co))
            comp = rec if (rec and (rec.get("status") or "").startswith("ended")) else st
            E = num(c.get("Eshape")); Ec = num(comp.get("Eshape")) if comp else None
            ratio = E/Ec if (E is not None and Ec) else None
            tS = num(c.get("stepMedian")); tT = num(st.get("stepMedian")) if st else None
            cost = tS/tT if (tS and tT) else None
            ma = num(c.get("minAlpha")); mr = num(rdf.get("minAlpha")) if rdf else None
            mass = num(c.get("massDevPct"))
            v = []
            v.append("P3-1 " + ("ok" if (c.get("status") or "").startswith("ended") else "FAIL (" + str(c.get("status")) + ")"))
            v.append("P3-2 " + ("ok" if (ma is not None and mr is not None and ma >= mr) else "FAIL" if (ma is not None and mr is not None) else "n/a"))
            v.append("P3-3 " + ("ok" if (ratio is not None and ratio <= 1.3) else "FAIL" if ratio is not None else "n/a"))
            v.append("P3-4 " + ("ok" if (mass is not None and abs(mass) < 1e-8) else "FAIL" if mass is not None else "n/a"))
            v.append("P4 " + ("ok" if (cost is not None and cost <= 0.7) else "FAIL" if cost is not None else "n/a"))
            out.append(dict(family=fam, N=N, Co=co, ofCo=c.get("ofCo"), np=c.get("np"), Eshape=f(E), EshapeComp=f(Ec), ratio=g(ratio), compRow=(comp or {}).get("row", "-"),
                            L1af=f(num(c.get("L1af"))), L1x=f(num(c.get("L1x"))), L1n=f(num(c.get("L1n"))), minAlpha=g(ma, 6), minAlphaRDF=g(mr, 6), excess=g(num(c.get("excess"))),
                            massDevPct=g(mass), stepMedian_s=g(tS), stepSignedTets_s=g(tT), costRatio=g(cost), wall_s=c.get("wallMedian", c.get("wall", "")), status=c.get("status"),
                            identical=c.get("identical", ""), verdict="; ".join(v)))
    table(out, ["family", "N", "Co", "ofCo", "np", "Eshape", "EshapeComp", "ratio", "compRow", "L1af", "L1x", "L1n", "minAlpha", "minAlphaRDF", "excess", "massDevPct", "stepMedian_s", "stepSignedTets_s", "costRatio", "status", "verdict"], "P3 / P4 the standing matrix, the onePlane row against the record")
    write_csv(os.path.join(args.out, "kiss_standing_verdicts.csv"), out)
    # P5 and the counters from the kiss logs
    recs = collect(args.studies_dir, "standing_*_kiss")
    rows = []
    for r in sorted(recs, key=lambda r: (r["study"], r["co"], r["rep"])):
        op = r["onePlane"]
        rows.append(dict(study=r["study"], N=r["N"], Co=r["co"], rep=r["rep"], steps=r["steps"], status=r["status"], seedsHalf=r.get("seedsHalf"), tHalf=g(r.get("tHalf")),
                         seedsMax=r["seedsMax"], seedsSteps=r["seedsSteps"], sliversHalf=r.get("sliversHalf"), sliversMax=r["sliversMax"], sliversSteps=r.get("sliversSteps"),
                         clipTrisPerStep=g((op.get("clipTrisSum") or 0)/max(r["steps"] or 1, 1)) if op else "",
                         donors="up %s dn %s ring %s rem %s bey %s out %s" % tuple(int(op.get(k + "Sum", 0)) for k in ("donorUp", "donorDown", "donorRing", "donorRemote", "donorBeyond", "donorOutside")) if op else "",
                         noPlane="full %s empty %s wisp %s wispPl %s" % tuple(int(op.get(k + "Sum", 0)) for k in ("full", "empty", "wisps", "wispPlanes")) if op else "",
                         bound="liq %s gas %s min %s skipped %s" % (int(op.get("boundLiqSum", 0)), int(op.get("boundGasSum", 0)), g(op.get("boundMin")), int(op.get("boundSkippedSum", 0))) if op else "",
                         clamped=int(op.get("clampedSum", 0)) if op else "", reset="0:%s 1:%s" % (r["rule2"].get("reset0Sum", ""), r["rule2"].get("reset1Sum", "")) if r["rule2"] else ""))
    table(rows, ["study", "N", "Co", "rep", "steps", "status", "seedsHalf", "tHalf", "seedsMax", "seedsSteps", "sliversHalf", "sliversMax", "clipTrisPerStep", "donors", "noPlane", "bound", "clamped", "reset"], "P5 the seeds and the one-plane counters of the kiss rows")
    write_csv(os.path.join(args.out, "kiss_standing_counters.csv"), rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=("plane", "identity", "rotation", "standing"))
    ap.add_argument("--studies-dir", default=os.path.join(REPO, "studies", "davof"))
    ap.add_argument("--out", default=os.path.join(REPO, "studies", "davof", "kiss_tables"))
    ap.add_argument("--standing-tables", default=None, help="the directory of make_standing_tables.py's CSVs (mode standing)")
    a = ap.parse_args()
    {"plane": mode_plane, "identity": mode_identity, "rotation": mode_rotation, "standing": mode_standing}[a.mode](a)


if __name__ == "__main__":
    main()
