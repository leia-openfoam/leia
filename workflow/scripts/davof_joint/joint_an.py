#!/usr/bin/env python3
"""
The clipping rules of the swept-face flux scored against the pieces (2026-11-05;
the runs of joint-run.sbatch under /work/scratch/tm83tomy/davof/agents/joint/runs).

Per case: the converter's advectorVoF.csv (the end norms L1(alphaf), L1(x),
E_shape, L1(n), the mass deviation, the run maxima of L1(alphaf) and L1(n)) and
the log (per step: the state before the bounding and its passes, the
conservative bounding's min / excess, the content bound, the half-space clips,
the faces by path of the joint / hybrid clip, timeStepTime). Writes
csv/joint_cases.csv (one row per case), csv/joint_cells.csv (one row per cell,
row and variant with the medians over the repetitions), and the LaTeX table
bodies of the note (tex2/): the plane gate, the timed cells (errors relative to
the pieces, the wall clocks and the per-step ratios), the boundedness before
and after the bounding, the clips and the faces by path, the verdict against
the registered rule. Prints the markdown tables for REPORT-swept.md.

    python3 joint_an.py [--runs DIR] [--out DIR]
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import statistics

nan = float("nan")


def logs(c):
    p = c + "/log.advectorVoF"
    if not os.path.isfile(p):
        return {}
    log = open(p).read()
    f = lambda pat: [float(x) for x in re.findall(pat, log)]
    before_min = f(r"before bounding: min\(alpha\) = (\S+),")
    before_max = f(r"before bounding: min\(alpha\) = \S+, max\(alpha\) = 1 \+ (\S+),")
    before_out = [int(x) for x in re.findall(r"cells out of \[0, 1\] by 1e-10: (\d+)", log)]
    passes = [int(x) for x in re.findall(r"bounding passes (\d+)", log)]
    mins = f(r"After conservative bounding: min\(alpha\) = (\S+),")
    exc = f(r"After conservative bounding: min\(alpha\) = \S+, max\(alpha\) = 1 \+ (\S+)")
    cb = [int(x) for x in re.findall(r"content bound on (\d+) faces", log)]
    clips = [int(x) for x in re.findall(r"half-space clips (\d+)", log)]
    faces = [int(x) for x in re.findall(r"swept(?:FaceCells)?[^\n]*?faces (\d+)", log)]
    st = f(r"timeStepTime = (\S+) s")
    d = dict(
        steps=len(re.findall(r"\nTime = ", log)),
        minBefore=min(before_min) if before_min else nan,
        excessBefore=max(before_max) if before_max else nan,
        outBeforeMax=max(before_out) if before_out else -1,
        outBeforeSum=sum(before_out) if before_out else -1,
        stepsWithOut=sum(1 for x in before_out if x > 0),
        passesSum=sum(passes) if passes else -1,
        stepsWithPass=sum(1 for x in passes if x > 0),
        minAlpha=min(mins) if mins else nan,
        excess=max(exc) if exc else nan,
        contentBoundSum=sum(cb) if cb else -1,
        clipsPerStep=statistics.median(clips) if clips else nan,
        clipsSum=sum(clips) if clips else -1,
        stepMed=statistics.median(st) if st else nan,
        clock=float(re.findall(r"ClockTime = (\S+) s", log)[-1]) if "ClockTime" in log else nan,
    )
    m = re.search(r"sweptClip (\w+), sweptHexClip (\d)", log)
    d["keys"] = m.group(0) if m else ""
    m = re.search(r"sweptHexClip (on|requested)", log)
    d["hexDetected"] = m.group(1) if m else ""
    # the faces by path, summed over the steps
    j = re.findall(r"joint clip: faces (\d+) \(by the maximum (\d+)\), pieces on a saddle (\d+), on the caps (\d+), on no plane (\d+)", log)
    if j:
        d.update(jointFaces=sum(int(x[0]) for x in j), jointMax=sum(int(x[1]) for x in j), jointSaddle=sum(int(x[2]) for x in j),
                 jointCaps=sum(int(x[3]) for x in j), jointNoPlane=sum(int(x[4]) for x in j))
    hb = re.findall(r"hybrid clip: faces inside the donor (\d+), hybrid (\d+) \(with the maximum (\d+), a saddle piece by its own plane (\d+)\), fallbacks to the pieces (\d+)", log)
    if hb:
        d.update(hybInside=sum(int(x[0]) for x in hb), hybFaces=sum(int(x[1]) for x in hb), hybMax=sum(int(x[2]) for x in hb),
                 hybSaddle=sum(int(x[3]) for x in hb), hybFallback=sum(int(x[4]) for x in hb))
    sw = re.findall(r"swept faces (\d+)", log)
    if sw:
        d["sweptFacesSum"] = sum(int(x) for x in sw)
    return d


def metrics(c):
    p = c + "/advectorVoF.csv"
    if not os.path.isfile(p):
        return {}
    R = list(csv.DictReader(open(p)))
    if not R:
        return {}
    f = lambda k: [float(r[k]) for r in R if r.get(k) not in ("", None)]
    last = R[-1]
    g = lambda k: float(last[k]) if last.get(k) not in (None, "") else nan
    d = dict(L1af=g("L1_alphaf"), L1x=g("LCentreDiff1"), Eshape=g("Eshape"), L1n=g("LNormalDiff1"), massDevPct=g("devPercent"), status=last.get("status", ""))
    d["L1afMax"] = max(f("L1_alphaf")) if f("L1_alphaf") else nan
    d["L1nMax"] = max(f("LNormalDiff1")) if f("LNormalDiff1") else nan
    d["EshapeRel"] = d["Eshape"]/float(last["sumAlphaV"]) if last.get("sumAlphaV") else nan
    return d


def tokens(c):
    src = open(c + "/.source").read().strip() if os.path.isfile(c + "/.source") else c
    return json.load(open(src + "/case_params.json"))["tokens"]


def fnum(x, digits=3):
    if x != x:
        return "--"
    if x == 0:
        return "0"
    if 1e-2 <= abs(x) < 1e3:
        return ("%%.%dg" % digits) % x
    m, e = (("%%.%de" % (digits - 1)) % x).split("e")
    return "%s\\,\\cdot\\,10^{%d}" % (m, int(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="/work/scratch/tm83tomy/davof/agents/joint/runs")
    ap.add_argument("--out", default="/work/scratch/tm83tomy/davof/agents/joint")
    args = ap.parse_args()
    os.makedirs(args.out + "/csv", exist_ok=True); os.makedirs(args.out + "/tex2", exist_ok=True)
    rel = lambda a, b: (a - b)/b if (b == b and b != 0) else nan

    # ---------------------------------------------------------------- the plane gate
    gate = lambda d: d.get("L1afMax", nan) < 1e-11 and d.get("L1nMax", nan) < 1e-11 and d.get("minAlpha", nan) >= -1e-11 and d.get("excess", nan) <= 1e-11
    prow = []
    for c in sorted(glob.glob(args.runs + "/plane*/*/planeTranslation3D_0*")):
        if not os.path.isfile(c + "/advectorVoF.csv"):
            continue
        v = open(c + "/.variant").read().strip(); t = tokens(c)
        d = {**metrics(c), **logs(c)}
        prow.append(dict(N=int(float(t["N_CELLS"])), Co=float(t["TARGET_ADVECTIVE_CO"]), arm=t["TPF_VELOCITY_ARM"], row=t["TPF_ROW"].replace("davof-planeSurf-", ""), variant=v, gate=gate(d), **d))
    prow.sort(key=lambda r: (r["N"], r["Co"], r["arm"], r["row"], r["variant"]))
    print("### The plane gate (N 20 / 40, Co 0.2 / 0.95, both arms and rows; the gate: run maxima of L1(alphaf), L1(n) < 1e-11, min alpha >= -1e-11, excess <= 1e-11)\n")
    print("| N | Co | arm | row | variant | gate | L1(alphaf) max | L1(n) max | min alpha before / after | excess before / after | out-of-[0,1] steps | bounding passes | clips/step | steps |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    with open(args.out + "/tex2/table_plane.tex", "w") as fh:
        for r in prow:
            print("| %d | %g | %s | %s | %s | %s | %.2e | %.2e | %.1e / %.1e | %.1e / %.1e | %d | %d | %.0f | %d |" % (
                r["N"], r["Co"], r["arm"], r["row"], r["variant"], "pass" if r["gate"] else "FAIL", r["L1afMax"], r["L1nMax"], r["minBefore"], r["minAlpha"], r["excessBefore"], r["excess"], r["stepsWithOut"], r["passesSum"], r["clipsPerStep"], r["steps"]))
            fh.write("%d & %g & %s & %s & %s & %s & $%s$ & $%s$ & $%s$ & $%s$ & %d & %d & %.0f \\\\\n" % (
                r["N"], r["Co"], r["arm"], r["row"], r["variant"], "pass" if r["gate"] else "fail", fnum(r["L1afMax"], 2), fnum(r["L1nMax"], 2), fnum(r["minBefore"], 2), fnum(r["excessBefore"], 2), r["stepsWithOut"], r["passesSum"], r["clipsPerStep"]))
    with open(args.out + "/csv/joint_plane.csv", "w", newline="") as fh:
        if prow:
            w = csv.DictWriter(fh, fieldnames=list(prow[0].keys())); w.writeheader(); w.writerows(prow)
    byv = {}
    for r in prow:
        byv.setdefault(r["variant"], []).append(r["gate"])
    print("\nplane: " + "; ".join("%s %d / %d pass" % (v, sum(g), len(g)) for v, g in sorted(byv.items())))

    # ---------------------------------------------------------------- the timed cells
    crow = []
    for d in sorted(glob.glob(args.runs + "/cmp/*_N*_co*")):
        m = re.match(r".*/(\w+)_N(\d+)_co(\d)$", d)
        if not m:
            continue
        flow, N, co = m.group(1), int(m.group(2)), int(m.group(3))
        T = [row for tf in sorted(glob.glob(d + "/timing*.csv")) for row in csv.DictReader(open(tf))]
        for c in sorted(glob.glob(d + "/*_r[0-9]")):
            if not os.path.isfile(c + "/log.advectorVoF"):
                continue
            b = os.path.basename(c); mm = re.match(r"(\w+?)_(\w+)_r(\d+)$", b)
            row, v, rep = mm.group(1), mm.group(2), int(mm.group(3))
            tr = [t for t in T if t["row"] == row and t["variant"] == v and int(t["rep"]) == rep]
            wall = float(tr[0]["wall_s"]) if tr and tr[0].get("wall_s") else nan
            np_ = int(tr[0]["np"]) if tr else 0
            crow.append(dict(flow=flow, N=N, co=co, np=np_, row=row, variant=v, rep=rep, wall=wall, **metrics(c), **logs(c)))
    with open(args.out + "/csv/joint_cases.csv", "w", newline="") as fh:
        if crow:
            keys = sorted(set(k for r in crow for k in r.keys()))
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(crow)
    # per cell, row, variant: the medians over the repetitions; the norms from rep 1
    cells = {}
    for r in crow:
        cells.setdefault((r["flow"], r["N"], r["co"], r["row"], r["variant"]), []).append(r)
    crows = []
    for key, L in sorted(cells.items()):
        L.sort(key=lambda r: r["rep"])
        r0 = L[0]
        walls = [r["wall"] for r in L if r["wall"] == r["wall"]]
        steps = [r["stepMed"] for r in L if r["stepMed"] == r["stepMed"]]
        e = dict(flow=key[0], N=key[1], co=key[2], row=key[3], variant=key[4], np=r0["np"], reps=len(L),
                 wallMed=statistics.median(walls) if walls else nan, wallMin=min(walls) if walls else nan, wallMax=max(walls) if walls else nan,
                 stepMed=statistics.median(steps) if steps else nan)
        for k in ("L1af", "L1x", "Eshape", "L1n", "massDevPct", "minAlpha", "excess", "minBefore", "excessBefore", "outBeforeMax", "outBeforeSum", "stepsWithOut", "passesSum", "stepsWithPass", "contentBoundSum", "clipsPerStep", "steps", "keys", "hexDetected", "status",
                  "jointFaces", "jointMax", "jointSaddle", "jointCaps", "jointNoPlane", "hybInside", "hybFaces", "hybMax", "hybSaddle", "hybFallback", "sweptFacesSum"):
            e[k] = r0.get(k, nan)
        crows.append(e)
    with open(args.out + "/csv/joint_cells.csv", "w", newline="") as fh:
        if crows:
            keys = sorted(set(k for r in crows for k in r.keys()))
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(crows)
    ref = {(e["flow"], e["N"], e["co"], e["row"]): e for e in crows if e["variant"] == "pieces"}
    flowName = {"Rot": "rotation", "LineDiag": "diagonal", "LineX": "x-line"}
    print("\n### The timed cells: the variants against the pieces (errors relative to the pieces; the wall clock median [min, max] of the repetitions; the step the median timeStepTime)\n")
    print("| cell | OF Co | row | variant | E_shape (rel.) | L1(alphaf) (rel.) | L1(x) (rel.) | L1(n) (rel.) | min alpha before / after | excess before / after | steps out of [0,1] / passes | clips/step (rel.) | wall s [min, max] | pieces / variant | step s |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    verdict = []
    with open(args.out + "/tex2/table_cmp.tex", "w") as fh, open(args.out + "/tex2/table_bound2.tex", "w") as fb, open(args.out + "/tex2/table_cost2.tex", "w") as fc:
        for e in crows:
            p = ref.get((e["flow"], e["N"], e["co"], e["row"]))
            if not p:
                continue
            rE, rA, rX, rN = rel(e["Eshape"], p["Eshape"]), rel(e["L1af"], p["L1af"]), rel(e["L1x"], p["L1x"]), rel(e["L1n"], p["L1n"])
            rW = p["wallMed"]/e["wallMed"] if e["wallMed"] == e["wallMed"] and e["wallMed"] > 0 else nan
            rS = p["stepMed"]/e["stepMed"] if e["stepMed"] == e["stepMed"] and e["stepMed"] > 0 else nan
            rC = e["clipsPerStep"]/p["clipsPerStep"] if p["clipsPerStep"] == p["clipsPerStep"] and p["clipsPerStep"] > 0 else nan
            cell = "%s N %d, %d rank%s" % (flowName.get(e["flow"], e["flow"]), e["N"], e["np"], "s" if e["np"] != 1 else "")
            print("| %s | %s | %s | %s | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.4e (%+.1e) | %.1e / %.1e | %.1e / %.1e | %d / %d | %.0f (%.2f) | %.1f [%.1f, %.1f] | %.2f | %.3f |" % (
                cell, "1" if e["co"] else "0.2", e["row"], e["variant"], e["Eshape"], rE, e["L1af"], rA, e["L1x"], rX, e["L1n"], rN, e["minBefore"], e["minAlpha"], e["excessBefore"], e["excess"], e["stepsWithOut"], e["passesSum"], e["clipsPerStep"], rC, e["wallMed"], e["wallMin"], e["wallMax"], rW, e["stepMed"]))
            lab = "%s & %d & %d & %s & %s & %s" % (flowName.get(e["flow"], e["flow"]), e["N"], e["np"], "1" if e["co"] else "0.2", e["row"], e["variant"])
            fh.write("%s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & %.2f & %.2f \\\\\n" % (lab, fnum(e["Eshape"], 4) + ("" if e["variant"] == "pieces" else " (%+.1e)" % rE), fnum(e["L1af"], 4) + ("" if e["variant"] == "pieces" else " (%+.1e)" % rA),
                                                                             fnum(e["L1x"], 4) + ("" if e["variant"] == "pieces" else " (%+.1e)" % rX), fnum(e["L1n"], 4) + ("" if e["variant"] == "pieces" else " (%+.1e)" % rN), fnum(e["massDevPct"], 2), e["wallMed"], rW))
            fb.write("%s & $%s$ & $%s$ & $%s$ & $%s$ & %d / %d & %d / %d & %s \\\\\n" % (lab, fnum(e["minBefore"], 2), fnum(e["excessBefore"], 2), fnum(e["minAlpha"], 2), fnum(e["excess"], 2), e["stepsWithOut"], e["steps"], e["passesSum"], e["stepsWithPass"], str(e["contentBoundSum"])))
            extra = ""
            if e["variant"] == "joint" and e.get("jointFaces", nan) == e.get("jointFaces", nan):
                extra = "%d / %d / %d / %d / %d" % (e["jointFaces"], e["jointMax"], e["jointSaddle"], e["jointCaps"], e["jointNoPlane"])
            if e["variant"].startswith("hybrid") and e.get("hybInside", nan) == e.get("hybInside", nan):
                extra = "%d / %d / %d / %d / %d" % (e["hybInside"], e["hybFaces"], e["hybMax"], e["hybSaddle"], e["hybFallback"])
            fc.write("%s & %.0f & %.2f & %.3f & %.2f & %.1f & %.2f & %s \\\\\n" % (lab, e["clipsPerStep"], rC, e["stepMed"], rS, e["wallMed"], rW, extra))
            if e["variant"] != "pieces" and e["row"] == "ring":
                okB = e["minBefore"] >= -1e-10 and e["excessBefore"] <= 1e-10
                okE = all(x == x and x <= 0.05 for x in (rE, rA, rX, rN)) or all(x == x and x <= 0.05 for x in (rE, rA, rX))
                okT = rS == rS and rS >= 1.5
                verdict.append((cell, "1" if e["co"] else "0.2", e["variant"], okB, max(rE, rA, rX, rN), okE, rS, okT))
    print("\n### The registered rule per cell (ring row): bounded to 1e-10 before the bounding; within 1.05 of the pieces' errors; at least 1.5x faster per step\n")
    print("| cell | OF Co | variant | bounded before | max rel. error | within 1.05 | step speed-up | >= 1.5 |")
    print("|---|---|---|---|---|---|---|---|")
    with open(args.out + "/tex2/verdict2.tex", "w") as fh:
        for cell, co, v, okB, mx, okE, rS, okT in verdict:
            print("| %s | %s | %s | %s | %+.2e | %s | %.2f | %s |" % (cell, co, v, "yes" if okB else "no", mx, "yes" if okE else "no", rS, "yes" if okT else "no"))
            fh.write("%s & %s & %s & %s & $%s$ & %s & %.2f & %s \\\\\n" % (cell, co, v, "yes" if okB else "no", fnum(mx, 2), "yes" if okE else "no", rS, "yes" if okT else "no"))
    # the hex check: the written fields are compared in the job's log (cmp); here the norms
    for e in crows:
        if e["variant"] in ("piecesHex", "hybridHex"):
            base = e["variant"][:-3]
            p = next((x for x in crows if x["variant"] == base and (x["flow"], x["N"], x["co"], x["row"]) == (e["flow"], e["N"], e["co"], e["row"])), None)
            if p:
                print("hex check %s N %d co %d %s: %s vs %s: dE_shape %+.1e dL1af %+.1e dL1x %+.1e; detected '%s'; step %.3f vs %.3f s" % (
                    e["flow"], e["N"], e["co"], e["row"], e["variant"], base, rel(e["Eshape"], p["Eshape"]), rel(e["L1af"], p["L1af"]), rel(e["L1x"], p["L1x"]), e["hexDetected"], e["stepMed"], p["stepMed"]))


if __name__ == "__main__":
    main()
