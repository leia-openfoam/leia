#!/usr/bin/env python3
"""The scoring of the decisions 1-3 matrix (2026-11-07): G1 (the eight-case identity:
every setting against closure, the counters asserted), G4 (the wall clock table), G5 (the
parallel identity: np k against np 1 after reconstruction). G2 / G3 are scored by
enright_diag_an.py with the tags matrix/enright32 and matrix/enright64.
    python3 enright_matrix_an.py G1 | G4 | G5
Writes csv/matrix_<group>.csv and prints the table. The runs' per-step counters are read
from their logs (rank fallbacks, signed tets, bounded faces, the rule's counter, the
validation maximum, L1(alphaf) of the plane cases from the converter).
"""
import csv, glob, os, re, statistics, subprocess, sys

W = os.environ.get("ENRIGHT_W", "/work/scratch/tm83tomy/davof/agents/enright")
OUT = W + "/csv"; os.makedirs(OUT, exist_ok=True)
CMP = os.path.dirname(os.path.abspath(__file__)) + "/enright_compare.py"
PY = sys.executable


def counters(c):
    txt = open(c + "/log.advectorVoF").read() if os.path.isfile(c + "/log.advectorVoF") else ""
    g = lambda r: [float(x) for x in re.findall(r, txt)]
    d = dict(steps=len(re.findall(r"\nTime = ", txt)), ended="\nEnd" in txt, abort=txt.count("ABORT"),
             rankFb=sum(g(r"rank fallbacks (\d+)")), tets=sum(g(r"signed tets (\d+)")), closurePieces=sum(g(r"piece volumes \([^)]*\): closure (\d+)")),
             validMax=max(g(r"max \|V_closure - V_tets\|/V_c ([0-9.eE+-]+)") or [float("nan")]),
             boundedFaces=sum(g(r"bounded face fractions \([^)]*\): faces (\d+)")), boundedZero=sum(g(r"s = 0 on (\d+)")), boundedPartial=sum(g(r", partial (\d+), max")),
             maxOneMinusS=max(g(r"max \|1 - s\| ([0-9.eE+-]+)") or [0.0]), recomputed=sum(g(r"recomputed (\d+) \(")), recompEvalMax=max(g(r"positioning evaluations mean [0-9.eE+-]+ max (\d+)") or [0.0]),
             ruleOutside=sum(g(r"alpha outside \(tol, 1 - tol\) (\d+)")), ruleOutsidePlane=sum(g(r"alpha outside \(tol, 1 - tol\) \d+, with a plane (\d+)")),
             ruleReset=sum(g(r"reset cells (\d+)")), ruleResetFaces=sum(g(r"reset cells \d+, faces (\d+)")),
             minAlpha=min(g(r"min\(alpha1\) = (\S+) at cell") or [float("nan")]), excess=max([-x for x in g(r"1-max\(alpha1\) = (\S+) at cell")] or [float("nan")]),
             minBefore=min(g(r"before the bounding: min\(alpha\) = (\S+),") or [float("nan")]), excBefore=max(g(r"before the bounding: min\(alpha\) = \S+, max\(alpha\) = 1 \+ (\S+)") or [float("nan")]),
             clock=(g(r"ClockTime = (\S+) s") or [float("nan")])[-1], stepMed=statistics.median(g(r"timeStepTime = (\S+) s")) if g(r"timeStepTime = (\S+) s") else float("nan"))
    csvf = c + "/advectorVoF.csv"
    if os.path.isfile(csvf):
        rows = [r for r in csv.DictReader(open(csvf)) if r.get("L1_alphaf")]
        d["L1af"] = float(rows[-1]["L1_alphaf"]) if rows else float("nan")
        rows2 = [r for r in csv.DictReader(open(csvf)) if r.get("Eshape")]
        d["Eshape"] = float(rows2[-1]["Eshape"]) if rows2 else float("nan")
    return d


def compare(ref, case):
    r = subprocess.run([PY, CMP, ref, case], capture_output=True, text=True)
    p = r.stdout.strip().split(",")
    return dict(maxDalpha=p[3], nDalpha=p[4], maxDalphaf=p[5], nDalphaf=p[6], identical=p[7]) if len(p) >= 8 else dict(maxDalpha="nan", nDalpha="nan", maxDalphaf="nan", nDalphaf="nan", identical="0")


def g1():
    D = W + "/matrix/G1_identity"
    rows = []
    cases = sorted(set(os.path.basename(c).split("__")[0] for c in glob.glob(D + "/*__closure")))
    for case in cases:
        ref = "%s/%s__closure" % (D, case)
        for c in sorted(glob.glob("%s/%s__*" % (D, case))):
            setting = os.path.basename(c).split("__")[1]
            d = counters(c); cm = compare(ref, c) if setting != "closure" else dict(maxDalpha="0", nDalpha="0", maxDalphaf="0", nDalphaf="0", identical="1")
            rows.append(dict(case=case, setting=setting, **d, **cm))
    with open(OUT + "/matrix_G1.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("case setting steps End | rankFb tets closure validMax | bounded(faces zero partial) max|1-s| recomp | rule outside/plane reset | min excess | L1af | maxDalpha nDalpha maxDalphaf nDalphaf identical")
    for r in rows:
        print("%-20s %-15s %3d %d | %5d %6d %7d %s | %5d %4d %4d %.1e %4d | %4d/%-4d %4d | %.1e %.1e | %s | %s %s %s %s %s" % (
            r["case"], r["setting"], r["steps"], r["ended"], r["rankFb"], r["tets"], r["closurePieces"], ("%.1e" % r["validMax"]) if r["validMax"] == r["validMax"] else "-",
            r["boundedFaces"], r["boundedZero"], r["boundedPartial"], r["maxOneMinusS"], r["recomputed"], r["ruleOutside"], r["ruleOutsidePlane"], r["ruleReset"],
            r["minAlpha"], r["excess"], ("%.3e" % r["L1af"]) if r.get("L1af", float("nan")) == r.get("L1af") else "-", r["maxDalpha"], r["nDalpha"], r["maxDalphaf"], r["nDalphaf"], r["identical"]))


def g4():
    D = W + "/matrix/G4_timing"
    t = list(csv.DictReader(open(D + "/timing.csv")))
    rows = []
    for x in t:
        m = re.match(r"(\w+?)_np(\d+)_r(\d+)", x["name"])
        c = D + "/" + x["name"]; d = counters(c)
        rows.append(dict(setting=m.group(1), np=int(m.group(2)), rep=int(m.group(3)), wall=float(x["wall_s"]), clock=float(x["clock_s"] or "nan"), stepMed=float(x["median_step_s"] or "nan"), steps=int(x["steps"]), rc=x["rc"], rankFb=d["rankFb"], tets=d["tets"], minAlpha=d["minAlpha"], Eshape=d.get("Eshape", float("nan"))))
    with open(OUT + "/matrix_G4.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("setting np | wall median [min, max] s | step median s | ratio to closure | rankFb tets | min alpha | Eshape(end)")
    for np_ in (1, 8):
        base = None
        for s in ("closure", "cTT", "sT"):
            rs = [r for r in rows if r["setting"] == s and r["np"] == np_]
            if not rs: continue
            med = statistics.median(r["wall"] for r in rs)
            if s == "closure": base = med
            print("%-8s %d | %.1f [%.1f, %.1f] | %.4f | %.3f | %d %d | %.1e | %.3e" % (s, np_, med, min(r["wall"] for r in rs), max(r["wall"] for r in rs), statistics.median(r["stepMed"] for r in rs), med/base if base else float("nan"), rs[0]["rankFb"], rs[0]["tets"], rs[0]["minAlpha"], rs[0]["Eshape"]))


def g5():
    D = W + "/matrix/G5_parallel"
    rows = []
    for c in sorted(glob.glob(D + "/*_np*")):
        m = re.match(r"(.+)_np(\d+)$", os.path.basename(c))
        if not m: continue
        np_ = int(m.group(2))
        if np_ > 1 and not os.path.isdir(c + "/0.") and glob.glob(c + "/processor0/[0-9]*"):
            # reconstruct the last time for the comparison (the fields only)
            tds = sorted((float(os.path.basename(d)), os.path.basename(d)) for d in glob.glob(c + "/processor0/[0-9]*"))
            subprocess.run("cd %s && reconstructPar -time %s -fields '(alpha.water alphaf.water)' > log.reconstructPar 2>&1" % (c, tds[-1][1]), shell=True)
        ref = "%s/%s_np1" % (D, m.group(1))
        d = counters(c); cm = compare(ref, c) if np_ > 1 else dict(maxDalpha="0", nDalpha="0", maxDalphaf="0", nDalphaf="0", identical="1")
        rows.append(dict(setting=m.group(1), np=np_, **d, **cm))
    with open(OUT + "/matrix_G5.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("setting np steps | rankFb tets bounded ruleOutside | min excess | maxDalpha nDalpha maxDalphaf nDalphaf identical (vs np 1)")
    for r in rows:
        print("%-16s %d %3d | %5d %6d %5d %4d | %.1e %.1e | %s %s %s %s %s" % (r["setting"], r["np"], r["steps"], r["rankFb"], r["tets"], r["boundedFaces"], r["ruleOutside"], r["minAlpha"], r["excess"], r["maxDalpha"], r["nDalpha"], r["maxDalphaf"], r["nDalphaf"], r["identical"]))


if __name__ == "__main__":
    {"G1": g1, "G4": g4, "G5": g5}[sys.argv[1]]()
