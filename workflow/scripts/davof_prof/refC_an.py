"""C (2026-11-05): the flux isolated. Per cell the rows swept-exact and swept-ring (the reference rerun,
ref/cmp), plane-clipped-exact and plane-clipped-ring (ref/cmpC: the preprint's loop, prism flux) and
plicRDF (the comparison, csv/cmp_table.csv): the norms, the bounds, the mass, the wall clock at equal
ranks (median [min, max] of the repetitions), and the flux's share of the exact-normal advantage,
(L1x_plane-exact - L1x_swept-exact) / (L1x_plicRDF - L1x_swept-exact). Writes csv/fluxC.csv.
    python3 refC_an.py"""
import csv, glob, os, re, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
W = "/work/scratch/tm83tomy/davof/agents/swept"
nan = float("nan")
def logs(c):
    log = open(c + "/log.advectorVoF").read()
    mins = [float(x) for x in re.findall(r"After conservative bounding: min\(alpha\) = (\S+),", log)]
    exc = [float(x) for x in re.findall(r"max\(alpha\) = 1 \+ (\S+)", log)] + [0.0]
    cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
    st = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)]
    ab = re.search(r"ABORT at step (\d+)", log)
    return dict(minAlpha=min(mins) if mins else nan, excess=max(exc), ofCo=max(cos) if cos else nan, steps=len(re.findall(r"\nTime = ", log)), step_med=statistics.median(st) if st else nan,
                abort=("diverged, the abort guard at step %s of %d" % (ab.group(1), len(re.findall(r"\nTime = ", log))) if ab else ""))
def metrics(c):
    p = c + "/advectorVoF.csv"
    if not os.path.isfile(p): return {}
    R = list(csv.DictReader(open(p)))
    if not R: return {}
    last = R[-1]; g = lambda k: float(last[k]) if last.get(k) not in (None, "") else nan
    return dict(L1af=g("L1_alphaf"), L1x=g("LCentreDiff1"), Eshape=g("Eshape"), L1n=g("LNormalDiff1"), massDevPct=g("devPercent"), status=last.get("status", ""))
def walls(d, meth):
    T = [row for tf in sorted(glob.glob(d + "/timing*.csv")) for row in csv.DictReader(open(tf))]
    w = [float(t["wall_s"]) for t in T if t["method"] == meth and t.get("wall_s")]
    reps = sorted(set(t["rep"] for t in T if t["method"] == meth))
    return w, reps
cmp = {}
for r in csv.DictReader(open(W + "/csv/cmp_table.csv")):
    cmp[(r["flow"], int(r["N"]), 1 if float(r["ofCo"]) > 0.5 else 0, r["method"])] = r
rows = []
print("| cell | OF Co | loop, normal | L1(x) | E_shape | L1(alphaf) | L1(n) | min / excess | mass % | wall median [min, max] s | step s |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for FLOW in ("Rot", "LineDiag"):
    for N in (32, 64, 128):
        for co in (0, 1):
            d = "%s/ref/cmp/%s_N%d_co%d" % (W, FLOW, N, co); dC = "%s/ref/cmpC/%s_N%d_co%d" % (W, FLOW, N, co)
            cell = {}
            for lab, dd, meth in (("swept, exact", d, "exact"), ("swept, ring", d, "ring"), ("plane-clipped, exact", dC, "prismExact"), ("plane-clipped, ring", dC, "prismRing")):
                w, reps = walls(dd, meth) if os.path.isdir(dd) else ([], [])
                if not reps: continue
                c = "%s/%s_r%s" % (dd, meth, reps[0])
                a = {**metrics(c), **logs(c)}; a["wall"] = w; cell[lab] = a
            ref = cmp.get((FLOW, N, co, "rdf"))
            if ref:
                cell["plicRDF + isoAdvector"] = dict(L1x=float(ref["L1x"]), Eshape=float(ref["Eshape"]), L1af=nan, L1n=float(ref["L1n"]), minAlpha=float(ref["minAlpha"]), excess=float(ref["excess"]),
                                                    massDevPct=float(ref["massDevPct"]), wall=[float(ref["wall_med"]), float(ref["wall_min"]), float(ref["wall_max"])], step_med=float(ref["step_med"]), ofCo=float(ref["ofCo"]), steps=int(ref["steps"]))
            if not cell: continue
            for lab, a in cell.items():
                w = a["wall"]
                ab = a.get("abort", "")
                fmt = lambda x: ("%.4e" % x) if x == x else "-"
                print("| %s N %d, %d ranks | %.3f | %s | %s | %s | %s | %s | %.1e / %.1e%s | %.1e | %.1f [%.1f, %.1f] | %.3f |" % (
                    "rotation" if FLOW == "Rot" else "diagonal", N, {32: 1, 64: 8, 128: 64}[N], a["ofCo"], lab, fmt(a["L1x"]), fmt(a["Eshape"]), fmt(a["L1af"]), fmt(a["L1n"]),
                    a["minAlpha"], a["excess"], (" (%s)" % ab) if ab else "", a["massDevPct"], statistics.median(w) if w else nan, min(w) if w else nan, max(w) if w else nan, a["step_med"]))
                a["abort_"] = ab
                rows.append(dict(flow=FLOW, N=N, co=co, ofCo=a["ofCo"], row=lab, L1x=a["L1x"], Eshape=a["Eshape"], L1af=a["L1af"], L1n=a["L1n"], minAlpha=a["minAlpha"], excess=a["excess"], massDevPct=a["massDevPct"],
                                 wall_med=statistics.median(w) if w else nan, wall_min=min(w) if w else nan, wall_max=max(w) if w else nan, reps=len(w), step_med=a["step_med"], steps=a["steps"], status=a.get("abort", "") or "completed"))
            se, pe, rdf = cell.get("swept, exact"), cell.get("plane-clipped, exact"), cell.get("plicRDF + isoAdvector")
            if se and pe and rdf:
                share = (pe["L1x"] - se["L1x"]) / (rdf["L1x"] - se["L1x"]) if rdf["L1x"] != se["L1x"] else nan
                shareE = (pe["Eshape"] - se["Eshape"]) / (rdf["Eshape"] - se["Eshape"]) if rdf["Eshape"] != se["Eshape"] else nan
                print("|  | | the flux's share of the exact-normal advantage: L1(x) %.3f, E_shape %.3f; plane-clipped / swept (exact): L1(x) %.3f, E_shape %.3f, wall %.2f | | | | | | | | |" % (
                    share, shareE, pe["L1x"]/se["L1x"], pe["Eshape"]/se["Eshape"], statistics.median(pe["wall"])/statistics.median(se["wall"]) if pe["wall"] and se["wall"] else nan))
                rows[-1]["share_L1x"] = share; rows[-1]["share_Eshape"] = shareE
            sr, pr = cell.get("swept, ring"), cell.get("plane-clipped, ring")
            if sr and pr:
                print("|  | | plane-clipped / swept (ring): L1(x) %.3f, E_shape %.3f, L1(alphaf) %.3f, wall %.2f; plane-clipped-ring / plicRDF wall %.2f | | | | | | | | |" % (
                    pr["L1x"]/sr["L1x"], pr["Eshape"]/sr["Eshape"], pr["L1af"]/sr["L1af"], statistics.median(pr["wall"])/statistics.median(sr["wall"]) if pr["wall"] and sr["wall"] else nan,
                    statistics.median(pr["wall"])/rdf["wall"][0] if (pr["wall"] and rdf) else nan))
keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open(W + "/csv/fluxC.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
print("wrote", W + "/csv/fluxC.csv", len(rows))
