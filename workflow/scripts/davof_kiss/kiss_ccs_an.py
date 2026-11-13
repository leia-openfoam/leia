#!/usr/bin/env python3
"""The ccsFail lines of kiss_ccsdiag_np1: per case the classes (maxIter with a plane, maxIter
without, no plane with another exit) and the quantiles of alpha, |m_c|/V^(2/3), |n|."""
import glob, json, os, re, math, sys
S = "/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof/" + (sys.argv[1] if len(sys.argv) > 1 else "kiss_ccsdiag_np1")
L = re.compile(r"ccsFail t (\d+) cell (\d+) alpha (\S+) mRel (\S+) nMag (\S+) exit (\d) plane (\d) evals (\d+) rFinal (\S+)")
def q(xs, ps=(0, 0.1, 0.5, 0.9, 1.0)):
    xs = sorted(xs)
    if not xs: return "-"
    return " ".join("%.2g" % xs[min(int(p*(len(xs)-1)+0.5), len(xs)-1)] for p in ps)
for c in sorted(glob.glob(S + "/planeTranslation3D_*")):
    t = json.load(open(c + "/case_params.json"))["tokens"]
    log = open(c + "/log.davofAdvectorFoam", errors="replace").read()
    rows = [(int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)), int(m.group(6)), int(m.group(7)), int(m.group(8)), float(m.group(9))) for m in L.finditer(log)]
    status = "ABORT" if "ABORT" in log else ("done" if os.path.isfile(c + "/advectorVoF.csv") else "?")
    print("== N %s surfCellTol %s wispTol %s: %s, ccsFail lines %d (cap 60000)" % (t["N_CELLS"], t["TPF_SURF_CELL_TOL"], t["TPF_RECON_WISP_TOL"], status, len(rows)))
    classes = {"maxIter, plane found": [r for r in rows if r[5] == 4 and r[6] == 1], "maxIter, no plane": [r for r in rows if r[5] == 4 and r[6] == 0],
               "no plane, exit 1 (tolerance met at the first cut)": [r for r in rows if r[5] == 1 and r[6] == 0], "no plane, exit 2/3": [r for r in rows if r[5] in (2, 3) and r[6] == 0],
               "no plane, exit 0 (zero normal)": [r for r in rows if r[5] == 0 and r[6] == 0]}
    for k, v in classes.items():
        if not v: continue
        print("  %-52s n %6d | alpha min/10/50/90/max %s | mRel %s | |n| %s | evals %s | |rFinal| %s" % (
            k, len(v), q([min(r[2], 1 - r[2]) if r[2] > 0.5 else r[2] for r in v]), q([r[3] for r in v]), q([r[4] for r in v]), q([float(r[7]) for r in v]), q([abs(r[8]) for r in v])))
    # the alpha classes of all failures: below 1e-14 (isoFaceTol), 1e-14 .. 1e-8, above
    al = [min(r[2], 1 - r[2]) for r in rows]
    print("  failures by alpha: <= 1e-14: %d; 1e-14 .. 1e-10: %d; 1e-10 .. 1e-8: %d; 1e-8 .. 1e-3: %d; > 1e-3: %d; alpha <= 0 (negative round-off): %d" % (
        sum(a <= 1e-14 for a in al), sum(1e-14 < a <= 1e-10 for a in al), sum(1e-10 < a <= 1e-8 for a in al), sum(1e-8 < a <= 1e-3 for a in al), sum(a > 1e-3 for a in al), sum(r[2] <= 0 for r in rows)))
    nm = [r[4] for r in rows]
    print("  failures with |n| = 0: %d; with |n| < 1e-10 (relative to the unit normal the model returns |n| ~ 1 or the raw m): %d" % (sum(x == 0 for x in nm), sum(0 < x < 1e-10 for x in nm)))
    nev = len(re.findall(r"floor exits", log))
