#!/usr/bin/env python3
"""The contact line study with the contact cell chosen by the largest |m_k| (analysis variant), next to the standing measure
(the outermost candidate). Reads the candidate dumps of studies/davof/cl_all_* (contactLineError key dumpCandidates).
Reference: M. Fricke, T. Maric, D. Bothe, Contact line advection using the geometrical Volume-of-Fluid method,
J. Comput. Phys. 407 (2020) 109221, doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.
    python3 contactline_an2.py <leia dir> <out.md>"""
import sys, os, glob, json, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contactline_an as a
R = sys.argv[1]; OUT = sys.argv[2]
COS = a.COS_IN
NS = a.NS
ROWS = a.ROWS
RS = {"plicRDF": "(1)", "pointIso": "(2)", "pointIso-davof": "(3)"}
def read(c):
    f = os.path.join(c, "postProcessing", "contactLineError", "dump", "candidates.dat")
    cand = collections.defaultdict(list)
    for l in open(f):
        w = l.split()
        if w[0] == "C":
            cand[(float(w[1]), int(w[2]))].append((int(float(w[4])), float(w[10]), float(w[11]), float(w[13])))   # cell, x, theta, |m|
    return cand
def outer(side):
    return (lambda l: min(l, key=lambda c: c[1])) if side == 0 else (lambda l: max(l, key=lambda c: c[1]))
def bigm(l):
    return max(l, key=lambda c: c[3])
def evalc(cand, arm, side, pick):
    eth = ecl = 0.0
    for (t, s), l in cand.items():
        if s != side or not l:
            continue
        c = pick(l); xr, thr = a.ref(arm, t, side)
        eth = max(eth, abs(c[2] - thr)); ecl = max(ecl, abs(c[1] - xr)/a.R0)
    return eth, ecl
def ndiff(cand, side):
    n = tot = 0
    for (t, s), l in cand.items():
        if s != side or not l:
            continue
        tot += 1
        if outer(side)(l)[0] != bigm(l)[0]:
            n += 1
    return n, tot
# the standing study values (committed tables' source)
standing = {(r["arm"], r["row"], r["N"], r["coN"]): r for r in a.load()}
D = {}
for c in sorted(glob.glob(os.path.join(R, "studies", "davof", "cl_all_*", "*_0????"))):
    tk = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
    arm, row, N = tk["TPF_VELOCITY_ARM"], tk["TPF_ROW"], int(float(tk["N_CELLS"]))
    co = min(COS, key=lambda x: abs(x - float(tk["TARGET_ADVECTIVE_CO"])))
    if not os.path.isfile(os.path.join(c, "postProcessing", "contactLineError", "dump", "candidates.dat")):
        continue
    cand = read(c)
    r = {"dir": c}
    for side in (0, 1):
        e1 = evalc(cand, arm, side, outer(side)); e2 = evalc(cand, arm, side, bigm)
        r["std_th%d" % side], r["std_cl%d" % side] = e1
        r["big_th%d" % side], r["big_cl%d" % side] = e2
        r["nd%d" % side], r["nt%d" % side] = ndiff(cand, side)
    D[(arm, row, N, co)] = r
L = ["\n## Contact cell by the largest |m_k| (analysis variant next to the standing measure: the outermost candidate); Co 0.2 and 0.5 (the method's range)\n"]
L.append("|m_k| is the magnitude of the reconstruction's normal (area normal of the davof reconstruction in row (3), the normal of plicRDF in rows (1) and (2)); the choice is made within a row. Runs with a dump: %d.\n" % len(D))
# reproduction of the committed standing values
worst = 0.0; nrep = 0
for k, r in D.items():
    s = standing.get(k)
    if not s:
        continue
    for side in (0, 1):
        for q, qs in (("std_th%d" % side, "Eth%d" % side), ("std_cl%d" % side, "Ecl%d" % side)):
            if s.get(qs):
                worst = max(worst, abs(r[q]/s[qs] - 1)); nrep += 1
L.append("Control: the standing (outermost) values of the dump runs against the committed tables: %d values compared, largest relative difference %.2e.\n" % (nrep, worst))
def order(pts):
    return a.order(pts)
def col(r, k):
    return r["big_" + k]
cols = (("th0", "E_theta L"), ("th1", "E_theta R"), ("cl0", "E_cl L"), ("cl1", "E_cl R"))
ords = {}
for arm, title in (("clB", "Case B"), ("clC", "Case C"), ("clD", "Case D")):
    for co in COS:
        L.append("\n### %s, OpenFOAM Co %s: E^inf with the largest-|m_k| choice (standing value in brackets)\n" % (title, co))
        L.append("| row | N | E_theta L | E_theta R | E_cl L | E_cl R | steps where the two choices differ L / R |")
        L.append("|---|---|---|---|---|---|---|")
        for row in ROWS:
            for N in NS:
                r = D.get((arm, row, N, co))
                if not r:
                    continue
                L.append("| %s | %d | %.3e (%.3e) | %.3e (%.3e) | %.3e (%.3e) | %.3e (%.3e) | %d of %d / %d of %d |" % (
                    RS[row], N, r["big_th0"], r["std_th0"], r["big_th1"], r["std_th1"], r["big_cl0"], r["std_cl0"], r["big_cl1"], r["std_cl1"],
                    r["nd0"], r["nt0"], r["nd1"], r["nt1"]))
        L.append("\nOrders over dx/R0 (least squares; largest |m_k| choice, the standing in brackets) and the ratios of E^inf:\n")
        L.append("| row | E_theta L | E_theta R | E_cl L | E_cl R |")
        L.append("|---|---|---|---|---|")
        for row in ROWS:
            cells = []
            for k, nm in cols:
                pts = [(1.0/N/a.R0, D[(arm, row, N, co)]["big_" + k]) for N in NS if (arm, row, N, co) in D]
                pts2 = [(1.0/N/a.R0, D[(arm, row, N, co)]["std_" + k]) for N in NS if (arm, row, N, co) in D]
                o, o2 = order(pts), order(pts2)
                ords[(arm, co, row, k)] = o
                cells.append("%s (%s)" % ("-" if o is None else "%.2f" % o, "-" if o2 is None else "%.2f" % o2))
            L.append("| %s | %s |" % (RS[row], " | ".join(cells)))
        L.append("\n| N | ratio | E_theta L | E_theta R | E_cl L | E_cl R |")
        L.append("|---|---|---|---|---|---|")
        for N in NS:
            for other, nm in (("pointIso-davof", "(1)/(3)"), ("pointIso", "(1)/(2)")):
                r1, r2 = D.get((arm, "plicRDF", N, co)), D.get((arm, other, N, co))
                if r1 and r2:
                    L.append("| %d | %s | %s |" % (N, nm, " | ".join("%.3g" % (r1["big_" + k]/r2["big_" + k]) if r2["big_" + k] > 0 else "inf" for k, _ in cols)))
# the prediction
L.append("\n## Data of the prediction (main session) and of PW1 to PW3 under the largest-|m_k| choice\n")
chg = []
for (arm, row, N, co), r in D.items():
    if row in ("plicRDF", "pointIso"):
        for k, _ in cols:
            chg.append((abs(r["big_" + k]/r["std_" + k] - 1) if r["std_" + k] > 0 else 0.0, arm, row, N, co, k))
chg.sort(reverse=True)
L.append("- Rows (1) and (2): relative change of the 96 values (3 cases x 2 Co x 4 N x 4 columns, 2 rows = 192 values counted: %d); largest %.2f %% (%s); values above 5 %%: %d." % (len(chg), 100*chg[0][0], "%s %s N %d Co %s %s" % chg[0][1:], sum(1 for c in chg if c[0] > 0.05)))
for c in chg[:5]:
    L.append("  - %.2f %%: %s %s N %d Co %s %s" % ((100*c[0],) + c[1:]))
npos = ntot = 0; neg = []
for (arm, co, row, k), o in ords.items():
    if row == "pointIso-davof" and o is not None:
        ntot += 1
        if o > 0:
            npos += 1
        else:
            neg.append((arm, co, k, o))
L.append("- Row (3) orders (3 cases x 2 Co x 4 columns = 24): positive in %d of %d; not positive: %s." % (npos, ntot, "; ".join("%s Co%s %s %.2f" % x for x in neg) or "none"))
viol = []
for arm in ("clB", "clD"):
    for co in COS:
        for N in NS:
            r1, r3 = D.get((arm, "plicRDF", N, co)), D.get((arm, "pointIso-davof", N, co))
            if r1 and r3 and not r3["big_cl1"] < r1["big_cl1"]:
                viol.append((arm, co, N, r3["big_cl1"], r1["big_cl1"]))
L.append("- E_cl^inf right of cases B and D, row (3) below row (1) at every N and Co (16 cells): violations %d: %s." % (len(viol), "; ".join("%s Co%s N%d %.3e vs %.3e" % v for v in viol) or "none"))
rat = []
for (arm, row, N, co), r3 in D.items():
    if row != "pointIso-davof":
        continue
    r1 = D.get((arm, "plicRDF", N, co))
    for k in ("th0", "th1"):
        if r1 and r1["big_" + k] > 0:
            rat.append((r3["big_" + k]/r1["big_" + k], arm, N, co, k))
inside = sum(1 for x in rat if x[0] <= 2.0)
L.append("- E_theta^inf of row (3) within 2 times row (1) (both sides, all cases, N, Co: %d values): %d, i.e. %.1f %% (the prediction asks >= 80 %%)." % (len(rat), inside, 100.0*inside/max(len(rat), 1)))
rat.sort(reverse=True)
L.append("  The 5 worst (case, N, Co, side): " + "; ".join("%s N%d Co%s %s ratio %.2f" % (x[1], x[2], x[3], "LR"[int(x[4][-1])], x[0]) for x in rat[:5]))
L.append("- PW1 (case D, E_cl^inf ratio row (1) / row (2) and row (1) / row (3), at least 2 required), left/right:")
for N in NS:
    for co in COS:
        r1, r2, r3 = (D.get(("clD", row, N, co)) for row in ROWS)
        if r1 and r2 and r3:
            L.append("  - N %d Co %s: (1)/(2) %.3g / %.3g; (1)/(3) %.3g / %.3g" % (N, co, r1["big_cl0"]/r2["big_cl0"], r1["big_cl1"]/r2["big_cl1"], r1["big_cl0"]/r3["big_cl0"], r1["big_cl1"]/r3["big_cl1"]))
L.append("- PW2 (cases B, C): orders of E_cl^inf of rows (2) and (3) (>= 1 required) and E_cl row (2), (3) <= row (1) per cell (counts):")
for arm in ("clB", "clC"):
    for co in COS:
        L.append("  - %s Co %s: order (2) L %s R %s; (3) L %s R %s; cells with E_cl <= row (1): (2) %d of 8, (3) %d of 8" % (
            arm, co, *["%.2f" % ords[(arm, co, row, k)] if ords.get((arm, co, row, k)) is not None else "-" for row, k in (("pointIso", "cl0"), ("pointIso", "cl1"), ("pointIso-davof", "cl0"), ("pointIso-davof", "cl1"))],
            sum(1 for N in NS for k in ("cl0", "cl1") if (arm, "pointIso", N, co) in D and D[(arm, "pointIso", N, co)]["big_" + k] <= D[(arm, "plicRDF", N, co)]["big_" + k]),
            sum(1 for N in NS for k in ("cl0", "cl1") if (arm, "pointIso-davof", N, co) in D and D[(arm, "pointIso-davof", N, co)]["big_" + k] <= D[(arm, "plicRDF", N, co)]["big_" + k])))
L.append("- PW3 (E_theta^inf row (3) <= row (1)): cells with row (3) <= row (1), both sides, all cases, N, Co: %d of %d; orders of E_theta (rows (1), (2), (3)) are in the tables above." % (
    sum(1 for x in rat if x[0] <= 1.0), len(rat)))
open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L)[:6000])
