#!/usr/bin/env python3
"""Row (3) outliers of the contact line study (docs/davof/REPORT-swept.md "Locating the row (3) outliers").
Reads the candidate dumps (contactLineError key dumpCandidates) of studies/davof/cl_dump_{r1,r3}_<case> and writes a markdown section.
Reference: M. Fricke, T. Maric, D. Bothe, J. Comput. Phys. 407 (2020) 109221, doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.
    python3 cl_outliers.py <leia dir> <out.md>"""
import sys, os, glob, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import contactline_an as a
R = sys.argv[1]
OUT = sys.argv[2]
CASES = (("C02N1024", "clC", 0.2, 1024), ("C05N128", "clC", 0.5, 128), ("B02N256", "clB", 0.2, 256))
def load(rowtag, tag):
    c = glob.glob(os.path.join(R, "studies/davof/cl_dump_%s_%s/*_0????" % (rowtag, tag)))[0]
    f = glob.glob(os.path.join(c, "postProcessing/contactLineError/dump/candidates.dat"))[0]
    cand = collections.defaultdict(list)   # (t, side) -> list of dict
    win = {}
    for l in open(f):
        w = l.split()
        if w[0] == "C":
            t = float(w[1]); side = int(w[2])
            cand[(t, side)].append(dict(chosen=int(w[3]), cell=int(float(w[4])), alpha=float(w[5]), mkrel=float(w[6]), n=(float(w[7]), float(w[8]), float(w[9])),
                                        x=float(w[10]), th=float(w[11]), af=float(w[12]), mk=float(w[13])))
        elif w[0] == "W":
            win[(float(w[1]), int(w[2]))] = [float(x) for x in w[3:]]
    return c, cand, win
def evaluate(cand, arm, pick):
    out = {}
    for side in (0, 1):
        eth, ecl, ts = 0.0, 0.0, []
        for (t, s), lst in sorted(cand.items()):
            if s != side or not lst:
                continue
            ch = pick(lst)
            xr, thr = a.ref(arm, t, side)
            e = abs(ch["th"] - thr); x = abs(ch["x"] - xr)/a.R0
            eth = max(eth, e); ecl = max(ecl, x)
        out[side] = (eth, ecl)
    return out
def outer(side):
    return lambda lst: (min(lst, key=lambda c: c["x"]) if side == 0 else max(lst, key=lambda c: c["x"]))
def bigm(lst):
    return max(lst, key=lambda c: c["mk"])
_cs = a.load()
BY = {(r["arm"], r["row"], r["N"], r["coN"]): r for r in _cs}
L = ["\n## row (3) outliers (print-only dump of every candidate contact cell; Co 0.2 and 0.5; facts only)\n"]
worst_c1024 = None
for tag, arm, co, N in CASES:
    d = {}
    for rowtag in ("r1", "r3"):
        d[rowtag] = load(rowtag, tag)
    L.append("\n### %s, %s, Co %s, N %d\n" % (tag, arm, co, N))
    # (c) the standing measure (the outermost cell), the largest-|m| choice, per row and side
    L.append("(c) E_theta^inf (rad) and E_cl^inf, standing choice (outermost) and the analysis variant (largest |m_k|):\n")
    L.append("| row | side | E_theta outermost (dump run) | E_theta standing study | E_theta largest |m| | E_cl outermost (dump run) | E_cl standing study | E_cl largest |m| |")
    L.append("|---|---|---|---|---|---|---|---|")
    for rowtag, nm in (("r1", "(1) plicRDF"), ("r3", "(3) pointIso-davof")):
        c, cand, win = d[rowtag]
        for side in (0, 1):
            e1 = evaluate(cand, arm, outer(side))[side]
            e2 = evaluate(cand, arm, bigm)[side]
            st = BY.get((arm, "plicRDF" if rowtag == "r1" else "pointIso-davof", N, co))
            sth = st["Eth%d" % side] if st else float("nan"); scl = st["Ecl%d" % side] if st else float("nan")
            L.append("| %s | %s | %.3e | %.3e | %.3e | %.3e | %.3e | %.3e |" % (nm, "LR"[side], e1[0], sth, e2[0], e1[1], scl, e2[1]))
    # (a) worst steps, (b) chosen alpha outside (1e-3, 1 - 1e-3)
    for rowtag, nm in (("r3", "(3) pointIso-davof"), ("r1", "(1) plicRDF")):
        c, cand, win = d[rowtag]
        c1, cand1, win1 = d["r1"]
        times = sorted(set(t for (t, s) in cand))
        for side in (0, 1):
            rows = []
            nout = 0
            nsteps = 0
            for t in times:
                lst = cand.get((t, side), [])
                if not lst:
                    continue
                nsteps += 1
                ch = outer(side)(lst)
                if ch["alpha"] <= 1e-3 or ch["alpha"] >= 1 - 1e-3:
                    nout += 1
                xr, thr = a.ref(arm, t, side)
                rows.append((abs(ch["th"] - thr), t, lst, ch, xr, thr))
            rows.sort(key=lambda r: -r[0])
            L.append("\n%s, side %s: steps with a candidate %d; steps in which the chosen cell has alpha_k outside (1e-3, 1 - 1e-3): %d; steps with more than one candidate: %d\n" % (
                nm, "LR"[side], nsteps, nout, sum(1 for t in times if len(cand.get((t, side), [])) > 1)))
            L.append("(a) the 5 worst steps by |theta_num - theta_ref| of the chosen cell (all candidates of the step; * = chosen):\n")
            L.append("| step | t | cell | alpha_k | |m|/V^(2/3) | x_num | theta_num (rad) | theta_ref | alphaf wall (before the step) | |dtheta| | row (1) |dtheta| at that step |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|")
            for e, t, lst, ch, xr, thr in rows[:5]:
                step = times.index(t)
                l1 = cand1.get((t, side), [])
                e1 = "-"
                if l1:
                    c1x = outer(side)(l1)
                    e1 = "%.3e" % abs(c1x["th"] - thr)
                for cc in sorted(lst, key=lambda c: c["x"]):
                    L.append("| %d | %.5f | %d%s | %.3e | %.3f | %.6f | %.5f | %.5f | %.3e | %.3e | %s |" % (step, t, cc["cell"], "*" if cc is ch else "", cc["alpha"], cc["mkrel"], cc["x"], cc["th"], thr, cc["af"], abs(cc["th"] - thr), e1))
            if tag == "C02N1024" and rowtag == "r3" and side == 0 and rows:
                worst_c1024 = (rows[0][1], side)
    if tag == "C02N1024":
        # (d) the wall cells around the chosen one at the worst step of row (3)
        for side in (0, 1):
            c3, cand3, win3 = d["r3"]
            times = sorted(set(t for (t, s) in cand3))
            rows = []
            for t in times:
                lst = cand3.get((t, side), [])
                if lst:
                    ch = outer(side)(lst); xr, thr = a.ref(arm, t, side)
                    rows.append((abs(ch["th"] - thr), t))
            rows.sort(reverse=True)
            if not rows:
                continue
            t = rows[0][1]
            L.append("\n(d) case C N 1024, row (3), side %s, the worst step t = %.5f: the 10 wall cells around the chosen cell (x, alpha_k before the step, alphaf of the wall face before the step), then row (1) at the same t:\n" % ("LR"[side], t))
            for rowtag, nm in (("r3", "row (3)"), ("r1", "row (1)")):
                w = d[rowtag][2].get((t, side))
                if w:
                    cells = ["(%.5f, %.3e, %.3e)" % (w[i], w[i + 1], w[i + 2]) for i in range(0, len(w) - 2, 3)]
                    L.append("- %s: %s" % (nm, "; ".join(cells)))
open(OUT, "w").write("\n".join(L) + "\n")
print("\n".join(L)[:12000])
