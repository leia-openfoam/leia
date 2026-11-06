"""The decision 2/3 matrix of the plane loop (iso/dec/<cell>/<bff>-<rule>/): per run the
status, steps, E_shape, L1(n), L1(alphaf), min alpha, max excess, interface cells (first/max),
the mass deviation, the wall clock, and the counters' sums (bounded faces, recomputed; rule
cells outside, reset). Markdown tables per cell. Claude, 2026-11-07."""
import csv, glob, os, re, sys
D = "/work/scratch/tm83tomy/davof/agents/swept/iso/dec"
RE_T = re.compile(r"t = (\S+),\s+sum\(alpha\*V\) = (\S+),\s+dev = (\S+)%,\s+1-max\(alpha1\) = (\S+) at cell \d+,\s+min\(alpha1\) = (\S+)")
RE_IC = re.compile(r"davof: (\d+) interface cells")
RE_B = re.compile(r"bounded face fractions \(\w+\): faces (\d+), s = 0 on (\d+), partial (\d+), max \|1 - s\| (\S+)(?:, recomputed (\d+) \(no plane (\d+), not bracketed (\d+)\))?")
RE_R = re.compile(r"interface-cell rule \w+: cells with \|m_c\| above the threshold and alpha outside \(tol, 1 - tol\) (\d+), with a plane (\d+)(?:, reset cells (\d+), faces (\d+), max \|dalphaf\| (\S+))?")
def f(x, d=3):
    try: v = float(x)
    except Exception: return "-"
    if v == 0: return "0"
    return ("%%.%dg" % d) % v if 1e-2 <= abs(v) < 1e3 else ("%%.%de" % (d - 1)) % v
order_b = ["none", "scale", "recompute"]; order_r = ["normalOnly", "normalAndAlpha", "normalAndAlphaReset"]
for cell in sorted(glob.glob(D + "/*")):
    print("\n### %s\n" % os.path.basename(cell))
    print("| bounded face fraction | interface-cell rule | status | steps | E_shape | L1(n) | L1(alphaf) | min alpha | overshoot | interface cells first / max | mass dev. % | bounded faces per step (s=0 / partial / recomputed) | rule: outside per step (reset cells) | wall s |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in order_b:
        for r in order_r:
            d = "%s/%s-%s" % (cell, b, r); log = d + "/log.advectorVoF"
            if not os.path.exists(log): print("| %s | %s | missing |" % (b, r)); continue
            t = open(log, errors="replace").read()
            steps = len(re.findall(r"\nTime = ", t))
            tl = RE_T.findall(t); ic = [int(x) for x in RE_IC.findall(t)]
            mina = min(float(x[4]) for x in tl) if tl else float("nan"); exc = max(-float(x[3]) for x in tl) if tl else float("nan")
            mdev = max(abs(float(x[2])) for x in tl) if tl else float("nan")
            ab = re.search(r"ABORT at step (\d+), t = (\S+?):", t)
            status = ("ABORT step %s (t %.2f)" % (ab.group(1), float(ab.group(2).rstrip(":")))) if ab else ("completed" if "^End" in t or "\nEnd\n" in t else "incomplete")
            bs = [m for m in RE_B.findall(t)]; nb = sum(int(m[0]) for m in bs); nz = sum(int(m[1]) for m in bs); npart = sum(int(m[2]) for m in bs); nrec = sum(int(m[4] or 0) for m in bs)
            rs = [m for m in RE_R.findall(t)]; nout = sum(int(m[0]) for m in rs); nrc = sum(int(m[2] or 0) for m in rs)
            es = l1n = l1af = "-"
            csvp = d + "/advectorVoF.csv"
            if os.path.exists(csvp) and not ab:
                rows = list(csv.DictReader(open(csvp)))
                if rows:
                    last = rows[-1]; es = f(last.get("Eshape", "")); l1n = f(last.get("LNormalDiff1", "")); l1af = f(last.get("L1_alphaf", ""))
            wall = ""
            mw = re.findall(r"ExecutionTime = (\S+) s", t); wall = f(mw[-1], 4) if mw else "-"
            ns = max(steps, 1)
            print("| %s | %s | %s | %d | %s | %s | %s | %s | %s | %s / %s | %s | %.1f (%.1f / %.1f / %.1f) | %.1f (%.1f) | %s |" % (
                b, r, status, steps, es, l1n, l1af, f(mina, 2), f(exc, 2), ic[0] if ic else "-", max(ic) if ic else "-", f(mdev, 2),
                nb/ns, nz/ns, npart/ns, nrec/ns, nout/ns, nrc/ns, wall))
