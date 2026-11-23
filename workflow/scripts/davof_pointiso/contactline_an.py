#!/usr/bin/env python3
"""Tables and verdict data of the 2D contact line study (docs/davof/REPORT-swept.md "Contact line advection").
Reference of the tests and the measures: M. Fricke, T. Maric, D. Bothe, Contact line advection using the geometrical
Volume-of-Fluid method, Journal of Computational Physics 407 (2020) 109221, doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.
Reads studies/davof/contactLine_N<N> (the rows plicRDF + isoAdvector and pointIso), the raw contact points and angles of every step
(postProcessing/contactLineError, the function object contactLineError, the same code for every row) and the references
(Eqs. 12, 14, 38 of the paper; case D the limit c1 -> 0 of Eq. 14 with v0 = 0).
    python3 workflow/scripts/davof_pointiso/contactline_an.py [--out FILE] [--series DIR]
"""
import argparse, csv, glob, json, math, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "davof_kiss"))
import kiss_an as ka
REPO = ka.REPO
R0 = 0.2
TH0 = math.pi/3
X0L = 0.4 - math.sqrt(0.2**2 - 0.1**2)
X0R = 0.4 + math.sqrt(0.2**2 - 0.1**2)
ARMS = {"clB": dict(v0=-0.2, c1=0.1, c2=-2.0, tau=0.0), "clC": dict(v0=-0.2, c1=0.1, c2=-2.0, tau=0.2), "clD": dict(v0=0.0, c1=0.0, c2=-2.0, tau=0.0)}
NS = (128, 256, 512, 1024)
COS = (0.2, 0.5, 0.95)
ROWS = ("plicRDF", "pointIso", "pointIso-davof")
def ref(arm, t, side):
    """x_cl and theta of the exact solution; side 0 left ('-'), 1 right ('+')"""
    a = ARMS[arm]
    sgn = -1.0 if side == 0 else 1.0
    s = a["tau"]*math.sin(math.pi*t/a["tau"])/math.pi if a["tau"] > 0 else t
    x0 = X0L if side == 0 else X0R
    cot = 1.0/math.tan(TH0)
    if a["c1"] != 0:
        e = math.exp(a["c1"]*s)
        x = x0*e + a["v0"]/a["c1"]*(e - 1)
        th = math.pi/2 + math.atan(-cot*e*e + sgn*a["c2"]*(e*e - 1)/(2*a["c1"]))
    else:
        x = x0 + a["v0"]*s
        th = math.pi/2 + math.atan(-cot + sgn*a["c2"]*s)
    return x, th
def read_dat(c):
    fs = sorted(glob.glob(os.path.join(c, "postProcessing", "contactLineError", "*", "contactLineError.dat")))
    if not fs:
        return []
    rows = []
    for l in open(fs[-1]):
        if l.startswith("#") or not l.strip():
            continue
        rows.append([float(x) for x in l.split()])
    return rows
def case_status(c):
    p = os.path.join(c, "log.davofAdvectorFoam")
    if not os.path.isfile(p):
        return "no log", ""
    log = open(p, errors="replace").read()
    if "ABORT at step" in log:
        return "aborted", re.search(r"ABORT at step[^\n]*", log).group(0)[:120]
    if "[stack trace]" in log:
        return "crashed", ""
    if "\nEnd" in log:
        return "completed", ""
    return "incomplete", ""
def devmax(c):
    p = os.path.join(c, "log.davofAdvectorFoam")
    if not os.path.isfile(p):
        return None
    d = [abs(float(x)) for x in re.findall(r"dev = (\S+)%", open(p, errors="replace").read())]
    return max(d)/100 if d else None
def planes_boundary(c):
    p = os.path.join(c, "log.davofAdvectorFoam")
    out = dict(bndDisp=None, bndDist=None, bndCells=None)
    if not os.path.isfile(p):
        return out
    log = open(p, errors="replace").read()
    m = re.findall(r"cells with boundary points (\d+), with more than three non-collinear boundary points (\d+), boundary points: max \|dx\| = (\S+) m, max distance to the fitted plane/V\^\(1/3\) = (\S+)", log)
    if m:
        out = dict(bndCells=max(int(x[0]) for x in m), bndDisp=max(float(x[2]) for x in m), bndDist=max(float(x[3]) for x in m))
    return out
def load():
    out = []
    for N in NS:
        for c in sorted(glob.glob(os.path.join(REPO, "studies", "davof", "contactLine_N%d" % N, "*_[0-9][0-9][0-9][0-9][0-9]"))
                        + glob.glob(os.path.join(REPO, "studies", "davof", "contactLineD_N%d" % N, "*_[0-9][0-9][0-9][0-9][0-9]"))):
            r = ka.read_case(c)
            if not r:
                continue
            tok = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
            r["arm"] = tok.get("TPF_VELOCITY_ARM"); r["row"] = tok.get("TPF_ROW"); r["N"] = N
            r["coN"] = min(COS, key=lambda x: abs(x - float(tok["TARGET_ADVECTIVE_CO"])))
            r["status2"], r["reason"] = case_status(c)
            r["dev"] = devmax(c)
            r.update(planes_boundary(c))
            rows = read_dat(c)
            r["dat"] = rows
            r["dir"] = c
            r["dx"] = 1.0/N
            for side in (0, 1):
                errs_th, errs_x, miss, errs_th1 = [], [], 0, []
                for row in rows:
                    t = row[0]
                    n_c, x, th = row[1 + 4*side], row[2 + 4*side], row[3 + 4*side]
                    if n_c == 0:
                        miss += 1
                        continue
                    xr, thr = ref(r["arm"], t, side)
                    errs_th.append(abs(th - thr)); errs_x.append(abs(x - xr)/R0)
                    if n_c == 1:
                        errs_th1.append(abs(th - thr))
                r["Eth%d" % side] = max(errs_th) if errs_th else None
                r["Eth1c%d" % side] = max(errs_th1) if errs_th1 else None
                r["Ecl%d" % side] = max(errs_x) if errs_x else None
                r["miss%d" % side] = miss
                r["multi%d" % side] = sum(1 for row in rows if row[1 + 4*side] > 1)
            r["E1"] = rows[-1][10] if rows else None
            r["minA"] = min((row[11] for row in rows), default=None)
            r["maxA"] = max((row[12] for row in rows), default=None)
            out.append(r)
    return out
def f(v, nd=3):
    return "-" if v is None else "%.*e" % (nd, v)
def order(pts):
    pts = [(x, y) for x, y in pts if y is not None and y > 0]
    if len(pts) < 2:
        return None
    xs = [math.log(x) for x, y in pts]; ys = [math.log(y) for x, y in pts]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    return sum((x - mx)*(y - my) for x, y in zip(xs, ys))/sum((x - mx)**2 for x in xs)
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(REPO, "studies", "davof", "contactline_tables.md"))
    ap.add_argument("--series", default=os.path.join(REPO, "studies", "davof", "contactline_series"))
    a = ap.parse_args()
    cs = load()
    by = {(r["arm"], r["row"], r["N"], r["coN"]): r for r in cs}
    L = []
    L.append("Contact line study, rows: plicRDF = plicRDF + isoAdvector (OpenFOAM), pointIso = pointIsoAdvection + plicRDF. Ranks: 1 for N 128 and 256, 8 for N 512 and 1024 (the same for both rows).")
    L.append("Errors: E_theta^inf in rad (deg in the second column), E_cl^inf = max |x_num - x_ref|/R0, per contact point (left L, right R), over the steps with a detected contact cell; 'miss' = steps without one; 'multi' = steps with more than one candidate cell on a side.")
    for arm, title in (("clB", "Case B (Eq. 36, steady)"), ("clC", "Case C (Eq. 37, periodic, T = 2 tau)"), ("clD", "Case D (no slip, rolling)")):
        for co in COS:
            L.append("\n### %s, OpenFOAM Co target %s\n" % (title, co))
            L.append("| row | N | ranks | OF Co | steps | status | E_theta L (rad / deg) | E_theta R (rad / deg) | E_theta L, one candidate (rad) | E_theta R, one candidate (rad) | E_cl L | E_cl R | miss L / R | multi L / R | E_1 | liquid volume |dev| | min alpha | max alpha - 1 | wall (s) |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for row in ROWS:
                for N in NS:
                    r = by.get((arm, row, N, co))
                    if r is None:
                        continue
                    L.append("| %s | %d | %s | %s | %s | %s | %s / %s | %s / %s | %s | %s | %s | %s | %d / %d | %d / %d | %s | %s | %s | %s | %s |" % (
                        row, N, r.get("np"), "%.3f" % r["ofCo"] if r.get("ofCo") else "-", r.get("steps"), r["status2"],
                        f(r["Eth0"]), "-" if r["Eth0"] is None else "%.3f" % math.degrees(r["Eth0"]),
                        f(r["Eth1"]), "-" if r["Eth1"] is None else "%.3f" % math.degrees(r["Eth1"]),
                        f(r["Eth1c0"]), f(r["Eth1c1"]), f(r["Ecl0"]), f(r["Ecl1"]), r["miss0"], r["miss1"], r["multi0"], r["multi1"],
                        f(r["E1"]) if arm == "clC" else "-", f(r["dev"]), f(r["minA"]), f(r["maxA"]), "%.0f" % r["wall"] if r.get("wall") else "-"))
            L.append("\nOrders over dx/R0 (least squares over the completed rungs):\n")
            L.append("| row | E_theta L | E_theta R | E_theta L (one candidate) | E_theta R (one candidate) | E_cl L | E_cl R | E_1 |")
            L.append("|---|---|---|---|---|---|---|---|")
            for row in ROWS:
                rs = [by[(arm, row, N, co)] for N in NS if (arm, row, N, co) in by and by[(arm, row, N, co)]["status2"] == "completed"]
                o = lambda k: order([(r["dx"]/R0, r[k]) for r in rs])
                L.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (row, *["-" if o(k) is None else "%.2f" % o(k) for k in ("Eth0", "Eth1", "Eth1c0", "Eth1c1", "Ecl0", "Ecl1")], "-" if (arm != "clC" or o("E1") is None) else "%.2f" % o("E1")))
    # Gw1: step 0
    L.append("\n### Gw1: the initial reconstruction (step 0), N 128, Co 0.5\n")
    L.append("| case | row | side | x_num | x_ref | theta_num (deg) | theta_ref (deg) | candidates |")
    L.append("|---|---|---|---|---|---|---|---|")
    for arm in ARMS:
        for row in ROWS:
            r = by.get((arm, row, 128, 0.5))
            if r and r["dat"]:
                d = r["dat"][0]
                for side in (0, 1):
                    xr, thr = ref(arm, d[0], side)
                    L.append("| %s | %s | %s | %.10f | %.10f | %s | %.6f | %d |" % (arm, row, "LR"[side], d[2 + 4*side], xr, "%.6f" % math.degrees(d[3 + 4*side]) if d[1 + 4*side] else "-", math.degrees(thr), d[1 + 4*side]))
    # Gw2: boundary points of the pointIso rows, case D
    L.append("\n### Gw2: boundary points of pointIso (all cases; max over the steps)\n")
    L.append("| case | N | Co | max displacement of a boundary point (m) | max distance of a boundary point to the fitted plane / V^(1/3) | max cells with boundary points |")
    L.append("|---|---|---|---|---|---|")
    for arm in ARMS:
        for N in NS:
            for co in COS:
                r = by.get((arm, "pointIso", N, co))
                if r and r.get("bndDisp") is not None:
                    L.append("| %s | %d | %s | %s | %s | %s |" % (arm, N, co, f(r["bndDisp"]), f(r["bndDist"]), r["bndCells"]))
    # series for N 256
    os.makedirs(a.series, exist_ok=True)
    for r in cs:
        if r["N"] == 256 and r["dat"]:
            fn = os.path.join(a.series, "series_N256_%s_Co%s_%s.csv" % (r["arm"], r["coN"], r["row"]))
            with open(fn, "w") as fo:
                fo.write("t,nCandL,xL_num,xL_ref,thetaL_num_deg,thetaL_ref_deg,nCandR,xR_num,xR_ref,thetaR_num_deg,thetaR_ref_deg,volume,L1toInitial,minAlpha,maxAlphaMinus1\n")
                for d in r["dat"]:
                    xl, tl = ref(r["arm"], d[0], 0); xr, tr = ref(r["arm"], d[0], 1)
                    fo.write("%.17g,%d,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n" % (
                        d[0], d[1], d[2], xl, math.degrees(d[3]), math.degrees(tl), d[5], d[6], xr, math.degrees(d[7]), math.degrees(tr), d[9], d[10], d[11], d[12]))
    # verdict data
    L.append("\n### Data for the pre-registered statements\n")
    L.append("- PW1 (case D, E_cl^inf ratio row (1) / row (2) = plicRDF / pointIso, and row (1) / row (3) = plicRDF / pointIso-davof, per point; at least 2 required): ")
    for N in NS:
        for co in COS:
            a3 = by.get(("clD", "pointIso-davof", N, co)); a1 = by.get(("clD", "plicRDF", N, co))
            if a1 and a3 and a1.get("Ecl0") is not None and a3.get("Ecl0") is not None:
                L.append("  - row (3) N %d Co %s: L %.3g, R %.3g (E_cl: row1 %s / %s, row3 %s / %s)" % (N, co, a1["Ecl0"]/a3["Ecl0"] if a3["Ecl0"] > 0 else float("inf"), a1["Ecl1"]/a3["Ecl1"] if a3["Ecl1"] > 0 else float("inf"), f(a1["Ecl0"]), f(a1["Ecl1"]), f(a3["Ecl0"]), f(a3["Ecl1"])))
            a1, a2 = by.get(("clD", "plicRDF", N, co)), by.get(("clD", "pointIso", N, co))
            if a1 and a2 and a1.get("Ecl0") is not None and a2.get("Ecl0") is not None:
                rl = a1["Ecl0"]/a2["Ecl0"] if a2["Ecl0"] > 0 else float("inf")
                rr = a1["Ecl1"]/a2["Ecl1"] if a2["Ecl1"] > 0 else float("inf")
                L.append("  - N %d Co %s: L %.3g, R %.3g (E_cl: row1 %s / %s, row2 %s / %s)" % (N, co, rl, rr, f(a1["Ecl0"]), f(a1["Ecl1"]), f(a2["Ecl0"]), f(a2["Ecl1"])))
    L.append("- PW2 (cases B, C: orders of E_cl^inf of pointIso >= 1 and E_cl row (2) <= row (1)): see the order tables; pointwise E_cl row2 <= row1:")
    for arm in ("clB", "clC"):
        for N in NS:
            for co in COS:
                a1, a2 = by.get((arm, "plicRDF", N, co)), by.get((arm, "pointIso", N, co))
                if a1 and a2 and a1.get("Ecl0") is not None and a2.get("Ecl0") is not None:
                    L.append("  - %s N %d Co %s: L %s <= %s %s, R %s <= %s %s" % (arm, N, co, f(a2["Ecl0"]), f(a1["Ecl0"]), a2["Ecl0"] <= a1["Ecl0"], f(a2["Ecl1"]), f(a1["Ecl1"]), a2["Ecl1"] <= a1["Ecl1"]))
    L.append("- PW3 (E_theta^inf row (3) <= row (1)): row (3) is not wired (see the report); E_theta row (2) <= row (1):")
    for arm in ARMS:
        for N in NS:
            for co in COS:
                a1, a2 = by.get((arm, "plicRDF", N, co)), by.get((arm, "pointIso", N, co))
                if a1 and a2 and a1.get("Eth0") is not None and a2.get("Eth0") is not None:
                    L.append("  - %s N %d Co %s: L %s vs %s (%s), R %s vs %s (%s)" % (arm, N, co, f(a2["Eth0"]), f(a1["Eth0"]), a2["Eth0"] <= a1["Eth0"], f(a2["Eth1"]), f(a1["Eth1"]), a2["Eth1"] <= a1["Eth1"]))
    L.append("- PW4 (liquid volume |dev| <= 1e-10, alpha within 1e-8 of [0, 1]), counterexamples:")
    bad = [r for r in cs if (r["dev"] is None or r["dev"] > 1e-10 or (r["minA"] is not None and r["minA"] < -1e-8) or (r["maxA"] is not None and r["maxA"] > 1e-8))]
    for r in bad:
        L.append("  - %s %s N %d Co %s: dev %s, min alpha %s, max alpha - 1 %s" % (r["arm"], r["row"], r["N"], r["coN"], f(r["dev"]), f(r["minA"]), f(r["maxA"])))
    if not bad:
        L.append("  - none")
    L.append("- PW5 (wall clock row (2) / row (1), at most 2):")
    for arm in ARMS:
        for N in NS:
            for co in COS:
                a1, a2 = by.get((arm, "plicRDF", N, co)), by.get((arm, "pointIso", N, co))
                if a1 and a2 and a1.get("wall") and a2.get("wall"):
                    L.append("  - %s N %d Co %s: %.2f (%s s / %s s)" % (arm, N, co, a2["wall"]/a1["wall"], "%.0f" % a2["wall"], "%.0f" % a1["wall"]))
    open(a.out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))
if __name__ == "__main__":
    main()
