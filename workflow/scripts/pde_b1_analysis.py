#!/usr/bin/env python3
"""pde_b1_analysis.py: the read-out of step B, round 1 (the marker test, PDE plan section 6).

    pde_b1_analysis.py --new <clone> [--old <clone for the identity check>] [--out <dir>]

Reads interfaceMarkers.csv of bMarker2Dtranslation, bMarker2Dstrain and bMarker2DstrainPar, takes the
rows at t = 0, T/2 and T (the row with the time nearest to each), prints E_Q and E_POS (L1, L2; never
L_inf) with the orders between consecutive N, and checks B1 to B8 of the config headers. With --old,
also the identity of a0Identity2Dvortex and aStatic2Dedge between the clones (the libleiaCore change
of 2026-11-18 must be inert). Writes <out>/b1_summary.txt. Stdlib only.
"""
import argparse, csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pde_a2_analysis import tok  # noqa: E402
from pde_a3_analysis import cmp  # noqa: E402

COLS = ("E_Q_L1", "E_Q_L2", "E_POS_L1", "E_POS_L2")


def rows(case):
    f = os.path.join(case, "interfaceMarkers.csv")
    if not (os.path.isfile(f) and os.path.getsize(f)):
        return []
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(f))]


def at(rs, t):
    return min(rs, key=lambda r: abs(r["TIME"] - t)) if rs else None


def orders(v):
    return [math.log(e1 / e2) / math.log(n2 / n1) if e1 > 0 and e2 > 0 else float("nan")
            for (n1, e1), (n2, e2) in zip(v, v[1:])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--old"); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "b1")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda root, st: os.path.join(root, "studies", "pdeLevelSet", st)
    verdict = {}

    def judge(tag, ok):
        verdict[tag] = verdict.get(tag, True) and bool(ok)

    if a.old:
        say("Inertness of the libleiaCore change (new clone against old clone, tolerance 0):")
        n = bad = 0
        for st, files in (("aStatic2Dedge", ("leiaTestVelocityExtension.csv", "leiaTestVelocityExtensionExact.csv")),
                          ("a0Identity2Dvortex", ("leiaLevelSetFoam.csv", "gradPsiError.csv"))):
            for c in sorted(glob.glob(os.path.join(S(a.new, st), "2Dvortex_[0-9]*"))):
                for f in files:
                    v, d = cmp(os.path.join(S(a.old, st), os.path.basename(c), f), os.path.join(c, f))
                    n += 1; bad += v != "PASS"
                    if v != "PASS":
                        say(f"  {st} {os.path.basename(c)} {f}: {v} {' | '.join(d[-2:])[:160]}")
        say(f"  -> {n - bad} of {n} identical")
        judge("INERT", n and not bad)
        say()

    data = {}
    for st in ("bMarker2Dtranslation", "bMarker2Dstrain", "bMarker2DstrainPar"):
        d = {}
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "markerTest2D_[0-9]*"))):
            t = tok(c)
            d[(t["VELOCITY_EXTENSION"], int(t["N_CELLS"]))] = (c, rows(c), float(t["END_TIME"]))
        data[st] = d
        nf = [(k, r["N_FOUND"]) for k, (_, rs, _) in d.items() for r in rs if r["N_FOUND"] != r["N_MARKERS"]]
        empty = [k for k, (_, rs, T) in d.items() if not rs or abs(rs[-1]["TIME"] - T) > 1e-9 * max(T, 1)]
        say(f"{st}: {len(d)} cases; rows with N_FOUND != N_MARKERS: {len(nf)}; cases not at T: {empty}")
        judge("B7" if st.endswith("Par") else "B6", not nf)
        judge("complete", not empty)

    for st in ("bMarker2Dtranslation", "bMarker2Dstrain"):
        say()
        say(f"===== {st}")
        d = data[st]
        for m in ("none", "projectionBlend", "closestPoint"):
            say(f"  == {m}")
            for lab, frac in (("t=0", 0.0), ("T/2", 0.5), ("T", 1.0)):
                for col in COLS:
                    v = []
                    for (mm, N), (_, rs, T) in sorted(d.items(), key=lambda kv: kv[0][1]):
                        r = at(rs, frac * T)
                        if mm == m and r:
                            v.append((N, r[col]))
                    o = orders(v)
                    say(f"    {lab:4s} {col:9s} " + "  ".join(f"N{n}:{e:.3e}" for n, e in v)
                        + "   orders " + " / ".join(f"{x:.2f}" for x in o))
                    if lab == "t=0" and col == "E_Q_L1":
                        judge("B5", o and all(x >= 1.8 for x in o))
                    if lab != "t=0" and col in ("E_Q_L1", "E_POS_L1"):
                        if st == "bMarker2Dtranslation":
                            judge("B2", o and all(x >= 1.5 for x in o))
                        elif col == "E_POS_L1":
                            judge("B2s", o and all(x >= 1.5 for x in o))
                        elif m == "none" and lab == "T":
                            judge("B3", o and all(x <= 0.3 for x in o) and all(e >= 0.05 for _, e in v))
                        elif m != "none":
                            judge("B4", o and all(x >= 1.5 for x in o))
        # B1 (translation): arms equal; B4 ratio (strain)
        for N in sorted({k[1] for k in d}):
            base = d.get(("none", N))
            for m in ("projectionBlend", "closestPoint"):
                arm = d.get((m, N))
                if not base or not arm:
                    continue
                for lab, frac in (("T/2", 0.5), ("T", 1.0)):
                    rb, ra = at(base[1], frac * base[2]), at(arm[1], frac * arm[2])
                    if not rb or not ra:
                        continue
                    if st == "bMarker2Dtranslation":
                        rel = max(abs(ra[c] - rb[c]) / abs(rb[c]) for c in COLS if rb[c] != 0)
                        judge("B1", rel <= 1e-6)
                        say(f"  B1 {m} N {N} {lab}: max rel. diff against none {rel:.2e}")
                    elif lab == "T":
                        ratio = ra["E_Q_L1"] / rb["E_Q_L1"] if rb["E_Q_L1"] > 0 else float("nan")
                        if N == max(k[1] for k in d):
                            judge("B4", ratio <= 0.1)
                        judge("active", ratio < 0.1 or N != max(k[1] for k in d))
                        say(f"  B4 {m} N {N} T: E_Q_L1 arm/none = {ratio:.3e}")

    say()
    say("===== bMarker2DstrainPar against bMarker2Dstrain (np 4 against serial):")
    for (m, N), (c, rs, T) in sorted(data["bMarker2DstrainPar"].items()):
        s = data["bMarker2Dstrain"].get((m, N))
        if not s or not rs or not s[1]:
            say(f"  {m} N {N}: missing"); judge("B8", False); continue
        lim = 0.05 if m == "closestPoint" else 0.01
        worst = 0.0
        for frac in (0.5, 1.0):
            rp, rr = at(rs, frac * T), at(s[1], frac * T)
            for col in ("E_Q_L1", "E_POS_L1"):
                worst = max(worst, abs(rp[col] - rr[col]) / abs(rr[col]))
        judge("B8", worst <= lim)
        say(f"  {m:16s} N {N}: max rel. diff of E_Q_L1, E_POS_L1 at T/2 and T {worst:.2e} (limit {lim})")

    say()
    say("VERDICTS: " + ", ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in sorted(verdict.items())))
    open(os.path.join(out, "b1_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
