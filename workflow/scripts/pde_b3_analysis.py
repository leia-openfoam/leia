#!/usr/bin/env python3
"""pde_b3_analysis.py: the read-out of the band-width sweep (bBand2Dstrain*, H1 and H2) and of the
raw-velocity control (bRaw2Dstrain, W1 and W2).

    pde_b3_analysis.py --new <clone of the sweep> [--raw <clone of bRaw2Dstrain>] [--out <dir>]

Growth rate per unit time: from the onset (the first row with E_Q_L1 > 2 x its t = 0 value) to the first
row with E_Q_L1 > 5 x the onset value (or the last row), rate = ln(e2/e1)/(t2 - t1). L1 only (never
L_inf). Writes <out>/b3_summary.txt. Stdlib only.
"""
import argparse, csv, glob, json, math, os, sys


def rows(case):
    f = os.path.join(case, "interfaceMarkers.csv")
    if not (os.path.isfile(f) and os.path.getsize(f)):
        return []
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(f))]


def rate(rs, col="E_Q_L1"):
    e0 = rs[0][col]
    i = next((k for k, r in enumerate(rs) if r[col] > 2 * e0), None)
    if i is None:
        return float("nan"), None
    j = next((k for k in range(i, len(rs)) if rs[k][col] > 5 * rs[i][col]), len(rs) - 1)
    if j <= i:
        return float("nan"), rs[i]["TIME"]
    return math.log(rs[j][col] / rs[i][col]) / (rs[j]["TIME"] - rs[i]["TIME"]), rs[i]["TIME"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--raw"); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "b3")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda root, st: os.path.join(root, "studies", "pdeLevelSet", st)
    if a.raw:
        say("bRaw2Dstrain (no extension; E_QRAW = q_h - |F^-T grad psi0|):")
        v = {}
        for c in sorted(glob.glob(os.path.join(S(a.raw, "bRaw2Dstrain"), "markerTest2D_*"))):
            N = int(json.load(open(os.path.join(c, "case_params.json")))["tokens"]["N_CELLS"])
            rs = rows(c)
            if rs:
                v[N] = rs
                say(f"  N {N}: t=0 E_QRAW_L1 {rs[0]['E_QRAW_L1']:.3e} (E_Q_L1 {rs[0]['E_Q_L1']:.3e}), "
                    f"T/2 {min(rs, key=lambda r: abs(r['TIME'] - 0.25))['E_QRAW_L1']:.3e}, T {rs[-1]['E_QRAW_L1']:.3e}")
        Ns = sorted(v)
        for lab, frac in (("T/2", 0.5), ("T", 1.0)):
            e = [min(v[N], key=lambda r: abs(r["TIME"] - frac * v[N][-1]["TIME"]))["E_QRAW_L1"] for N in Ns]
            o = [math.log(e1 / e2) / math.log(n2 / n1) for (n1, e1), (n2, e2) in zip(zip(Ns, e), list(zip(Ns, e))[1:])]
            say(f"  W2 {lab}: orders " + " / ".join(f"{x:.2f}" for x in o) + f"  [{'PASS' if o and all(x >= 1.5 for x in o) else 'FAIL'}]")
        w1 = all(abs(v[N][0]["E_QRAW_L1"] - v[N][0]["E_Q_L1"]) <= 1e-12 * v[N][0]["E_Q_L1"] for N in Ns)
        say(f"  W1 t = 0 row: E_QRAW = E_Q  [{'PASS' if w1 else 'FAIL'}]")
        say()
    res = {}
    for st in ("bBand2DstrainPB2", "bBand2DstrainPB4", "bBand2DstrainPB8", "bBand2DstrainCP"):
        for c in sorted(glob.glob(os.path.join(S(a.new, st), "markerTest2D_*"))):
            t = json.load(open(os.path.join(c, "case_params.json")))["tokens"]
            rs = rows(c)
            if not rs:
                say(f"{st} {os.path.basename(c)}: no CSV"); continue
            key = (f"projectionBlend w {t['PB_W1']}/{t['PB_W2']}" if t["VELOCITY_EXTENSION"] == "projectionBlend"
                   else f"closestPoint nLayers {t['N_LAYERS']}")
            r, ton = rate(rs)
            res[key] = r
            say(f"{key:30s} rate {r:7.2f} per unit time (onset t {ton})  E_Q_L1(T) {rs[-1]['E_Q_L1']:.3e}  "
                f"E_POS_L1(T) {rs[-1]['E_POS_L1']:.3e}  ({len(rs) - 1} steps)")
    say()
    g = lambda k: res.get(k, float("nan"))
    p2, p4, p8 = g("projectionBlend w 2/4"), g("projectionBlend w 4/8"), g("projectionBlend w 8/16")
    if p4 == p4 and p4 > 0:
        r8, r2 = p8 / p4, p2 / p4
        h1 = 1.6 <= r8 <= 2.5 and 0.4 <= r2 <= 0.65
        fals = max(p2, p4, p8) <= 1.2 * min(p2, p4, p8)
        say(f"H1 rate(8,16)/rate(4,8) = {r8:.2f}, rate(2,4)/rate(4,8) = {r2:.2f}  [{'PASS' if h1 else 'FAIL'}]"
            + ("  FALSIFIED (all three within 20 %)" if fals else ""))
    c2, c3, c5 = g("closestPoint nLayers 2"), g("closestPoint nLayers 3"), g("closestPoint nLayers 5")
    h2 = c5 > 1.15 * c3 and c3 > 1.15 * c2
    say(f"H2 closestPoint rates nLayers 2 / 3 / 5: {c2:.2f} / {c3:.2f} / {c5:.2f}  [{'PASS' if h2 else 'FAIL'}]")
    open(os.path.join(out, "b3_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
