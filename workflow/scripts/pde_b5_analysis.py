#!/usr/bin/env python3
"""pde_b5_analysis.py: the read-out of the medial-axis test (bSkel2DstrainA{03,1}N{128,256,512}, S1, S2).

    pde_b5_analysis.py --new <clone> [--out <dir>]

L1 only. Stdlib only."""
import argparse, glob, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pde_b3_analysis import rows, rate  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "b5")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    oT = {}
    for tag in ("03", "1"):
        v = {}
        for N in (128, 256, 512):
            cs = glob.glob(os.path.join(a.new, "studies", "pdeLevelSet", f"bSkel2DstrainA{tag}N{N}", "markerTest2D_*"))
            if cs and rows(cs[0]):
                v[N] = rows(cs[0])
                r, ton = rate(v[N])
                say(f"a = {'0.3' if tag == '03' else '1'} N {N}: E_Q_L1 T/2 {min(v[N], key=lambda x: abs(x['TIME'] - 0.25))['E_Q_L1']:.3e} "
                    f"T {v[N][-1]['E_Q_L1']:.3e}  E_POS_L1 T {v[N][-1]['E_POS_L1']:.3e}  onset t {ton}  rate {r:.2f}")
        Ns = sorted(v)
        for lab, tt in (("T/2", 0.25), ("T", 0.5)):
            for col in ("E_Q_L1", "E_POS_L1"):
                e = [min(v[N], key=lambda x: abs(x["TIME"] - tt))[col] for N in Ns]
                o = [math.log(e[i] / e[i + 1]) / math.log(Ns[i + 1] / Ns[i]) for i in range(len(Ns) - 1)]
                say(f"   {lab} {col}: orders " + " / ".join(f"{x:.2f}" for x in o))
                if col == "E_Q_L1":
                    oT[(tag, lab)] = o
    s1 = all(x >= 1.5 for lab in ("T/2", "T") for x in oT.get(("03", lab), [0]))
    s2 = all(x03 - x1 >= 0.5 for x03, x1 in zip(oT.get(("03", "T"), []), oT.get(("1", "T"), []))) and bool(oT.get(("1", "T")))
    fal = any(x < 1.2 for lab in ("T/2", "T") for x in oT.get(("03", lab), [0]))
    say(f"S1 [{'PASS' if s1 else 'FAIL'}]  S2 [{'PASS' if s2 else 'FAIL'}]" + ("  FALSIFIED (a = 0.3 below order 1.2)" if fal else ""))
    open(os.path.join(out, "b5_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
