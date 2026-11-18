#!/usr/bin/env python3
"""pde_b4_analysis.py: the read-out of the band-edge test (bEdge2Dstrain*, E1 and E2).

    pde_b4_analysis.py --new <clone> [--out <dir>]

E_Q_L1 and E_POS_L1 at T/2 and T with the orders over N 128-256-512 for the blend of fixed physical
width, and the domain-wide arm against the (4, 8) arm of bBand2DstrainPB4. L1 only. Stdlib only.
"""
import argparse, glob, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pde_b3_analysis import rows, rate  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    out = a.out or os.path.join(a.new, "studies", "pdeLevelSet", "b4")
    os.makedirs(out, exist_ok=True)
    L = []
    say = lambda s="": (L.append(s), print(s))
    S = lambda st: glob.glob(os.path.join(a.new, "studies", "pdeLevelSet", st, "markerTest2D_*"))
    v = {}
    for N in (128, 256, 512):
        cs = S(f"bEdge2DstrainN{N}")
        if cs and rows(cs[0]):
            v[N] = rows(cs[0])
            r, ton = rate(v[N])
            say(f"fixed physical blend N {N}: rate {r:.2f}, onset t {ton}, E_Q_L1 T/2 "
                f"{min(v[N], key=lambda x: abs(x['TIME'] - 0.25))['E_Q_L1']:.3e}, T {v[N][-1]['E_Q_L1']:.3e}, "
                f"E_POS_L1 T {v[N][-1]['E_POS_L1']:.3e}")
        else:
            say(f"fixed physical blend N {N}: no CSV")
    Ns = sorted(v)
    ok = len(Ns) == 3
    for lab, tt in (("T/2", 0.25), ("T", 0.5)):
        for col in ("E_Q_L1", "E_POS_L1"):
            e = [min(v[N], key=lambda x: abs(x["TIME"] - tt))[col] for N in Ns]
            o = [math.log(e[i] / e[i + 1]) / math.log(Ns[i + 1] / Ns[i]) for i in range(len(Ns) - 1)]
            if col == "E_Q_L1":
                ok = ok and all(x >= 1.5 for x in o)
            say(f"  {lab} {col}: orders " + " / ".join(f"{x:.2f}" for x in o))
    say(f"E1 [{'PASS' if ok else 'FAIL'}]")
    cs, ref = S("bEdge2DstrainWide"), S("bBand2DstrainPB4")
    if cs and rows(cs[0]):
        w = rows(cs[0]); r, ton = rate(w)
        eref = rows(ref[0])[-1]["E_Q_L1"] if ref and rows(ref[0]) else 2.194e-2
        e2 = w[-1]["E_Q_L1"] <= 0.1 * eref and (ton is None or ton >= 0.3)
        say(f"domain-wide arm N 256: E_Q_L1(T) {w[-1]['E_Q_L1']:.3e} (ratio to (4, 8): {w[-1]['E_Q_L1'] / eref:.3f}), "
            f"onset t {ton}, rate {r:.2f}  E2 [{'PASS' if e2 else 'FAIL'}]")
    open(os.path.join(out, "b4_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
