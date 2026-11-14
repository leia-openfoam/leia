#!/usr/bin/env python3
"""pde_a4_curvstencil.py: the read-out of a4CurvStencil2D (PDE level-set plan, step A0.4).

    pde_a4_curvstencil.py [--study studies/pdeLevelSet/a4CurvStencil2D]

Per gradient scheme: (1) the depth -- the smallest k with N_CHANGED_CELL = 0 (and FACE) at every N;
(2) E_KAPPA_CELL_L1 and E_KAPPA_FACE_L1 per (type, k) over N with the observed orders (L1, L2 only).
Writes <study>/a4_summary.txt. Stdlib only.
"""
import argparse, csv, glob, json, math, os, sys


def orders(v):
    return " / ".join(f"{math.log(a / b) / math.log(n2 / n1):.2f}" if a > 0 and b > 0 else "nan"
                      for (n1, a), (n2, b) in zip(v, v[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="studies/pdeLevelSet/a4CurvStencil2D")
    a = ap.parse_args()
    rows = []
    for case in sorted(glob.glob(os.path.join(a.study, "*_[0-9]*"))):
        N = int(json.load(open(os.path.join(case, "case_params.json")))["tokens"]["N_CELLS"])
        f = os.path.join(case, "leiaTestCurvatureStencil.csv")
        if os.path.isfile(f) and os.path.getsize(f):
            for r in csv.DictReader(open(f)):
                r["N"] = N; rows.append(r)
    L = []
    say = lambda s="": (L.append(s), print(s))
    Ns = sorted({r["N"] for r in rows})
    say(f"{len(rows)} rows, N = {Ns}")
    for sch in sorted({r["GRAD_SCHEME"] for r in rows}):
        rs = [r for r in rows if r["GRAD_SCHEME"] == sch]
        say(f"== {sch}")
        for typ in ("random", "step"):
            for col, lab in (("N_CHANGED_CELL", "cell"), ("N_CHANGED_FACE", "face")):
                ks = sorted({int(r["K"]) for r in rs if r["TYPE"] == typ})
                clean = [k for k in ks if all(int(r[col]) == 0 for r in rs if r["TYPE"] == typ and int(r["K"]) == k)]
                depth = min([k for k in ks if all(kk in clean for kk in ks if kk >= k)], default=None)
                say(f"   {typ:6s} {lab}: unchanged for k >= {depth}  (changed counts at the finest N: "
                    + ", ".join(f"k{int(r['K'])}:{r[col]}" for r in rs if r["TYPE"] == typ and r["N"] == Ns[-1]) + ")")
        for col in ("E_KAPPA_CELL_L1", "E_KAPPA_CELL_L2", "E_KAPPA_FACE_L1", "E_KAPPA_FACE_L2"):
            say(f"   {col}")
            for typ, k in [("none", 0)] + [(t, k) for t in ("random", "step") for k in range(1, 7)]:
                v = sorted((r["N"], float(r[col])) for r in rs if r["TYPE"] == typ and int(r["K"]) == k)
                if v:
                    say(f"     {typ:6s} k {k}: " + "  ".join(f"N{n}:{e:.2e}" for n, e in v) + "   orders " + orders(v))
    open(os.path.join(a.study, "a4_summary.txt"), "w").write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
