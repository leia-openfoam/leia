#!/usr/bin/env python3
"""pde_a_static_analysis.py: step A of the PDE level-set plan: each extension against u(Pi(x)).

    pde_a_static_analysis.py [--study studies/pdeLevelSet/aStatic2D]

Reads leiaTestVelocityExtensionExact.csv of every case, prints per model the errors and the
observed orders (L1 and L2 only) in the shells |d| <= 2h and 4h, the band |d| <= W and the face
flux, against the raw velocity (what no extension scores), and writes <study>/aStatic_table.csv.
Stdlib only.
"""
import argparse, csv, glob, json, math, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="studies/pdeLevelSet/aStatic2D")
    a = ap.parse_args()
    rows = []
    for case in sorted(glob.glob(os.path.join(a.study, "*_[0-9]*"))):
        f = os.path.join(case, "leiaTestVelocityExtensionExact.csv")
        t = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
        r = list(csv.DictReader(open(f)))[0] if os.path.isfile(f) and os.path.getsize(f) else {}
        rows.append(dict(model=str(t["VELOCITY_EXTENSION"]), N=int(t["N_CELLS"]), case=os.path.basename(case), **r))
    out = os.path.join(a.study, "aStatic_table.csv")
    keys = ["model", "N", "case"] + [k for k in (rows[0].keys() if rows else []) if k not in ("model", "N", "case")]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    cols = ("E_EXT_L1_2H", "E_EXT_L2_2H", "E_EXT_L1_4H", "E_EXT_L1_W", "E_EXT_L2_W", "E_FLUX_L1_W")
    print(f"{len(rows)} cases; columns: errors at the finest rung, then orders between rungs")
    for m in sorted({r["model"] for r in rows}):
        ms = sorted([r for r in rows if r["model"] == m], key=lambda r: r["N"])
        missing = [r["N"] for r in ms if "E_EXT_L1_W" not in r]
        print(f"== {m}" + (f"   (no exact CSV at N {missing})" if missing else ""))
        for c in cols + ("E_RAW_L1_W",):
            v = [(r["N"], float(r[c])) for r in ms if c in r and r[c] not in ("", None)]
            od = []
            for (n1, e1), (n2, e2) in zip(v, v[1:]):
                od.append(f"{math.log(e1 / e2) / math.log(n2 / n1):.2f}" if e1 > 0 and e2 > 0 else "nan")
            print(f"   {c:14s} " + "  ".join(f"N{n}:{e:.2e}" for n, e in v) + "   orders " + " / ".join(od))
    return 0


if __name__ == "__main__":
    sys.exit(main())
