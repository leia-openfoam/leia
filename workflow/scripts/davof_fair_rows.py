#!/usr/bin/env python3
"""The fair rows of the oblique comparison (2026-11-03; predictions O-b, O-c,
A-b, A-c, V-1 in the header of config/davof/sphereRotation3D_hexFineLineDiag
RingSwept.yaml): per row a diagonal and an x-line study, both N, read as
davof_static_faces.measure() reads them (the app staticFaceFractions over the
written times, the per-step .dat, the converter's csv), and

    Delta = s_Ec,D - s_Ec,X              (the oblique excess of Ec's late slope)
    phi   = (s_S,D - s_S,X)/(s_L,D - s_L,X)   (rows with a loop alphaf)

with the fixed E_shape at the end (the app's Evol at T: sum |alpha -
alpha_exact| V against the exact sphere, which the fixed volumeFractionError
writes) and, where the .dat is of the fixed object (header "exactState
exactSphere"), V-1: max over the written times of |Eshape - Evol| / Evol.

    python3 workflow/scripts/davof_fair_rows.py --row LABEL DIAG_STUDY X_STUDY [--row ...] [--ref LABEL] [--csv out.csv]
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import davof_static_faces as dsf  # noqa: E402


def cases(study, rowtok=None):
    """{N: (case, tokens)} of a study, only the cases of TPF_ROW rowtok when
    given (dsf.cases keys by N alone)"""
    if rowtok is None:
        return dsf.cases(study)
    import glob
    import json
    out = {}
    for c in sorted(glob.glob(os.path.join(study, "*_0*"))):
        p = os.path.join(c, "case_params.json")
        if "." in os.path.basename(c) or not os.path.isfile(p):
            continue
        tok = json.load(open(p))["tokens"]
        if tok.get("TPF_ROW") == rowtok:
            out[int(float(tok["N_CELLS"]))] = (c, tok)
    return out


def v1(case):
    """max over the written times of |Eshape(.dat) - Evol(app)| / Evol, and
    whether the .dat is of the fixed object"""
    dat = os.path.join(case, "postProcessing", "volumeFractionError", "0", "volumeFractionError.dat")
    fixed = "exactState exactSphere" in open(dat).readline()
    vn, Vf = dsf.read_table(dat)
    sn, S = dsf.read_table(os.path.join(case, "postProcessing", "staticFaceFractions", "staticFaceFractions.dat"))
    tv, Es = Vf[:, 0], Vf[:, vn.index("Eshape")]
    ts, Ev = S[:, sn.index("Time")], S[:, sn.index("Evol")]
    rel, n = 0.0, 0
    for t, e in zip(ts, Ev):
        j = np.where(np.abs(tv - t) < 1e-9)[0]
        if len(j) and e > 0:
            rel = max(rel, abs(Es[j[0]] - e) / e)
            n += 1
    return fixed, rel, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", nargs=3, action="append", metavar=("LABEL", "DIAG", "XLINE"), required=True)
    ap.add_argument("--ref", default="", help="the row whose Delta is the reference of O-b / O-c")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    out, delta = [], {}
    for label, sd, sx in a.row:
        # LABEL=TPF_ROW selects that row's cases of a study holding several
        # rows (2026-11-03, the paraboloid studies)
        rowtok = label.split("=", 1)[1] if "=" in label else None
        label = label.split("=", 1)[0]
        D = {N: dsf.measure(c, tok) for N, (c, tok) in cases(sd, rowtok).items()}
        X = {N: dsf.measure(c, tok) for N, (c, tok) in cases(sx, rowtok).items()}
        for N in sorted(D):
            if N not in X:
                continue
            d, x = D[N], X[N]
            r = dict(row=label, N=N)
            for arm, m in (("D", d), ("X", x)):
                for k in ("s_L1loop", "s_L1static", "s_Ec", "L1loop_end", "L1static_end", "Ec_end", "Eshape_end", "L1x_end"):
                    r[f"{k}_{arm}"] = m.get(k)
                # the fixed E_shape at the end: the app's Evol at T
                sn, S = dsf.read_table(os.path.join(m["case"], "postProcessing", "staticFaceFractions", "staticFaceFractions.dat"))
                r[f"EshapeFixed_end_{arm}"] = float(S[-1, sn.index("Evol")])
                r[f"tEnd_{arm}"] = float(S[-1, 0])
                fixed, rel, n = v1(m["case"])
                r[f"V1_{arm}"] = rel if fixed else float("nan")
                r[f"V1n_{arm}"] = n if fixed else 0
            r["Delta"] = r["s_Ec_D"] - r["s_Ec_X"]
            r["ratioEc"] = r["s_Ec_D"] / r["s_Ec_X"] if r["s_Ec_X"] else float("nan")
            den = (r["s_L1loop_D"] or 0) - (r["s_L1loop_X"] or 0)
            r["phi"] = ((r["s_L1static_D"] - r["s_L1static_X"]) / den
                        if r["L1loop_end_D"] is not None and r["L1loop_end_D"] >= 0 and den else float("nan"))
            delta[(label, N)] = r["Delta"]
            out.append(r)

    keys = ["row", "N"] + [f"{k}_{arm}" for arm in ("D", "X") for k in
                           ("tEnd", "L1loop_end", "s_L1loop", "L1static_end", "s_L1static", "Ec_end", "s_Ec",
                            "EshapeFixed_end", "Eshape_end", "L1x_end", "V1")] + ["Delta", "ratioEc", "phi"]
    fmt = lambda v: ("%.4e" % v) if isinstance(v, float) else str(v)
    print("\t".join(keys))
    for r in out:
        print("\t".join(fmt(r.get(k)) for k in keys))
    if a.ref:
        print("== the excess against the reference row %s" % a.ref)
        for r in out:
            ref = delta.get((a.ref, r["N"]))
            if ref:
                print("  %s N %d: Delta %.4e = %.3f of the reference's %.4e; s_Ec,D/s_Ec,X %.3f; phi %s"
                      % (r["row"], r["N"], r["Delta"], r["Delta"] / ref, ref, r["ratioEc"], fmt(r["phi"])))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
