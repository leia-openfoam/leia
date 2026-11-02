#!/usr/bin/env python3
"""Is the oblique build-up of L1(alphaf) a property of the transported alpha
or of the DAVOF face update? (2026-11-02; measures and the pre-registered
predictions K-0, H-alpha, H-face, K-1, E-1, R-1 in the header of
config/davof/sphereRotation3D_hexFineLineDiagW.yaml.)

Per case: postProcessing/staticFaceFractions/staticFaceFractions.dat (the
fork's app staticFaceFractions over the written times), the per-step
areaFractionError.dat and volumeFractionError.dat, and the converter's
advectorVoF.csv (L1(x) = LCentreDiff1). Late slope: the least-squares slope
over the written times in [0.25 T, T] (and, for the loop, over every step,
as davof_line_confound.py). Then, per N, with s_L the loop's, s_S the static,
s_0 the floor's late slope of L1(alphaf), D the diagonal, X the x-line,

    phi   = (s_S,D - s_S,X) / (s_L,D - s_L,X)       (ours)
    phi_R = (s_R,D - s_R,X) / (s_L,D - s_L,X)       (plicRDF's static)

    python3 workflow/scripts/davof_static_faces.py \
        --ours-d studies/davof/sphereRotation3D_hexFineLineDiagW --ours-x studies/davof/sphereRotation3D_hexFineLineXW \
        --rdf-d studies/davof/sphereRotation3D_hexFineLineDiagRDF --rdf-x studies/davof/sphereRotation3D_hexFineLineXRDF \
        --orig-d studies/davof/sphereRotation3D_hexFineLineDiag --orig-x studies/davof/sphereRotation3D_hexFineLineX [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os

import numpy as np


def read_table(path):
    names, rows = None, []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            q = s[1:].split()
            if q and q[0] == "Time":
                names = q
            continue
        rows.append([float(x) for x in s.split()[:len(names)]])
    return names, np.array(rows)


def cases(study):
    out = {}
    for c in sorted(glob.glob(os.path.join(study, "*_0*"))):
        p = os.path.join(c, "case_params.json")
        if "." in os.path.basename(c) or not os.path.isfile(p):
            continue
        tok = json.load(open(p))["tokens"]
        out[int(float(tok["N_CELLS"]))] = (c, tok)
    return out


def slope(t, y, T):
    m = (t >= 0.25 * T - 1e-12) & (t <= T * (1 + 1e-9))
    return float(np.polyfit(t[m], y[m], 1)[0]), int(m.sum())


def measure(c, tok):
    T = float(tok["END_TIME"])
    sn, S = read_table(os.path.join(c, "postProcessing", "staticFaceFractions", "staticFaceFractions.dat"))
    col = {k: S[:, i] for i, k in enumerate(sn)}
    t = col["Time"]
    r = dict(case=os.path.relpath(c), N=int(float(tok["N_CELLS"])), samples=len(t))
    for k in ("L1loop", "L1static", "L1staticMean", "L1staticExact", "Ec", "Eout"):
        r["s_" + k], r["nLate"] = slope(t, col[k], T) if (col[k] >= 0).all() else (float("nan"), 0)
        r[k + "_end"] = float(col[k][-1])
    r["nPlaneFail"] = int(col["nPlaneFail"].sum())
    r["nPlanes_end"] = int(col["nPlanes"][-1])
    an, A = read_table(os.path.join(c, "postProcessing", "areaFractionError", "0", "areaFractionError.dat"))
    ta, La = A[:, 0], A[:, 1]
    if (La >= 0).all():
        r["s_loopStep"], _ = slope(ta, La, T)
        # the app's L1 of the loop's alphaf against the .dat at the written times
        rel = []
        for ti, li in zip(t, col["L1loop"]):
            j = np.where(np.abs(ta - ti) < 1e-9)[0]
            if len(j) and li >= 0:
                rel.append(abs(li - La[j[0]]) / max(abs(La[j[0]]), 1e-300))
        r["loopCheckMaxRel"] = max(rel) if rel else float("nan")
        r["loopCheckN"] = len(rel)
    vn, Vf = read_table(os.path.join(c, "postProcessing", "volumeFractionError", "0", "volumeFractionError.dat"))
    tv, Es = Vf[:, 0], Vf[:, vn.index("Eshape")]
    r["Eshape_end"] = float(Es[-1])
    r["s_Eshape"], _ = slope(tv, Es, T)
    # Evol of the app at the written times against volumeFractionError's Eshape
    ev = []
    for ti, e in zip(t, col["Evol"]):
        j = np.where(np.abs(tv - ti) < 1e-9)[0]
        if len(j) and Es[j[0]] > 0:
            ev.append(e / Es[j[0]])
    r["EvolOverEshape"] = (min(ev), max(ev)) if ev else (float("nan"), float("nan"))
    rows = list(csv.DictReader(open(os.path.join(c, "advectorVoF.csv"))))
    try:
        r["L1x_end"] = float(rows[-1]["LCentreDiff1"])
    except (KeyError, ValueError):
        r["L1x_end"] = float("nan")
    return r


def same_bytes(a, b):
    return os.path.isfile(a) and os.path.isfile(b) and open(a, "rb").read() == open(b, "rb").read()


def main():
    ap = argparse.ArgumentParser()
    for k in ("ours-d", "ours-x", "rdf-d", "rdf-x", "orig-d", "orig-x"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    S = {k: cases(getattr(a, k.replace("-", "_"))) for k in ("ours-d", "ours-x", "rdf-d", "rdf-x", "orig-d", "orig-x")}
    M = {k: {N: measure(c, tok) for N, (c, tok) in S[k].items()} for k in ("ours-d", "ours-x", "rdf-d", "rdf-x")}

    print("== K-0: the .dat files of the W re-runs against the originals")
    for w, o in (("ours-d", "orig-d"), ("ours-x", "orig-x")):
        for N in sorted(S[w]):
            for fo in ("areaFractionError", "volumeFractionError"):
                p = os.path.join("postProcessing", fo, "0", fo + ".dat")
                print("  %s N %d %s: %s" % (w, N, fo, "IDENTICAL" if same_bytes(os.path.join(S[w][N][0], p), os.path.join(S[o][N][0], p)) else "DIFFERENT"))
    for k in M:
        for N, r in sorted(M[k].items()):
            print("  %s N %d: samples %d (late %d), loop check max rel %s over %s times, failed planes %d, Evol/Eshape %.6g..%.6g"
                  % (k, N, r["samples"], r["nLate"], "%.2e" % r["loopCheckMaxRel"] if "loopCheckMaxRel" in r else "-",
                     r.get("loopCheckN", "-"), r["nPlaneFail"], r["EvolOverEshape"][0], r["EvolOverEshape"][1]))

    keys = ["set", "N", "s_loopStep", "s_L1loop", "s_L1static", "s_L1staticMean", "s_L1staticExact", "s_Ec", "s_Eout",
            "s_Eshape", "L1loop_end", "L1static_end", "L1staticExact_end", "Ec_end", "Eout_end", "Eshape_end", "L1x_end",
            "nPlanes_end"]
    print("== per case (late slopes over the written times in [0.25 T, T]; s_loopStep over every step)")
    print("\t".join(keys))
    out = []
    for k in ("ours-d", "ours-x", "rdf-d", "rdf-x"):
        for N, r in sorted(M[k].items()):
            r["set"] = k
            out.append(r)
            print("\t".join(("%.4e" % r[q]) if isinstance(r.get(q), float) else str(r.get(q, "")) for q in keys))

    print("== verdicts")
    for N in sorted(M["ours-d"]):
        D, X, RD, RX = M["ours-d"][N], M["ours-x"][N], M["rdf-d"].get(N), M["rdf-x"].get(N)
        den = D["s_L1loop"] - X["s_L1loop"]
        phi = (D["s_L1static"] - X["s_L1static"]) / den
        phiC = ((D["s_L1static"] - D["s_L1staticExact"]) - (X["s_L1static"] - X["s_L1staticExact"])) / den
        k1 = max(abs(D["s_L1staticExact"]), abs(X["s_L1staticExact"])) <= 0.1 * X["s_L1loop"]
        print("  N %d: loop s_D %.4e s_X %.4e (rho %.3f); static s_D %.4e s_X %.4e (rho %.3f); phi %.3f; floor s_0 D %.3e X %.3e -> K-1 %s; phi from s_S - s_0 %.3f; mean rule phi %.3f"
              % (N, D["s_L1loop"], X["s_L1loop"], D["s_L1loop"] / X["s_L1loop"], D["s_L1static"], X["s_L1static"],
                 D["s_L1static"] / X["s_L1static"] if X["s_L1static"] else float("nan"), phi,
                 D["s_L1staticExact"], X["s_L1staticExact"], "HELD" if k1 else "NOT HELD", phiC,
                 (D["s_L1staticMean"] - X["s_L1staticMean"]) / den))
        print("        cell level: s_Ec D %.4e X %.4e (ratio %.3f); s_Eout D %.4e X %.4e"
              % (D["s_Ec"], X["s_Ec"], D["s_Ec"] / X["s_Ec"] if X["s_Ec"] else float("nan"), D["s_Eout"], X["s_Eout"]))
        if RD and RX:
            phiR = (RD["s_L1static"] - RX["s_L1static"]) / den
            print("        plicRDF: static s_D %.4e s_X %.4e; phi_R %.3f; s_Ec D %.4e X %.4e (ratio %.3f); Eshape end D %.4e X %.4e; L1(x) end D %.4e X %.4e"
                  % (RD["s_L1static"], RX["s_L1static"], phiR, RD["s_Ec"], RX["s_Ec"],
                     RD["s_Ec"] / RX["s_Ec"] if RX["s_Ec"] else float("nan"), RD["Eshape_end"], RX["Eshape_end"],
                     RD["L1x_end"], RX["L1x_end"]))

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
