#!/usr/bin/env python3
"""Where the face-fraction error of a straight translation builds up
(config/davof/sphereRotation3D_hexFineLine*Faces.yaml, 2026-11-02).

Input per case: postProcessing/areaFractionError/faceError/<timeIndex>_<rank>.csv
(fork areaFractionError, key faceErrorCsv): per write the faces counted in
L1(alphaf) with centre x_f, area vector S_f, alphaf and alphaf_exact, and a
comment line `# t <t> origin <x> <y> <z>` (the exact sphere's centre c(t)).

Per face e = |alphaf - alphaf_exact|; L1 = sum e / nFaces, nFaces the rows of
the step (checked against the L1 column of areaFractionError.dat). Bins:
  dir     the face direction, the largest component of |S_f| (x, y, z);
  side    the exact outward normal n = (x_f - c)/|x_f - c| against the travel
          direction u = U/|U|, theta = acos(n.u): front [0, 60), tangent
          [60, 120), rear [120, 180] degrees;
  theta   theta in six bins of 30 degrees;
  axis    the angle of n to the nearest mesh axis, acos(max_k |n_k|): [0, 15),
          [15, 30), [30, 45), [45, 54.8] degrees;
  dirside dir x side.
Per bin and step the contribution C_b = sum_{bin} e / nFaces (the C_b of a
binning add up to L1) and the count n_b. Per bin: C_b and n_b/nFaces at the
first step and at the end, the late slope of C_b (least squares against t
over [0.25 T, T], T = END_TIME, the measure of davof_line_confound.py; the
slopes of a binning add up to the slope of L1) and its share of the total
late slope, and the mean error per face sum_{bin} e / n_b at the first step
and at the end.

    python3 workflow/scripts/davof_face_bins.py studies/davof/<study> [...] [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os
import re

import numpy as np

SIDE = [("front", 0, 60), ("tangent", 60, 120), ("rear", 120, 180.0001)]
THETA = [(f"{a}-{a + 30}", a, a + 30 + (1e-4 if a == 150 else 0)) for a in range(0, 180, 30)]
AXIS = [("0-15", 0, 15), ("15-30", 15, 30), ("30-45", 30, 45), ("45-55", 45, 55)]


def pick(bins, v):
    for name, lo, hi in bins:
        if lo <= v < hi:
            return name
    return bins[-1][0]


def read_step(files):
    t = o = None
    rows = []
    for fn in files:
        with open(fn) as fh:
            for line in fh:
                if line.startswith("#"):
                    m = re.match(r"#\s*t\s+(\S+)\s+origin\s+(\S+)\s+(\S+)\s+(\S+)", line)
                    if m:
                        t = float(m.group(1))
                        o = np.array([float(m.group(i)) for i in (2, 3, 4)])
                    continue
                if line.startswith("x,"):
                    continue
                rows.append([float(v) for v in line.split(",")])
    return t, o, np.array(rows)


def analyse(case):
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    U = np.array([float(tok.get(k) or 0) for k in ("TPF_UX", "TPF_UY", "TPF_UZ")])
    u = U / np.linalg.norm(U)
    T = float(tok["END_TIME"])
    d = os.path.join(case, "postProcessing", "areaFractionError", "faceError")
    by = {}
    for fn in glob.glob(os.path.join(d, "*_*.csv")):
        by.setdefault(int(os.path.basename(fn).split("_")[0]), []).append(fn)
    dat = {}
    p = glob.glob(os.path.join(case, "postProcessing", "areaFractionError", "*", "areaFractionError.dat"))
    if p:
        for line in open(p[0]):
            s = line.split()
            if s and not s[0].startswith("#"):
                dat[f"{float(s[0]):.9g}"] = float(s[1])
    binnings = {"dir": None, "side": SIDE, "theta": THETA, "axis": AXIS, "dirside": None}
    hist = {b: {} for b in binnings}
    ts, L1s, nF, checkMax = [], [], [], 0.0
    for k in sorted(by):
        t, o, R = read_step(by[k])
        if R.size == 0:
            continue
        x, S, af, ae = R[:, 0:3], R[:, 3:6], R[:, 6], R[:, 7]
        e = np.abs(af - ae)
        n = len(e)
        L1 = e.sum() / n
        ref = dat.get(f"{t:.9g}")
        if ref is not None and ref > 0:
            checkMax = max(checkMax, abs(L1 - ref) / ref)
        r = x - o
        nn = r / np.linalg.norm(r, axis=1)[:, None]
        theta = np.degrees(np.arccos(np.clip(nn @ u, -1, 1)))
        axis = np.degrees(np.arccos(np.clip(np.abs(nn).max(axis=1), -1, 1)))
        dirs = np.array(["xyz"[i] for i in np.abs(S).argmax(axis=1)])
        side = np.array([pick(SIDE, v) for v in theta])
        labels = {"dir": dirs, "side": side,
                  "theta": np.array([pick(THETA, v) for v in theta]),
                  "axis": np.array([pick(AXIS, v) for v in axis]),
                  "dirside": np.char.add(np.char.add(dirs, "/"), side)}
        ts.append(t); L1s.append(L1); nF.append(n)
        for b, lab in labels.items():
            for name in np.unique(lab):
                m = lab == name
                hist[b].setdefault(name, {})[len(ts) - 1] = (e[m].sum() / n, int(m.sum()), e[m].sum())
    ts, L1s = np.array(ts), np.array(L1s)
    late = (ts >= 0.25 * T) & (ts <= T * (1 + 1e-9))
    sTot = float(np.polyfit(ts[late], L1s[late], 1)[0])
    out = []
    for b, bins in hist.items():
        for name in sorted(bins):
            series = bins[name]
            C = np.array([series.get(i, (0.0, 0, 0.0))[0] for i in range(len(ts))])
            cnt = np.array([series.get(i, (0.0, 0, 0.0))[1] for i in range(len(ts))])
            esum = np.array([series.get(i, (0.0, 0, 0.0))[2] for i in range(len(ts))])
            s = float(np.polyfit(ts[late], C[late], 1)[0])
            mean = lambda i: esum[i] / cnt[i] if cnt[i] else float("nan")
            out.append(dict(case=case, N=int(float(tok["N_CELLS"])), binning=b, bin=name,
                            C_first=C[0], C_end=C[-1], frac_first=cnt[0] / nF[0], frac_end=cnt[-1] / nF[-1],
                            lateSlope=s, slopeShare=s / sTot if sTot else float("nan"),
                            meanFirst=mean(0), meanEnd=mean(len(ts) - 1)))
    info = dict(case=case, steps=len(ts), L1_first=L1s[0], L1_end=L1s[-1], lateSlope=sTot,
                nFaces_first=nF[0], nFaces_end=nF[-1], L1_check_maxrel=checkMax)
    return info, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    allrows = []
    for st in a.studies:
        for case in sorted(glob.glob(os.path.join(st, "*_0*"))):
            if "." in os.path.basename(case) or not os.path.isdir(
                    os.path.join(case, "postProcessing", "areaFractionError", "faceError")):
                continue
            info, rows = analyse(case)
            print("== %s steps %d L1 first %.4e end %.4e late slope %.4e nFaces %d -> %d; L1 from the rows vs .dat max rel %.2e"
                  % (info["case"], info["steps"], info["L1_first"], info["L1_end"], info["lateSlope"],
                     info["nFaces_first"], info["nFaces_end"], info["L1_check_maxrel"]))
            print("binning\tbin\tC_first\tC_end\tfrac_first\tfrac_end\tlateSlope\tslopeShare\tmeanFirst\tmeanEnd")
            for r in rows:
                print("%s\t%s\t%.3e\t%.3e\t%.3f\t%.3f\t%.3e\t%.3f\t%.3e\t%.3e" % (
                    r["binning"], r["bin"], r["C_first"], r["C_end"], r["frac_first"], r["frac_end"],
                    r["lateSlope"], r["slopeShare"], r["meanFirst"], r["meanEnd"]))
            allrows += rows
    if a.csv and allrows:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
            w.writeheader()
            w.writerows(allrows)


if __name__ == "__main__":
    main()
