#!/usr/bin/env python3
"""The consistency of the swept flux under oblique travel
(config/davof/sphereRotation3D_hexFineLine*Cons*.yaml, 2026-11-02; fork
davof-swept c635611, keys sweptConsistencyLog and sweptConsistent).

Per case:
  film     the cells of the written alpha at the end by band, [1e-12, 1e-6),
           [1e-6, 1e-3), [1e-3, 1e-2), [1e-2, 0.99), [0.99, 1 - 1e-6); the
           film = the cells with alpha in [1e-12, 1e-3), with their volume
           (in cell volumes); the interface cells of the log's last step;
  diag     (postProcessing/sweptConsistency/summary.dat, sweptConsistencyLog
           true) summed over the steps: the algebraic faces probed, those
           whose swept volume reaches a cell in another state (a surface, a
           wisp, a pure cell of the other state), the probe fallbacks, the
           mismatch of the consistent rule minus the algebraic flux (sum, sum
           of |.|, largest |.|; in cell volumes h^3), the wisp pieces of the
           computed faces and their mismatch (0/1 rule minus alpha V_piece),
           the residue (sum over the steps of sum_c |R_c|, largest |R|/V);
           per step means of the counts; sweptConsistent's added candidates,
           faces without a surface and their fallbacks;
  residue  (residue_<step>.csv: per step the cells with |R| > 1e-15 V, R in
           cell volumes) the cells that received a residue in any step, their
           overlap with the film at the end (the share of the film cells that
           received one, the share of the residue cells that end in the film),
           the cumulative signed residue (sum over the steps of R/V) in the
           film cells against the film's liquid (sum of alpha); the share of
           sum |R| by the side of the sphere, theta = the angle between
           x_c - c(t^n) and the travel direction (front [0, 60), tangent
           [60, 120), rear [120, 180] degrees; c(t^n) = x0 + U0 t^n, t^n =
           (step - 1) dt, the state the step starts from), over all steps and
           over the late window t^n in [0.25 T, T]; and the share of sum |R|
           by the distance of the cell centre from the exact sphere,
           |x_c - c| - R0 in cells h: [0, 1), [1, 2), [2, inf) and inside.

    python3 workflow/scripts/davof_consistency.py studies/davof/<study> [...] [--csv out.csv]
"""
import argparse
import csv
import glob
import json
import math
import os
import re

import numpy as np

BANDS = [(1e-12, 1e-6), (1e-6, 1e-3), (1e-3, 1e-2), (1e-2, 0.99), (0.99, 1 - 1e-6)]
SIDE = [("front", 0, 60), ("tangent", 60, 120), ("rear", 120, 180.0001)]
DIST = [("inside", -1e9, 0), ("0-1h", 0, 1), ("1-2h", 1, 2), (">2h", 2, 1e9)]
SUMMARY = ["time", "step", "nProbed", "nReach", "nReachInterface", "nReachWisp", "nReachPure",
           "nProbeFallback", "algMisSum", "algMisAbs", "algMisMax", "nWispPieces", "wispMisSum",
           "wispMisAbs", "wispMisMax", "nResidueCells", "residueAbs", "residueMaxRel", "nCandAdded",
           "nNoSurface", "nConsFallback"]


def f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def cases(study):
    for c in sorted(glob.glob(os.path.join(study, "*_0*"))):
        if "." in os.path.basename(c) or not os.path.isfile(os.path.join(c, "case_params.json")):
            continue
        yield c


def read_internal(path):
    text = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\(", text)
    if m:
        return np.array(text[m.end():text.index(")", m.end())].split(), dtype=float)
    m = re.search(r"internalField\s+uniform\s+(\S+);", text)
    return float(m.group(1))


def measure(case):
    tok = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    N = int(float(tok["N_CELLS"]))
    h = 1.0 / N
    T = f(tok["END_TIME"])
    dt = f(tok["MAX_DELTA_T"])
    U = np.array([f(tok.get(k)) or 0.0 for k in ("TPF_UX", "TPF_UY", "TPF_UZ")])
    x0 = np.array([f(tok.get(k)) for k in ("TPF_SPHERE_CX", "TPF_SPHERE_CY", "TPF_SPHERE_CZ")])
    R0 = f(tok.get("TPF_SPHERE_RADIUS")) or 0.15
    u = U / max(np.linalg.norm(U), 1e-300)
    r = dict(study=os.path.basename(os.path.dirname(case)), case=os.path.basename(case), N=N,
             arm=tok.get("TPF_VELOCITY_ARM"), consistent=tok.get("TPF_SWEPT_CONSISTENT", "false"))

    # the film at the end
    ts = sorted([d for d in os.listdir(case) if re.match(r"^[0-9.eE+-]+$", d) and f(d) and f(d) > 0
                 and os.path.isfile(os.path.join(case, d, "alpha.water"))], key=float)
    a = read_internal(os.path.join(case, ts[-1], "alpha.water"))
    r["tEnd"] = ts[-1]
    for lo, hi in BANDS:
        r["n[%g,%g)" % (lo, hi)] = int(((a >= lo) & (a < hi)).sum())
    film = (a >= 1e-12) & (a < 1e-3)
    r["nFilm"] = int(film.sum())
    r["filmVolume"] = float(a[film].sum())
    r["minAlphaEnd"] = float(a.min())
    log = os.path.join(case, "log.advectorVoF")
    if os.path.isfile(log):
        last = None
        for line in open(log, errors="replace"):
            m = re.match(r"davof: (\d+) interface cells", line)
            if m:
                last = int(m.group(1))
        r["interfaceCellsEnd"] = last

    # the diagnostic
    sp = os.path.join(case, "postProcessing", "sweptConsistency", "summary.dat")
    if not os.path.isfile(sp):
        return r
    S = np.loadtxt(sp, comments="#", ndmin=2)
    col = {k: S[:, i] for i, k in enumerate(SUMMARY)}
    h3 = h ** 3
    r["steps"] = len(S)
    for k in ("nProbed", "nReach", "nReachInterface", "nReachWisp", "nReachPure", "nWispPieces",
              "nResidueCells", "nCandAdded", "nNoSurface"):
        r[k + "PerStep"] = float(col[k].mean())
    for k in ("nProbeFallback", "nConsFallback"):
        r[k] = int(col[k].sum())
    for k in ("algMisSum", "algMisAbs", "wispMisSum", "wispMisAbs", "residueAbs"):
        r[k] = float(col[k].sum() / h3)
    for k in ("algMisMax", "wispMisMax"):
        r[k] = float(col[k].max() / h3)
    r["residueMaxRel"] = float(col["residueMaxRel"].max())

    # the residues per cell and by side and distance
    Rcum = {}
    Rabs = {}
    side = {s[0]: [0.0, 0.0] for s in SIDE}
    dist = {s[0]: 0.0 for s in DIST}
    total = [0.0, 0.0]
    for fn in glob.glob(os.path.join(case, "postProcessing", "sweptConsistency", "residue_*.csv")):
        step = int(re.search(r"residue_(\d+)\.csv", fn).group(1))
        tn = (step - 1) * dt
        late = tn >= 0.25 * T - 1e-12
        d = np.genfromtxt(fn, delimiter=",", names=True, ndmin=1)
        if d.size == 0:
            continue
        c = x0 + U * tn
        X = np.stack([d["x"], d["y"], d["z"]], axis=1) - c
        rr = np.linalg.norm(X, axis=1)
        cosT = np.clip((X @ u) / np.maximum(rr, 1e-300), -1, 1)
        th = np.degrees(np.arccos(cosT))
        Ra = np.abs(d["R"])
        for name, lo, hi in SIDE:
            m = (th >= lo) & (th < hi)
            side[name][0] += float(Ra[m].sum())
            if late:
                side[name][1] += float(Ra[m].sum())
        dd = (rr - R0) / h
        for name, lo, hi in DIST:
            m = (dd >= lo) & (dd < hi)
            dist[name] += float(Ra[m].sum())
        total[0] += float(Ra.sum())
        if late:
            total[1] += float(Ra.sum())
        for cell, R in zip(d["cell"].astype(int), d["R"]):
            Rcum[cell] = Rcum.get(cell, 0.0) + R
            Rabs[cell] = Rabs.get(cell, 0.0) + abs(R)
    cellsR = np.array(sorted(Rabs), dtype=int)
    r["nResidueCellsAny"] = int(cellsR.size)
    if cellsR.size:
        inFilm = film[cellsR]
        r["filmShareWithResidue"] = float(inFilm.sum() / max(r["nFilm"], 1))
        r["residueCellsShareInFilm"] = float(inFilm.mean())
        rc = np.array([Rcum[c] for c in cellsR])
        r["filmNetResidue"] = float(rc[inFilm].sum())
        r["filmNetResiduePos"] = float(np.maximum(rc[inFilm], 0).sum())
        r["netResidueAll"] = float(rc.sum())
    for name, _, _ in SIDE:
        r["share_" + name] = side[name][0] / total[0] if total[0] > 0 else None
        r["shareLate_" + name] = side[name][1] / total[1] if total[1] > 0 else None
    for name, _, _ in DIST:
        r["shareDist_" + name] = dist[name] / total[0] if total[0] > 0 else None
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    out = []
    for st in a.studies:
        for c in cases(st):
            r = measure(c)
            out.append(r)
            print("== %s/%s N %d arm %s consistent %s" % (r["study"], r["case"], r["N"], r["arm"], r["consistent"]))
            for k, v in r.items():
                if k in ("study", "case", "N", "arm", "consistent"):
                    continue
                print("  %-26s %s" % (k, ("%.4g" % v) if isinstance(v, float) else v))
    if a.csv and out:
        keys = []
        for r in out:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(out)


if __name__ == "__main__":
    main()
