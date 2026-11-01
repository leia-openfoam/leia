#!/usr/bin/env python3
"""The decomposition comparison of the parallel swept flux (2026-10-30, fork
davof-swept-par): the volume fractions and the per-face results of a case
against a reference case of the same settings, step by step.

Per case and per common time step:
  - maxAlphaDiff: max over the cells of |alpha - alpha_ref| (the fields at the
    case root; a parallel case after reconstructPar);
  - from the per-face logs postProcessing/sweptFaceLog/<step>_<rank>.csv of
    both cases (fork key sweptFaceLog true; one row per face computed, with
    the face centre, the face fraction, the flux in the face's canonical
    orientation, the remote flag and the processor-face flag): the faces are
    matched by their centres (rounded to 1e-9), maxAlphafDiff and
    maxDVfDiff (the flux difference in units of the cell volume h^3, h = 1/N
    from case_params.json), the faces of the reference missing in the case
    (nMissing: no value), those only in the case (nExtra), those computed
    twice (nDup: two ranks), and the crossing faces (remote or processor-face
    flag in the case) with a difference above --tol in alphaf or dVf/h^3
    (nCrossDiff) against the other faces above --tol (nOtherDiff, the
    control).

    python3 davof_decomposition_compare.py --ref <case> --case <case> [...] --out cmp.csv
"""
import argparse
import csv
import glob
import json
import os
import re
import sys


def read_internal(path):
    text = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\(", text)
    if m:
        n = int(m.group(1))
        end = text.index(")", m.end())
        vals = [float(v) for v in text[m.end():end].split()]
        assert len(vals) == n, path
        return vals
    m = re.search(r"internalField\s+uniform\s+(\S+);", text)
    if m:
        return float(m.group(1))
    raise ValueError(f"no internalField in {path}")


def times(case):
    out = {}
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0 and os.path.isfile(os.path.join(case, d, "alpha.water")):
            out[f"{t:.9g}"] = d
    return out


def face_logs(case):
    """step -> {key: (alphaf, dVf, remote, proc, count)}"""
    steps = {}
    for f in glob.glob(os.path.join(case, "postProcessing", "sweptFaceLog", "*_*.csv")):
        step = int(os.path.basename(f).split("_")[0])
        d = steps.setdefault(step, {})
        with open(f) as fh:
            for row in csv.DictReader(fh):
                key = (round(float(row["x"]), 9), round(float(row["y"]), 9), round(float(row["z"]), 9))
                if key in d:
                    a, v, r, p, c = d[key]
                    d[key] = (a, v, r, p, c + 1)
                else:
                    d[key] = (float(row["alphaf"]), float(row["dVf"]),
                              int(row["remote"]), int(row["procFace"]), 1)
    return steps


def case_n(case):
    try:
        with open(os.path.join(case, "case_params.json")) as fh:
            p = json.load(fh)
        for k in ("N_CELLS",):
            if k in p:
                return int(float(p[k]))
            if "params" in p and k in p["params"]:
                return int(float(p["params"][k]))
    except (OSError, ValueError):
        pass
    return None


def nprocs(case):
    return len(glob.glob(os.path.join(case, "processor*"))) or 1


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--case", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--N", type=int, default=None)
    a = ap.parse_args(argv[1:])

    ref_t = times(a.ref)
    ref_f = face_logs(a.ref)
    N = a.N or case_n(a.ref)
    V = (1.0/N)**3 if N else 1.0
    rows = []
    for case in a.case:
        ct = times(case)
        cf = face_logs(case)
        common = sorted(set(ref_t) & set(ct), key=float)
        for k, tk in enumerate(common, start=1):
            ar = read_internal(os.path.join(a.ref, ref_t[tk], "alpha.water"))
            ac = read_internal(os.path.join(case, ct[tk], "alpha.water"))
            if isinstance(ar, float) or isinstance(ac, float):
                dmax = float("nan")
            else:
                dmax = max(abs(x - y) for x, y in zip(ar, ac))
            row = {"case": os.path.basename(os.path.normpath(case)), "np": nprocs(case),
                   "step": k, "time": tk, "maxAlphaDiff": f"{dmax:.3e}"}
            fr, fc = ref_f.get(k), cf.get(k)
            if fr is not None and fc is not None:
                mx_af = mx_dv = 0.0
                n_cross = n_cross_diff = n_other_diff = n_remote = 0
                missing = sum(1 for key in fr if key not in fc)
                extra = sum(1 for key in fc if key not in fr)
                dup = sum(1 for v in fc.values() if v[4] > 1)
                for key, (af, dv, r, p, c) in fc.items():
                    if key not in fr:
                        continue
                    afr, dvr = fr[key][0], fr[key][1]
                    daf, ddv = abs(af - afr), abs(dv - dvr)/V
                    mx_af, mx_dv = max(mx_af, daf), max(mx_dv, ddv)
                    if r:
                        n_remote += 1
                    if r or p:
                        n_cross += 1
                        if max(daf, ddv) > a.tol:
                            n_cross_diff += 1
                    elif max(daf, ddv) > a.tol:
                        n_other_diff += 1
                row.update({"nFacesRef": len(fr), "nFaces": len(fc), "nMissing": missing,
                            "nExtra": extra, "nDup": dup, "maxAlphafDiff": f"{mx_af:.3e}",
                            "maxDVfDiff": f"{mx_dv:.3e}", "nRemote": n_remote,
                            "nCross": n_cross, "nCrossDiff": n_cross_diff,
                            "nOtherDiff": n_other_diff})
            rows.append(row)
    cols = ["case", "np", "step", "time", "maxAlphaDiff", "maxAlphafDiff", "maxDVfDiff",
            "nFacesRef", "nFaces", "nMissing", "nExtra", "nDup", "nRemote", "nCross",
            "nCrossDiff", "nOtherDiff"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    for r in rows:
        print(",".join(str(r.get(c, "")) for c in cols))


if __name__ == "__main__":
    main(sys.argv)
