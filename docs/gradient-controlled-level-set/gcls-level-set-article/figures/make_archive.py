#!/usr/bin/env python3
"""The secondary-data archive of the pre-print, one folder per software version.

The pre-print reads its numbers only from data/archive/<version>/. This script writes that
folder from the raw output: the method-gate studies (on the cluster), the preserved laptop runs
(runs/gcls-laptop-20260927), and the pre-fix gate record. Nothing in the archive is typed by hand.

  make_archive.py gate   --studies-dir STUDIES --summary-dir STUDIES/methodGate2D_summary \
                         --out ARCHIVE/gate
  make_archive.py prefix --studies-dir STUDIES --summary-dir STUDIES/methodGate2D_summary_pre-... \
                         --out ARCHIVE/prefix
  make_archive.py laptop --runs runs/gcls-laptop-20260927 --out ARCHIVE/laptop
  make_archive.py manifest --out ARCHIVE      # MANIFEST.csv: every file, rows, bytes, sha256

Only the Python standard library is used, so the same file runs on the laptop and on the
cluster login node. Reduced time histories keep about 300 rows per case plus the last row.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
CLASSIFIER = os.path.join(REPO, "workflow", "scripts", "foam_log_state.sh")

DROPLET_COLS = ["TIME", "maxMagU", "meanMagUPrime", "l2MagUPrime", "phaseVolumeRelError",
                "centroidX", "centroidY", "zeroSetRadialL2", "m2Amplitude", "pLaplace",
                "kErrL2Band", "gradPsiL2ErrorBand", "minGradPsiBand", "maxGradPsiBand"]
KIN_COLS = {"leiaSemiLagrangeLevelSetFoam.csv": ["TIME", "E_GEOM_ALPHA_REL", "E_VOL_ALPHA_REL"],
            "gradPsiError.csv": ["TIME", "E_NARROW_L2_GRAD_PSI", "E_NARROW_L1_GRAD_PSI"]}
TWO_PHASE_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"


def rows(path):
    # A run that stopped on a floating-point exception can leave its last line incomplete;
    # a row with a missing column is dropped.
    with open(path) as f:
        return [r for r in csv.DictReader(f) if None not in r.values() and None not in r]


def reduce_rows(R, n=300):
    if not R:
        return []
    k = max(1, len(R) // n)
    out = R[::k]
    if out[-1] is not R[-1]:
        out.append(R[-1])
    return out


def state(case_dir):
    logs = sorted(glob.glob(os.path.join(case_dir, "log.leia*Foam")) +
                  glob.glob(os.path.join(case_dir, "log.solver")))
    if not logs:
        return "NOLOG", ""
    r = subprocess.run(["bash", CLASSIFIER, logs[0]], capture_output=True, text=True)
    parts = r.stdout.split()
    steps = next((p.split("=", 1)[1] for p in parts if p.startswith("steps=")), "")
    return (parts[0] if parts else "UNKNOWN"), steps


def write(path, header, table):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(table)
    print(f"[archive] {os.path.relpath(path)} ({len(table)} rows)")


def n_cells(case_dir):
    try:
        tok = json.load(open(os.path.join(case_dir, "case_params.json"))).get("tokens", {})
        return int(float(tok.get("N_CELLS")))
    except (OSError, ValueError, TypeError):
        m = re.search(r"^n_cells\s+(\d+);", open(os.path.join(case_dir, "system", "blockMeshDict")).read(), re.M)
        return int(m.group(1)) if m else -1


def lib_stamp(case_dir):
    p = os.path.join(case_dir, "leia.version")
    if not os.path.isfile(p):
        return ""
    for line in open(p):
        if line.startswith("libleiaCore"):
            return line.split()[1]
    return ""


# --------------------------------------------------------------------------- gate / prefix
def study_dir(studies_dir, study, suffix_glob):
    """The study folder; with suffix_glob the renamed twin (the latest match), because the
    PRESERVE step renamed the pre-fix studies with a dated suffix and their manifests still
    name the original study."""
    if not suffix_glob:
        return os.path.join(studies_dir, study)
    m = sorted(glob.glob(os.path.join(studies_dir, study + suffix_glob)))
    return m[-1] if m else os.path.join(studies_dir, study + "__missing__")


def cmd_gate(a, copy_summaries=True, suffix_glob=""):
    os.makedirs(a.out, exist_ok=True)
    cands = sorted(d for d in os.listdir(a.summary_dir)
                   if os.path.isfile(os.path.join(a.summary_dir, d, "manifest.json")))
    hist = {}
    cases = []
    for c in cands:
        sd = os.path.join(a.summary_dir, c)
        if copy_summaries:
            for f in ("summary.csv", "orders.csv", "vsBaseline.csv", "seam.csv", "verdict.txt",
                      "candidate.json", "exact1D.csv"):
                if os.path.isfile(os.path.join(sd, f)):
                    os.makedirs(os.path.join(a.out, "summaries", c), exist_ok=True)
                    shutil.copy2(os.path.join(sd, f), os.path.join(a.out, "summaries", c, f))
        man = json.load(open(os.path.join(sd, "manifest.json")))
        for arm in man["arms"]:
            study = study_dir(a.studies_dir, arm["study"], suffix_glob)
            for case in sorted(glob.glob(os.path.join(study, "*_0*"))):
                if not os.path.isdir(case):
                    continue
                st, steps = state(case)
                N = n_cells(case)
                cases.append([c, arm["arm"], N, arm["np"], st, steps, lib_stamp(case),
                              os.path.relpath(case, a.studies_dir)])
                if arm["kind"] == "twoPhase":
                    p = os.path.join(case, TWO_PHASE_CSV)
                    if os.path.isfile(p):
                        for r in reduce_rows(rows(p)):
                            hist.setdefault(arm["arm"], []).append(
                                [c, N, arm["np"], st] + [r.get(k, "") for k in DROPLET_COLS])
                elif arm["arm"] in ("shear", "seamNp1", "seamNp8"):
                    merged = {}
                    for fname, cols in KIN_COLS.items():
                        p = os.path.join(case, fname)
                        if os.path.isfile(p):
                            for r in rows(p):
                                merged.setdefault(r["TIME"], {}).update({k: r.get(k, "") for k in cols})
                    R = [{**v, "TIME": t} for t, v in sorted(merged.items(), key=lambda kv: float(kv[0]))]
                    for r in reduce_rows(R):
                        hist.setdefault(arm["arm"], []).append(
                            [c, N, arm["np"], st] + [r.get(k, "") for fname in KIN_COLS
                                                      for k in KIN_COLS[fname] if k != "TIME"]
                            + [r["TIME"]])
    write(os.path.join(a.out, "cases.csv"),
          ["candidate", "arm", "N", "np", "state", "steps", "libleiaCore", "case"], cases)
    for arm, table in hist.items():
        if arm in ("shear", "seamNp1", "seamNp8"):
            head = ["candidate", "N", "np", "state"] + [k for f in KIN_COLS for k in KIN_COLS[f] if k != "TIME"] + ["TIME"]
        else:
            head = ["candidate", "N", "np", "state"] + DROPLET_COLS
        write(os.path.join(a.out, f"histories_{arm}.csv"), head, table)


def cmd_prefix(a):
    cmd_gate(a, copy_summaries=True, suffix_glob="_pre-*")
    # The 0.1 s translating runs of the pre-fix gate, kept as *_translating_endTime0p1_20260927.
    table = []
    for study in sorted(glob.glob(os.path.join(a.studies_dir, "methodGate2D_*_translating_endTime0p1_20260927"))):
        cand = os.path.basename(study).split("_")[1]
        for case in sorted(glob.glob(os.path.join(study, "*_0*"))):
            p = os.path.join(case, TWO_PHASE_CSV)
            if not os.path.isfile(p):
                continue
            R = rows(p)
            st, steps = state(case)
            table.append([cand, n_cells(case), st, steps, R[-1]["TIME"] if R else "", lib_stamp(case)])
    write(os.path.join(a.out, "translating_endTime0p1.csv"),
          ["candidate", "N", "state", "steps", "lastTime", "libleiaCore"], table)


# --------------------------------------------------------------------------- laptop
LAPTOP_TRANSLATING = [
    # run folder, box length [mm], N, np, droplet start x [mm], variant, binaries
    ("seamcheck/trans100_serial", 10, 100, 1, 2.5, "reference", "before the fixes"),
    ("seamcheck/trans100_np4", 10, 100, 4, 2.5, "reference", "before the fixes"),
    ("seamcheck/trans100_np4fix", 10, 100, 4, 2.5, "reference", "coupled-face fix"),
    ("latecheck/std142", 10, 142, 8, 2.5, "reference", "both fixes"),
    ("latecheck/none", 10, 100, 4, 2.5, "curvatureExtension none", "both fixes"),
    ("latecheck/rk2", 10, 100, 4, 2.5, "footIntegrator rk2", "both fixes"),
    ("latecheck/midpoint", 10, 100, 4, 2.5, "capillaryForceCentring midpoint", "both fixes"),
    ("latecheck/ratio1", 10, 100, 4, 2.5, "density ratio 1", "both fixes"),
    ("latecheck/longbox", 20, 100, 8, 2.5, "reference", "both fixes"),
    ("latecheck/long142", 20, 142, 8, 2.5, "reference", "both fixes"),
    ("latecheck/long200", 20, 200, 16, 2.5, "reference", "both fixes"),
    ("latecheck/long100_T25", 20, 100, 4, 2.5, "reference", "both fixes"),
    ("latecheck/long142_T25", 20, 142, 8, 2.5, "reference", "both fixes"),
    ("latecheck/long100_x5", 20, 100, 4, 5.0, "reference", "both fixes"),
    ("latecheck/box40", 40, 100, 8, 2.5, "reference", "both fixes"),
    ("latecheck/box40_142", 40, 142, 16, 2.5, "reference", "both fixes"),
]
U0 = 0.05


def jump_time(R, factor=10.0, t0=0.03):
    base = None
    for r in R:
        t, v = float(r["TIME"]), float(r["l2MagUPrime"])
        if t < t0:
            continue
        if base is None:
            base = v
        if v > factor * base:
            return t, float(r["centroidX"])
    return None, None


def colscaled(A, B, col, n):
    va = [float(r[col]) for r in A[:n]]
    vb = [float(r[col]) for r in B[:n]]
    sc = max(max(abs(x) for x in va), max(abs(x) for x in vb), 1e-300)
    return max(abs(x - y) for x, y in zip(va, vb)) / sc


def proc_face_values(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(r"(procBoundary(\d+)to(\d+))\s*\{", txt):
        i = m.end(); depth = 1; j = i
        while depth:
            depth += {"{": 1, "}": -1}.get(txt[j], 0)
            j += 1
        blk = txt[i:j]
        mv = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", blk)
        mu = re.search(r"value\s+uniform\s+([-+0-9.eE]+)", blk)
        out[(int(m.group(2)), int(m.group(3)))] = ([float(x) for x in mv.group(2).split()] if mv
                                                  else ([float(mu.group(1))] if mu else None))
    return out


def internal_field(path, vector=False):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*(\d+)\s*\(", txt)
    if not m:
        mu = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
        return float(mu.group(1)) if mu and not vector else None
    n = int(m.group(2)); body = txt[m.end():]
    if vector:
        v = re.findall(r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", body)[:n]
        return [tuple(float(x) for x in t) for t in v]
    return [float(x) for x in body.split(")")[0].split()[:n]]


def cmd_laptop(a):
    R0 = a.runs
    # 1. the translating runs: state, divergence, the growth onset, and reduced histories
    runs, hist = [], []
    for folder, L, N, npr, x0, variant, binaries in LAPTOP_TRANSLATING:
        d = os.path.join(R0, folder)
        p = os.path.join(d, TWO_PHASE_CSV)
        if not os.path.isfile(p):
            continue
        R = rows(p)
        st, steps = state(d)
        tj, cj = jump_time(R)
        endT = re.search(r"^endTime\s+([0-9.eE+-]+);", open(os.path.join(d, "system", "controlDict")).read(), re.M).group(1)
        last = R[-1]
        runs.append([os.path.basename(folder), L, N, npr, x0, variant, binaries, endT, st, steps,
                     last["TIME"], "" if tj is None else f"{tj:.6g}",
                     "" if cj is None else f"{cj * 1e3:.4g}", lib_stamp(d)])
        for r in reduce_rows(R):
            lead = (float(r["centroidX"]) - (x0 * 1e-3 + U0 * float(r["TIME"]))) * 1e3
            hist.append([os.path.basename(folder), L, N, npr, x0, variant] +
                        [r.get(k, "") for k in DROPLET_COLS] + [f"{lead:.6g}"])
    write(os.path.join(a.out, "translating_runs.csv"),
          ["run", "boxLength_mm", "N", "np", "startX_mm", "variant", "binaries", "endTime", "state",
           "steps", "lastTime", "growthOnsetTime", "centroidAtOnset_mm", "libleiaCore"], runs)
    write(os.path.join(a.out, "translating_histories.csv"),
          ["run", "boxLength_mm", "N", "np", "startX_mm", "variant"] + DROPLET_COLS + ["centroidLead_mm"], hist)

    # 2. np 4 against serial, column-scaled maximum difference over the first 0.05 s
    S = rows(os.path.join(R0, "seamcheck/trans100_serial", TWO_PHASE_CSV))
    O = rows(os.path.join(R0, "seamcheck/trans100_np4", TWO_PHASE_CSV))
    F = rows(os.path.join(R0, "seamcheck/trans100_np4fix", TWO_PHASE_CSV))
    M = rows(os.path.join(R0, "seamcheck/m_np4", TWO_PHASE_CSV))
    n = min(4604, len(S), len(O), len(F), len(M))
    table = []
    for col in S[0]:
        if col in ("TIME", "ELAPSED_CPU_TIME", "ELAPSED_CLOCK_TIME") or not col:
            continue
        try:
            table.append([col, f"{colscaled(S, O, col, n):.3e}", f"{colscaled(S, F, col, n):.3e}",
                          f"{colscaled(S, M, col, n):.3e}"])
        except (ValueError, KeyError):
            pass
    write(os.path.join(a.out, "seam_np4_vs_serial.csv"),
          ["column", "beforeTheFixes", "coupledFaceFix", "bothFixes"], table)

    # 3. the face density on the two sides of every processor face, one step and five steps
    table = []
    for run, label in (("seamcheck/trans64", "before"), ("seamcheck/trans64fix", "after")):
        d = os.path.join(R0, run)
        procs = sorted(int(x[9:]) for x in os.listdir(d) if x.startswith("processor"))
        for t in ("2.12129e-05", "0.0001060645"):
            for fld in ("phi", "alphaf", "rhof"):
                per = {p: proc_face_values(os.path.join(d, f"processor{p}", t, fld)) for p in procs}
                nface = ndiff = 0; worst = 0.0
                for p in procs:
                    for (i, j), va in per[p].items():
                        if i > j or va is None:
                            continue
                        vb = per[j].get((j, i))
                        if vb is None:
                            continue
                        if len(va) == 1 and len(vb) > 1: va = va * len(vb)
                        if len(vb) == 1 and len(va) > 1: vb = vb * len(va)
                        for x, y in zip(va, vb):
                            dd = abs(x + y) if fld == "phi" else abs(x - y)
                            s = max(abs(x), abs(y), 1e-300)
                            nface += 1; ndiff += dd > 1e-12 * s; worst = max(worst, dd / s)
                table.append([label, t, fld, nface, ndiff, f"{worst:.3e}"])
    write(os.path.join(a.out, "seam_face_density.csv"),
          ["code", "time", "field", "facesCompared", "facesDiffering", "worstRelative"], table)

    # 4. the Eulerian two-phase solver with and without the shared mass flux (1000 steps)
    table = []
    for run, binary, ratio, npr in (("trans100_new", "shared rhoLENT mass flux", 838.8, 4),
                                    ("trans100_new_serial", "shared rhoLENT mass flux", 838.8, 1),
                                    ("trans100_old", "frozen density", 838.8, 4),
                                    ("trans100r1_new", "shared rhoLENT mass flux", 1.0, 4),
                                    ("trans100r1_old", "frozen density", 1.0, 4)):
        d = os.path.join(R0, "eulgate", run, "0.010861")
        al = internal_field(os.path.join(d, "alpha.water"))
        U = internal_field(os.path.join(d, "U"), vector=True)
        C = internal_field(os.path.join(d, "C"), vector=True)
        V = internal_field(os.path.join(d, "V"))
        if not isinstance(V, list):
            V = [V] * len(al)
        vol = sum(x * v for x, v in zip(al, V))
        xc = sum(x * v * c[0] for x, v, c in zip(al, V, C)) / vol
        mag = [math.sqrt((u[0] - U0) ** 2 + u[1] ** 2 + u[2] ** 2) for u in U]
        tot = sum(V)
        L2 = math.sqrt(sum(m * m * v for m, v in zip(mag, V)) / tot)
        L1 = sum(m * v for m, v in zip(mag, V)) / tot
        table.append([run, binary, ratio, npr, 1000, 0.010861,
                      f"{(xc - 0.0025) / (U0 * 0.010861):.6f}", f"{L2:.4e}", f"{L1:.4e}"])
    write(os.path.join(a.out, "eulerian_mass_flux.csv"),
          ["run", "solver", "densityRatio", "np", "steps", "time", "travelledFraction",
           "L2_U_minus_U0", "L1_U_minus_U0"], table)


def cmd_manifest(a):
    import hashlib
    table = []
    for dp, _, fns in os.walk(a.out):
        for fn in sorted(fns):
            if fn == "MANIFEST.csv":
                continue
            p = os.path.join(dp, fn)
            data = open(p, "rb").read()
            nrows = data.count(b"\n") - 1 if fn.endswith(".csv") else ""
            table.append([os.path.relpath(p, a.out), nrows, len(data), hashlib.sha256(data).hexdigest()])
    table.sort()
    write(os.path.join(a.out, "MANIFEST.csv"), ["file", "rows", "bytes", "sha256"], table)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("gate", "prefix"):
        s = sub.add_parser(name)
        s.add_argument("--studies-dir", required=True)
        s.add_argument("--summary-dir", required=True)
        s.add_argument("--out", required=True)
    s = sub.add_parser("laptop")
    s.add_argument("--runs", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("manifest")
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    {"gate": cmd_gate, "prefix": cmd_prefix, "laptop": cmd_laptop, "manifest": cmd_manifest}[a.cmd](a)


if __name__ == "__main__":
    main()
