#!/usr/bin/env python3
"""Convert the output of the TwoPhaseFlow advection benchmark solver
(advectorVoF of the fork, branch davof) into the per-case CSV the leia
workflow aggregates: one row per time step, joined on TIME.

Sources, all optional (a missing one leaves its columns blank):

    postProcessing/areaFractionError/<t>/areaFractionError.dat
        Time L1 L2 Linf nFaces nInterfaceCells consistency
        (the transported face fractions against the exact DAVOF state of the
        translated function, over the faces of the interface cells; -1 for
        the schemes without face fractions)
    postProcessing/reconstructionError/<t>/error.dat
        Time ... LNormalDiff1 LNormalDiff2 LNormalDiffInf areaSum areaExact
        nInterface LCentreDiff1 LCentreDiff2 LCentreDiffInf
        (every step with reconstructionErrorEveryStep true, else the end time)
    postProcessing/volumeFractionError/<t>/volumeFractionError.dat
        Time Eshape Emass Ebound ...
    log.advectorVoF
        per step: the FIRST `davof: N interface cells (...), normalModel X
        (...): S slivers, R replaced[, K kept m][, F fit fallbacks]` line after
        `Time = ` (the advection step's reconstruction; the function objects
        reconstruct again), `isoAdvection: After conservative bounding:
        min(alpha) = ..., max(alpha) = 1 + ...`, `t = ..., sum(alpha*V) = ...,
        dev = ...%, 1-max(alpha1) = ..., min(alpha1) = ...`, `ExecutionTime =
        ... ClockTime = ... timeStepTime = ...`; and the abort guard's
        `advectorVoF: ABORT at step <i>, t = <t>: <reason> (...)` before `End`.

The latest time directory under each postProcessing/<fo>/ is read (a restart
writes a new one); rows are de-duplicated on the rounded Time, the last one
kept. The join is an OUTER join on TIME (9 significant digits): a row of any
source makes a row of the CSV, blank where a source has no row; `step` is the
rank of the time in the union.

Columns, in this order:

    TIME step L1_alphaf L2_alphaf Linf_alphaf nFaces nInterfaceCells consistency
    LNormalDiff1 LNormalDiff2 LNormalDiffInf LCentreDiff1 LCentreDiff2 LCentreDiffInf
    areaSum areaExact nInterface Eshape Emass Ebound
    nInterfaceAdv nSliver nReplaced nKept nFitFallback minAlphaBound maxAlphaBound
    sumAlphaV devPercent oneMinusMaxAlpha minAlpha executionTime clockTime timeStepTime
    status abortStep abortReason

status (constant over the rows): `diverged` when the ABORT line exists,
`completed` on `End` without it, `crashed` when .leia_solver_failed exists
(the Snakefile's record of a non-zero exit), else `incomplete`. abortStep and
abortReason are blank unless the guard fired.

Exit 0 on partial data (a diverged or crashed run is a result); exit 1 only
when neither a .dat file nor the log exists.

    python3 tpf_advection_metrics.py --case . --out advectorVoF.csv
"""
import argparse
import csv
import glob
import os
import re
import sys

AF_COLS = {"L1": "L1_alphaf", "L2": "L2_alphaf", "Linf": "Linf_alphaf",
           "nFaces": "nFaces", "nInterfaceCells": "nInterfaceCells",
           "consistency": "consistency"}
REC_COLS = {k: k for k in ("LNormalDiff1", "LNormalDiff2", "LNormalDiffInf",
                           "LCentreDiff1", "LCentreDiff2", "LCentreDiffInf",
                           "areaSum", "areaExact", "nInterface")}
VOL_COLS = {k: k for k in ("Eshape", "Emass", "Ebound")}
LOG_COLS = ("nInterfaceAdv", "nSliver", "nReplaced", "nKept", "nFitFallback",
            "minAlphaBound", "maxAlphaBound", "sumAlphaV", "devPercent",
            "oneMinusMaxAlpha", "minAlpha", "executionTime", "clockTime",
            "timeStepTime")
COLUMNS = (["TIME", "step"] + list(AF_COLS.values()) + list(REC_COLS.values())
           + list(VOL_COLS.values()) + list(LOG_COLS)
           + ["status", "abortStep", "abortReason"])

_RE_DAVOF = re.compile(
    r"^davof: (\d+) interface cells \(wispTol [^)]*\), normalModel \S+ \([^)]*\): "
    r"(\d+) slivers, (\d+) replaced(?:, (\d+) kept m)?(?:, (\d+) fit fallbacks)?")
_RE_BOUND = re.compile(
    r"^isoAdvection: After conservative bounding: min\(alpha\) = (\S+), "
    r"max\(alpha\) = 1 \+ (\S+)")
_RE_T = re.compile(
    r"^t = ([^,]+),\s*sum\(alpha\*V\) = ([^,]+),\s*dev = ([^%]+)%,\s*"
    r"1-max\(alpha1\) = (\S+) at cell \S+,\s*min\(alpha1\) = (\S+) at cell")
_RE_EXEC = re.compile(
    r"^ExecutionTime = (\S+) s\s+ClockTime = (\S+) s\s+timeStepTime = (\S+) s")
_RE_ABORT = re.compile(
    r"^advectorVoF: ABORT at step (\d+), t = (\S+): (.*?) \(min\(alpha\)")


def tkey(t):
    """The join key of a time: 9 significant digits, so that the .dat files
    (writePrecision 12 or 16) and the log (the same precision) agree."""
    return f"{float(t):.9g}"


def read_dat(fo_dir, name):
    """The rows of <fo_dir>/<latest time>/<name> as {column: string}, keyed by
    the time key, the last row of a repeated time kept."""
    tdirs = [d for d in glob.glob(os.path.join(fo_dir, "*")) if os.path.isdir(d)]
    out = {}
    if not tdirs:
        return out, False
    try:
        latest = max(tdirs, key=lambda d: float(os.path.basename(d)))
    except ValueError:
        latest = sorted(tdirs)[-1]
    path = os.path.join(latest, name)
    if not os.path.isfile(path):
        return out, False
    names = None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                parts = s[1:].split()
                if parts and parts[0] == "Time":
                    names = parts
                continue
            if names is None:
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                k = tkey(parts[0])
            except ValueError:
                continue
            out[k] = dict(zip(names, parts))
    return out, True


def read_log(path):
    """Per-step log values keyed by the time key, plus the run status."""
    steps = {}      # tkey -> dict of LOG_COLS
    cur = None      # the values collected since the last 'Time = ' line
    abort = None
    ended = False
    if not os.path.isfile(path):
        return steps, abort, ended, False
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                cur = {}
                continue
            # The guard's line and the final "End" come outside a step's block
            # (after the ExecutionTime line of the last step, or in place of it).
            m = _RE_ABORT.match(line)
            if m:
                abort = (m.group(1), m.group(2), m.group(3))
                if cur is not None and "_t" in cur:
                    steps[tkey(cur["_t"])] = cur
                continue
            if line.startswith("End"):
                ended = True
                continue
            if cur is None:
                continue
            m = _RE_DAVOF.match(line)
            if m and "nInterfaceAdv" not in cur:
                cur["nInterfaceAdv"] = m.group(1)
                cur["nSliver"] = m.group(2)
                cur["nReplaced"] = m.group(3)
                cur["nKept"] = m.group(4) or "0"
                cur["nFitFallback"] = m.group(5) or "0"
                continue
            m = _RE_BOUND.match(line)
            if m:
                cur["minAlphaBound"] = m.group(1)
                cur["maxAlphaBound"] = m.group(2)
                continue
            m = _RE_T.match(line)
            if m:
                cur["sumAlphaV"] = m.group(2)
                cur["devPercent"] = m.group(3).strip()
                cur["oneMinusMaxAlpha"] = m.group(4)
                cur["minAlpha"] = m.group(5)
                cur["_t"] = m.group(1).strip()
                continue
            m = _RE_EXEC.match(line)
            if m:
                cur["executionTime"] = m.group(1)
                cur["clockTime"] = m.group(2)
                cur["timeStepTime"] = m.group(3)
                if "_t" in cur:
                    steps[tkey(cur["_t"])] = cur
                cur = None
                continue
    return steps, abort, ended, True


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--case", default=".", help="case directory (default .)")
    ap.add_argument("--out", default="advectorVoF.csv", help="output CSV, relative to --case")
    ap.add_argument("--log", default="log.advectorVoF", help="solver log, relative to --case")
    a = ap.parse_args(argv)
    case = a.case
    pp = os.path.join(case, "postProcessing")
    af, have_af = read_dat(os.path.join(pp, "areaFractionError"), "areaFractionError.dat")
    rec, have_rec = read_dat(os.path.join(pp, "reconstructionError"), "error.dat")
    vol, have_vol = read_dat(os.path.join(pp, "volumeFractionError"), "volumeFractionError.dat")
    log, abort, ended, have_log = read_log(os.path.join(case, a.log))
    if not (have_af or have_rec or have_vol or have_log):
        print(f"[tpf_advection_metrics] no .dat file under {pp} and no {a.log}: nothing to convert",
              file=sys.stderr)
        return 1

    if abort is not None:
        status = "diverged"
    elif ended:
        status = "completed"
    elif os.path.isfile(os.path.join(case, ".leia_solver_failed")):
        status = "crashed"
    else:
        status = "incomplete"
    abort_step = abort[0] if abort else ""
    abort_reason = abort[2] if abort else ""

    keys = set(af) | set(rec) | set(vol) | set(log)
    keys = sorted(keys, key=float)
    rows = []
    for step, k in enumerate(keys, start=1):
        row = {c: "" for c in COLUMNS}
        row["TIME"] = k
        row["step"] = step
        for src, cols in ((af, AF_COLS), (rec, REC_COLS), (vol, VOL_COLS)):
            r = src.get(k)
            if r:
                for src_col, out_col in cols.items():
                    row[out_col] = r.get(src_col, "")
        r = log.get(k)
        if r:
            for c in LOG_COLS:
                row[c] = r.get(c, "")
        row["status"] = status
        row["abortStep"] = abort_step
        row["abortReason"] = abort_reason
        rows.append(row)

    out = os.path.join(case, a.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"[tpf_advection_metrics] wrote {out}: {len(rows)} steps, status {status}"
          + (f" (abort at step {abort_step}: {abort_reason})" if abort else "")
          + f"; sources: areaFractionError {len(af)}, reconstructionError {len(rec)}, "
            f"volumeFractionError {len(vol)}, log {len(log)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
