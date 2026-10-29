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
        (...): S slivers, R replaced[, K kept m][, F fit fallbacks][, G ring
        fallbacks]` line after
        `Time = ` (the advection step's reconstruction; the function objects
        reconstruct again), the FIRST `davof: curvature tensor on N cells, F
        tensor fallbacks` line (the planeSagitta face clipping, 2026-10-15),
        the face update's `isoAdvection: alphaf set by a plane on ...` line
        with its optional clause `; sagitta corrected F faces, C capped, K
        clamped, S skipped (cos), M skipped (chord)` (the same model),
        `isoAdvection: After conservative bounding:
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
    nInterfaceAdv nSliver nReplaced nKept nFitFallback nRingFallback
    nTensorCells nTensorFallback nSagittaCorrected nSagittaCapped nSagittaClamped
    nSagittaSkipCos nSagittaSkipChord minAlphaBound maxAlphaBound
    sumAlphaV devPercent oneMinusMaxAlpha minAlpha executionTime clockTime timeStepTime
    status abortStep abortReason tpfCommit tpfVersion

(the tensor and sagitta counters are blank for the rows whose face clipping
is plane: the lines do not exist in their logs)

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
import subprocess
import sys

AF_COLS = {"L1": "L1_alphaf", "L2": "L2_alphaf", "Linf": "Linf_alphaf",
           "nFaces": "nFaces", "nInterfaceCells": "nInterfaceCells",
           "consistency": "consistency"}
REC_COLS = {k: k for k in ("LNormalDiff1", "LNormalDiff2", "LNormalDiffInf",
                           "LCentreDiff1", "LCentreDiff2", "LCentreDiffInf",
                           "areaSum", "areaExact", "nInterface")}
VOL_COLS = {k: k for k in ("Eshape", "Emass", "Ebound")}
LOG_COLS = ("nInterfaceAdv", "nSliver", "nReplaced", "nKept", "nFitFallback",
            "nRingFallback",
            "nTensorCells", "nTensorFallback", "nSagittaCorrected", "nSagittaCapped",
            "nSagittaClamped", "nSagittaSkipCos", "nSagittaSkipChord",
            "clippingSurface", "nCurvedFaces", "nCurvedSecondRing", "nCurvedSkipped",
            "nCurvedClips", "nCurvedLinearCuts", "nCurvedLeaves", "nFitCells",
            "nFitConverged", "nFitStalled", "nFitSingular", "nFitNoBracket",
            "fitResMax", "fitResMean", "nAnchorFail", "anchorShiftMax", "anchorShiftMean",
            "nAnchorEvals", "nAnchorEdge",
            "defectMax", "defectL1", "volumeRule", "nRuleFallback",
            "nSweptFaces", "nSweptFlips", "nSweptOutOfDonor", "nSweptClipTris",
            "sweptClosureMax",
            "nSweptPieces", "nSweptBeyond", "sweptBeyondMax", "nSweptRankFallback",
            "nSweptPartialCells", "nSweptDegenerate", "sweptPieceClosureMax",
            "nClassDonor", "nClassInterface", "nClassFull", "nClassEmpty", "nClassOutside",
            "massAdvRel", "massBoundRel", "massClipRel", "massBoundaryRel",
            "minAlphaBound", "maxAlphaBound", "sumAlphaV", "devPercent",
            "oneMinusMaxAlpha", "minAlpha", "executionTime", "clockTime",
            "timeStepTime")
COLUMNS = (["TIME", "step"] + list(AF_COLS.values()) + list(REC_COLS.values())
           + list(VOL_COLS.values()) + list(LOG_COLS)
           + ["status", "abortStep", "abortReason", "tpfCommit", "tpfVersion"])

_RE_DAVOF = re.compile(
    r"^davof: (\d+) interface cells \(wispTol [^)]*\), normalModel \S+ \([^)]*\): "
    r"(\d+) slivers, (\d+) replaced(?:, (\d+) kept m)?(?:, (\d+) fit fallbacks)?"
    r"(?:, (\d+) ring fallbacks)?")
_RE_TENSOR = re.compile(
    r"^davof: curvature tensor on (\d+) cells, (\d+) tensor fallbacks")
_RE_SAGITTA = re.compile(
    r"; sagitta corrected (\d+) faces, (\d+) capped, (\d+) clamped, "
    r"(\d+) skipped \(cos\), (\d+) skipped \(chord\)")
_RE_CURVED = re.compile(
    r"^isoAdvection: curved clipping \((\w+)\): alphaf and dVf set on (\d+) faces "
    r"\((\d+) second ring, (\d+) boundary\), (\d+) skipped, (\d+) clips, "
    r"(\d+) linear cuts, (\d+) quadrature leaves"
    # the residual and the face shift are printed only with the fork's
    # anchorDiagnostics on (davof-perf, 2026-10-25), the anchor's evaluation
    # count since the Newton anchor; both clauses optional, the groups kept
    r"(?:; fit: cells (\d+), converged (\d+), stalled (\d+), singular (\d+), "
    r"no bracket (\d+), no frame (\d+)(?:, residual max ([^\s;]+) mean ([^\s;]+))?)?"
    r"(?:; volume anchor: failed (\d+)(?:, face shift max ([^\s;]+) mean ([^\s;]+))?, leaves (\d+)"
    r"(?:, evaluations (\d+) \(per cell [^)]*\))?(?:, edge (\d+))?)?"
    r"(?:; surface defect max ([^\s;]+) L1 ([^\s;]+) over (\d+) cells)?"
    r"(?:; volume rule (\w+) fallbacks (\d+))?"
    # the swept-face flux's diagnostics (fork branch davof-swept, 2026-10-28,
    # fluxIntegration sweptFace): the faces computed, the sign flips of the
    # triangle list, the faces whose pre-image leaves the donor cell, the
    # clipped triangles, the largest closure residual of the lists
    r"(?:; swept faces (\d+), flips (\d+), out of donor (\d+), "
    r"clipped triangles (\d+), closure max ([^\s;]+))?"
    # the pre-image classification (sweptClassify, 2026-10-29): the faces by
    # the class of the cells their pre-image reaches beyond the donor
    # the piecewise clip (fluxIntegration sweptFaceCells, 2026-10-29)
    r"(?:; cells: pieces (\d+), beyond layer (\d+) \(max ([^)]+)\), rank fallbacks (\d+), "
    r"partial cells (\d+)(?:, degenerate (\d+))?, piece closure max ([^\s;]+))?"
    r"(?:; pre-image classes: donor (\d+), interface (\d+), full (\d+), "
    r"empty (\d+), outside (\d+))?")
# the mass through the stages of the step (the curved loop, 2026-10-29; the
# boundary outflow since fork 78464c0)
_RE_MASS = re.compile(
    r"^isoAdvection: mass change relative to sum\(alpha V\) at t\^n: advection ([^,\s]+), "
    r"bounding ([^,\s]+), clipping ([^,\s]+)(?:, boundary outflow (\S+))?")
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
                cur["nRingFallback"] = m.group(6) or "0"
                continue
            m = _RE_TENSOR.match(line)
            if m and "nTensorCells" not in cur:
                cur["nTensorCells"] = m.group(1)
                cur["nTensorFallback"] = m.group(2)
                continue
            if line.startswith("isoAdvection: alphaf set by a plane"):
                m = _RE_SAGITTA.search(line)
                if m:
                    cur["nSagittaCorrected"] = m.group(1)
                    cur["nSagittaCapped"] = m.group(2)
                    cur["nSagittaClamped"] = m.group(3)
                    cur["nSagittaSkipCos"] = m.group(4)
                    cur["nSagittaSkipChord"] = m.group(5)
                continue
            m = _RE_CURVED.match(line)
            if m:
                cur["clippingSurface"] = m.group(1)
                cur["nCurvedFaces"] = m.group(2)
                cur["nCurvedSecondRing"] = m.group(3)
                cur["nCurvedSkipped"] = m.group(5)
                cur["nCurvedClips"] = m.group(6)
                cur["nCurvedLinearCuts"] = m.group(7)
                cur["nCurvedLeaves"] = m.group(8)
                if m.group(9) is not None:
                    cur["nFitCells"] = m.group(9)
                    cur["nFitConverged"] = m.group(10)
                    cur["nFitStalled"] = m.group(11)
                    cur["nFitSingular"] = m.group(12)
                    cur["nFitNoBracket"] = m.group(13)
                    cur["fitResMax"] = m.group(15) or ""
                    cur["fitResMean"] = m.group(16) or ""
                if m.group(17) is not None:
                    cur["nAnchorFail"] = m.group(17)
                    cur["anchorShiftMax"] = m.group(18) or ""
                    cur["anchorShiftMean"] = m.group(19) or ""
                    cur["nAnchorEvals"] = m.group(21) or ""
                    cur["nAnchorEdge"] = m.group(22) or ""
                if m.group(23) is not None:
                    cur["defectMax"] = m.group(23)
                    cur["defectL1"] = m.group(24)
                if m.group(26) is not None:
                    cur["volumeRule"] = m.group(26)
                    cur["nRuleFallback"] = m.group(27)
                if m.group(28) is not None:
                    cur["nSweptFaces"] = m.group(28)
                    cur["nSweptFlips"] = m.group(29)
                    cur["nSweptOutOfDonor"] = m.group(30)
                    cur["nSweptClipTris"] = m.group(31)
                    cur["sweptClosureMax"] = m.group(32)
                if m.group(33) is not None:
                    cur["nSweptPieces"] = m.group(33)
                    cur["nSweptBeyond"] = m.group(34)
                    cur["sweptBeyondMax"] = m.group(35)
                    cur["nSweptRankFallback"] = m.group(36)
                    cur["nSweptPartialCells"] = m.group(37)
                    cur["nSweptDegenerate"] = m.group(38) or ""
                    cur["sweptPieceClosureMax"] = m.group(39)
                if m.group(40) is not None:
                    cur["nClassDonor"] = m.group(40)
                    cur["nClassInterface"] = m.group(41)
                    cur["nClassFull"] = m.group(42)
                    cur["nClassEmpty"] = m.group(43)
                    cur["nClassOutside"] = m.group(44)
                continue
            m = _RE_MASS.match(line)
            if m:
                cur["massAdvRel"] = m.group(1)
                cur["massBoundRel"] = m.group(2)
                cur["massClipRel"] = m.group(3)
                cur["massBoundaryRel"] = m.group(4) or ""
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


def tpf_provenance():
    """(short git commit of the TwoPhaseFlow fork, its platforms version stamp), from
    LEIA_TPF_DIR (exported by etc/leia-env.sh); "unknown" where not available. The solver
    binary the case ran is the fork's build at that commit (the clone-local platforms/)."""
    tpf = os.environ.get("LEIA_TPF_DIR", "")
    commit, version = "unknown", "unknown"
    if tpf and os.path.isdir(tpf):
        try:
            commit = subprocess.run(["git", "-C", tpf, "rev-parse", "--short", "HEAD"],
                                    capture_output=True, text=True, timeout=20).stdout.strip() or "unknown"
            dirty = subprocess.run(["git", "-C", tpf, "status", "--porcelain", "--untracked-files=no",
                                    "--", "src", "apps"], capture_output=True, text=True, timeout=20).stdout.strip()
            if dirty:
                commit += "-dirty"
        except Exception:  # noqa: BLE001 - provenance must never fail the converter
            pass
        opts = os.environ.get("WM_OPTIONS", "")
        stamp = os.path.join(tpf, "platforms", opts, ".openfoam-version")
        if opts and os.path.isfile(stamp):
            version = open(stamp).read().strip()
    return commit, version


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
    tpf_commit, tpf_version = tpf_provenance()
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
        row["tpfCommit"] = tpf_commit
        row["tpfVersion"] = tpf_version
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
