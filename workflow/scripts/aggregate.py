#!/usr/bin/env python3
"""Join per-case CSVs + ``case_params.json`` into one study database CSV.

For each generated case directory we read its parameter vector from
``case_params.json`` and the final-time row of every per-case CSV the solver /
function objects already write (``leiaLevelSetFoam.csv``, ``leiaSetFields.csv``,
``gradPsiError.csv``), prefixing those columns by their source.  One row per
variation; missing CSVs leave empty cells (the case is still reported).
"""
import csv
import json
import os

import fvschemes    # discretization actually used, read from the rendered case dicts
import method_label  # the single definition of a method label

PER_CASE_CSVS = [
    "leiaLevelSetFoam.csv",
    "leiaRedistancedLevelSetFoam.csv",   # geometric-redistancing solver (same columns)
    "leiaSemiLagrangeLevelSetFoam.csv",  # semi-Lagrangian solver (same columns)
    "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv",  # two-phase droplet metrics
                                         # (TIME,maxMagU,meanMagU,pLaplace)
    "leiaSetFields.csv",
    "gradPsiError.csv",
    "leiaTestVelocityExtension.csv",   # static t=0 extension verification
    "leiaTestMeanCurvature.csv",       # static curvature-accuracy test
    "leiaTestRedistance.csv",          # static redistancing gate (circle)
]


def _find_csv(case_dir, name):
    for cand in (os.path.join(case_dir, name),
                 os.path.join(case_dir, "processor0", name)):
        if os.path.isfile(cand) and os.path.getsize(cand) > 0:
            return cand
    return None


def _final_row(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows[-1] if rows else {}


def _last_time(path):
    """TIME on the last row, or None if the CSV is empty or has no TIME."""
    row = _final_row(path)
    return _num(row.get("TIME")) if row else None


# A RUNNING CASE IS NOT A FINISHED CASE.
# The per-step CSV is created when the solver STARTS, so its presence proves
# nothing, and _final_row() returns whatever last row it happens to contain. For
# a case still in flight -- or one killed by a walltime limit, which leaves no
# .leia_solver_failed marker because nothing exited non-zero -- that row is a
# mid-trajectory state published under the t=T column heading.
#
# Measured: aggregating sdplsOrderAblation while its N=256 rungs were still
# running gave the sourceless arm a "final" band gradient error of 1.92 at
# N=256 against 2.84 at N=181, read off a case sitting at t=0.72 of 1.0. That
# reads as the error turning over at the fine end -- a convergence claim -- when
# it is only an earlier time. Blank those metrics and record why.
_END_TIME_TOL = 0.99


def _reached(t_last, t_target):
    return (t_last is not None and t_target is not None
            and t_last >= _END_TIME_TOL * t_target)


def _half_time_row(path, t_half):
    """Row whose TIME is nearest t_half (the per-step CSVs log every step).
    At maximal deformation (t = T/2) no reversal cancellation has happened yet,
    so these values measure the forward-deformation error the final-time row
    hides. Returns {} if the CSV has no TIME column."""
    best, best_d = {}, None
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            t = _num(row.get("TIME"))
            if t is None:
                return {}
            d = abs(t - t_half)
            if best_d is None or d < best_d:
                best, best_d = row, d
    return best


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _write_error_table(records, database_path):
    """Flat, simple convergence table alongside the full database: one row per
    case holding only the study values — period T, mesh spacing h, and the
    gradient / shape / volume-conservation errors (final time). Plain columns
    (no source prefixes), sorted by (T, h). Blank where a source is absent.
    """
    err_path = database_path[:-len("_database.csv")] + "_errors.csv" \
        if database_path.endswith("_database.csv") else database_path + ".errors.csv"
    # Curated CSVs are copied verbatim into docs/**/data/tables, where the study
    # name survives only in the filename -- and filenames get renamed. Carry it
    # in the data.
    _base = os.path.basename(err_path)
    study = _base[:-len("_errors.csv")] if _base.endswith("_errors.csv") else _base
    # Shape/volume errors come from whichever advection solver ran; both write
    # identical column names, so pick the CSV that is present.
    def solver_of(rec):
        if rec.get("leiaSemiLagrangeLevelSetFoam.E_GEOM_ALPHA", "") != "":
            return "leiaSemiLagrangeLevelSetFoam"
        if rec.get("leiaRedistancedLevelSetFoam.E_GEOM_ALPHA", "") != "":
            return "leiaRedistancedLevelSetFoam"
        if rec.get("leiaTestRedistance.E_LINF_BAND_PSI", "") != "":
            return "leiaTestRedistance"
        return "leiaLevelSetFoam"

    def sget(rec, key):
        s = solver_of(rec)
        return rec.get(f"{s}.{key}", "")

    # Method label: the single definition lives in method_label.py, shared
    # with make_bench_fields_fig.py so the CSV label and the figure filenames
    # can never diverge again.
    def method_of(rec):
        return method_label.method_label(rec)

    def shget(rec, key):
        s = solver_of(rec)
        return rec.get(f"half.{s}.{key}", "")

    cols = ["method", "velocityExtension", "reconstruction", "redistancer",
            "redistTrigger", "solver", "phaseIndicator",
            "T", "h", "maxCellSize", "cfl",
            # Wall-clock cost: total seconds and seconds per unit simulated
            # time (row-count-independent -- reporting may be gated to write
            # times, so per-CSV-row is NOT per-time-step).
            "totalClockTime", "clockPerSimTime",
            "gradientError", "shapeError", "volumeError",
            # Band-restricted gradient error + the same metrics at maximal
            # deformation t = T/2 (before any reversal cancellation) + the
            # narrow-band min |grad psi| flattening diagnostic.
            "gradientErrorBand", "gradientErrorHalf", "gradientErrorBandHalf",
            # GLOBAL max ||grad psi| - 1| -- the far field, which every band
            # metric is blind to. The SDPLS source acts wherever psi is large
            # and the strain does not cancel, amplifying psi in proportion to
            # itself: measured 119.6 at N=256 on the 2D steady vortex, GROWING
            # under refinement (order -1.44). This is the column the mollifier
            # (arXiv:2208.01269 eq. 24) is supposed to fix, so a study of it
            # that does not carry this column cannot show its own result.
            "gradientErrorMax", "gradientErrorMaxHalf",
            "volumeErrorHalf", "minGradPsiBand", "minGradPsiBandHalf",
            # The band MEAN and MAX of |grad psi|. The min alone is a
            # flattening sentinel; a source that relaxes toward a wrong TARGET
            # (sdplsBeta drives |grad psi| -> beta - a) shifts the CENTRE of the
            # band distribution, which the min barely registers. The spread
            # (max - min) is the second half of that test: a target offset moves
            # the mean without widening the spread.
            "meanGradPsiBand", "meanGradPsiBandHalf",
            "maxGradPsiBand", "maxGradPsiBandHalf",
            # Band statistics of the SDPLS normal strain a = n.grad(U).n, so
            # `beta - strain` can be compared against the measured band
            # |grad psi| directly. Blank when no SDPLS source is active.
            "minStrainBand", "meanStrainBand", "maxStrainBand",
            "minStrainBandHalf", "meanStrainBandHalf", "maxStrainBandHalf",
            # The sdplsBeta target. Without it every beta value collapses to the
            # same `method` string at the same h, and make_convergence_table.py
            # would fit ONE regression straight through a beta sweep and report
            # it as an order.
            "sdplsBeta",
            # Static t=0 extension verification (leiaTestVelocityExtension;
            # blank for advection studies).
            "anchorLayers", "uextDiv", "eNormalL2", "eNormalLinf",
            "eNormalRawL2", "ratioL2", "eNormalL2In", "eNormalL2Out",
            # Discretization provenance -- what the solver actually read, not
            # what a token requested. Without these, two runs that differ only
            # in div(phi,psi) are indistinguishable in the curated CSV.
            # Non-empty = the solver exited non-zero for this case (its exit
            # code). The metric columns are then blank BY MEASUREMENT, not by
            # omission -- e.g. `beta` legitimately blows up on the vortex.
            "solverFailed",
            # Provenance. Without these a curated CSV cannot be traced back to
            # the run that produced it -- and that has already gone wrong here:
            # sdplsStability_errors.csv reports the POST-sign-fix band min while
            # the case directories beside it are pre-fix, and
            # benchVortexEulerT2's 30-row errors CSV sits next to a 15-row
            # database from a different run. Nothing in either file said so.
            "study", "caseDir", "np", "gitCommit", "runDate",
            # `on` = reversed flow (cos(pi t/tau)); the interface returns to its
            # initial shape at t=T and the sourceless drift cancels by symmetry,
            # which flatters doing nothing. `off` = steady vortex, sustained
            # strain, no cancellation -- the non-reversing stress test. Which one
            # a row came from decides whether its t=T column means anything, so
            # it belongs in the curated table.
            "oscillation",
            # Did the run actually get to END_TIME? A per-step CSV exists from
            # the first time step, so its presence never proved completion.
            "endTimeReached", "lastTime",
            ] + fvschemes.COLUMNS
    rows = []
    for rec in records:
        n = _num(rec.get("N_CELLS"))
        mcs = _num(rec.get("MAX_CELL_SIZE"))
        # Mesh spacing h for the convergence fit. Hexahedral meshes have a
        # uniform N per direction -> h = 1/N (domain length 1). Polyhedral
        # (cfMesh) meshes have no uniform N (N_CELLS is a dummy pin), so the
        # characteristic length IS the target maxCellSize -> h = maxCellSize.
        if rec.get("mesh") == "poly":
            h = mcs if (mcs and mcs > 0) else ""
        else:
            h = (1.0 / n) if n else ""
        clk = _num(sget(rec, "ELAPSED_CLOCK_TIME"))
        tfin = _num(sget(rec, "TIME"))
        rows.append({
            "method": method_of(rec),
            "totalClockTime": clk if clk is not None else "",
            "clockPerSimTime":
                (clk/tfin) if (clk is not None and tfin) else "",
            "velocityExtension": rec.get("VELOCITY_EXTENSION", "none"),
            # semi-Lagrangian reconstruction (blank for the Eulerian solver).
            "reconstruction": rec.get("SL_RECONSTRUCTION", "")
            if solver_of(rec) == "leiaSemiLagrangeLevelSetFoam" else "",
            # Geometric-redistancing solver / static gate: which redistancer
            # model and trigger produced the row (blank for other solvers).
            "redistancer": rec.get("REDISTANCER", "")
            if solver_of(rec) in ("leiaRedistancedLevelSetFoam",
                                  "leiaTestRedistance") else "",
            "redistTrigger": rec.get("REDIST_TRIGGER", "")
            if solver_of(rec) in ("leiaRedistancedLevelSetFoam",
                                  "leiaTestRedistance") else "",
            "solver": solver_of(rec),
            "phaseIndicator": rec.get("PHASE_INDICATOR", ""),
            "T": rec.get("END_TIME", ""),                 # oscillation period / end time
            "h": h,                                        # 1/N (hex) or maxCellSize (poly)
            "maxCellSize": (mcs if (mcs and mcs > 0) else ""),  # cfMesh target size (poly only)
            "cfl": rec.get("CFL", ""),                     # max Courant (SL sweeps 0.5/1.0)
            "gradientError": rec.get("gradPsiError.E_L2_GRAD_PSI", ""),
            # SHAPE ERROR IS ONLY AN ERROR ON A REVERSED FLOW. E_GEOM_ALPHA
            # compares alpha against alpha0, which errorCalculation.H sets to the
            # INITIAL field (unless a case supplies an explicit alphaEnd). With
            # `oscillation off` the vortex is steady, the interface never returns,
            # and that difference measures how far the interface travelled -- not
            # an error. Publishing it would let make_convergence_table fit and
            # print a meaningless "shape order". The raw value is untouched in
            # <study>_database.csv; only the curated table blanks it, and the
            # `oscillation` column below says why.
            # Judged on the RENDERED value for the same reason as `oscillation`
            # below: a case whose template hardcodes the setting has no token.
            "shapeError": ("" if (rec.get("oscillationRendered")
                                  or rec.get("OSCILLATION", "on")) == "off"
                           else sget(rec, "E_GEOM_ALPHA")),
            # Volume error stays valid either way: it measures conservation
            # against the initial volume, which a divergence-free velocity
            # preserves regardless of reversal. So does the gradient error, which
            # is intrinsic and needs no reference solution at all.
            # PREFER THE RENDERED DICTIONARY over the token. The 3D cases
            # hardcode `oscillation on` with no @!OSCILLATION!@ placeholder, so
            # their config's OSCILLATION is dropped as an unreferenced token and
            # this column would be blank -- which defeats the reversed-metric
            # suppression downstream. fvschemes reads what the solver parsed.
            "oscillation": (rec.get("oscillationRendered")
                            or rec.get("OSCILLATION", "")),
            # 1 = reached END_TIME, 0 = still running or cut short (see
            # _reached()); blank when there is no per-step CSV to judge from.
            "endTimeReached": rec.get("endTimeReached", ""),
            "lastTime": rec.get("lastTime", ""),
            "volumeError": sget(rec, "E_VOL_ALPHA_REL"),
            "gradientErrorBand": rec.get("gradPsiError.E_NARROW_L2_GRAD_PSI", ""),
            "gradientErrorHalf": rec.get("half.gradPsiError.E_L2_GRAD_PSI", ""),
            "gradientErrorMax": rec.get("gradPsiError.E_MAX_GRAD_PSI", ""),
            "gradientErrorMaxHalf": rec.get("half.gradPsiError.E_MAX_GRAD_PSI", ""),
            "gradientErrorBandHalf": rec.get("half.gradPsiError.E_NARROW_L2_GRAD_PSI", ""),
            "volumeErrorHalf": shget(rec, "E_VOL_ALPHA_REL"),
            "minGradPsiBand": rec.get("gradPsiError.NARROW_MIN_MAG_GRAD_PSI", ""),
            "minGradPsiBandHalf": rec.get("half.gradPsiError.NARROW_MIN_MAG_GRAD_PSI", ""),
            "meanGradPsiBand": rec.get("gradPsiError.NARROW_MEAN_MAG_GRAD_PSI", ""),
            "meanGradPsiBandHalf": rec.get("half.gradPsiError.NARROW_MEAN_MAG_GRAD_PSI", ""),
            "maxGradPsiBand": rec.get("gradPsiError.NARROW_MAX_MAG_GRAD_PSI", ""),
            "maxGradPsiBandHalf": rec.get("half.gradPsiError.NARROW_MAX_MAG_GRAD_PSI", ""),
            "minStrainBand": rec.get("gradPsiError.NARROW_MIN_R", ""),
            "meanStrainBand": rec.get("gradPsiError.NARROW_MEAN_R", ""),
            "maxStrainBand": rec.get("gradPsiError.NARROW_MAX_R", ""),
            "minStrainBandHalf": rec.get("half.gradPsiError.NARROW_MIN_R", ""),
            "meanStrainBandHalf": rec.get("half.gradPsiError.NARROW_MEAN_R", ""),
            "maxStrainBandHalf": rec.get("half.gradPsiError.NARROW_MAX_R", ""),
            # Only meaningful for the sdplsBeta source; blank elsewhere so the
            # column does not imply a target the other arms never had.
            "sdplsBeta": rec.get("SDPLS_BETA", "")
            if rec.get("SDPLS_SOURCE", "noSource") == "beta" else "",
            "anchorLayers": rec.get("ANCHOR_LAYERS", ""),
            "uextDiv": rec.get("UEXT_DIV", ""),
            "eNormalL2": rec.get("leiaTestVelocityExtension.E_NORMAL_L2", ""),
            "eNormalLinf": rec.get("leiaTestVelocityExtension.E_NORMAL_LINF", ""),
            "eNormalRawL2": rec.get("leiaTestVelocityExtension.E_NORMAL_RAW_L2", ""),
            "ratioL2": rec.get("leiaTestVelocityExtension.RATIO_L2", ""),
            "eNormalL2In": rec.get("leiaTestVelocityExtension.E_NORMAL_L2_IN", ""),
            "eNormalL2Out": rec.get("leiaTestVelocityExtension.E_NORMAL_L2_OUT", ""),
            "solverFailed": rec.get("solverFailed", ""),
            "study": study,
            "caseDir": rec.get("case_dir", ""),
            "np": rec.get("np", ""),
            "gitCommit": rec.get("gitCommit", ""),
            "runDate": rec.get("runDate", ""),
            **{c: rec.get(c, "") for c in fvschemes.COLUMNS},
        })
    rows.sort(key=lambda r: (r["velocityExtension"], r["phaseIndicator"],
                             _num(r["T"]) or 0.0,
                             r["h"] if isinstance(r["h"], float) else 0.0))
    with open(err_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[aggregate] wrote {err_path}: {len(rows)} rows (T, h, gradient/shape/volume)")
    return err_path


def build_database(case_dirs, out_path):
    # AN EMPTY INPUT MUST NOT PRODUCE AN EMPTY TABLE.
    # Writing zero rows over an existing curated CSV destroys it silently, and
    # the caller most likely mis-globbed rather than genuinely having no cases.
    # Measured: re-aggregating studies/sdplsBetaSweep, whose case directories
    # had never been pulled locally (only its tables had), replaced a complete
    # 28-row table with a header. Refuse instead; a study with no cases is a
    # bug in the caller, not a result.
    case_dirs = list(case_dirs)
    if not case_dirs:
        raise ValueError(
            f"aggregate: no case directories given for {out_path!r}. Refusing "
            f"to write an empty table over an existing one -- check the glob, "
            f"or that the study was pulled with its case directories.")

    records, columns = [], []

    def add_col(c):
        if c not in columns:
            columns.append(c)

    for case_dir in case_dirs:
        rec = {}
        meta_path = os.path.join(case_dir, "case_params.json")
        if os.path.isfile(meta_path):
            with open(meta_path) as fh:
                meta = json.load(fh)
            for key in ("case", "index", "mesh", "mode", "np",
                        "gitCommit", "runDate"):
                rec[key] = meta.get(key, "")
            for k, v in meta.get("tokens", {}).items():
                rec[k] = v
        else:
            rec["index"] = os.path.basename(case_dir)
        rec["case_dir"] = os.path.relpath(case_dir, os.path.dirname(out_path) or ".")

        # The discretization the solver ACTUALLY read, parsed from the rendered
        # system/fvSchemes + system/controlDict. Study tokens are not a reliable
        # record: 2Dvortex hardcodes div(phi,psi) with no DIV token, a DIV token
        # may hold an unexpanded $alias, and gradPsiSdpls/gradUSdpls are not
        # tokens at all. See workflow/scripts/fvschemes.py.
        rec.update(fvschemes.read_discretization(case_dir))

        # A diverged/crashed solve is a RESULT: the Snakefile records the exit
        # code here rather than stalling the study, and the row survives with
        # blank metrics. Blank = the solver exited 0.
        _failed = os.path.join(case_dir, ".leia_solver_failed")
        rec["solverFailed"] = ""
        if os.path.isfile(_failed):
            try:
                rec["solverFailed"] = open(_failed).read().strip() or "1"
            except OSError:
                rec["solverFailed"] = "1"

        t_end = _num(rec.get("END_TIME"))
        rec["endTimeReached"] = ""
        rec["lastTime"] = ""
        for csv_name in PER_CASE_CSVS:
            path = _find_csv(case_dir, csv_name)
            if not path:
                continue
            prefix = csv_name[:-len(".csv")]
            # See _reached(): a per-step CSV exists from the first time step, so
            # only a row that actually reached t_end may fill the t=T columns.
            t_last = _last_time(path)
            if t_last is not None:
                rec["lastTime"] = f"{t_last:.6g}"
                rec["endTimeReached"] = ("1" if _reached(t_last, t_end)
                                         else "0" if t_end else "")
            # Fill the final-time columns only for a run that got there. Still
            # running, or cut short without a non-zero exit, leaves them blank
            # rather than publishing an earlier state as the final one. The
            # half-time block below is judged separately: a run that reached T/2
            # but not T still has a valid T/2 row.
            if not t_end or _reached(t_last, t_end):
                for k, v in _final_row(path).items():
                    if k is not None:
                        rec[f"{prefix}.{k.strip()}"] = (v or "").strip()
            # Also the row at maximal deformation (t = T/2): the reversal has
            # not cancelled anything yet, so these columns carry the honest
            # forward-deformation error of the advection studies. Same rule --
            # a run that never got to T/2 has no T/2 row to report.
            if t_end and _reached(t_last, 0.5*t_end):
                for k, v in _half_time_row(path, 0.5*t_end).items():
                    if k is not None:
                        rec[f"half.{prefix}.{k.strip()}"] = (v or "").strip()

        for c in rec:
            add_col(c)
        records.append(rec)

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for rec in records:
            writer.writerow({c: rec.get(c, "") for c in columns})
    print(f"[aggregate] wrote {out_path}: {len(records)} rows x {len(columns)} cols")
    _write_error_table(records, out_path)
    return out_path
