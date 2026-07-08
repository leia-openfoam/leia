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

PER_CASE_CSVS = [
    "leiaLevelSetFoam.csv",
    "leiaSemiLagrangeLevelSetFoam.csv",  # semi-Lagrangian solver (same columns)
    "leiaSetFields.csv",
    "gradPsiError.csv",
    "leiaTestVelocityExtension.csv",   # static t=0 extension verification
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
    # Shape/volume errors come from whichever advection solver ran; both write
    # identical column names, so pick the CSV that is present.
    def solver_of(rec):
        if rec.get("leiaSemiLagrangeLevelSetFoam.E_GEOM_ALPHA", "") != "":
            return "leiaSemiLagrangeLevelSetFoam"
        return "leiaLevelSetFoam"

    def sget(rec, key):
        s = solver_of(rec)
        return rec.get(f"{s}.{key}", "")

    def shget(rec, key):
        s = solver_of(rec)
        return rec.get(f"half.{s}.{key}", "")

    cols = ["velocityExtension", "reconstruction", "solver", "phaseIndicator",
            "T", "h", "cfl",
            "gradientError", "shapeError", "volumeError",
            # Band-restricted gradient error + the same metrics at maximal
            # deformation t = T/2 (before any reversal cancellation) + the
            # narrow-band min |grad psi| flattening diagnostic.
            "gradientErrorBand", "gradientErrorHalf", "gradientErrorBandHalf",
            "volumeErrorHalf", "minGradPsiBand", "minGradPsiBandHalf",
            # Static t=0 extension verification (leiaTestVelocityExtension;
            # blank for advection studies).
            "anchorLayers", "uextDiv", "eNormalL2", "eNormalLinf",
            "eNormalRawL2", "ratioL2", "eNormalL2In", "eNormalL2Out"]
    rows = []
    for rec in records:
        n = _num(rec.get("N_CELLS"))
        rows.append({
            "velocityExtension": rec.get("VELOCITY_EXTENSION", "none"),
            # semi-Lagrangian reconstruction (blank for the Eulerian solver).
            "reconstruction": rec.get("SL_RECONSTRUCTION", "")
            if solver_of(rec) == "leiaSemiLagrangeLevelSetFoam" else "",
            "solver": solver_of(rec),
            "phaseIndicator": rec.get("PHASE_INDICATOR", ""),
            "T": rec.get("END_TIME", ""),                 # oscillation period / end time
            "h": (1.0 / n) if n else "",                  # domain length 1 / cells
            "cfl": rec.get("CFL", ""),                     # max Courant (SL sweeps 0.5/1.0)
            "gradientError": rec.get("gradPsiError.E_L2_GRAD_PSI", ""),
            "shapeError": sget(rec, "E_GEOM_ALPHA"),
            "volumeError": sget(rec, "E_VOL_ALPHA_REL"),
            "gradientErrorBand": rec.get("gradPsiError.E_NARROW_L2_GRAD_PSI", ""),
            "gradientErrorHalf": rec.get("half.gradPsiError.E_L2_GRAD_PSI", ""),
            "gradientErrorBandHalf": rec.get("half.gradPsiError.E_NARROW_L2_GRAD_PSI", ""),
            "volumeErrorHalf": shget(rec, "E_VOL_ALPHA_REL"),
            "minGradPsiBand": rec.get("gradPsiError.NARROW_MIN_MAG_GRAD_PSI", ""),
            "minGradPsiBandHalf": rec.get("half.gradPsiError.NARROW_MIN_MAG_GRAD_PSI", ""),
            "anchorLayers": rec.get("ANCHOR_LAYERS", ""),
            "uextDiv": rec.get("UEXT_DIV", ""),
            "eNormalL2": rec.get("leiaTestVelocityExtension.E_NORMAL_L2", ""),
            "eNormalLinf": rec.get("leiaTestVelocityExtension.E_NORMAL_LINF", ""),
            "eNormalRawL2": rec.get("leiaTestVelocityExtension.E_NORMAL_RAW_L2", ""),
            "ratioL2": rec.get("leiaTestVelocityExtension.RATIO_L2", ""),
            "eNormalL2In": rec.get("leiaTestVelocityExtension.E_NORMAL_L2_IN", ""),
            "eNormalL2Out": rec.get("leiaTestVelocityExtension.E_NORMAL_L2_OUT", ""),
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
            for key in ("case", "index", "mesh", "mode", "np"):
                rec[key] = meta.get(key, "")
            for k, v in meta.get("tokens", {}).items():
                rec[k] = v
        else:
            rec["index"] = os.path.basename(case_dir)
        rec["case_dir"] = os.path.relpath(case_dir, os.path.dirname(out_path) or ".")

        t_end = _num(rec.get("END_TIME"))
        for csv_name in PER_CASE_CSVS:
            path = _find_csv(case_dir, csv_name)
            if not path:
                continue
            prefix = csv_name[:-len(".csv")]
            for k, v in _final_row(path).items():
                if k is not None:
                    rec[f"{prefix}.{k.strip()}"] = (v or "").strip()
            # Also the row at maximal deformation (t = T/2): the reversal has
            # not cancelled anything yet, so these columns carry the honest
            # forward-deformation error of the advection studies.
            if t_end:
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
