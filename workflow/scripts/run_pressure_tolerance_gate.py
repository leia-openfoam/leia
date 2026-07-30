#!/usr/bin/env python3
"""Sweep only pressure algebra tolerance in the frozen exact-circle gate."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import sys

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import run_pressure_compatibility_gate as pressure_gate  # noqa: E402


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        nargs="?",
        default=str(REPO / "config" / "staticISTPressureCompatibilityGate.yaml"),
    )
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()

    pressure_gate.require_commands()
    config_path = Path(args.config).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    gate = config["pressure_tolerance_gate"]
    study_dir = pressure_gate.checked_study_dir(gate["study_name"])
    if args.fresh and study_dir.exists():
        shutil.rmtree(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)

    correction_count = int(gate["nonorthogonal_correctors"])
    if correction_count < pressure_gate.MIN_PERTURBED_NONORTHOGONAL_CORRECTORS:
        raise SystemExit(
            "pressure-tolerance gate requires at least "
            f"{pressure_gate.MIN_PERTURBED_NONORTHOGONAL_CORRECTORS} "
            "non-orthogonal correctors on perturbed meshes"
        )

    tokens = pressure_gate.base_tokens(config)
    rows: list[dict] = []
    case_index = 0
    for resolution in gate["resolutions"]:
        delta_t = float(tokens["CAPILLARY_DT_COEFF"])/float(resolution)**1.5
        for tolerance in gate["pressure_tolerances"]:
            tolerance_value = float(tolerance)
            tolerance_id = f"{tolerance_value:.0e}".replace("-", "m")
            case_dir = study_dir / f"N{resolution}_tol{tolerance_id}"
            variant_tokens = dict(tokens)
            variant_tokens.update(
                {
                    "N_CELLS": str(resolution),
                    "END_TIME": f"{delta_t:.12g}",
                    "WRITE_INTERVAL": f"{delta_t:.12g}",
                    "SURFACE_TENSION_FORCE": gate["surface_tension_force"],
                    "N_NON_ORTHOGONAL_CORRECTORS": str(correction_count),
                    "PRESSURE_NONORTHOGONAL_SCHEME": gate[
                        "pressure_nonorthogonal_scheme"
                    ],
                    "PRESSURE_SOLVER": gate["pressure_solver"],
                    "PRESSURE_TOLERANCE": f"{tolerance_value:.12g}",
                    "PRESSURE_REL_TOL": str(gate["pressure_rel_tol"]),
                }
            )
            pressure_gate.make_case(case_dir, variant_tokens, f"{case_index:05d}")
            case_index += 1
            shutil.copytree(case_dir / "0.org", case_dir / "0")

            mesh_status = pressure_gate.build_mesh(case_dir, "perturbed", gate)
            set_status = (
                pressure_gate.run_command(
                    ["leiaSetFields", "-alphaName", "alpha.water"],
                    case_dir,
                    "log.setFields",
                )
                if mesh_status == 0
                else 125
            )
            env = dict(os.environ)
            env["SL_FREEZE_INTERFACE"] = "1"
            env.pop("LEIA_DIAGNOSTIC_CONSTANT_RAUF", None)
            solve_status = (
                pressure_gate.run_command(
                    [pressure_gate.SOLVER],
                    case_dir,
                    "log.solve",
                    env=env,
                    timeout=900,
                )
                if set_status == 0
                else 125
            )
            if solve_status == 0:
                time_value, vectors = pressure_gate.latest_velocity(case_dir)
                max_u = pressure_gate.max_magnitude(vectors)
            else:
                time_value = math.nan
                max_u = math.nan
            rows.append(
                {
                    "geometry": "exactCircle",
                    "mesh_variant": "perturbed",
                    "N_CELLS": resolution,
                    "delta_t_s": delta_t,
                    "pressure_nonorthogonal_scheme": gate[
                        "pressure_nonorthogonal_scheme"
                    ],
                    "nNonOrthogonalCorrectors": correction_count,
                    "pressure_solver": gate["pressure_solver"],
                    "pressure_tolerance": tolerance_value,
                    "pressure_rel_tol": gate["pressure_rel_tol"],
                    "solver_status": solve_status,
                    "max_cell_velocity_m_per_s": max_u,
                    "written_time_s": time_value,
                }
            )
            print(
                f"[pressure-tolerance] N={resolution} "
                f"tol={tolerance_value:.0e}: status={solve_status}, "
                f"max|U|={max_u:.12g} m/s"
            )

    result = study_dir / "pressure_tolerance_gate.csv"
    publication = (REPO / gate["publication_table"]).resolve()
    write_rows(result, rows)
    write_rows(publication, rows)
    print(f"[pressure-tolerance] result: {result}")
    print(f"[pressure-tolerance] publication: {publication}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
