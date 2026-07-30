#!/usr/bin/env python3
"""Converge the corrected pressure flux on perturbed exact-circle meshes.

The sweep changes only PIMPLE.nNonOrthogonalCorrectors.  Its decision
observable is the maximum written cell-centred velocity after one frozen
capillary step; force norms are deliberately excluded.
"""

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
    sweep = config["nonorthogonal_correction_sweep"]
    study_dir = pressure_gate.checked_study_dir(sweep["study_name"])
    if args.fresh and study_dir.exists():
        shutil.rmtree(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)

    tokens = pressure_gate.base_tokens(config)
    rows: list[dict] = []
    case_index = 0
    for resolution in sweep["resolutions"]:
        delta_t = float(tokens["CAPILLARY_DT_COEFF"])/float(resolution)**1.5
        for correction_count in sweep["correction_counts"]:
            case_dir = study_dir / f"N{resolution}_nNonOrth{correction_count}"
            variant_tokens = dict(tokens)
            variant_tokens.update(
                {
                    "N_CELLS": str(resolution),
                    "END_TIME": f"{delta_t:.12g}",
                    "WRITE_INTERVAL": f"{delta_t:.12g}",
                    "SURFACE_TENSION_FORCE": sweep["surface_tension_force"],
                    "N_NON_ORTHOGONAL_CORRECTORS": str(correction_count),
                }
            )
            pressure_gate.make_case(
                case_dir, variant_tokens, f"{case_index:05d}"
            )
            case_index += 1
            shutil.copytree(case_dir / "0.org", case_dir / "0")

            mesh_status = pressure_gate.build_mesh(case_dir, "perturbed", sweep)
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
            solve_status = (
                pressure_gate.run_command(
                    [pressure_gate.SOLVER], case_dir, "log.solve", env=env,
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
                    "nNonOrthogonalCorrectors": correction_count,
                    "delta_t_s": delta_t,
                    "solver_status": solve_status,
                    "max_cell_velocity_m_per_s": max_u,
                    "written_time_s": time_value,
                }
            )
            print(
                f"[nonorth-sweep] N={resolution} "
                f"nNonOrth={correction_count}: status={solve_status}, "
                f"max|U|={max_u:.12g} m/s"
            )

    result = study_dir / "nonorthogonal_correction_sweep.csv"
    publication = (REPO / sweep["publication_table"]).resolve()
    write_rows(result, rows)
    write_rows(publication, rows)
    print(f"[nonorth-sweep] result: {result}")
    print(f"[nonorth-sweep] publication: {publication}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
