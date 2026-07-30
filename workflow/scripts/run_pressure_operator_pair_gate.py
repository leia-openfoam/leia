#!/usr/bin/env python3
"""Compare internally paired corrected and uncorrected pressure operators.

Both the capillary pressure-potential flux and pressure Laplacian use the same
snGrad correction keyword in each variant.  The sole decision observable is
the maximum written cell-centred velocity after one frozen capillary step.
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
    gate = config["pressure_operator_pair_gate"]
    study_dir = pressure_gate.checked_study_dir(gate["study_name"])
    if args.fresh and study_dir.exists():
        shutil.rmtree(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)

    tokens = pressure_gate.base_tokens(config)
    rows: list[dict] = []
    case_index = 0
    for resolution in gate["resolutions"]:
        delta_t = float(tokens["CAPILLARY_DT_COEFF"])/float(resolution)**1.5
        fields: dict[str, list[tuple[float, float, float]]] = {}
        metadata: dict[str, dict] = {}

        for variant in gate["variants"]:
            correction_count = int(variant["nonorthogonal_correctors"])
            if correction_count < pressure_gate.MIN_PERTURBED_NONORTHOGONAL_CORRECTORS:
                raise SystemExit(
                    "pressure-operator pair requires at least "
                    f"{pressure_gate.MIN_PERTURBED_NONORTHOGONAL_CORRECTORS} "
                    "non-orthogonal correctors on perturbed meshes"
                )
            case_dir = study_dir / f"N{resolution}_{variant['id']}"
            variant_tokens = dict(tokens)
            variant_tokens.update(
                {
                    "N_CELLS": str(resolution),
                    "END_TIME": f"{delta_t:.12g}",
                    "WRITE_INTERVAL": f"{delta_t:.12g}",
                    "SURFACE_TENSION_FORCE": gate["surface_tension_force"],
                    "N_NON_ORTHOGONAL_CORRECTORS": str(correction_count),
                    "PRESSURE_NONORTHOGONAL_SCHEME": variant[
                        "pressure_nonorthogonal_scheme"
                    ],
                }
            )
            pressure_gate.make_case(
                case_dir, variant_tokens, f"{case_index:05d}"
            )
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
                fields[variant["id"]] = vectors
                max_u = pressure_gate.max_magnitude(vectors)
            else:
                time_value = math.nan
                max_u = math.nan
            metadata[variant["id"]] = {
                "status": solve_status,
                "time": time_value,
                "max_u": max_u,
                "scheme": variant["pressure_nonorthogonal_scheme"],
                "correctors": correction_count,
            }
            print(
                f"[pressure-operator-pair] N={resolution} {variant['id']}: "
                f"status={solve_status}, max|U|={max_u:.12g} m/s"
            )

        corrected_id = "corrected"
        uncorrected_id = "uncorrected"
        if corrected_id not in fields or uncorrected_id not in fields:
            raise SystemExit("a paired pressure-operator solve failed")
        max_difference, difference_cell = pressure_gate.max_velocity_difference(
            fields[corrected_id], fields[uncorrected_id]
        )
        rows.append(
            {
                "geometry": "exactCircle",
                "mesh_variant": "perturbed",
                "N_CELLS": resolution,
                "delta_t_s": delta_t,
                "surface_tension_force": gate["surface_tension_force"],
                "corrected_nNonOrthogonalCorrectors": metadata[corrected_id][
                    "correctors"
                ],
                "uncorrected_nNonOrthogonalCorrectors": metadata[uncorrected_id][
                    "correctors"
                ],
                "corrected_solver_status": metadata[corrected_id]["status"],
                "uncorrected_solver_status": metadata[uncorrected_id]["status"],
                "corrected_max_cell_velocity_m_per_s": metadata[corrected_id][
                    "max_u"
                ],
                "uncorrected_max_cell_velocity_m_per_s": metadata[uncorrected_id][
                    "max_u"
                ],
                "max_cell_velocity_difference_m_per_s": max_difference,
                "max_difference_cell_index": difference_cell,
            }
        )
        print(
            f"[pressure-operator-pair] N={resolution} "
            f"max|U_corrected-U_uncorrected|={max_difference:.12g} m/s"
        )

    result = study_dir / "pressure_operator_pair_gate.csv"
    publication = (REPO / gate["publication_table"]).resolve()
    write_rows(result, rows)
    write_rows(publication, rows)
    print(f"[pressure-operator-pair] result: {result}")
    print(f"[pressure-operator-pair] publication: {publication}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
