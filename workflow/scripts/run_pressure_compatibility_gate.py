#!/usr/bin/env python3
"""Run the frozen exact-circle capillary pressure-compatibility gate.

The sole decision observable is the maximum cell-centred velocity produced by
the complete PIMPLE solve.  No surface-force norm is used as a pass/fail proxy.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import foam_param  # noqa: E402
import materialize  # noqa: E402


SOLVER = "leiaSemiLagrangianLevelSetTwoPhaseFoam"
MIN_PERTURBED_NONORTHOGONAL_CORRECTORS = 8


def checked_study_dir(name: str) -> Path:
    studies = (REPO / "studies").resolve()
    target = (studies / name).resolve()
    if target.parent != studies or not name.startswith("staticISTPressure"):
        raise SystemExit(f"refusing unsafe pressure-gate study path: {target}")
    return target


def run_command(
    command: list[str],
    cwd: Path,
    log_name: str,
    env: dict[str, str] | None = None,
    timeout: int = 360,
) -> int:
    with (cwd / log_name).open("w", encoding="utf-8") as log:
        try:
            result = subprocess.run(
                command,
                cwd=cwd,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            return result.returncode
        except subprocess.TimeoutExpired:
            log.write(f"\npressure-gate timeout after {timeout} s\n")
            return 124


def require_commands() -> None:
    required = ("blockMesh", "leiaPerturbMesh", "leiaSetFields", SOLVER)
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        raise SystemExit(
            "OpenFOAM environment is not active; missing: " + ", ".join(missing)
        )


def base_tokens(config: dict) -> dict[str, str]:
    case_name = config["case"]
    base_case = REPO / "cases" / case_name
    axes, constants, ignored = foam_param.build_token_grid(
        str(REPO / "cases" / f"{case_name}.parameter"),
        str(REPO / "cases" / "default.parameter"),
        str(base_case),
        axes_override=config.get("axes_override"),
        collapse_other_axes=True,
    )
    variations = foam_param.expand(axes, constants)
    if ignored:
        print(f"[pressure-gate] ignored parameters: {ignored}")
    if len(variations) != 1:
        raise SystemExit(
            f"pressure gate must materialize one base case, got {len(variations)}"
        )
    return dict(variations[0])


def make_case(case_dir: Path, tokens: dict[str, str], index: str) -> None:
    materialize.materialize(
        str(REPO / "cases" / "oscISTDroplet2D"),
        tokens,
        str(case_dir),
        1,
        "hex",
        "serial",
        2,
        "oscISTDroplet2D",
        index,
    )


def build_mesh(case_dir: Path, mesh_variant: str, gate: dict) -> int:
    status = run_command(["blockMesh"], case_dir, "log.blockMesh")
    if status or mesh_variant == "uniform":
        return status

    status = run_command(
        [
            "leiaPerturbMesh",
            "-alpha",
            str(gate["perturbation_alpha"]),
            "-seed",
            str(gate["perturbation_seed"]),
        ],
        case_dir,
        "log.perturbMesh",
    )
    if status:
        return status

    time_mesh = case_dir / "0" / "polyMesh"
    constant_mesh = case_dir / "constant" / "polyMesh"
    if not time_mesh.is_dir():
        raise SystemExit(f"perturbed mesh output missing: {time_mesh}")
    for source in time_mesh.iterdir():
        destination = constant_mesh / source.name
        if source.is_dir():
            if destination.exists():
                shutil.rmtree(destination)
            shutil.copytree(source, destination)
        else:
            shutil.copy2(source, destination)
    shutil.rmtree(time_mesh)
    return 0


def numeric_time_dirs(case_dir: Path) -> list[tuple[float, Path]]:
    found: list[tuple[float, Path]] = []
    for path in case_dir.iterdir():
        if path.is_dir():
            try:
                found.append((float(path.name), path))
            except ValueError:
                pass
    return sorted(found)


def read_internal_vectors(path: Path) -> list[tuple[float, float, float]]:
    text = path.read_text(encoding="utf-8", errors="strict")
    match = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s+"
        r"(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise ValueError(f"cannot parse nonuniform internal vector field: {path}")
    expected = int(match.group(1))
    values = [
        tuple(float(item) for item in vector.split())
        for vector in re.findall(r"\(([^()]+)\)", match.group(2))
    ]
    if len(values) != expected or any(len(value) != 3 for value in values):
        raise ValueError(
            f"parsed {len(values)} vectors from {path}, expected {expected}"
        )
    return values


def latest_velocity(case_dir: Path) -> tuple[float, list[tuple[float, float, float]]]:
    candidates = [item for item in numeric_time_dirs(case_dir) if item[0] > 0]
    if not candidates:
        raise ValueError(f"no written positive-time field in {case_dir}")
    time_value, time_dir = candidates[-1]
    return time_value, read_internal_vectors(time_dir / "U")


def max_magnitude(vectors: list[tuple[float, float, float]]) -> float:
    return max(
        math.sqrt(sum(component*component for component in value))
        for value in vectors
    )


def max_velocity_difference(
    left: list[tuple[float, float, float]],
    right: list[tuple[float, float, float]],
) -> tuple[float, int]:
    if len(left) != len(right):
        raise ValueError("pressure-gate velocity fields have different sizes")
    differences = [
        math.sqrt(sum((a - b)**2 for a, b in zip(lhs, rhs)))
        for lhs, rhs in zip(left, right)
    ]
    index = max(range(len(differences)), key=differences.__getitem__)
    return differences[index], index


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

    require_commands()
    config_path = Path(args.config).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    study_dir = checked_study_dir(config["study_name"])
    if args.fresh and study_dir.exists():
        shutil.rmtree(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)

    tokens = base_tokens(config)
    gate = config["pressure_compatibility_gate"]
    perturbed_nonorth = int(
        gate["nonorthogonal_correctors"]["perturbed"]
    )
    if perturbed_nonorth < MIN_PERTURBED_NONORTHOGONAL_CORRECTORS:
        raise SystemExit(
            "perturbed pressure gate requires at least "
            f"{MIN_PERTURBED_NONORTHOGONAL_CORRECTORS} non-orthogonal "
            "correctors (velocity-converged exact-circle gate)"
        )
    rows: list[dict] = []
    case_index = 0

    for resolution in gate["resolutions"]:
        delta_t = float(tokens["CAPILLARY_DT_COEFF"])/float(resolution)**1.5
        for mesh_variant in gate["mesh_variants"]:
            fields: dict[str, list[tuple[float, float, float]]] = {}
            metadata: dict[str, dict] = {}

            for variant in gate["variants"]:
                case_dir = (
                    study_dir
                    / f"N{resolution}_{mesh_variant}"
                    / variant["id"]
                )
                variant_tokens = dict(tokens)
                variant_tokens.update(
                    {
                        "N_CELLS": str(resolution),
                        "END_TIME": f"{delta_t:.12g}",
                        "WRITE_INTERVAL": f"{delta_t:.12g}",
                        "N_NON_ORTHOGONAL_CORRECTORS": str(
                            gate["nonorthogonal_correctors"][mesh_variant]
                        ),
                        "SURFACE_TENSION_FORCE": variant[
                            "surface_tension_force"
                        ],
                    }
                )
                make_case(case_dir, variant_tokens, f"{case_index:05d}")
                case_index += 1
                shutil.copytree(case_dir / "0.org", case_dir / "0")

                mesh_status = build_mesh(case_dir, mesh_variant, gate)
                set_status = (
                    run_command(
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
                    run_command([SOLVER], case_dir, "log.solve", env=env)
                    if set_status == 0
                    else 125
                )

                if solve_status == 0:
                    time_value, vectors = latest_velocity(case_dir)
                    fields[variant["id"]] = vectors
                    max_u = max_magnitude(vectors)
                else:
                    time_value = math.nan
                    max_u = math.nan
                metadata[variant["id"]] = {
                    "status": solve_status,
                    "time": time_value,
                    "max_u": max_u,
                }
                print(
                    f"[pressure-gate] N={resolution} {mesh_variant} "
                    f"{variant['id']}: status={solve_status}, max|U|={max_u:.12g}"
                )

            csf_id = "constantCurvatureCSF"
            potential_id = "capillaryPressurePotential"
            if csf_id not in fields or potential_id not in fields:
                raise SystemExit("a pressure-compatibility solve failed")
            max_difference, difference_cell = max_velocity_difference(
                fields[csf_id], fields[potential_id]
            )
            rows.append(
                {
                    "geometry": "exactCircle",
                    "radius_m": 1.0e-3,
                    "curvature_per_m": 1000,
                    "N_CELLS": resolution,
                    "mesh_variant": mesh_variant,
                    "delta_t_s": delta_t,
                    "nNonOrthogonalCorrectors": gate[
                        "nonorthogonal_correctors"
                    ][mesh_variant],
                    "interface_frozen": 1,
                    "phase_indicator_geometry": "analyticImplicitSurface",
                    "mass_flux": "rhoLENT-central",
                    "csf_flux": "sigma*kappa*snGrad(alpha)*magSf",
                    "pressure_potential_flux": "snGrad(sigma*kappa*alpha)*magSf",
                    "csf_solver_status": metadata[csf_id]["status"],
                    "pressure_potential_solver_status": metadata[potential_id][
                        "status"
                    ],
                    "csf_max_cell_velocity_m_per_s": metadata[csf_id]["max_u"],
                    "pressure_potential_max_cell_velocity_m_per_s": metadata[
                        potential_id
                    ]["max_u"],
                    "max_cell_velocity_difference_m_per_s": max_difference,
                    "max_difference_cell_index": difference_cell,
                }
            )
            print(
                f"[pressure-gate] N={resolution} {mesh_variant} "
                f"max|U_csf-U_potential|={max_difference:.12g} m/s"
            )

    result = study_dir / "pressure_compatibility_gate.csv"
    publication = (REPO / config["publication_table"]).resolve()
    write_rows(result, rows)
    write_rows(publication, rows)
    print(f"[pressure-gate] result: {result}")
    print(f"[pressure-gate] publication: {publication}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
