#!/usr/bin/env python3
"""Separate curvature, geometry/localization and pressure-coupling velocity errors."""

from __future__ import annotations

import argparse
import csv
import json
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
    if (
        target.parent != studies
        or not name.startswith("oscISTVariableCurvatureOracle")
    ):
        raise SystemExit(f"refusing unsafe oracle study path: {target}")
    return target


def run_command(
    command: list[str], cwd: Path, log_name: str,
    env: dict[str, str] | None = None, timeout: int = 360,
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
            log.write(f"\noracle runner timeout after {timeout} s\n")
            return 124


def require_commands() -> None:
    missing = [
        name for name in ("blockMesh", "leiaPerturbMesh", "leiaSetFields", SOLVER)
        if shutil.which(name) is None
    ]
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
        print(f"[velocity-oracle] ignored parameters: {ignored}")
    if len(variations) != 1:
        raise SystemExit(f"oracle gate must materialize one base case, got {len(variations)}")
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


def build_mesh(case_dir: Path, variant: str, config: dict) -> int:
    status = run_command(["blockMesh"], case_dir, "log.blockMesh")
    if status or variant == "uniform":
        return status
    status = run_command(
        [
            "leiaPerturbMesh",
            "-alpha",
            str(config["perturbation_alpha"]),
            "-seed",
            str(config["perturbation_seed"]),
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
        if not path.is_dir():
            continue
        try:
            found.append((float(path.name), path))
        except ValueError:
            continue
    return sorted(found)


def read_internal_vectors(path: Path) -> list[tuple[float, float, float]]:
    text = path.read_text(encoding="utf-8", errors="strict")
    uniform = re.search(
        r"internalField\s+uniform\s+\(([^)]+)\)\s*;", text
    )
    if uniform:
        value = tuple(float(item) for item in uniform.group(1).split())
        raise ValueError(
            f"uniform U in {path} has no stored cell count: {value}"
        )
    match = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s+"
        r"(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise ValueError(f"cannot parse internal vector field: {path}")
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
    return max(math.sqrt(sum(component*component for component in value)) for value in vectors)


def first_solver_metrics(case_dir: Path) -> dict[str, float]:
    path = case_dir / f"{SOLVER}.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    return {
        "min_grad_psi": float(row["minGradPsiBand"]),
        "max_grad_psi": float(row["maxGradPsiBand"]),
        "rho_lent_reset_l1": float(row["rhoLENTResetL1"]),
    }


def velocity_difference(
    candidate: list[tuple[float, float, float]],
    oracle: list[tuple[float, float, float]],
) -> tuple[float, int]:
    if len(candidate) != len(oracle):
        raise ValueError("candidate and oracle velocity fields have different sizes")
    differences = [
        math.sqrt(sum((a - b)**2 for a, b in zip(cand, exact)))
        for cand, exact in zip(candidate, oracle)
    ]
    index = max(range(len(differences)), key=differences.__getitem__)
    return differences[index], index


def topology(case_dir: Path) -> dict[str, int | bool]:
    pattern = re.compile(
        r"connectedInterface curvature: (\d+) elements, (\d+) component\(s\), "
        r"(\d+) open endpoint\(s\), (\d+) fits, (\d+) fallback"
    )
    matches = pattern.findall(
        (case_dir / "log.solve").read_text(encoding="utf-8", errors="replace")
    )
    if not matches:
        return {"valid": False, "elements": 0, "components": 0, "open": -1}
    elements, components, open_ends, fits, fallbacks = map(int, matches[-1])
    return {
        "valid": elements > 0 and components == 1 and open_ends == 0
        and fits == elements and fallbacks == 0,
        "elements": elements,
        "components": components,
        "open": open_ends,
    }


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
        default=str(REPO / "config" / "oscISTVariableCurvatureOracle.yaml"),
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
    gate = config["oracle_gate"]
    perturbed_nonorth = int(
        gate["nonorthogonal_correctors"]["perturbed"]
    )
    if perturbed_nonorth < MIN_PERTURBED_NONORTHOGONAL_CORRECTORS:
        raise SystemExit(
            "perturbed velocity oracle requires at least "
            f"{MIN_PERTURBED_NONORTHOGONAL_CORRECTORS} non-orthogonal "
            "correctors (velocity-converged exact-circle gate)"
        )
    rows: list[dict] = []
    case_index = 0
    for resolution in gate["resolutions"]:
        dt = float(tokens["CAPILLARY_DT_COEFF"])/float(resolution)**1.5
        for mesh_variant in gate["mesh_variants"]:
            fields: dict[str, list[tuple[float, float, float]]] = {}
            metadata: dict[str, dict] = {}
            for variant in gate["variants"]:
                case_dir = (
                    study_dir / f"N{resolution}_{mesh_variant}" / variant["id"]
                )
                variant_tokens = dict(tokens)
                variant_tokens.update(
                    {
                        "N_CELLS": str(resolution),
                        "END_TIME": f"{dt:.12g}",
                        "WRITE_INTERVAL": f"{dt:.12g}",
                        "N_NON_ORTHOGONAL_CORRECTORS": str(
                            gate["nonorthogonal_correctors"][mesh_variant]
                        ),
                        "CURVATURE_ESTIMATOR": variant["curvature_estimator"],
                        "CURVATURE_GEOMETRY_SOURCE": variant[
                            "curvature_geometry_source"
                        ],
                        "PHASE_INDICATOR_GEOMETRY_SOURCE": variant[
                            "phase_indicator_geometry_source"
                        ],
                    }
                )
                make_case(case_dir, variant_tokens, f"{case_index:05d}")
                case_index += 1
                shutil.copytree(
                    case_dir / "0.org", case_dir / "0", dirs_exist_ok=True
                )
                mesh_status = build_mesh(case_dir, mesh_variant, gate)
                set_status = run_command(
                    ["leiaSetFields", "-alphaName", "alpha.water"],
                    case_dir,
                    "log.setFields",
                ) if mesh_status == 0 else 125
                env = dict(os.environ)
                env["SL_FREEZE_INTERFACE"] = "1"
                solve_status = run_command(
                    [SOLVER], case_dir, "log.solve", env=env
                ) if set_status == 0 else 125
                if solve_status == 0:
                    time_value, vectors = latest_velocity(case_dir)
                    fields[variant["id"]] = vectors
                else:
                    time_value = math.nan
                    vectors = []
                metadata[variant["id"]] = {
                    "status": solve_status,
                    "time": time_value,
                    "max_u": max_magnitude(vectors) if vectors else math.nan,
                    "topology": (
                        topology(case_dir)
                        if solve_status == 0 else {"valid": False}
                    ),
                    "metrics": (
                        first_solver_metrics(case_dir)
                        if solve_status == 0 else {}
                    ),
                }
                print(
                    f"[velocity-oracle] N={resolution} {mesh_variant} "
                    f"{variant['id']}: status={solve_status}, "
                    f"max|U|={metadata[variant['id']]['max_u']:.10g}"
                )

            curvature_oracle_id = "analyticCurvatureNumericalGeometry"
            geometry_oracle_id = "analyticGeometryLocalizationOracle"
            candidate_id = "connectedModePreservingCandidate"
            required = (curvature_oracle_id, geometry_oracle_id, candidate_id)
            if any(variant_id not in fields for variant_id in required):
                raise SystemExit("candidate or oracle solve failed")
            candidate_curvature_difference, candidate_curvature_cell = (
                velocity_difference
                (
                    fields[candidate_id],
                    fields[curvature_oracle_id],
                )
            )
            geometry_effect, geometry_effect_cell = velocity_difference(
                fields[geometry_oracle_id],
                fields[curvature_oracle_id],
            )
            candidate_geometry_difference, candidate_geometry_cell = (
                velocity_difference
                (
                    fields[candidate_id],
                    fields[geometry_oracle_id],
                )
            )
            curvature_oracle_max = metadata[curvature_oracle_id]["max_u"]
            geometry_oracle_max = metadata[geometry_oracle_id]["max_u"]
            rows.append(
                {
                    "geometry": "signedDistanceEllipse",
                    "N_CELLS": resolution,
                    "mesh_variant": mesh_variant,
                    "delta_t_s": dt,
                    "curvature_oracle_estimator": "analyticImplicitSurface",
                    "curvature_oracle_geometry": "levelSetField",
                    "geometry_oracle_estimator": "analyticImplicitSurface",
                    "geometry_oracle_geometry": "analyticImplicitSurface",
                    "candidate_estimator": "connectedFit+helmholtzPreserveModes",
                    "curvature_oracle_solver_status": metadata[curvature_oracle_id]["status"],
                    "geometry_oracle_solver_status": metadata[geometry_oracle_id]["status"],
                    "candidate_solver_status": metadata[candidate_id]["status"],
                    "curvature_oracle_topology_valid": int(
                        metadata[curvature_oracle_id]["topology"]["valid"]
                    ),
                    "geometry_oracle_topology_valid": int(
                        metadata[geometry_oracle_id]["topology"]["valid"]
                    ),
                    "candidate_topology_valid": int(
                        metadata[candidate_id]["topology"]["valid"]
                    ),
                    "curvature_oracle_interface_elements": metadata[
                        curvature_oracle_id
                    ]["topology"].get(
                        "elements", ""
                    ),
                    "geometry_oracle_interface_elements": metadata[
                        geometry_oracle_id
                    ]["topology"].get("elements", ""),
                    "curvature_oracle_max_cell_velocity_m_per_s": curvature_oracle_max,
                    "geometry_oracle_max_cell_velocity_m_per_s": geometry_oracle_max,
                    "candidate_max_cell_velocity_m_per_s": metadata[candidate_id]["max_u"],
                    "candidate_vs_curvature_oracle_max_difference_m_per_s": candidate_curvature_difference,
                    "candidate_vs_curvature_oracle_relative": (
                        candidate_curvature_difference/curvature_oracle_max
                    ),
                    "geometry_localization_effect_max_difference_m_per_s": geometry_effect,
                    "geometry_localization_effect_relative_to_curvature_oracle": (
                        geometry_effect/curvature_oracle_max
                    ),
                    "candidate_vs_geometry_oracle_max_difference_m_per_s": candidate_geometry_difference,
                    "candidate_vs_geometry_oracle_relative": (
                        candidate_geometry_difference/geometry_oracle_max
                    ),
                    "curvature_oracle_min_grad_psi_band": metadata[curvature_oracle_id]["metrics"][
                        "min_grad_psi"
                    ],
                    "curvature_oracle_max_grad_psi_band": metadata[curvature_oracle_id]["metrics"][
                        "max_grad_psi"
                    ],
                    "curvature_oracle_rho_lent_reset_l1": metadata[curvature_oracle_id]["metrics"][
                        "rho_lent_reset_l1"
                    ],
                    "candidate_curvature_difference_cell_index": candidate_curvature_cell,
                    "geometry_effect_cell_index": geometry_effect_cell,
                    "candidate_geometry_difference_cell_index": candidate_geometry_cell,
                }
            )
            print(
                f"[velocity-oracle] N={resolution} {mesh_variant} "
                f"candidate-curvature={candidate_curvature_difference:.10g}, "
                f"geometry-effect={geometry_effect:.10g}, "
                f"candidate-geometry={candidate_geometry_difference:.10g} m/s"
            )

    result = study_dir / "variable_curvature_velocity_oracle.csv"
    write_rows(result, rows)
    publication = (REPO / config["publication_table"]).resolve()
    write_rows(publication, rows)
    print(f"[velocity-oracle] result: {result}")
    print(f"[velocity-oracle] publication: {publication}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
