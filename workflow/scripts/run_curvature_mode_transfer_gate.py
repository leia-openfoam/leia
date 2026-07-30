#!/usr/bin/env python3
"""Run the manufactured connected-interface curvature Fourier-mode gate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import foam_param  # noqa: E402
import materialize  # noqa: E402


UTILITY = "leiaTestConnectedCurvatureModes"


def checked_study_dir(name: str) -> Path:
    studies = (REPO / "studies").resolve()
    target = (studies / name).resolve()
    if target.parent != studies or not name.startswith(
        "curvatureModeTransferGate"
    ):
        raise SystemExit(f"refusing unsafe mode-gate study path: {target}")
    return target


def run_command(
    command: list[str], cwd: Path, log_name: str, timeout: int = 180
) -> int:
    with (cwd / log_name).open("w", encoding="utf-8") as log:
        try:
            result = subprocess.run(
                command,
                cwd=cwd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            return result.returncode
        except subprocess.TimeoutExpired:
            log.write(f"\nmode-gate timeout after {timeout} s\n")
            return 124


def require_commands(mesh_variants: list[str]) -> None:
    required = ["blockMesh", UTILITY]
    if "perturbed" in mesh_variants:
        required.append("leiaPerturbMesh")
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        raise SystemExit(
            "OpenFOAM environment/build is not active; missing: "
            + ", ".join(missing)
        )


def base_tokens(config: dict) -> dict[str, str]:
    case_name = config["case"]
    base_case = REPO / "cases" / case_name
    axes, constants, ignored = foam_param.build_token_grid(
        str(REPO / "cases" / f"{case_name}.parameter"),
        str(REPO / "cases" / "default.parameter"),
        str(base_case),
        axes_override={"N_CELLS": [config["resolutions"][0]]},
        collapse_other_axes=True,
    )
    variations = foam_param.expand(axes, constants)
    if ignored:
        print(f"[mode-gate] ignored parameters: {ignored}")
    if len(variations) != 1:
        raise SystemExit(
            f"mode gate base must materialize one case, got {len(variations)}"
        )
    return dict(variations[0])


def make_case(case_dir: Path, tokens: dict[str, str], index: str) -> None:
    materialize.materialize(
        str(REPO / "cases" / "transISTDroplet2D"),
        tokens,
        str(case_dir),
        1,
        "hex",
        "serial",
        2,
        "transISTDroplet2D",
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
    # leiaPerturbMesh writes only the moved points into 0/polyMesh. Overlay
    # those files on the complete blockMesh topology in constant/polyMesh;
    # replacing the directory would discard faces/owner/neighbour/boundary.
    constant_mesh.mkdir(parents=True, exist_ok=True)
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


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise SystemExit(f"expected one data row in {path}, found {len(rows)}")
    return rows[0]


def quantitative_pass(row: dict[str, str], gate: dict) -> tuple[bool, str]:
    if row["topology_valid"] not in ("1", "true", "True"):
        return False, "topology"
    if int(row["N_CELLS"]) < int(gate["minimum_resolution"]):
        return True, "topology-only coarse point"
    transfer = float(row["regularized_amplitude_transfer"])
    phase = float(row["regularized_phase_error_deg"])
    relative_l2 = float(row["regularized_relative_L2"])
    failures = []
    if not (
        float(gate["amplitude_transfer_min"])
        <= transfer
        <= float(gate["amplitude_transfer_max"])
    ):
        failures.append("amplitude")
    if phase > float(gate["maximum_phase_error_deg"]):
        failures.append("phase")
    if relative_l2 > float(gate["maximum_relative_L2"]):
        failures.append("L2")
    return not failures, "ok" if not failures else "+".join(failures)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def candidate_summary(rows: list[dict], config: dict) -> list[dict]:
    gate = config["quantitative_gate"]
    summaries = []
    for candidate in config["candidates"]:
        cid = candidate["id"]
        subset = [row for row in rows if row["candidate"] == cid]
        row_pass = all(int(row["row_pass"]) == 1 for row in subset)
        refinement_pass = True
        refinement_reason = "ok"
        if gate.get("require_refinement_nonincrease", False):
            indexed = {
                (row["mesh_variant"], row["mode"], int(row["N_CELLS"])): row
                for row in subset
            }
            for mesh_variant in config["mesh_variants"]:
                for mode in config["modes"]:
                    row64 = indexed.get((mesh_variant, str(mode), 64))
                    row128 = indexed.get((mesh_variant, str(mode), 128))
                    if row64 and row128 and (
                        float(row128["regularized_relative_L2"])
                        > float(row64["regularized_relative_L2"])
                    ):
                        refinement_pass = False
                        refinement_reason = "relative L2 grows on refinement"
        passed = row_pass and refinement_pass
        quantitative = [
            row for row in subset
            if int(row["N_CELLS"]) >= int(gate["minimum_resolution"])
        ]
        failed_reasons = sorted(
            {
                row["row_reason"]
                for row in subset
                if int(row["row_pass"]) != 1
            }
        )
        if not refinement_pass:
            failed_reasons.append(refinement_reason)
        summaries.append(
            {
                "candidate": cid,
                "fit_half_width": candidate["fit_half_width"],
                "regularization_operator": candidate["regularization_operator"],
                "tangential_regularization": candidate[
                    "tangential_regularization"
                ],
                "regularization_iterations": candidate[
                    "regularization_iterations"
                ],
                "preserve_modes_through": candidate.get(
                    "preserve_modes_through", ""
                ),
                "gate_pass": int(passed),
                "row_limits_pass": int(row_pass),
                "refinement_pass": int(refinement_pass),
                "reason": "ok" if passed else "; ".join(failed_reasons),
                "minimum_amplitude_transfer_Nge64": min(
                    float(row["regularized_amplitude_transfer"])
                    for row in quantitative
                ),
                "maximum_amplitude_transfer_Nge64": max(
                    float(row["regularized_amplitude_transfer"])
                    for row in quantitative
                ),
                "maximum_phase_error_deg_Nge64": max(
                    float(row["regularized_phase_error_deg"])
                    for row in quantitative
                ),
                "maximum_relative_L2_Nge64": max(
                    float(row["regularized_relative_L2"])
                    for row in quantitative
                ),
            }
        )
    return summaries


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        nargs="?",
        default=str(REPO / "config" / "curvatureModeTransferGate.yaml"),
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="remove and regenerate only studies/curvatureModeTransferGate",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    require_commands(config["mesh_variants"])

    study_dir = checked_study_dir(config["study_name"])
    if args.fresh and study_dir.exists():
        shutil.rmtree(study_dir)
    study_dir.mkdir(parents=True, exist_ok=True)

    tokens0 = base_tokens(config)
    rows: list[dict] = []
    index = 0
    gate = config["quantitative_gate"]
    centre = "(" + " ".join(str(v) for v in config["center"]) + ")"

    for candidate in config["candidates"]:
        for resolution in config["resolutions"]:
            for mesh_variant in config["mesh_variants"]:
                case_dir = (
                    study_dir
                    / candidate["id"]
                    / f"N{resolution}_{mesh_variant}"
                )
                tokens = dict(tokens0)
                tokens.update(
                    {
                        "N_CELLS": str(resolution),
                        "CURVATURE_EXTENSION": "connectedInterface",
                        "CURVATURE_FACE_INTERPOLATION": "connectedInterface",
                        "CURVATURE_FACE_SOURCE": "registered",
                        "CURVATURE_FIT_HALF_WIDTH": str(
                            candidate["fit_half_width"]
                        ),
                        "CURVATURE_REGULARIZATION_OPERATOR": candidate[
                            "regularization_operator"
                        ],
                        "CURVATURE_TANGENTIAL_REGULARIZATION": str(
                            candidate["tangential_regularization"]
                        ),
                        "CURVATURE_REGULARIZATION_ITERATIONS": str(
                            candidate["regularization_iterations"]
                        ),
                        "CURVATURE_PRESERVE_MODES_THROUGH": str(
                            candidate.get("preserve_modes_through", 4)
                        ),
                        "CURVATURE_ESTIMATOR": candidate.get(
                            "estimator", "connectedFit"
                        ),
                    }
                )
                make_case(case_dir, tokens, f"{index:05d}")
                index += 1
                mesh_status = build_mesh(case_dir, mesh_variant, config)
                if mesh_status:
                    raise SystemExit(
                        f"mesh failed ({mesh_status}): {case_dir}"
                    )

                for mode in config["modes"]:
                    output = case_dir / f"mode_m{mode}.csv"
                    status = run_command(
                        [
                            UTILITY,
                            "-mode",
                            str(mode),
                            "-epsilon",
                            str(config["epsilon"]),
                            "-radius",
                            str(config["radius"]),
                            "-center",
                            centre,
                            "-meshVariant",
                            mesh_variant,
                            "-output",
                            output.name,
                        ],
                        case_dir,
                        f"log.mode_m{mode}",
                    )
                    if status not in (0, 2) or not output.exists():
                        raise SystemExit(
                            f"mode utility failed ({status}): {case_dir}, m={mode}"
                        )
                    element_file = case_dir / "curvatureModeElements.csv"
                    if element_file.exists():
                        element_file.replace(case_dir / f"elements_m{mode}.csv")
                    row = read_single_row(output)
                    passed, reason = quantitative_pass(row, gate)
                    row = {
                        "candidate": candidate["id"],
                        **row,
                        "utility_status": status,
                        "row_pass": int(passed),
                        "row_reason": reason,
                    }
                    rows.append(row)
                    print(
                        "[mode-gate] "
                        f"{candidate['id']} N={resolution} {mesh_variant} m={mode}: "
                        f"H={float(row['regularized_amplitude_transfer']):.4f}, "
                        f"phase={float(row['regularized_phase_error_deg']):.3f} deg, "
                        f"relL2={float(row['regularized_relative_L2']):.4f}, "
                        f"pass={passed} ({reason})"
                    )

    detailed_fields = list(rows[0].keys())
    results = study_dir / "curvature_mode_transfer_gate.csv"
    write_csv(results, rows, detailed_fields)

    summaries = candidate_summary(rows, config)
    summary_fields = list(summaries[0].keys())
    summary_path = study_dir / "curvature_mode_transfer_gate_summary.csv"
    write_csv(summary_path, summaries, summary_fields)

    publication = (REPO / config["publication_table"]).resolve()
    publication_summary = (REPO / config["publication_summary"]).resolve()
    write_csv(publication, rows, detailed_fields)
    write_csv(publication_summary, summaries, summary_fields)
    (study_dir / "gate_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )

    print(f"[mode-gate] detailed results: {results}")
    print(f"[mode-gate] candidate summary: {summary_path}")
    for row in summaries:
        print(
            f"[mode-gate] candidate {row['candidate']}: "
            f"{'PASS' if row['gate_pass'] else 'FAIL'}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
