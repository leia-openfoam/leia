#!/usr/bin/env python3
"""Run the N=64 curvature-delivery snapshot replay gate.

The source trajectory is generated once with exact constant curvature and the
full uncached-quadratic SL/rhoLENT solver.  Selected transported-interface
snapshots are then copied into fresh cases, U is reset by retaining the case's
uniform Utrans initial field, SL advection is frozen, and the complete solver is
advanced exactly one timestep.  The output table contains only the coupled
cell-centred max|U-Utrans| decision metric (plus provenance/status columns).
"""

from __future__ import annotations

import argparse
import csv
import json
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


def checked_study_dir(name: str) -> Path:
    studies = (REPO / "studies").resolve()
    target = (studies / name).resolve()
    if target.parent != studies or not name.startswith("transISTN64"):
        raise SystemExit(f"refusing unsafe replay study path: {target}")
    return target


def run_command(command: list[str], cwd: Path, log_name: str,
                env: dict[str, str] | None = None, timeout: int = 360) -> int:
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
            log.write(f"\nreplay runner timeout after {timeout} s\n")
            return 124


def require_commands() -> None:
    missing = [name for name in ("blockMesh", "leiaSetFields", SOLVER)
               if shutil.which(name) is None]
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
        print(f"[replay] ignored parameters: {ignored}")
    if len(variations) != 1:
        raise SystemExit(
            f"curvature replay base must materialize one case, got {len(variations)}"
        )
    base = dict(variations[0])
    override = os.environ.get("LEIA_MOMENTUM_PREDICTOR")
    if override:
        if override not in ("yes", "no"):
            raise SystemExit(
                "LEIA_MOMENTUM_PREDICTOR must be 'yes' or 'no', got "
                f"'{override}'"
            )
        base["MOMENTUM_PREDICTOR"] = override
        print(f"[replay] MOMENTUM_PREDICTOR override -> {override}")
    return base



def make_case(case_dir: Path, tokens: dict[str, str], index: str) -> dict:
    return materialize.materialize(
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


def numeric_time_dirs(case_dir: Path) -> list[tuple[float, Path]]:
    found = []
    for path in case_dir.iterdir():
        if not path.is_dir():
            continue
        try:
            found.append((float(path.name), path))
        except ValueError:
            pass
    return sorted(found)


def nearest_snapshot(case_dir: Path, target: float) -> tuple[float, Path]:
    candidates = numeric_time_dirs(case_dir)
    if not candidates:
        raise SystemExit(f"no numeric snapshot directories in {case_dir}")
    return min(candidates, key=lambda item: abs(item[0] - target))


def replace_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def prepare_replay_fields(replay_case: Path, snapshot: Path) -> None:
    replace_tree(replay_case / "0.org", replay_case / "0")
    for field in ("psi", "alpha.water"):
        src = snapshot / field
        if not src.exists():
            raise SystemExit(f"snapshot is missing {field}: {snapshot}")
        shutil.copy2(src, replay_case / "0" / field)


def read_metric(case_dir: Path) -> tuple[int, float | None]:
    path = case_dir / f"{SOLVER}.csv"
    if not path.exists():
        return 0, None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return 0, None
    return len(rows), max(float(row["maxMagUPrime"]) for row in rows)


def read_connected_topology(case_dir: Path) -> dict[str, int | bool | str]:
    """Audit that a quiet replay was produced by a real closed reconstruction."""
    log_path = case_dir / "log.solve"
    if not log_path.exists():
        return {"topology_valid": False, "topology_reason": "missing log"}
    pattern = re.compile(
        r"connectedInterface curvature: (\d+) elements, (\d+) component\(s\), "
        r"(\d+) open endpoint\(s\), (\d+) fits, (\d+) fallback"
    )
    matches = pattern.findall(log_path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        return {"topology_valid": False, "topology_reason": "missing audit line"}
    elements, components, open_ends, fits, fallback = map(int, matches[-1])
    valid = (
        elements > 0
        and components == 1
        and open_ends == 0
        and fits == elements
        and fallback == 0
    )
    return {
        "topology_valid": valid,
        "topology_reason": "ok" if valid else "non-closed or incomplete chain",
        "interface_elements": elements,
        "interface_components": components,
        "open_endpoints": open_ends,
        "curvature_fits": fits,
        "curvature_fallbacks": fallback,
    }


def write_results(path: Path, rows: list[dict]) -> None:
    fields = [
        "target_snapshot_time_s",
        "actual_snapshot_time_s",
        "variant",
        "surface_tension_force",
        "curvature_extension",
        "face_interpolation",
        "face_curvature_source",
        "fit_half_width",
        "regularization_operator",
        "tangential_regularization",
        "regularization_iterations",
        "preserve_modes_through",
        "topology_valid",
        "topology_reason",
        "interface_elements",
        "interface_components",
        "open_endpoints",
        "curvature_fits",
        "curvature_fallbacks",
        "solver_status",
        "metric_samples",
        "max_U_minus_Utrans_m_per_s",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        nargs="?",
        default=str(REPO / "config" / "transISTN64CurvatureDeliveryReplay.yaml"),
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="remove and regenerate only this validated replay study directory",
    )
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
    baseline = study_dir / "baseline_exact_constant"
    if not (baseline / f"{SOLVER}.csv").exists():
        print("[replay] generating quiet transported-interface baseline")
        make_case(baseline, tokens, "baseline")
        replace_tree(baseline / "0.org", baseline / "0")
        status = run_command(["blockMesh"], baseline, "log.blockMesh")
        if status:
            raise SystemExit(f"baseline blockMesh failed with status {status}")
        status = run_command(
            ["leiaSetFields", "-alphaName", "alpha.water"],
            baseline,
            "log.setFields",
        )
        if status:
            raise SystemExit(f"baseline leiaSetFields failed with status {status}")
        status = run_command([SOLVER], baseline, "log.solve", timeout=600)
        (baseline / "run.status").write_text(f"{status}\n", encoding="utf-8")
        if status:
            raise SystemExit(f"baseline solver failed with status {status}")
    else:
        print("[replay] reusing existing quiet baseline")

    with (baseline / "case_params.json").open(encoding="utf-8") as handle:
        baseline_meta = json.load(handle)
    dt = float(baseline_meta["tokens"]["MAX_DELTA_T"])

    replay_cfg = config["replay"]
    rows: list[dict] = []
    case_index = 0
    for target in replay_cfg["snapshots"]:
        actual, snapshot = nearest_snapshot(baseline, float(target))
        print(f"[replay] target t={target:g}: source snapshot {snapshot.name}")
        for variant in replay_cfg["variants"]:
            variant_id = variant["id"]
            replay_case = (
                study_dir
                / f"snapshot_{float(target):.5f}"
                / variant_id
            )
            variant_tokens = dict(tokens)
            variant_tokens.update(
                {
                    "END_TIME": f"{dt:.10g}",
                    "WRITE_INTERVAL": f"{dt:.10g}",
                    "SURFACE_TENSION_FORCE": variant["surface_tension_force"],
                    "CURVATURE_EXTENSION": variant["curvature_extension"],
                    "CURVATURE_ESTIMATOR": variant.get(
                        "estimator", "connectedFit"
                    ),
                    "CURVATURE_FACE_INTERPOLATION": variant["face_interpolation"],
                    "CURVATURE_FACE_SOURCE": variant.get
                    (
                        "face_curvature_source", "model"
                    ),
                    "CURVATURE_FIT_HALF_WIDTH": str(
                        variant.get("fit_half_width", 3)
                    ),
                    "CURVATURE_REGULARIZATION_OPERATOR": variant.get(
                        "regularization_operator", "helmholtz"
                    ),
                    "CURVATURE_TANGENTIAL_REGULARIZATION": str(
                        variant.get("tangential_regularization", 1)
                    ),
                    "CURVATURE_REGULARIZATION_ITERATIONS": str(
                        variant.get("regularization_iterations", 30)
                    ),
                    "CURVATURE_PRESERVE_MODES_THROUGH": str(
                        variant.get("preserve_modes_through", 4)
                    ),
                }
            )
            make_case(replay_case, variant_tokens, f"{case_index:05d}")
            case_index += 1
            replace_tree(baseline / "constant" / "polyMesh",
                         replay_case / "constant" / "polyMesh")
            prepare_replay_fields(replay_case, snapshot)

            env = dict(os.environ)
            env["SL_FREEZE_INTERFACE"] = "1"
            status = run_command([SOLVER], replay_case, "log.solve", env=env)
            (replay_case / "run.status").write_text(
                f"{status}\n", encoding="utf-8"
            )
            samples, metric = read_metric(replay_case)
            topology = (
                read_connected_topology(replay_case)
                if variant["curvature_extension"] == "connectedInterface"
                else {
                    "topology_valid": "",
                    "topology_reason": "not applicable",
                    "interface_elements": "",
                    "interface_components": "",
                    "open_endpoints": "",
                    "curvature_fits": "",
                    "curvature_fallbacks": "",
                }
            )
            print(
                f"[replay]   {variant_id}: status={status}, "
                f"samples={samples}, max|U-Utrans|={metric}"
            )
            rows.append(
                {
                    "target_snapshot_time_s": f"{float(target):.8g}",
                    "actual_snapshot_time_s": f"{actual:.10g}",
                    "variant": variant_id,
                    "surface_tension_force": variant["surface_tension_force"],
                    "curvature_extension": variant["curvature_extension"],
                    "face_interpolation": variant["face_interpolation"],
                    "face_curvature_source": variant.get
                    (
                        "face_curvature_source", "model"
                    ),
                    "fit_half_width": variant.get("fit_half_width", 3),
                    "regularization_operator": variant.get(
                        "regularization_operator", "helmholtz"
                    ),
                    "tangential_regularization": variant.get(
                        "tangential_regularization", 1
                    ),
                    "regularization_iterations": variant.get(
                        "regularization_iterations", 30
                    ),
                    "preserve_modes_through": variant.get(
                        "preserve_modes_through", ""
                    ),
                    **topology,
                    "solver_status": status,
                    "metric_samples": samples,
                    "max_U_minus_Utrans_m_per_s": (
                        "" if metric is None else f"{metric:.10g}"
                    ),
                }
            )

    results = study_dir / "curvature_delivery_replay.csv"
    write_results(results, rows)
    print(f"[replay] results: {results}")
    if config.get("publication_table"):
        publication = (REPO / config["publication_table"]).resolve()
        publication.parent.mkdir(parents=True, exist_ok=True)
        write_results(publication, rows)
        print(f"[replay] publication table: {publication}")
    valid = all(
        row["solver_status"] == 0 and row["topology_valid"] is not False
        for row in rows
    )
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
