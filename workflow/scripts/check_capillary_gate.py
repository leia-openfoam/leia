#!/usr/bin/env python3
"""Check one capillary-droplet metrics CSV against explicit research gates.

The solver writes this CSV at every time step.  This checker deliberately uses
only the standard library so it can run on a login node before a study is
accepted or its plots are generated.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


REQUIRED = {
    "TIME",
    "maxMagU",
    "maxMagUPrime",
    "phaseVolumeRelError",
    "centroidError",
    "pLaplace",
    "minGradPsiBand",
    "maxGradPsiBand",
    "rhoMassResidualRelL1",
}


def finite_rows(path: Path) -> tuple[list[dict[str, float]], list[str]]:
    errors: list[str] = []
    rows: list[dict[str, float]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED.difference(reader.fieldnames or [])
        if missing:
            return [], ["missing columns: " + ", ".join(sorted(missing))]
        for line_no, raw in enumerate(reader, start=2):
            try:
                row = {key: float(value) for key, value in raw.items() if key}
            except (TypeError, ValueError):
                errors.append(f"line {line_no}: incomplete or non-numeric row")
                continue
            bad = [key for key, value in row.items() if not math.isfinite(value)]
            if bad:
                errors.append(
                    f"line {line_no}: non-finite values in {', '.join(sorted(bad))}"
                )
                continue
            rows.append(row)
    if not rows:
        errors.append("no complete finite data rows")
    return rows, errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--target-time", type=float)
    parser.add_argument("--max-velocity", type=float)
    parser.add_argument(
        "--max-disturbance",
        type=float,
        help="maximum allowed max|U-U_ref| over all rows",
    )
    parser.add_argument("--max-phase-volume-relative-error", type=float)
    parser.add_argument("--max-centroid-error", type=float)
    parser.add_argument("--min-grad", type=float)
    parser.add_argument("--max-grad", type=float)
    parser.add_argument("--max-relative-mass-residual", type=float, default=1e-10)
    parser.add_argument("--expected-pressure-jump", type=float)
    parser.add_argument("--pressure-jump-tolerance", type=float)
    args = parser.parse_args()

    rows, errors = finite_rows(args.csv)
    if not rows:
        for error in errors:
            print(f"FAIL: {error}")
        return 1

    last = rows[-1]
    maxima = {key: max(row[key] for row in rows) for key in REQUIRED}
    minima = {key: min(row[key] for row in rows) for key in REQUIRED}

    if args.target_time is not None and last["TIME"] + 1e-12 < args.target_time:
        errors.append(
            f"ended at t={last['TIME']:.9g}, before target {args.target_time:.9g}"
        )
    if args.max_velocity is not None and maxima["maxMagU"] > args.max_velocity:
        errors.append(
            f"maxMagU={maxima['maxMagU']:.9g} exceeds {args.max_velocity:.9g}"
        )
    if (
        args.max_disturbance is not None
        and maxima["maxMagUPrime"] > args.max_disturbance
    ):
        errors.append(
            "maxMagUPrime="
            f"{maxima['maxMagUPrime']:.9g} exceeds {args.max_disturbance:.9g}"
        )
    if (
        args.max_phase_volume_relative_error is not None
        and maxima["phaseVolumeRelError"]
        > args.max_phase_volume_relative_error
    ):
        errors.append(
            "phaseVolumeRelError="
            f"{maxima['phaseVolumeRelError']:.9g} exceeds "
            f"{args.max_phase_volume_relative_error:.9g}"
        )
    if (
        args.max_centroid_error is not None
        and maxima["centroidError"] > args.max_centroid_error
    ):
        errors.append(
            f"centroidError={maxima['centroidError']:.9g} exceeds "
            f"{args.max_centroid_error:.9g}"
        )
    if args.min_grad is not None and minima["minGradPsiBand"] < args.min_grad:
        errors.append(
            "minGradPsiBand="
            f"{minima['minGradPsiBand']:.9g} is below {args.min_grad:.9g}"
        )
    if args.max_grad is not None and maxima["maxGradPsiBand"] > args.max_grad:
        errors.append(
            "maxGradPsiBand="
            f"{maxima['maxGradPsiBand']:.9g} exceeds {args.max_grad:.9g}"
        )
    if maxima["rhoMassResidualRelL1"] > args.max_relative_mass_residual:
        errors.append(
            "rhoMassResidualRelL1="
            f"{maxima['rhoMassResidualRelL1']:.9g} exceeds "
            f"{args.max_relative_mass_residual:.9g}"
        )
    if args.expected_pressure_jump is not None:
        if args.pressure_jump_tolerance is None:
            parser.error("--expected-pressure-jump requires --pressure-jump-tolerance")
        pressure_error = abs(last["pLaplace"] - args.expected_pressure_jump)
        if pressure_error > args.pressure_jump_tolerance:
            errors.append(
                f"final pressure-jump error={pressure_error:.9g} exceeds "
                f"{args.pressure_jump_tolerance:.9g}"
            )

    print(
        f"{args.csv}: rows={len(rows)}, finalTime={last['TIME']:.9g}, "
        f"maxMagU={maxima['maxMagU']:.9g}, "
        f"maxMagUPrime={maxima['maxMagUPrime']:.9g}, "
        f"maxPhaseVolumeRelError={maxima['phaseVolumeRelError']:.9g}, "
        f"maxCentroidError={maxima['centroidError']:.9g}, "
        f"gradPsi=[{minima['minGradPsiBand']:.9g},"
        f"{maxima['maxGradPsiBand']:.9g}], "
        f"maxRelMassResidual={maxima['rhoMassResidualRelL1']:.9g}"
    )
    for error in errors:
        print(f"FAIL: {error}")
    if errors:
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
