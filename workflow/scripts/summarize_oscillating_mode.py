#!/usr/bin/env python3
"""Summarize the N=64 mode-2 droplet trace, tolerating a partial crash row."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


TRACE_FIELDS = [
    "TIME",
    "maxMagU",
    "phaseVolumeRelError",
    "m2CosCoefficient",
    "m2SinCoefficient",
    "m2Amplitude",
    "m2Phase",
    "minGradPsiBand",
    "maxGradPsiBand",
    "pLaplace",
]


def complete_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return [
        row for row in rows
        if all(row.get(field) not in (None, "") for field in TRACE_FIELDS)
    ]


def zero_crossings(rows: list[dict[str, str]]) -> list[float]:
    crossings: list[float] = []
    for previous, current in zip(rows, rows[1:]):
        a = float(previous["m2CosCoefficient"])
        b = float(current["m2CosCoefficient"])
        if not ((a <= 0 < b) or (a >= 0 > b)):
            continue
        t0 = float(previous["TIME"])
        t1 = float(current["TIME"])
        crossings.append(t0 - a*(t1 - t0)/(b - a))
    return crossings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("trace", type=Path)
    args = parser.parse_args()

    rows = complete_rows(args.input)
    if not rows:
        raise SystemExit("no complete oscillating-droplet rows")

    crossings = zero_crossings(rows)
    periods = [crossings[i] - crossings[i - 2] for i in range(2, len(crossings))]
    measured_period = sum(periods)/len(periods) if periods else math.nan
    sigma = 0.07274
    rho_sum = 998.2 + 1.19
    radius = 1e-3
    exact_period = 2*math.pi/math.sqrt(6*sigma/(rho_sum*radius**3))

    initial = rows[0]
    final = rows[-1]
    max_u = max(rows, key=lambda row: float(row["maxMagU"]))
    max_volume = max(rows, key=lambda row: float(row["phaseVolumeRelError"]))
    max_amplitude = max(rows, key=lambda row: float(row["m2Amplitude"]))

    summary = {
        "candidate": "h4_helmholtz_l16_preserve_m2_m4",
        "N_CELLS": 64,
        "solver_status": "SIGFPE",
        "requested_end_time_s": 0.1,
        "last_complete_time_s": final["TIME"],
        "complete_samples": len(rows),
        "initial_m2_amplitude_m": initial["m2Amplitude"],
        "maximum_m2_amplitude_m": max_amplitude["m2Amplitude"],
        "maximum_m2_amplitude_time_s": max_amplitude["TIME"],
        "final_m2_amplitude_m": final["m2Amplitude"],
        "maximum_velocity_m_per_s": max_u["maxMagU"],
        "maximum_velocity_time_s": max_u["TIME"],
        "maximum_phase_volume_relative_error": max_volume[
            "phaseVolumeRelError"
        ],
        "measured_mode2_period_s": measured_period,
        "inviscid_mode2_period_s": exact_period,
        "mode2_period_relative_error": measured_period/exact_period - 1,
        "zero_crossings_used": len(crossings),
        "initial_min_grad_psi_band": initial["minGradPsiBand"],
        "initial_max_grad_psi_band": initial["maxGradPsiBand"],
        "final_min_grad_psi_band": final["minGradPsiBand"],
        "final_max_grad_psi_band": final["maxGradPsiBand"],
        "final_p_laplace_Pa": final["pLaplace"],
    }

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)

    args.trace.parent.mkdir(parents=True, exist_ok=True)
    with args.trace.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRACE_FIELDS)
        writer.writeheader()
        writer.writerows({field: row[field] for field in TRACE_FIELDS} for row in rows)

    print(f"[oscillation] summary: {args.summary}")
    print(f"[oscillation] trace: {args.trace}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
