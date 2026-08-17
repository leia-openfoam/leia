#!/usr/bin/env python3
"""r_m(drift): capillary-mode growth rate as a function of the psi-profile
drift, on the shipped solver.

THE HYPOTHESIS UNDER TEST. On the CLEAN base state every resolvable interface
mode is damped (measured at N=64: r_2 = -11 1/s, r_4 = -131 1/s, omega within
4-8% of the capillary dispersion relation), while the coupled run nevertheless
blows up -- and the delivered-force roughness is slaved to the profile drift
(driverAcross ~ gradPsiErrorBand^1.00, r^2 = 0.98). The missing link is whether
the drift RAISES THE LOOP GAIN PAST NEUTRAL: is the instability a mode not of
the clean state but of the DRIFTED state? This utility measures r_m on base
states taken from along a coupled run's own drift trajectory. A zero crossing
r_m(drift*) = 0 demonstrates the source and turns "SDF maintenance" into a
BUDGET: a maintenance operation pays off only if the perturbation it injects is
small against the margin to drift*.

DESIGN.
  * Base states: snapshots of a coupled production run, selected by their own
    minGradPsiBand (the run's CSV provides the mapping time -> drift).
  * U and p_rgh are ZEROED in the restart states. This isolates the
    psi-profile's effect on the loop gain from the accumulated velocity field;
    the solver re-establishes the Young-Laplace balance at t=0 and the mode
    evolves on the drifted profile alone. (A full-state arm -- snapshot
    velocity kept -- can be added later; it answers a different question.)
  * phi and old-time fields are removed from the restart so the solver
    regenerates them consistent with U = 0.
  * The mode measurement itself is interface_mode_trajectory.py, unchanged:
    perturb by eps*h*cos(m theta), fit a_m(t) = A0 exp(r t) cos(w t + phi).
    VALIDITY per (state, m): fit residual, and omega against the analytic
    capillary dispersion -- modes are only quoted where both hold. At N=128,
    m <= 6 keeps >= 13 cells per wavelength.

Usage:
  mode_rate_vs_drift.py <generation_case> [--targets 0.985,0.97,0.93,0.88,0.80]
      [--modes 2,4,6] [--steps 700] [--sample 10] [--eps 1e-2] [--out csv]
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

import numpy as np
from foamlib import FoamCase

HERE = os.path.dirname(os.path.abspath(__file__))
TRAJ = os.path.join(HERE, "interface_mode_trajectory.py")
SOLVER_CSV = "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv"

# fields the restart keeps; everything else in the snapshot is regenerated
KEEP = {"psi", "U", "p_rgh", "alpha.water", "NarrowBand"}


def _time_dirs(d):
    out = []
    for x in os.listdir(d):
        if not os.path.isdir(os.path.join(d, x)):
            continue
        try:
            out.append((float(x), x))
        except ValueError:
            pass
    return sorted(out)


def _read_run_csv(gen):
    rows = []
    with open(os.path.join(gen, SOLVER_CSV)) as fh:
        for r in csv.DictReader(fh):
            try:
                rows.append({k: float(v) for k, v in r.items()})
            except (TypeError, ValueError):
                continue
    return rows


def _zero_field(case_dir, name, vector):
    f = FoamCase(case_dir)[0][name]
    a = np.asarray(f.internal_field, dtype=float)
    f.internal_field = (np.zeros_like(a)).tolist() if not vector else \
        [[0.0, 0.0, 0.0]] * len(a)


def _make_restart(gen, tdir, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    for d in ("constant", "system"):
        shutil.copytree(os.path.join(gen, d), os.path.join(dst, d))
    for f in ("case_params.json",):
        p = os.path.join(gen, f)
        if os.path.exists(p):
            shutil.copy(p, dst)
    os.makedirs(os.path.join(dst, "0"))
    src = os.path.join(gen, tdir)
    for f in os.listdir(src):
        if f in KEEP and os.path.isfile(os.path.join(src, f)):
            shutil.copy(os.path.join(src, f), os.path.join(dst, "0", f))
    _zero_field(dst, "U", vector=True)
    _zero_field(dst, "p_rgh", vector=False)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("generation_case")
    ap.add_argument("--targets", default="0.985,0.97,0.93,0.88,0.80",
                    help="minGradPsiBand values to select base states at")
    ap.add_argument("--modes", default="2,4,6")
    ap.add_argument("--steps", type=int, default=700)
    ap.add_argument("--sample", type=int, default=10)
    ap.add_argument("--eps", type=float, default=1e-2)
    ap.add_argument("--out", default=None)
    ap.add_argument("--work", default=None)
    a = ap.parse_args(argv)

    gen = os.path.abspath(a.generation_case)
    work = os.path.abspath(a.work or (gen + "_rates"))
    os.makedirs(work, exist_ok=True)

    rows = _read_run_csv(gen)
    tds = [(t, d) for t, d in _time_dirs(gen) if t > 0]
    if not rows or not tds:
        print("generation run has no CSV rows / no snapshots", file=sys.stderr)
        return 1

    # drift level of each SNAPSHOT: the CSV row nearest its write time
    def drift_at(t):
        r = min(rows, key=lambda r: abs(r["TIME"] - t))
        return r["minGradPsiBand"], r["gradPsiL2ErrorBand"], r["maxMagU"]

    picks = []
    for tgt in [float(x) for x in a.targets.split(",")]:
        t, d = min(tds, key=lambda td: abs(drift_at(td[0])[0] - tgt))
        mg, ge, mu = drift_at(t)
        if any(p[1] == d for p in picks):
            continue
        picks.append((t, d, mg, ge, mu))
    print("[drift] selected base states:")
    for t, d, mg, ge, mu in picks:
        print(f"    t={t:<10.5g} minGradPsiBand={mg:.4f} "
              f"gradPsiL2ErrBand={ge:.4g} maxMagU(base)={mu:.3g}")

    allrows = []
    for t, d, mg, ge, mu in picks:
        tag = f"s{mg:.3f}".replace(".", "p")
        case = os.path.join(work, f"restart_{tag}")
        print(f"\n[drift] === base state t={t:.5g} (minGradPsi={mg:.4f}) ===")
        _make_restart(gen, d, case)
        out_csv = os.path.join(work, f"rates_{tag}.csv")
        rc = subprocess.call(
            [sys.executable, TRAJ, case, "--modes", a.modes,
             "--steps", str(a.steps), "--sample", str(a.sample),
             "--eps", str(a.eps), "--out", out_csv,
             "--work", case + "_traj"])
        if rc != 0:
            print(f"[drift] harness failed on {tag} -- recorded, continuing")
            continue
        with open(out_csv) as fh:
            for r in csv.DictReader(fh):
                r["snapshot_t"] = t
                r["minGradPsiBand"] = mg
                r["gradPsiL2ErrorBand"] = ge
                allrows.append(r)

    if a.out and allrows:
        cols = list(allrows[0].keys())
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(allrows)
        print(f"\n[drift] wrote {a.out} ({len(allrows)} rows)")

    # summary table
    print(f"\n{'minGradPsi':>11} {'t_snap':>9} " +
          " ".join(f"{'r_m'+m:>10}" for m in a.modes.split(",")))
    for t, d, mg, ge, mu in picks:
        line = f"{mg:>11.4f} {t:>9.5g} "
        for m in a.modes.split(","):
            r = [x for x in allrows
                 if x["snapshot_t"] == t and x["mode"] == m]
            line += f"{float(r[0]['r_m']):>10.1f} " if r else f"{'--':>10} "
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
