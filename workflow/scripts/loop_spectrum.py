#!/usr/bin/env python3
"""Dominant eigenvalue of the quasi-static psi-loop, by power iteration on the
SHIPPED solver step.

WHY THIS EXISTS. Every lever in the curvature-stabilization campaign has been
decided by running the coupled droplet to blow-up and comparing t_blow or a
fitted growth rate. That ranks amplifiers; it does not say what the instability
IS. Two facts the programme measured make the distinction decisive (plan
sec. 16.1): the growth rate splits as r = r0(N) + c(N)*dt, r0 does NOT vanish
and RISES toward fine grids, and the c*dt mechanism is ABSENT at N=64 and
present at N=128.

A purely explicit-coupling instability must vanish as dt -> 0, because the
per-step amplification tends to 1. A growth rate that survives dt -> 0 says the
SEMI-DISCRETE system already has a positive eigenvalue -- an ill-posedness of
the spatial operator, not a coupling defect -- and then no time scheme, no
coupling depth and no dt policy can ever remove it. That is consistent with
every exoneration in sec. 16.2 (backward/CN, adaptivity, GAMG, psiOuterCorrectors,
six Picard passes, the departure-point kernel: all flat). This utility measures
the eigenvalue directly instead of inferring it from a blow-up.

THE OPERATOR MEASURED. Perturb ONLY psi, restart the solver, take one step,
and read the resulting psi perturbation back:

    M(v) = [ Step(psi0 + eps*v) - Step(psi0) ]_psi / eps

with every other field reset to its base value before each application. This is
deliberately the QUASI-STATIC closure -- the loop sec. 16.3 names as "the one
mechanism left standing": psi -> fit -> kappa_f -> quasi-static velocity
response -> displacement -> psi. It is NOT the full state-space spectrum (which
would carry U and p_rgh perturbations too); it is the spectrum of exactly the
mechanism the programme claims is operating, which is what makes it a test of
that claim rather than another ranking.

The solver is driven by RESTART rather than reimplemented, so the operator
measured is the shipped pipeline including every guard, fallback and seam rule.
Cost is one solver startup per application; at N=64 serial a 40-iteration run is
a couple of minutes.

VALIDITY GATE. M is linear only to O(eps). The utility sweeps eps and reports
the eigenvalue at each: if lambda moves with eps the map is not linear there and
the number is meaningless. Report the plateau, not a single eps.

Usage:
  loop_spectrum.py <prepared_case_dir> [--iters 40] [--eps 1e-6,1e-5,1e-4]
                   [--band 3.0] [--steps 1] [--dt-divisor 1] [--out out.csv]

The case must already be preprocessed (blockMesh + setFields done, 0/psi is the
initialised field).
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
from foamlib import FoamCase

SOLVER = "leiaSemiLagrangianLevelSetTwoPhaseFoam"
BASHRC = os.path.expanduser("~/OpenFOAM/OpenFOAM-v2512/etc/bashrc")


def _set_entries(path, entries):
    """Rewrite whole-line dictionary entries in an OpenFOAM file."""
    with open(path) as fh:
        lines = fh.readlines()
    out = []
    for ln in lines:
        s = ln.strip()
        hit = None
        for key in entries:
            if s.startswith(key) and (len(s) == len(key) or s[len(key)] in " \t;"):
                hit = key
                break
        out.append(f"{hit:<16}{entries[hit]};\n" if hit else ln)
    with open(path, "w") as fh:
        fh.writelines(out)


def _run_step(case_dir, env):
    r = subprocess.run(
        ["bash", "-lc", f"source {BASHRC} >/dev/null 2>&1 && cd {case_dir} && {SOLVER}"],
        capture_output=True, text=True, env=env,
    )
    return r.returncode, r.stdout + r.stderr


def _time_dirs(case_dir):
    out = []
    for d in os.listdir(case_dir):
        p = os.path.join(case_dir, d)
        if not os.path.isdir(p):
            continue
        try:
            t = float(d)
        except ValueError:
            continue
        out.append((t, d))
    return sorted(out)


def _read_psi(case_dir, tdir):
    return np.asarray(
        FoamCase(case_dir)[float(tdir)]["psi"].internal_field, dtype=float
    )


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--iters", type=int, default=40)
    ap.add_argument("--eps", default="1e-6,1e-5,1e-4")
    ap.add_argument("--band", type=float, default=3.0,
                    help="perturbation support: |psi0| < band*h")
    ap.add_argument("--steps", type=int, default=4,
                    help="solver steps per operator application. MUST be >= 2: "
                         "the step advects psi with U^n, so a ONE-step map from "
                         "rest is exactly the identity for psi (measured: the "
                         "response equals the perturbation to 4e-18, machine "
                         "roundoff, at every amplitude over four decades). The "
                         "quasi-static loop needs two steps to close -- psi -> "
                         "force -> u, then u -> psi.")
    ap.add_argument("--dt-divisor", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--work", default=None)
    a = ap.parse_args(argv)

    src = os.path.abspath(a.case)
    work = os.path.abspath(a.work or (src + "_spectrum"))
    if os.path.exists(work):
        shutil.rmtree(work)
    shutil.copytree(src, work)
    # strip any previous output times, keep 0
    for t, d in _time_dirs(work):
        if t > 0:
            shutil.rmtree(os.path.join(work, d))

    cd = os.path.join(work, "system", "controlDict")
    with open(cd) as fh:
        txt = fh.read()
    dt0 = None
    for ln in txt.splitlines():
        s = ln.strip()
        if s.startswith("deltaT") and ";" in s:
            dt0 = float(s.split()[1].rstrip(";"))
            break
    if dt0 is None:
        print("could not read deltaT", file=sys.stderr)
        return 1
    dt = dt0 / a.dt_divisor
    horizon = dt * a.steps

    _set_entries(cd, {
        "deltaT": repr(dt),
        "endTime": repr(horizon),
        "startTime": "0",
        "startFrom": "startTime",
        "stopAt": "endTime",
        "writeControl": "timeStep",
        "writeInterval": str(a.steps),
        "adjustTimeStep": "no",
        "maxDeltaT": repr(dt),
        "purgeWrite": "0",
        "writeFormat": "ascii",
        "writePrecision": "15",
        "timeFormat": "general",
        "timePrecision": "12",
    })

    env = dict(os.environ)
    # pristine 0/ to restore before every application
    zero0 = os.path.join(work, "0")
    zeroBak = os.path.join(work, "0.spectrumBase")
    shutil.copytree(zero0, zeroBak)

    psi0 = _read_psi(work, "0")
    ncells = psi0.size

    # cell size from the mesh bounding box (2D hex cases: nz = 1)
    h = None
    try:
        with open(os.path.join(work, "case_params.json")) as fh:
            N = int(json.load(fh)["tokens"]["N_CELLS"])
        h = 0.01 / N
    except Exception:
        h = math.sqrt(1e-4 / ncells)
    support = np.abs(psi0) < a.band * h
    nsup = int(support.sum())
    print(f"[spectrum] cells={ncells} h={h:.6g} dt={dt:.6g} steps/apply={a.steps} "
          f"support={nsup} cells (|psi0| < {a.band}h)")

    def restore():
        shutil.rmtree(zero0)
        shutil.copytree(zeroBak, zero0)
        for t, d in _time_dirs(work):
            if t > 0:
                shutil.rmtree(os.path.join(work, d))

    def apply_step(psi_field):
        restore()
        f = FoamCase(work)[0]["psi"]
        f.internal_field = psi_field.tolist()
        rc, log = _run_step(work, env)
        ts = [x for x in _time_dirs(work) if x[0] > 0]
        if rc != 0 or not ts:
            print(log[-2500:], file=sys.stderr)
            raise RuntimeError(f"solver failed rc={rc}")
        return _read_psi(work, ts[-1][1])

    print("[spectrum] base trajectory ...")
    psi_base_after = apply_step(psi0)

    rows = []
    for eps_s in a.eps.split(","):
        eps = float(eps_s)
        rng = np.random.default_rng(20260815)
        v = np.zeros(ncells)
        v[support] = rng.uniform(-1.0, 1.0, nsup)
        v /= np.linalg.norm(v)

        lam = float("nan")
        hist = []
        for it in range(a.iters):
            psi_p = psi0 + eps * h * v          # eps is RELATIVE to h
            psi_a = apply_step(psi_p)
            Mv = (psi_a - psi_base_after) / (eps * h)
            Mv[~support] = 0.0
            nrm = np.linalg.norm(Mv)
            if nrm == 0.0:
                print("[spectrum] zero response -- eps below solver noise")
                break
            # For a NON-symmetric operator the power method converges v onto the
            # dominant eigenvector and ||Mv|| onto |lambda|; the Rayleigh
            # quotient converges too but is not a bound on the way. Both are
            # reported so a complex dominant pair (oscillating estimates) is
            # visible rather than averaged away.
            ray = float(v @ Mv)
            lam = nrm
            hist.append(lam)
            dev = float(np.linalg.norm(Mv - v))
            print(f"    it{it:>3d} |Mv|={nrm:.8f} rayleigh={ray:>11.8f} "
                  f"|Mv-v|={dev:.3e}")
            v = Mv / nrm
            if it >= 8 and abs(hist[-1] - hist[-2]) < 1e-7 * abs(hist[-1]):
                break
        horizon_t = a.steps * dt
        r = math.log(abs(lam)) / horizon_t if lam and abs(lam) > 0 else float("nan")
        print(f"[spectrum] eps={eps:>8g}*h  applications={len(hist):>3d}  "
              f"lambda^{a.steps}={lam:.8f}  r=ln(lambda)/({a.steps}*dt)"
              f"={r:>12.4f} 1/s")
        rows.append({
            "eps_over_h": eps, "lambda_k": lam, "rayleigh": ray, "r": r,
            "dt": dt, "steps_per_apply": a.steps, "applications": len(hist),
            "n_cells": ncells, "n_support": nsup, "h": h,
        })
        # keep the last eigenvector for the finest converged eps
        np.save(os.path.join(work, f"eigvec_eps{eps:g}.npy"), v)

    if a.out:
        import csv
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"[spectrum] wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
