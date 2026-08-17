#!/usr/bin/env python3
"""Per-mode GROWTH RATE and FREQUENCY of the coupled loop, from the trajectory
of a single interface-displacement mode on the shipped solver.

WHY A TRAJECTORY AND NOT AN AMPLITUDE RATIO. A single-horizon ratio
a_m(T)/a_m(0) conflates growth with oscillation. A mode-m interface
displacement on a droplet is a capillary wave: it oscillates at
omega_m = sqrt(sigma*m*(m^2-1)/(rho*R^3)) whether or not it grows. Measured at
N=64 over 8 steps the m=8 ratio read 0.548, which looks like strong damping and
is almost entirely cos(omega_8 T) = 0.515 -- the mode had simply rotated 0.16 of
a period. Separating the two requires sampling the trajectory and fitting

    a_m(t) = A0 * exp(r_m * t) * cos(omega_m * t + phi)

so r_m (the growth rate, the quantity that decides stability) and omega_m (the
physics, which must reproduce the capillary dispersion relation) come out as
separate numbers. omega_m agreeing with the analytic dispersion is the harness's
own validity check.

PERTURBATION SPACE. Interface displacement, not raw delta-psi: perturbing psi by
eps*h*cos(m theta) (a function of theta alone, so |grad psi| = 1 makes it a
radial displacement of the zero set). This avoids the level set's GAUGE
subspace, every member of which -- delta-psi vanishing on Gamma -- leaves alpha,
the force and the whole trajectory exactly unchanged and is therefore neutral by
construction. A power iteration on raw delta-psi converges onto that subspace
(measured: lambda -> 0.9996 at N=64) and never sees the instability.

Usage:
  interface_mode_trajectory.py <prepared_case> [--modes 2,4,6,8,12]
       [--steps 200] [--sample 4] [--eps 1e-2] [--dt-divisor 1] [--out t.csv]
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np
from foamlib import FoamCase

SOLVER = "leiaSemiLagrangianLevelSetTwoPhaseFoam"
BASHRC = os.path.expanduser("~/OpenFOAM/OpenFOAM-v2512/etc/bashrc")


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


def _set_controls(cd, dt, nsteps, sample):
    t = open(cd).read()
    subs = {"deltaT": repr(dt), "endTime": repr(dt * nsteps),
            "writeControl": "timeStep", "writeInterval": str(sample),
            "adjustTimeStep": "no", "maxDeltaT": repr(dt),
            "writeFormat": "ascii", "writePrecision": "15",
            "purgeWrite": "0", "startFrom": "startTime", "startTime": "0"}
    for k, v in subs.items():
        if re.search(rf"^\s*{k}\b", t, flags=re.M):
            t = re.sub(rf"^\s*{k}\b.*$", f"{k:<16}{v};", t, flags=re.M)
        else:
            t = t.rstrip() + f"\n{k:<16}{v};\n"
    open(cd, "w").write(t)


def _fit_decaying_cosine(ts, ys):
    """Fit y = A0 exp(r t) cos(w t + p). Grid-search w, linear-solve the rest."""
    ts = np.asarray(ts, float)
    ys = np.asarray(ys, float)
    span = ts[-1] - ts[0]
    best = None
    # frequency grid: up to ~6 full periods over the sampled window
    for w in np.linspace(0.0, 6 * 2 * math.pi / max(span, 1e-30), 4000):
        # y = e^{rt}(a cos wt + b sin wt): nonlinear in r, so scan r too via
        # a coarse-to-fine on the log-envelope after projecting out (a,b).
        Cc, Ss = np.cos(w * ts), np.sin(w * ts)
        # solve min ||y - e^{rt}(a Cc + b Ss)|| for r by 1-D golden search
        def resid(r):
            E = np.exp(r * ts)
            A = np.column_stack([E * Cc, E * Ss])
            coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
            return float(np.linalg.norm(ys - A @ coef)), coef
        lo, hi = -2.0e5, 2.0e5
        for _ in range(60):
            m1 = lo + (hi - lo) * 0.382
            m2 = lo + (hi - lo) * 0.618
            if resid(m1)[0] < resid(m2)[0]:
                hi = m2
            else:
                lo = m1
        r = 0.5 * (lo + hi)
        res, coef = resid(r)
        if best is None or res < best[0]:
            best = (res, r, w, coef)
    res, r, w, coef = best
    A0 = math.hypot(coef[0], coef[1])
    rel = res / max(np.linalg.norm(ys), 1e-300)
    return r, w, A0, rel


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--modes", default="2,4,6,8,12")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--sample", type=int, default=4)
    ap.add_argument("--eps", type=float, default=1e-2)
    ap.add_argument("--dt-divisor", type=float, default=1.0)
    ap.add_argument("--cut", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--work", default=None)
    a = ap.parse_args(argv)

    src = os.path.abspath(a.case)
    work = os.path.abspath(a.work or (src + "_traj"))
    if os.path.exists(work):
        shutil.rmtree(work)
    shutil.copytree(src, work)
    for t, d in _time_dirs(work):
        if t > 0:
            shutil.rmtree(os.path.join(work, d))

    cd = os.path.join(work, "system", "controlDict")
    dt0 = next(float(l.split()[1].rstrip(";")) for l in open(cd)
               if l.strip().startswith("deltaT") and ";" in l)
    dt = dt0 / a.dt_divisor
    _set_controls(cd, dt, a.steps, a.sample)

    if not os.path.exists(os.path.join(work, "0", "C")):
        subprocess.run(["bash", "-lc",
                        f"source {BASHRC} >/dev/null 2>&1 && cd {work} && "
                        f"postProcess -func writeCellCentres -time 0 >/dev/null 2>&1"],
                       capture_output=True, text=True)
    C = np.asarray(FoamCase(work)[0]["C"].internal_field, dtype=float)
    psi0 = np.asarray(FoamCase(work)[0]["psi"].internal_field, dtype=float)
    try:
        N = int(json.load(open(os.path.join(work, "case_params.json")))
                ["tokens"]["N_CELLS"])
        h = 0.01 / N
    except Exception:
        N, h = 0, math.sqrt(1e-4 / psi0.size)

    ctr = C[psi0 < 0].mean(axis=0)
    rel = C - ctr
    th_all = np.arctan2(rel[:, 1], rel[:, 0])
    cutm = np.abs(psi0) < a.cut * h
    th = th_all[cutm]
    R = float(np.mean(np.hypot(rel[cutm, 0], rel[cutm, 1])))

    zero0, zbak = os.path.join(work, "0"), os.path.join(work, "0.trajBase")
    shutil.copytree(zero0, zbak)

    def run(delta):
        shutil.rmtree(zero0)
        shutil.copytree(zbak, zero0)
        for t, d in _time_dirs(work):
            if t > 0:
                shutil.rmtree(os.path.join(work, d))
        FoamCase(work)[0]["psi"].internal_field = (psi0 + delta).tolist()
        r = subprocess.run(["bash", "-lc",
                            f"source {BASHRC} >/dev/null 2>&1 && cd {work} && {SOLVER}"],
                           capture_output=True, text=True)
        ts = [x for x in _time_dirs(work) if x[0] > 0]
        if not ts:
            print((r.stdout + r.stderr)[-1500:], file=sys.stderr)
            raise RuntimeError("solver produced no output")
        c = FoamCase(work)
        return [(t, np.asarray(c[float(d)]["psi"].internal_field, dtype=float))
                for t, d in ts]

    print(f"[traj] N={N} h={h:.6g} dt={dt:.6g} steps={a.steps} "
          f"sample every {a.sample}  R={R:.6g} m  cut cells={int(cutm.sum())}")
    print("[traj] base trajectory ...")
    base = dict(run(np.zeros_like(psi0)))

    rho, sigma = 1000.0, 0.0728
    rows = []
    print(f"\n{'m':>3} {'r_m [1/s]':>12} {'w_meas':>11} {'w_theory':>11} "
          f"{'w ratio':>8} {'fit rel':>8}")
    for ms in a.modes.split(","):
        m = int(ms)
        traj = run(a.eps * h * np.cos(m * th_all))
        ts, amps = [], []
        for t, psit in traj:
            d = (psit - base[t])[cutm]
            # projection onto cos(m theta) over the cut cells
            cc = np.cos(m * th)
            amps.append(float(2.0 * (d @ cc) / len(th)))
            ts.append(t)
        r, w, A0, rel_res = _fit_decaying_cosine(ts, amps)
        w_th = math.sqrt(sigma * m * (m * m - 1) / (rho * R ** 3))
        print(f"{m:>3} {r:>12.2f} {w:>11.1f} {w_th:>11.1f} "
              f"{w/w_th if w_th else float('nan'):>8.3f} {rel_res:>8.4f}")
        rows.append({"N": N, "h": h, "dt": dt, "steps": a.steps, "mode": m,
                     "r_m": r, "omega_measured": w, "omega_theory": w_th,
                     "omega_ratio": w / w_th if w_th else float("nan"),
                     "fit_rel_residual": rel_res, "A0": A0,
                     "eps_over_h": a.eps, "R": R})

    if a.out:
        with open(a.out, "w", newline="") as fh:
            wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wtr.writeheader()
            wtr.writerows(rows)
        print(f"\n[traj] wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
