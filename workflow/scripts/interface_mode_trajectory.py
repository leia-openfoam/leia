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
    """Fit y = A0 exp(r t) cos(w t + p).

    Done in SCALED variables tau = (t - t0)/span and yhat = y/max|y|, so the
    exponent search bracket is dimensionless: rho = r*span with |rho| <= 30
    (exp(30) ~ 1e13, comfortably inside double range) and Omega = w*span. An
    ABSOLUTE bracket on r overflows exp() as soon as the window is lengthened
    -- which is exactly what happened when the audit went from 0.55 to 2.5
    periods: LAPACK received non-finite entries and the whole fit was lost.
    Scaling makes the fit independent of how long the window is.
    """
    ts = np.asarray(ts, float)
    ys = np.asarray(ys, float)
    span = float(ts[-1] - ts[0])
    yscale = float(np.max(np.abs(ys)))
    if span <= 0.0 or yscale <= 0.0 or not np.all(np.isfinite(ys)):
        return float("nan"), float("nan"), float("nan"), float("nan")
    tau = (ts - ts[0]) / span
    yh = ys / yscale

    def resid_at(rho, Cc, Ss):
        E = np.exp(rho * tau)
        A = np.column_stack([E * Cc, E * Ss])
        if not np.all(np.isfinite(A)):
            return float("inf"), np.zeros(2)
        coef, *_ = np.linalg.lstsq(A, yh, rcond=None)
        return float(np.linalg.norm(yh - A @ coef)), coef

    best = None
    # up to ~6 full periods across the window
    for Om in np.linspace(0.0, 6 * 2 * math.pi, 3000):
        Cc, Ss = np.cos(Om * tau), np.sin(Om * tau)
        lo, hi = -30.0, 30.0
        for _ in range(70):
            m1 = lo + (hi - lo) * 0.382
            m2 = lo + (hi - lo) * 0.618
            if resid_at(m1, Cc, Ss)[0] < resid_at(m2, Cc, Ss)[0]:
                hi = m2
            else:
                lo = m1
        rho = 0.5 * (lo + hi)
        res, coef = resid_at(rho, Cc, Ss)
        if best is None or res < best[0]:
            best = (res, rho, Om, coef)

    res, rho, Om, coef = best
    r = rho / span
    w = Om / span
    A0 = math.hypot(coef[0], coef[1]) * yscale
    rel = res / max(float(np.linalg.norm(yh)), 1e-300)
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
        w_th = math.sqrt(sigma * m * (m * m - 1) / (rho * R ** 3))

        # ALWAYS dump the raw amplitude trajectory next to the fit. A fitted r
        # is not inspectable and a failed optimisation is not obviously
        # different from a physical result (a viscosity sweep produced
        # omega = 1.3 against a theoretical 653 and an r of +3000 before this
        # was here). One file per mode; look at it before quoting anything.
        if a.out:
            tj = a.out.replace(".csv", "") + f"_traj_m{m}.csv"
            with open(tj, "w", newline="") as fh:
                wtr = csv.writer(fh)
                wtr.writerow(["t", "a_m", "periods_theory"])
                for tt, aa in zip(ts, amps):
                    wtr.writerow([f"{tt:.10g}", f"{aa:.10e}",
                                  f"{tt*w_th/(2*math.pi):.4f}"])
            print(f"    m={m} trajectory -> {os.path.basename(tj)}")

        r, w, A0, rel_res = _fit_decaying_cosine(ts, amps)

        # WINDOW CONVERGENCE. A 4-parameter fit A0 exp(r t) cos(w t + phi) over
        # LESS THAN one period has a near-degenerate direction: envelope growth
        # trades against phase, so a tiny residual does NOT mean r is
        # determined. (Learned the hard way -- an r0 extrapolation from
        # 0.55-period windows gave +0.03 by 2-point Richardson and +3.3 by a
        # 4-point fit of the same data, against an expected -5.) Refit over
        # growing prefixes: r must settle, and the window must reach >= 1.5
        # periods, before the number is quotable.
        T_ref = 2 * math.pi / w_th if w_th > 0 else float("inf")
        sweep = []
        for frac in (0.4, 0.6, 0.8, 1.0):
            k = max(8, int(frac * len(ts)))
            rr, ww, _, res = _fit_decaying_cosine(ts[:k], amps[:k])
            sweep.append((ts[k - 1] / T_ref, rr, res))
        print(f"    m={m} window sweep (periods, r_m, fit rel):  " +
              "  ".join(f"[{p:.2f}P {rr:+.1f} {res:.3f}]" for p, rr, res in sweep))
        r_spread = max(s[1] for s in sweep) - min(s[1] for s in sweep)
        periods = sweep[-1][0]
        om_ok = (w_th > 0) and (0.75 <= w / w_th <= 1.25)
        quotable = (periods >= 1.5) and om_ok and \
            (r_spread <= max(1.0, 0.15 * abs(r)))
        print(f"    m={m} {'QUOTABLE' if quotable else 'NOT QUOTABLE'}: "
              f"{periods:.2f} periods, r spread {r_spread:.2f}")
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
