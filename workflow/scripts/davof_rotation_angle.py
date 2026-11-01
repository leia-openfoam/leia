#!/usr/bin/env python3
"""The face-fraction error against the direction of travel relative to the
mesh axes (config/davof/sphereRotation3D_hexFine*.yaml, 2026-11-02).

Per case (sphereRotation3D; the rotation, the circular translation without
spin at omega and at 2 omega, the straight translations), per step:
  t, the turned angle omega t, the direction of travel of the sphere's
  centre phi = atan2(U_z, U_x) (degrees), its angle to the nearest mesh axis
  psi in [0, 45] degrees, L1(alphaf), L1(x), Eshape
(written with --steps out.csv). The centre's velocity follows from the
case's tokens: U = omega (-d_z, 0, d_x) + (TPF_UX, 0, TPF_UZ) with
d(t) = d0 cos(omega t) + (k x d0) sin(omega t), d0 = c0 - axisPoint,
k = (0 -1 0) (the rotation about the axis through axisPoint with the angular
velocity (0 -omega 0); Rodrigues' formula, O. Rodrigues, J. Math. Pures
Appl. 5 (1840) 380-440, no DOI); on the straight translations the axis is
the sphere's centre, d = 0, and U is the arm's.

Summary per case (the pre-registered measures of the configs):
  - the smooth part S(t): the centred moving average of L1(alphaf) over a
    quarter of the turn, 2 pi/(4 omega) (an integer number of steps), which
    removes every harmonic of period T/4 and T/8 from S;
  - the ripple R = L1 - S fitted over the window where S is defined by
    c + sum_{m = 4, 8} (a_m cos(m omega t) + b_m sin(m omega t)) (least
    squares): r_m = sqrt(a_m^2 + b_m^2) / mean(L1) and the phase delta_m of
    A_m cos(m omega t - delta_m) in degrees (delta_4 = 0: the maxima at
    psi = 0, the travel along a mesh axis; +-180: at psi = 45 degrees);
  - the hump: the time of the maximum of S over the window, and per turn of
    a multi-turn case the maxima of S in each half of the window with the
    minimum of S over [0.4, 0.6] of the end time between them;
  - straight translations: the mean of L1(alphaf) over [0.1, 0.45] and the
    end value against the maximum.

Usage (on the cluster, from the leia clone):
    python3 workflow/scripts/davof_rotation_angle.py studies/davof/<study> [...] [--steps steps.csv] [--csv summary.csv]
"""
import argparse
import csv
import glob
import json
import math
import os

import numpy as np


def f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load(case):
    t = json.load(open(os.path.join(case, "case_params.json")))["tokens"]
    path = os.path.join(case, "advectorVoF.csv")
    rows = list(csv.DictReader(open(path))) if os.path.isfile(path) and os.path.getsize(path) else []
    return t, rows


def direction(tok, t):
    w = f(tok.get("TPF_OMEGA")) or 0.0
    ax, az = f(tok.get("TPF_AXIS_X")), f(tok.get("TPF_AXIS_Z"))
    cx, cz = f(tok.get("TPF_SPHERE_CX")), f(tok.get("TPF_SPHERE_CZ"))
    ux, uz = f(tok.get("TPF_UX")) or 0.0, f(tok.get("TPF_UZ")) or 0.0
    d0x, d0z = cx - ax, cz - az
    # k x d0 = (-d0z, 0, d0x)
    c, s = math.cos(w * t), math.sin(w * t)
    dx, dz = d0x * c - d0z * s, d0z * c + d0x * s
    Ux, Uz = w * (-dz) + ux, w * dx + uz
    phi = math.degrees(math.atan2(Uz, Ux))
    psi = abs(((phi + 45.0) % 90.0) - 45.0)
    return phi, psi, math.hypot(Ux, Uz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("studies", nargs="+")
    ap.add_argument("--steps", default="")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()
    steps_out, summ = [], []
    for st in a.studies:
        for case in sorted(glob.glob(os.path.join(st, "*_0*"))):
            if not os.path.isfile(os.path.join(case, "case_params.json")) or "." in os.path.basename(case):
                continue
            tok, rows = load(case)
            if not rows:
                continue
            N = int(float(tok["N_CELLS"]))
            w = f(tok.get("TPF_OMEGA")) or 0.0
            line = (abs(f(tok["TPF_AXIS_X"]) - f(tok["TPF_SPHERE_CX"])) < 1e-12
                    and abs(f(tok["TPF_AXIS_Z"]) - f(tok["TPF_SPHERE_CZ"])) < 1e-12)
            flow = "line" if line else tok.get("TPF_SETFLOW_MODE", "rotation")
            flow = {"function": "orbit"}.get(flow, flow)
            T = f(tok.get("END_TIME"))
            tt, L1, Lx, Es, ph, ps = [], [], [], [], [], []
            for r in rows:
                t, v = f(r.get("TIME")), f(r.get("L1_alphaf"))
                if t is None or v is None or v < 0:
                    continue
                phi, psi, sp = direction(tok, t)
                tt.append(t); L1.append(v); ph.append(phi); ps.append(psi)
                Lx.append(f(r.get("LCentreDiff1"))); Es.append(f(r.get("Eshape")))
                steps_out.append(dict(study=os.path.basename(st), case=os.path.basename(case), N=N, flow=flow,
                                      omega=w, t=t, omega_t_deg=math.degrees(w * t) if not line else 0.0,
                                      phi=phi, psi=psi, speed=sp, L1_alphaf=v, L1_x=Lx[-1], Eshape=Es[-1]))
            tt, L1 = np.array(tt), np.array(L1)
            dt = float(np.median(np.diff(tt)))
            rec = dict(study=os.path.basename(st), case=os.path.basename(case), N=N, flow=flow, omega=w,
                       arm=tok.get("TPF_VELOCITY_ARM"), steps=len(tt), T=T,
                       L1_first=float(L1[0]), L1_end=float(L1[-1]), L1_max=float(L1.max()),
                       t_max=float(tt[np.argmax(L1)]),
                       Lx_end=Lx[-1], Es_end=Es[-1], Es_25=None)
            if line:
                m = (tt >= 0.1 - 1e-12) & (tt <= 0.45 + 1e-12)
                rec.update(mean_01_045=float(L1[m].mean()) if m.any() else None,
                           end_over_max=float(L1[-1] / L1.max()), psi_mean=float(np.mean(ps)))
            elif w > 0:
                Tq = 2 * math.pi / (4 * w)          # a quarter of the turn
                n = int(round(Tq / dt))
                if n >= 4 and len(L1) > n + 4:
                    # the trapezoidal mean over one quarter turn, centred on
                    # the sample (n + 1 samples, half weights at the ends;
                    # exact for the harmonics of period T/4 and T/8); an odd
                    # n takes the plain mean of n samples (half a step off)
                    if n % 2 == 0:
                        ker = np.ones(n + 1); ker[0] = ker[-1] = 0.5; ker /= n
                    else:
                        ker = np.ones(n) / n
                    S = np.convolve(L1, ker, mode="valid")
                    h = (len(ker) - 1) // 2
                    ts = tt[h:len(tt) - h]
                    Lw = L1[h:len(tt) - h]
                    R = Lw - S
                    cols = [np.ones_like(ts)]
                    for mm in (4, 8):
                        cols += [np.cos(mm * w * ts), np.sin(mm * w * ts)]
                    X = np.vstack(cols).T
                    coef, *_ = np.linalg.lstsq(X, R, rcond=None)
                    mean = float(Lw.mean())
                    for j, mm in enumerate((4, 8)):
                        am, bm = coef[1 + 2 * j], coef[2 + 2 * j]
                        rec[f"r{mm}"] = float(math.hypot(am, bm) / mean)
                        rec[f"delta{mm}_deg"] = float(math.degrees(math.atan2(bm, am)))
                    resid = R - X @ coef
                    rec["ripple_resid_rel"] = float(np.sqrt(np.mean(resid ** 2)) / mean)
                    rec["S_tpeak"] = float(ts[np.argmax(S)])
                    rec["S_max"] = float(S.max()); rec["S_mean"] = float(S.mean())
                    half = 0.5 * T
                    m1, m2 = ts <= half, ts > half
                    mid = (ts >= 0.4 * T) & (ts <= 0.6 * T)
                    if m1.any() and m2.any() and mid.any():
                        rec["S_tpeak_first_half"] = float(ts[m1][np.argmax(S[m1])])
                        rec["S_tpeak_second_half"] = float(ts[m2][np.argmax(S[m2])])
                        rec["S_peak1_minus_midmin_rel"] = float((S[m1].max() - S[mid].min()) / mean)
                        rec["S_peak2_minus_midmin_rel"] = float((S[m2].max() - S[mid].min()) / mean)
                    rec["window"] = f"[{ts[0]:.4f}, {ts[-1]:.4f}]"
            summ.append(rec)
    keys = ["study", "N", "flow", "omega", "steps", "L1_first", "L1_end", "L1_max", "t_max", "S_tpeak",
            "S_tpeak_first_half", "S_tpeak_second_half", "S_peak1_minus_midmin_rel", "S_peak2_minus_midmin_rel",
            "r4", "delta4_deg", "r8", "delta8_deg", "ripple_resid_rel", "mean_01_045", "end_over_max",
            "Lx_end", "Es_end"]
    fmt = lambda v: ("%.4g" % v) if isinstance(v, float) else str(v)
    print("\t".join(keys))
    for r in summ:
        print("\t".join(fmt(r.get(k)) for k in keys))
    if a.csv and summ:
        allk = sorted({k for r in summ for k in r})
        with open(a.csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=allk); wr.writeheader(); wr.writerows(summ)
    if a.steps and steps_out:
        with open(a.steps, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(steps_out[0].keys())); wr.writeheader(); wr.writerows(steps_out)


if __name__ == "__main__":
    main()
