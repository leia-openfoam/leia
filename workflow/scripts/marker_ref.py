#!/usr/bin/env python3
"""Marker-traced quasi-exact interface reference for the 2D shear vortex.

The prescribed velocity is ANALYTIC, so the exact interface at any time t is
the image of the initial circle under the flow map. The classical reversed
benchmark dodges this by measuring only at t = T (where the exact solution is
the initial circle again) -- at the price of reversal cancellation. Here we
integrate marker points through the analytic field with a high-order ODE
solver (DOP853, rtol 1e-10), which gives a reference at ANY time, reversal or
not.

Shape error per case and snapshot time:

    E_shape(t) = sum_c |alpha_c(t) - f_ref,c(t)| * V_c  /  A_ref

where f_ref,c is the reference indicator function block-averaged onto the case
mesh from a fine raster (even-odd scanline fill of the marker polygon at
RASTER x RASTER pixel resolution >> h) and A_ref = pi R^2 is the conserved
reference area. The reference raster depends only on (T, oscillation, t), so it
is computed once and shared by every model and resolution.

Matches leia's shear2D velocity model exactly
(src/leiaLevelSet/velocityModel/prescribedVelocityModels.C):
    u =  sin(2 pi y) sin^2(pi x)
    v = -sin(2 pi x) sin^2(pi y)
scaled by cos(pi t / tau) when oscillation is on (tau = END_TIME), and the
2Dvortex initial circle: R = 0.15 centred at (0.5, 0.75).

Caveat: NMARK markers resolve the stretched filament down to a segment spacing
of roughly stretch * 2 pi R / NMARK; for heavy winding (T = 8 reversed at
t = T/2) treat the value as a diagnostic, not a convergence-grade metric.
"""
import csv
import glob
import json
import os

import numpy as np

CENTER = (0.5, 0.75)
RADIUS = 0.15
NMARK = 8192          # marker points on the initial circle
RASTER = 1024         # reference raster resolution (pixels across [0,1])


def _vel(t, xy, tau, oscillating):
    x = xy[0::2]
    y = xy[1::2]
    g = np.cos(np.pi*t/tau) if oscillating else 1.0
    out = np.empty_like(xy)
    out[0::2] = g*np.sin(2*np.pi*y)*np.sin(np.pi*x)**2
    out[1::2] = -g*np.sin(2*np.pi*x)*np.sin(np.pi*y)**2
    return out


def trace(times, tau, oscillating, n=NMARK):
    """Marker polygon (n, 2) at each requested time. times must include only
    t >= 0; t = 0 returns the initial circle exactly."""
    from scipy.integrate import solve_ivp
    th = 2*np.pi*np.arange(n)/n
    xy0 = np.empty(2*n)
    xy0[0::2] = CENTER[0] + RADIUS*np.cos(th)
    xy0[1::2] = CENTER[1] + RADIUS*np.sin(th)
    want = sorted({float(t) for t in times})
    out = {}
    if want and want[0] == 0.0:
        out[0.0] = xy0.reshape(-1, 2).copy()
        want = want[1:]
    if want:
        sol = solve_ivp(_vel, (0.0, want[-1]), xy0, t_eval=want,
                        args=(tau, oscillating), method="DOP853",
                        rtol=1e-10, atol=1e-12)
        for k, t in enumerate(sol.t):
            out[float(t)] = sol.y[:, k].reshape(-1, 2)
    return out


def rasterize(poly, M=RASTER):
    """Even-odd scanline fill of the closed marker polygon on an M x M
    cell-centre raster of the unit box. Returns a float mask (M, M)."""
    x = poly[:, 0]
    y = poly[:, 1]
    x2 = np.roll(x, -1)
    y2 = np.roll(y, -1)
    keep = y != y2                      # horizontal edges never cross a row
    xa, ya = x[keep], y[keep]
    xb, yb = x2[keep], y2[keep]
    mask = np.zeros((M, M), bool)
    xc = (np.arange(M) + 0.5)/M
    for j in range(M):
        yy = (j + 0.5)/M
        cross = (ya <= yy) != (yb <= yy)
        if not cross.any():
            continue
        xs = np.sort(xa[cross] + (yy - ya[cross])
                     * (xb[cross] - xa[cross])/(yb[cross] - ya[cross]))
        for k in range(0, len(xs) - 1, 2):
            i0 = np.searchsorted(xc, xs[k])
            i1 = np.searchsorted(xc, xs[k+1])
            mask[j, i0:i1] = True
    return mask.astype(float)


def block_average(mask, N):
    """Reference cell volume fractions (N, N) from the fine raster."""
    M = mask.shape[0]
    if M % N:
        raise ValueError(f"raster {M} not divisible by N={N}")
    b = M // N
    return mask.reshape(N, b, N, b).mean(axis=(1, 3))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _cases(study_dir):
    out = []
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*",
                                                   "case_params.json"))):
        try:
            meta = json.load(open(meta_path))
        except Exception:
            continue
        if meta.get("dims") != 2:
            continue
        tok = meta.get("tokens", {})
        T, N = _f(tok.get("END_TIME")), _f(tok.get("N_CELLS"))
        if T is None or N is None:
            continue
        out.append({
            "dir": os.path.dirname(meta_path),
            "model": tok.get("VELOCITY_EXTENSION", "none"),
            "osc": str(tok.get("OSCILLATION", "on")).lower(),
            "T": T, "N": int(N),
        })
    return out


def _alpha_at(case_dir, targets):
    """alpha (1D arrays) at the snapshot times nearest each target."""
    from foamlib import FoamCase
    case = FoamCase(case_dir)
    times = [float(t.time) for t in case]
    if len(times) < 2:
        return None
    ncells = np.asarray(case[0]["psi"].internal_field, float).size
    out = {}
    for tgt in targets:
        k = min(range(len(times)), key=lambda i: abs(times[i] - tgt))
        a = np.asarray(case[k]["alpha"].internal_field, float)
        if a.size == 1:                 # fully-lost phase: uniform field
            a = np.full(ncells, float(a.reshape(-1)[0]))
        out[tgt] = a
    return out


def shape_report(study_dir, figdir, prefix=""):
    """Shape error vs the marker reference for every case at t = 0, T/2, T.
    Writes <study>/<prefix>shape_reference.csv and a convergence figure.
    Returns the list of figure paths."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = _cases(study_dir)
    if not cases:
        return []
    aref = np.pi*RADIUS**2

    # Reference rasters, shared across models and resolutions.
    refs = {}                            # (T, osc, t) -> RASTER mask
    for T, osc in sorted({(c["T"], c["osc"]) for c in cases}):
        targets = [0.0, 0.5*T, T]
        polys = trace(targets, T, osc == "on")
        for t, poly in polys.items():
            refs[(T, osc, t)] = rasterize(poly)
        print(f"[marker_ref] traced+rasterized reference: T={T:g} osc={osc}")

    rows = []
    for c in sorted(cases, key=lambda c: (c["model"], c["N"])):
        targets = [0.0, 0.5*c["T"], c["T"]]
        try:
            alphas = _alpha_at(c["dir"], targets)
        except Exception as exc:  # noqa: BLE001
            print(f"[marker_ref] {c['dir']}: {type(exc).__name__}: {exc}")
            continue
        if alphas is None:
            continue
        N = c["N"]
        Vc = (1.0/N)**2
        for t in targets:
            key = (c["T"], c["osc"], t)
            if key not in refs:
                continue
            fref = block_average(refs[key], N)
            a = alphas[t]
            if a.size != N*N:
                continue
            e = float(np.abs(a.reshape(N, N) - fref).sum()*Vc/aref)
            rows.append({"model": c["model"], "N": N, "T": c["T"],
                         "t": t, "shapeErrRel": e})

    if not rows:
        return []
    csv_path = os.path.join(study_dir, f"{prefix}shape_reference.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "N", "T", "t",
                                           "shapeErrRel"])
        w.writeheader()
        w.writerows(rows)
    print(f"[marker_ref] wrote {csv_path}: {len(rows)} rows")

    # Convergence figure: shape error vs h, panels t = T/2 and t = T.
    models = sorted({r["model"] for r in rows})
    Ts = sorted({r["T"] for r in rows})
    cmap = plt.get_cmap("tab10")
    mcol = {m: cmap(i % 10) for i, m in enumerate(models)}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), squeeze=False)
    for ci, frac in enumerate((0.5, 1.0)):
        ax = axes[0][ci]
        allh, alle = [], []
        for m in models:
            for T in Ts:
                pts = sorted((1.0/r["N"], r["shapeErrRel"]) for r in rows
                             if r["model"] == m and r["T"] == T
                             and abs(r["t"] - frac*T) < 1e-9
                             and r["shapeErrRel"] > 0)
                if not pts:
                    continue
                h, e = zip(*pts)
                ax.loglog(h, e, marker="o", ms=4, lw=1.4, color=mcol[m],
                          ls=["-", "--", ":", "-."][Ts.index(T) % 4],
                          label=(m + (f", T={T:g}" if len(Ts) > 1 else ""))
                          if ci == 0 else None)
                allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.7, alpha=0.5)
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.7, alpha=0.5)
        ax.set_xlabel(r"$h$")
        ax.set_title("t = T/2" if frac == 0.5 else "t = T")
        ax.grid(True, which="both", alpha=0.3)
    axes[0][0].set_ylabel(r"$\|\alpha - \chi_{\mathrm{ref}}\|_{L_1} / A_{\mathrm{ref}}$")
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center",
                   ncol=min(6, len(labels)), fontsize=7, frameon=False)
    fig.suptitle("Shape error vs marker-traced reference (no reversal needed)",
                 fontsize=12)
    p = os.path.join(figdir, f"{prefix}shape_reference.png")
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return [p]


if __name__ == "__main__":
    import sys
    shape_report(sys.argv[1],
                 sys.argv[2] if len(sys.argv) > 2
                 else os.path.join(sys.argv[1], "figures"))
