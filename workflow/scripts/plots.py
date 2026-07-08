#!/usr/bin/env python3
"""Post-processing figures for a leia study, using foamlib to read the cases.

Produces, into ``<study>/figures/`` and (copied) into the reveal slide folder:
  * ``convergence_<model>.png`` -- one triptych per velocity-extension model
    (gradient / shape / volume error vs h, one line per T, O(h)/O(h^2) guides),
    with the velocity-extension model in the figure title;
  * ``interface_grid.png`` -- reversed-interface montage: initial (dashed) vs
    final (solid) alpha=0.5 contour for each (T, h) case  [2D structured cases,
    a single representative velocity-extension model];
  * ``alpha_evo_<model>_T<T>.png`` -- volume-fraction field at t = 0, T/2
    (max deformation) and T (reversed), rows = resolution;
  * ``graderr_<model>_T<T>.png`` -- signed-distance error |grad psi|-1 at the
    same snapshots (where the SDF property breaks);
  * ``volume_history_T8.png`` -- relative phase-volume error vs t/T (T=8, none).

Driven by the curated ``*_errors.csv`` (columns velocityExtension,T,h,
gradientError,shapeError,volumeError) written by aggregate.py, plus the per-case
``case_params.json`` for the montage geometry.
"""
import csv
import glob
import json
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from foamlib import FoamCase

# (column, subplot title) for the three convergence measures.
SERIES = [
    ("gradientError", r"$\||\nabla\psi|-1\|_{L_2}$"),
    ("shapeError",    "geometric shape error"),
    ("volumeError",   "volume-conservation error"),
]


def _rows(errors_csv):
    with open(errors_csv, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def convergence_triptych(rows, model, figdir, prefix=""):
    """One figure per velocity-extension model: gradient/shape/volume vs h,
    one line per period T, with O(h)/O(h^2) guides. Model name in the title."""
    sub = [r for r in rows if r.get("velocityExtension", "none") == model]
    if not sub:
        return None
    Ts = sorted({_f(r["T"]) for r in sub if _f(r["T"]) is not None})
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    for ax, (col, title) in zip(axes, SERIES):
        allh, alle = [], []
        for T in Ts:
            pts = sorted(
                (_f(r["h"]), _f(r[col])) for r in sub
                if _f(r["T"]) == T and _f(r.get(col)) not in (None, 0.0)
            )
            pts = [(h, e) for h, e in pts if h and e]
            if not pts:
                continue
            h, e = zip(*pts)
            ax.loglog(h, e, marker="o", ms=5, label=f"T = {T:g}")
            allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.8, alpha=0.6, label=r"$O(h)$")
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.8, alpha=0.6, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$"); ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"Velocity extension: {model}", fontsize=14)
    p = os.path.join(figdir, f"{prefix}convergence_{model}.png")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def _cases_2d(study_dir):
    """Discover materialized 2D cases from their case_params.json.
    Returns list of {dir, model, T, N}."""
    out = []
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
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
        # For SL studies the "model" that distinguishes cases is the
        # reconstruction, not the (constant) velocity extension.
        recon = tok.get("SL_RECONSTRUCTION", "")
        out.append({
            "dir": os.path.dirname(meta_path),
            "model": recon if recon else tok.get("VELOCITY_EXTENSION", "none"),
            "recon": recon,
            "cfl": _f(tok.get("CFL")),
            "T": T, "N": int(N),
        })
    return out


def interface_montage(study_dir, figdir, model="none", prefix="", steady=False):
    """Interface montage (initial dashed vs final solid alpha=0.5) over the
    (T, h) grid, for a single representative velocity-extension model.
    Domain geometry is inferred per case (square cells, x in [0,1])."""
    cases = [c for c in _cases_2d(study_dir) if c["model"] == model]
    if not cases:
        allc = _cases_2d(study_dir)
        if not allc:
            return None
        model = sorted({c["model"] for c in allc})[0]
        cases = [c for c in allc if c["model"] == model]
    grid = {(c["T"], c["N"]): c["dir"] for c in cases}
    Ts = sorted({c["T"] for c in cases})
    Ns = sorted({c["N"] for c in cases})
    fig, axes = plt.subplots(
        len(Ts), len(Ns),
        figsize=(3.0*len(Ns), 2.0*len(Ts)), squeeze=False,
    )
    any_ok = False
    for ir, T in enumerate(Ts):
        for ic, N in enumerate(Ns):
            ax = axes[ir][ic]
            ax.set_xticks([]); ax.set_yticks([])
            cdir = grid.get((T, N))
            if not cdir:
                ax.axis("off"); continue
            try:
                case = FoamCase(cdir)
                a0 = np.asarray(case[0]["alpha"].internal_field, float)
                aT = np.asarray(case[-1]["alpha"].internal_field, float)
                if aT.size == 1:      # fully-lost phase: uniform field
                    aT = np.full(a0.size, float(aT.reshape(-1)[0]))
                Nx = N; Ny = a0.size // Nx           # square cells, x in [0,1]
                ymax = Ny / Nx
                a0 = a0.reshape(Ny, Nx); aT = aT.reshape(Ny, Nx)
                xc = (np.arange(Nx)+0.5)/Nx
                yc = (np.arange(Ny)+0.5)/Ny*ymax
                X, Y = np.meshgrid(xc, yc)
                ax.contour(X, Y, aT, [0.5], colors="#D85A30", linewidths=2.2)
                ax.contour(X, Y, a0, [0.5], colors="k", linestyles="--", linewidths=1.0)
                ax.set_xlim(0, 1); ax.set_ylim(0, ymax); ax.set_aspect("equal")
                ax.set_title(f"T={T:g}, h={1.0/Nx:.3g}", fontsize=8)
                any_ok = True
            except Exception as exc:
                ax.text(0.5, 0.5, "n/a", ha="center", va="center", fontsize=8)
                ax.axis("off")
                print(f"[plots] {cdir}: {type(exc).__name__}: {exc}")
    if not any_ok:
        plt.close(fig); return None
    flow = "Steady vortex" if steady else "Reversed vortex"
    fig.suptitle(
        f"{flow} ({model}): initial (dashed) vs final (solid) interface",
        fontsize=11)
    p = os.path.join(figdir, f"{prefix}interface_grid.png")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def indicator_comparison(rows, figdir):
    """Overlay convergence of the phase indicators (e.g. geometric vs
    detrixheAslam) on one triptych: gradient/shape/volume vs h, colour per T,
    line style + marker per indicator. Coinciding curves ⇒ the analytic
    (tolerance-free) indicator matches the geometric one. Also prints the max
    relative indicator disagreement per error (they integrate the same LLS plane,
    so this should be tiny)."""
    inds = sorted({r.get("phaseIndicator", "") for r in rows if r.get("phaseIndicator")})
    if len(inds) < 2:
        return None
    Ts = sorted({_f(r["T"]) for r in rows if _f(r["T"]) is not None})
    lstyles = ["-", "--", ":", "-."]
    markers = ["o", "x", "s", "^"]
    istyle = {ind: (lstyles[k % 4], markers[k % 4]) for k, ind in enumerate(inds)}
    cmap = plt.get_cmap("viridis")
    tcol = {T: cmap(i/max(1, len(Ts)-1)) for i, T in enumerate(Ts)}

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    for ax, (col, title) in zip(axes, SERIES):
        allh, alle = [], []
        for T in Ts:
            for ind in inds:
                pts = sorted(
                    (_f(r["h"]), _f(r[col])) for r in rows
                    if r.get("phaseIndicator") == ind and _f(r["T"]) == T
                    and _f(r.get(col)) not in (None, 0.0)
                )
                pts = [(h, e) for h, e in pts if h and e]
                if not pts:
                    continue
                h, e = zip(*pts)
                ls, mk = istyle[ind]
                ax.loglog(h, e, color=tcol[T], linestyle=ls, marker=mk, ms=5,
                          label=f"{ind}, T={T:g}")
                allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.7, alpha=0.5, label=r"$O(h)$")
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.7, alpha=0.5, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$"); ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle("Phase indicator: " + " vs ".join(inds)
                 + " (bulk shear vortex)", fontsize=14)
    p = os.path.join(figdir, "convergence_indicators.png")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(p, dpi=130); plt.close(fig)

    # Quantify indicator disagreement (same T,h; relative to the first indicator).
    ref = inds[0]
    for col, _ in SERIES:
        worst = 0.0
        for T in Ts:
            for r in rows:
                if r.get("phaseIndicator") != ref or _f(r["T"]) != T:
                    continue
                a = _f(r.get(col))
                b = next((_f(o.get(col)) for o in rows
                          if o.get("phaseIndicator") == inds[1]
                          and _f(o["T"]) == T and o["h"] == r["h"]), None)
                if a not in (None, 0.0) and b is not None:
                    worst = max(worst, abs(b-a)/abs(a))
        print(f"[plots] indicator max rel. diff ({inds[1]} vs {ref}) {col}: {worst:.3e}")
    return p


# --- field-evolution atlas: alpha and |grad psi|-1 at t = 0, T/2, T ----------
_SNAP_LABELS = ("t = 0", "t = T/2  (max deformation)", "t = T  (reversed)")
_SNAP_LABELS_STEADY = ("t = 0", "t = T/2", "t = T  (max winding)")


def _snapshot_indices(case, T):
    """Time indices of the snapshots nearest {0, T/2, T}. Returns (idxs, times)
    or None if the case has too few writes."""
    times = [float(t.time) for t in case]
    if len(times) < 2:
        return None
    idxs = [min(range(len(times)), key=lambda k: abs(times[k] - tgt))
            for tgt in (0.0, 0.5*float(T), float(T))]
    return idxs, times


def _read2d(case, idx, name, N, ncells=None):
    """Read an internal field and reshape to (Ny, Nx) for a structured hex mesh
    with square cells and x in [0,1] (same convention as interface_montage).
    A totally-lost phase writes alpha as a *uniform* field (size 1); broadcast it
    to the mesh size (from ``ncells``) so the panel still renders (e.g. empty)."""
    arr = np.asarray(case[idx][name].internal_field, float)
    if arr.size == 1 and ncells:
        arr = np.full(ncells, float(arr.reshape(-1)[0]))
    Nx = N
    Ny = arr.size // Nx
    return arr.reshape(Ny, Nx), Ny


def _grid_by_model_T(cases):
    g = {}
    for c in cases:
        g.setdefault((c["model"], c["T"]), {})[c["N"]] = c["dir"]
    return g


def _field_montage(model, T, ngrid, figdir, kind, prefix="", steady=False):
    """One montage per (model, T): rows = resolution N, cols = {t=0, T/2, T}.
    kind='alpha' -> filled volume fraction (+ 0.5 contour, dashed initial on the
    reversed column); kind='graderr' -> |grad psi|-1 field (+ psi=0 interface)."""
    Ns = sorted(ngrid)
    data = {}
    emax = []
    for r, N in enumerate(Ns):
        try:
            case = FoamCase(ngrid[N])
            snap = _snapshot_indices(case, T)
            if snap is None:
                continue
            idxs, times = snap
            xc = (np.arange(N) + 0.5) / N
            # psi is a signed distance -> always non-uniform, so its length is a
            # reliable cell count to broadcast any uniform (fully-lost) alpha to.
            ncells = np.asarray(case[idxs[0]]["psi"].internal_field, float).size
            a0 = _read2d(case, idxs[0], "alpha", N, ncells)[0] if kind == "alpha" else None
            for c, ti in enumerate(idxs):
                if kind == "alpha":
                    f2d, Ny = _read2d(case, ti, "alpha", N, ncells)
                    psi2d = None
                else:
                    psi2d, Ny = _read2d(case, ti, "psi", N, ncells)
                    h = 1.0 / N
                    gy, gx = np.gradient(psi2d, h, h)
                    f2d = np.abs(np.sqrt(gx*gx + gy*gy) - 1.0)
                    emax.append(float(np.nanpercentile(f2d, 99)))
                ymax = Ny / N
                yc = (np.arange(Ny) + 0.5) / Ny * ymax
                X, Y = np.meshgrid(xc, yc)
                data[(r, c)] = dict(X=X, Y=Y, f=f2d, psi=psi2d, a0=a0, ymax=ymax)
        except Exception as exc:  # noqa: BLE001
            print(f"[plots] {kind} {model} T={T:g} N={N}: {type(exc).__name__}: {exc}")
    if not data:
        return None
    vmax = float(np.clip(max(emax), 0.5, 3.0)) if (kind == "graderr" and emax) else 1.0
    cmap = "Blues" if kind == "alpha" else "inferno"
    fig, axes = plt.subplots(
        len(Ns), 3, figsize=(9.8, 2.7*len(Ns)), squeeze=False, layout="constrained")
    im = None
    for r, N in enumerate(Ns):
        for c in range(3):
            ax = axes[r][c]
            ax.set_xticks([]); ax.set_yticks([])
            d = data.get((r, c))
            if d is None:
                ax.axis("off"); continue
            im = ax.imshow(d["f"], origin="lower", extent=[0, 1, 0, d["ymax"]],
                           vmin=0.0, vmax=(1.0 if kind == "alpha" else vmax),
                           cmap=cmap, aspect="equal", interpolation="nearest")
            if kind == "alpha":
                ax.contour(d["X"], d["Y"], d["f"], [0.5], colors="#D85A30", linewidths=1.6)
                if c == 2 and d["a0"] is not None:
                    ax.contour(d["X"], d["Y"], d["a0"], [0.5],
                               colors="k", linestyles="--", linewidths=1.0)
            elif d["psi"] is not None:
                ax.contour(d["X"], d["Y"], d["psi"], [0.0], colors="#4fc3f7", linewidths=1.2)
            if r == 0:
                ax.set_title((_SNAP_LABELS_STEADY if steady else _SNAP_LABELS)[c],
                             fontsize=9)
            if c == 0:
                ax.set_ylabel(f"h = 1/{N}", fontsize=9)
    field_name = ("Volume fraction  " + r"$\alpha$") if kind == "alpha" \
        else r"Signed-distance error  $|\nabla\psi|-1$"
    flow = "steady vortex" if steady else "reversed vortex"
    fig.suptitle(f"{field_name} — {model}, {flow}, T = {T:g}", fontsize=12)
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
    stem = "alpha_evo" if kind == "alpha" else "graderr"
    p = os.path.join(figdir, f"{prefix}{stem}_{model}_T{T:g}.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def alpha_evolution(study_dir, figdir, prefix="", steady=False):
    """One alpha montage per (model, T) over the vortex."""
    grid = _grid_by_model_T(_cases_2d(study_dir))
    out = []
    for key in sorted(grid):
        p = _field_montage(key[0], key[1], grid[key], figdir, "alpha",
                           prefix, steady)
        if p:
            out.append(p)
    return out


def gradient_error_field(study_dir, figdir, prefix="", steady=False):
    """One |grad psi|-1 montage per (model, T): where the SDF property breaks."""
    grid = _grid_by_model_T(_cases_2d(study_dir))
    out = []
    for key in sorted(grid):
        p = _field_montage(key[0], key[1], grid[key], figdir, "graderr",
                           prefix, steady)
        if p:
            out.append(p)
    return out


def volume_history_T8(study_dir, figdir, model="none", T=8.0):
    """Relative phase-volume error vs t/T for the three resolutions (T=8, model
    none): the loss completes at t/T=0.5 (max deformation) and never recovers."""
    cases = [c for c in _cases_2d(study_dir)
             if c["model"] == model and abs(c["T"] - T) < 1e-9]
    if not cases:
        return None
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    plotted = False
    for c in sorted(cases, key=lambda c: c["N"]):
        csvp = os.path.join(c["dir"], "leiaLevelSetFoam.csv")
        if not os.path.isfile(csvp):
            continue
        ts, ev = [], []
        with open(csvp, newline="") as fh:
            for row in csv.DictReader(fh):
                t, e = _f(row.get("TIME")), _f(row.get("E_VOL_ALPHA_REL"))
                if t is not None and e is not None:
                    ts.append(t/T); ev.append(e)
        if ts:
            ax.plot(ts, ev, lw=1.8, label=f"h = 1/{c['N']}")
            plotted = True
    if not plotted:
        plt.close(fig); return None
    ax.axvline(0.5, color="k", ls=":", lw=1, alpha=0.7)
    ax.text(0.49, 0.03, "max deformation", rotation=90, va="bottom", ha="right",
            fontsize=8, alpha=0.7)
    ax.set_xlabel("t / T"); ax.set_ylabel(r"relative volume error  $E_{\mathrm{vol}}$")
    ax.set_title(f"Phase-volume loss vs time — {model}, T = {T:g}")
    ax.set_ylim(-0.03, 1.05); ax.grid(alpha=0.3); ax.legend()
    p = os.path.join(figdir, "volume_history_T8.png")
    fig.tight_layout(); fig.savefig(p, dpi=130); plt.close(fig)
    return p


# --- static (t=0) extension verification: e = |nHat . grad(Uext)| ------------
def _static_cases(study_dir):
    """Discover static-study cases: {dir, model, N, anchor}."""
    out = []
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*", "case_params.json"))):
        try:
            meta = json.load(open(meta_path))
        except Exception:
            continue
        tok = meta.get("tokens", {})
        N, A = _f(tok.get("N_CELLS")), _f(tok.get("ANCHOR_LAYERS"))
        if N is None or A is None:
            continue
        out.append({
            "dir": os.path.dirname(meta_path),
            "model": tok.get("VELOCITY_EXTENSION", "none"),
            "N": int(N), "anchor": int(A),
            "div": tok.get("UEXT_DIV", ""),
        })
    return out


def static_convergence(rows, figdir):
    """Convergence of the static normal-derivative norm: L2 (and Linf) of
    e = |nHat . grad(Uext)| over |psi| <= nLayers*h vs h. One series per
    (model [, uextDiv if swept] [, anchor if swept]); dashed raw-velocity
    reference; O(h)/O(h^2) guides."""
    sub = [r for r in rows if _f(r.get("eNormalL2")) is not None]
    if not sub:
        return None
    # steadyUpwind-family models are the only ones sensitive to UEXT_DIV;
    # collapse the sweep for the others to avoid duplicate series.
    DIV_SENSITIVE = ("steadyUpwind", "steadyUpwindLinear", "closestPoint")
    anchors = sorted({r.get("anchorLayers", "") for r in sub})
    divs = sorted({r.get("uextDiv", "") for r in sub})

    def series_key(r):
        m = r["velocityExtension"]
        d = r.get("uextDiv", "") if (len(divs) > 1 and m in DIV_SENSITIVE) else ""
        a = r.get("anchorLayers", "") if len(anchors) > 1 else ""
        return (m, d, a)

    keys = sorted({series_key(r) for r in sub if r["velocityExtension"] != "none"})
    cmap = plt.get_cmap("tab10")
    style = {}
    for i, k in enumerate(keys):
        style[k] = dict(color=cmap(i % 10),
                        linestyle=["-", "--", "-.", ":"][i % 4],
                        marker="osD^vP*X"[i % 8])

    def label_of(k):
        m, d, a = k
        lab = m
        if d:
            lab += ", " + ("upwind" if "Upwind" == d[:6] and "Linear" not in d
                           else "linUpwind" if "Linear" in d else d)
        if a != "":
            lab += f", anchor {a}"
        return lab

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    for ax, col, title in zip(
            axes, ("eNormalL2", "eNormalLinf"),
            (r"$\| \hat{n}\cdot\nabla U_{ext} \|_{L_2}$",
             r"$\| \hat{n}\cdot\nabla U_{ext} \|_{L_\infty}$")):
        allh, alle = [], []
        raw = sorted({(_f(r["h"]),
                       _f(r["eNormalRawL2"] if col == "eNormalL2"
                          else r.get("eNormalLinf")))
                      for r in sub if r["velocityExtension"] == "none"})
        raw = [(h, e) for h, e in raw if h and e]
        if raw:
            h, e = zip(*raw)
            ax.loglog(h, e, color="grey", ls="--", marker="s", ms=4,
                      label="raw U (= none)")
            allh += list(h); alle += list(e)
        for k in keys:
            pts = sorted({(_f(r["h"]), _f(r[col])) for r in sub
                          if series_key(r) == k
                          and r["velocityExtension"] != "none"})
            pts = [(h, e) for h, e in pts if h and e]
            if not pts:
                continue
            h, e = zip(*pts)
            ax.loglog(h, e, ms=4, label=label_of(k), **style[k])
            allh += list(h); alle += list(e)
        if allh:
            h0 = max(allh)
            e0 = max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.8, alpha=0.5, label=r"$O(h)$")
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.8, alpha=0.5, label=r"$O(h^2)$")
        ax.set_xlabel(r"$h$"); ax.set_title(title)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle(r"Static (t=0) extension verification: "
                 r"$e=|\hat{n}\cdot\nabla U_{ext}|$ over $|\psi|\leq n_{layers} h$",
                 fontsize=13)
    p = os.path.join(figdir, "static_convergence.png")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def static_e_montage(study_dir, figdir):
    """One montage per ANCHOR_LAYERS value: rows = resolution, cols = (raw U,
    anisotropicDiffusion, pseudoTime) showing e at t=0 with the psi=0 contour.
    Shared per-row vmax from the raw field so improvements read as darker."""
    cases = _static_cases(study_dir)
    if not cases:
        return []
    # If the UEXT_DIV scheme is swept, montage the (second-order) linearUpwind
    # variant; the convergence plot carries both.
    divs = sorted({c["div"] for c in cases})
    if len(divs) > 1:
        pref = [d for d in divs if "Linear" in d] or divs[:1]
        cases = [c for c in cases if c["div"] == pref[0]]
    anchors = sorted({c["anchor"] for c in cases})
    Ns = sorted({c["N"] for c in cases})
    order = ["none", "anisotropicDiffusion", "pseudoTime", "steadyUpwind",
             "steadyUpwindLinear", "closestPoint", "meshWave"]
    present = {c["model"] for c in cases}
    models = [m for m in order if m in present]
    titles = ["raw U" if m == "none" else m for m in models]
    ncols = len(models)
    out = []
    for A in anchors:
        grid = {(c["model"], c["N"]): c["dir"] for c in cases if c["anchor"] == A}
        fig, axes = plt.subplots(len(Ns), ncols,
                                 figsize=(3.3*ncols, 2.9*len(Ns)),
                                 squeeze=False, layout="constrained")
        im = None
        for r, N in enumerate(Ns):
            vmax = None
            for ci, m in enumerate(models):
                ax = axes[r][ci]
                ax.set_xticks([]); ax.set_yticks([])
                cdir = grid.get((m, N))
                if not cdir:
                    ax.axis("off"); continue
                try:
                    case = FoamCase(cdir)
                    name = "eNormalU" if m == "none" else "eNormalUext"
                    f2d, Ny = _read2d(case, 0, name, N)
                    psi2d, _ = _read2d(case, 0, "psi", N)
                    if vmax is None:
                        vmax = float(np.nanpercentile(f2d, 99))  # raw column first
                    ymax = Ny/N
                    im = ax.imshow(f2d, origin="lower", extent=[0, 1, 0, ymax],
                                   vmin=0, vmax=vmax, cmap="inferno",
                                   aspect="equal", interpolation="nearest")
                    xc = (np.arange(N)+0.5)/N
                    yc = (np.arange(Ny)+0.5)/Ny*ymax
                    X, Y = np.meshgrid(xc, yc)
                    ax.contour(X, Y, psi2d, [0.0], colors="#4fc3f7", linewidths=1.1)
                    if r == 0:
                        ax.set_title(titles[ci], fontsize=10)
                    if ci == 0:
                        ax.set_ylabel(f"h = 1/{N}", fontsize=9)
                except Exception as exc:  # noqa: BLE001
                    print(f"[plots] static {m} N={N} A={A}: {type(exc).__name__}: {exc}")
                    ax.axis("off")
        fig.suptitle(r"$e=|\hat{n}\cdot\nabla U|$ at $t=0$"
                     + f" — nAnchorLayers = {A}", fontsize=12)
        if im is not None:
            fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
        p = os.path.join(figdir, f"static_e_anchor{A}.png")
        fig.savefig(p, dpi=130); plt.close(fig)
        out.append(p)
    return out


# --- non-reversing analysis ---------------------------------------------------
def _study_oscillation(study_dir):
    """'on' (reversed benchmark) or 'off' (steady non-reversing vortex), read
    from the first case's tokens. Defaults to 'on' (all legacy studies)."""
    for meta_path in sorted(glob.glob(os.path.join(study_dir, "*",
                                                   "case_params.json")))[:1]:
        try:
            tok = json.load(open(meta_path)).get("tokens", {})
            return str(tok.get("OSCILLATION", "on")).lower()
        except Exception:
            pass
    return "on"


def maxdef_convergence(rows, figdir, prefix=""):
    """Band gradient error vs h per model at t = T/2 (top row: maximal
    deformation, BEFORE any reversal cancellation) and at t = T (bottom row:
    final). Columns: T. The drop from the top row to the bottom row for model
    'none' is the reversibility bias of the classical benchmark — errors that
    cancelled on the way back."""
    sub = [r for r in rows if _f(r.get("gradientErrorBandHalf")) is not None]
    if not sub:
        return None
    Ts = sorted({_f(r["T"]) for r in sub if _f(r["T"]) is not None})
    models = sorted({r.get("velocityExtension", "none") for r in sub})
    cmap = plt.get_cmap("tab10")
    mcol = {m: cmap(i % 10) for i, m in enumerate(models)}
    fig, axes = plt.subplots(2, len(Ts), figsize=(3.4*len(Ts)+1.6, 7.4),
                             squeeze=False)
    for ci, T in enumerate(Ts):
        for ri, col in enumerate(("gradientErrorBandHalf", "gradientErrorBand")):
            ax = axes[ri][ci]
            allh, alle = [], []
            for m in models:
                pts = sorted((_f(r["h"]), _f(r[col])) for r in sub
                             if r.get("velocityExtension") == m
                             and _f(r["T"]) == T
                             and _f(r.get(col)) not in (None, 0.0))
                pts = [(h, e) for h, e in pts if h and e]
                if not pts:
                    continue
                h, e = zip(*pts)
                ax.loglog(h, e, marker="o", ms=4, lw=1.4, color=mcol[m],
                          label=m if (ri == 0 and ci == 0) else None)
                allh += list(h); alle += list(e)
            if allh:
                h0, e0 = max(allh), max(alle)
                hr = np.array([min(allh), max(allh)])
                ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.7, alpha=0.5)
                ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.7, alpha=0.5)
            if ri == 0:
                ax.set_title(f"T = {T:g}", fontsize=10)
            if ci == 0:
                ax.set_ylabel(("t = T/2 (max deformation)" if ri == 0
                               else "t = T (final)") + "\n"
                              + r"$\||\nabla\psi|-1\|_{L_2}$ (band)", fontsize=9)
            ax.set_xlabel(r"$h$")
            ax.grid(True, which="both", alpha=0.3)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center",
                   ncol=min(6, len(models)), fontsize=8, frameon=False)
    fig.suptitle("Band gradient error at maximal deformation (t = T/2, "
                 "no cancellation) vs final time", fontsize=12)
    p = os.path.join(figdir, f"{prefix}maxdef_convergence.png")
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def timeseries_dashboard(study_dir, figdir, prefix=""):
    """Steady (non-reversing) study headline: per-model time series at each
    resolution. Rows: band gradient error (log), band min |grad psi|
    (flattening / conditioning err_x ~ err_psi/|grad psi|), relative volume
    error (log). The crossover time t* where 'none' stops winning is read
    directly off the top row."""
    cases = _cases_2d(study_dir)
    if not cases:
        return None
    Ns = sorted({c["N"] for c in cases})
    models = sorted({c["model"] for c in cases})
    cmap = plt.get_cmap("tab10")
    mcol = {m: cmap(i % 10) for i, m in enumerate(models)}
    specs = [
        ("gradPsiError.csv", "E_NARROW_L2_GRAD_PSI",
         r"$\||\nabla\psi|-1\|_{L_2}$ (band)", "log"),
        ("gradPsiError.csv", "NARROW_MIN_MAG_GRAD_PSI",
         r"$\min_{\mathrm{band}}|\nabla\psi|$", "linear"),
        ("leiaLevelSetFoam.csv", "E_VOL_ALPHA_REL",
         r"$E_{\mathrm{vol}}$ (relative)", "log"),
    ]
    fig, axes = plt.subplots(len(specs), len(Ns),
                             figsize=(5.2*len(Ns), 3.0*len(specs)),
                             squeeze=False)
    plotted = False
    for ci, N in enumerate(Ns):
        for c in sorted((c for c in cases if c["N"] == N),
                        key=lambda c: c["model"]):
            for ri, (csvname, col, label, yscale) in enumerate(specs):
                csvp = os.path.join(c["dir"], csvname)
                if not os.path.isfile(csvp):
                    continue
                ts, ys = [], []
                with open(csvp, newline="") as fh:
                    for row in csv.DictReader(fh):
                        t, y = _f(row.get("TIME")), _f(row.get(col))
                        # |y| >= 1e30: pre-guard CSVs wrote +-VGREAT for an
                        # empty narrow band (annihilated phase) -- drop those.
                        if t is not None and y is not None and abs(y) < 1e30:
                            ts.append(t); ys.append(y)
                if not ts:
                    continue
                ax = axes[ri][ci]
                ax.plot(ts, ys, lw=1.5, color=mcol[c["model"]],
                        label=c["model"] if (ri == 0 and ci == 0) else None)
                plotted = True
    if not plotted:
        plt.close(fig); return None
    for ci, N in enumerate(Ns):
        for ri, (_, _, label, yscale) in enumerate(specs):
            ax = axes[ri][ci]
            ax.set_yscale(yscale)
            ax.grid(alpha=0.3)
            if ri == 0:
                ax.set_title(f"h = 1/{N}", fontsize=10)
            if ri == len(specs) - 1:
                ax.set_xlabel("t")
            if ci == 0:
                ax.set_ylabel(label, fontsize=9)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center",
                   ncol=min(6, len(models)), fontsize=8, frameon=False)
    fig.suptitle("Steady (non-reversing) vortex: error evolution per model",
                 fontsize=12)
    p = os.path.join(figdir, f"{prefix}timeseries.png")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


# --- semi-Lagrangian: reconstruction convergence + cross-study head-to-head --
def reconstruction_convergence(rows, figdir, Tpick=2.0):
    """SL order figure: band gradient error and shape error vs h, one line per
    reconstruction, at a representative period T, in two columns (CFL 0.5, 1.0).
    Shows the order jump linearTaylor O(h) -> nested/quadratic O(h^2) and its
    robustness to the Courant number."""
    sub = [r for r in rows if r.get("reconstruction", "")
           and _f(r.get("gradientErrorBand")) is not None]
    if not sub:
        return None
    recons = sorted({r["reconstruction"] for r in sub})
    cfls = sorted({_f(r["cfl"]) for r in sub if _f(r.get("cfl")) is not None})
    if not cfls:
        cfls = [None]
    cmap = plt.get_cmap("tab10")
    rcol = {m: cmap(i % 10) for i, m in enumerate(recons)}
    series = [("gradientErrorBand", r"$\||\nabla\psi|-1\|_{L_2}$ (band)"),
              ("shapeError", "geometric shape error")]

    fig, axes = plt.subplots(len(series), len(cfls),
                             figsize=(5.0*len(cfls)+1, 4.0*len(series)),
                             squeeze=False)
    for cj, cfl in enumerate(cfls):
        for ri, (col, ylabel) in enumerate(series):
            ax = axes[ri][cj]
            allh, alle = [], []
            for m in recons:
                pts = sorted((_f(r["h"]), _f(r[col])) for r in sub
                             if r["reconstruction"] == m
                             and (cfl is None or _f(r.get("cfl")) == cfl)
                             and abs((_f(r["T"]) or -1) - Tpick) < 1e-9
                             and _f(r.get(col)) not in (None, 0.0))
                pts = [(h, e) for h, e in pts if h and e]
                if not pts:
                    continue
                h, e = zip(*pts)
                ax.loglog(h, e, marker="o", ms=5, lw=1.6, color=rcol[m],
                          label=(m if (ri == 0 and cj == 0) else None))
                allh += list(h); alle += list(e)
            if allh:
                h0, e0 = max(allh), max(alle)
                hr = np.array([min(allh), max(allh)])
                ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.8, alpha=0.5, label=r"$O(h)$"
                          if (ri == 0 and cj == 0) else None)
                ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.8, alpha=0.5, label=r"$O(h^2)$"
                          if (ri == 0 and cj == 0) else None)
            if ri == 0:
                ax.set_title(f"CFL = {cfl:g}" if cfl is not None else "", fontsize=11)
            ax.set_xlabel(r"$h$")
            if cj == 0:
                ax.set_ylabel(ylabel, fontsize=10)
            ax.grid(True, which="both", alpha=0.3)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center",
                   ncol=min(6, len(labels)), fontsize=8, frameon=False)
    fig.suptitle(f"Semi-Lagrangian reconstruction convergence "
                 f"(reversed vortex, T = {Tpick:g})", fontsize=13)
    p = os.path.join(figdir, "sl_reconstruction_convergence.png")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def sl_vs_extension(study_dir, figdir, Tpick=2.0):
    """Cross-study head-to-head at fixed T: the SL quadraticWLSQ solver against
    the velocity-extension solver (none + meshWave), band gradient error and
    shape error vs h on one axis. Skips gracefully if a study dir is absent."""
    studies_root = os.path.dirname(os.path.abspath(study_dir))
    sources = [
        ("bulkVortexSLHighRes", "reconstruction", "quadraticWLSQ",
         "SL quadraticWLSQ"),
        ("bulkVortexSL", "reconstruction", "quadraticWLSQ", "SL quadraticWLSQ"),
        ("bulkVortexHighRes", "velocityExtension", "none", "Eulerian none"),
        ("bulkVortexHighRes", "velocityExtension", "meshWave",
         "Eulerian meshWave"),
        ("bulkVortexConvergence", "velocityExtension", "none", "Eulerian none"),
    ]
    # Collect one series per (label): first matching study wins.
    picked = {}
    for study, key, val, label in sources:
        if label in picked:
            continue
        hits = glob.glob(os.path.join(studies_root, study, "*_errors.csv"))
        if not hits:
            continue
        rows = _rows(hits[0])
        picked[label] = [r for r in rows if r.get(key, "") == val]
    if len(picked) < 2:
        return None

    series = [("gradientErrorBand", r"$\||\nabla\psi|-1\|_{L_2}$ (band)"),
              ("shapeError", "geometric shape error")]
    cmap = plt.get_cmap("tab10")
    lcol = {lab: cmap(i % 10) for i, lab in enumerate(picked)}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), squeeze=False)
    for ci, (col, ylabel) in enumerate(series):
        ax = axes[0][ci]
        allh, alle = [], []
        for label, rws in picked.items():
            pts = sorted((_f(r["h"]), _f(r[col])) for r in rws
                         if abs((_f(r.get("T")) or -1) - Tpick) < 1e-9
                         and _f(r.get(col)) not in (None, 0.0)
                         and (_f(r.get("cfl")) in (None, 0.5)
                              or r.get("cfl", "") == ""))
            pts = [(h, e) for h, e in pts if h and e]
            if not pts:
                continue
            h, e = zip(*pts)
            ax.loglog(h, e, marker="o", ms=5, lw=1.7, color=lcol[label],
                      label=(label if ci == 0 else None))
            allh += list(h); alle += list(e)
        if allh:
            h0, e0 = max(allh), max(alle)
            hr = np.array([min(allh), max(allh)])
            ax.loglog(hr, e0*(hr/h0)**1, "k--", lw=0.8, alpha=0.5)
            ax.loglog(hr, e0*(hr/h0)**2, "k:", lw=0.8, alpha=0.5)
        ax.set_xlabel(r"$h$"); ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(True, which="both", alpha=0.3)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center",
                   ncol=min(4, len(labels)), fontsize=8, frameon=False)
    fig.suptitle(f"Semi-Lagrangian vs velocity extension "
                 f"(reversed vortex, T = {Tpick:g}, CFL 0.5)", fontsize=13)
    p = os.path.join(figdir, "sl_vs_extension.png")
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def sl_alpha_montage(study_dir, figdir):
    """alpha montage per reconstruction (CFL 0.5) over the reversed vortex."""
    cases = [c for c in _cases_2d(study_dir)
             if c.get("recon") and (c.get("cfl") in (None, 0.5))]
    out = []
    grid = _grid_by_model_T(cases)
    for key in sorted(grid):
        p = _field_montage(key[0], key[1], grid[key], figdir, "alpha",
                           prefix="sl_", steady=False)
        if p:
            out.append(p)
    return out


def main(study_dir, slides_dir):
    hits = glob.glob(os.path.join(study_dir, "*_errors.csv"))
    if not hits:
        print(f"[plots] no *_errors.csv in {study_dir}; nothing to plot")
        return []
    rows = _rows(hits[0])
    figdir = os.path.join(study_dir, "figures")
    os.makedirs(figdir, exist_ok=True)
    if slides_dir:
        os.makedirs(slides_dir, exist_ok=True)

    # Static (t=0) extension-verification study -> its own figure family.
    if any(_f(r.get("eNormalL2")) is not None for r in rows):
        outputs = [p for p in [static_convergence(rows, figdir)] if p]
        outputs += static_e_montage(study_dir, figdir)
        if slides_dir:
            for p in outputs:
                shutil.copy(p, os.path.join(slides_dir, os.path.basename(p)))
        print(f"[plots] wrote {len(outputs)} static-verification figures to {figdir}"
              + (f" and {slides_dir}" if slides_dir else ""))
        return outputs

    # Semi-Lagrangian study (reconstruction axis) -> its own figure family.
    if any(r.get("reconstruction", "") for r in rows):
        outputs = []
        rc = reconstruction_convergence(rows, figdir)
        if rc:
            outputs.append(rc)
        ov = sl_vs_extension(study_dir, figdir)
        if ov:
            outputs.append(ov)
        try:
            outputs += sl_alpha_montage(study_dir, figdir)
        except Exception as exc:  # noqa: BLE001
            print(f"[plots] SL alpha montage skipped: {type(exc).__name__}: {exc}")
        if slides_dir:
            for p in outputs:
                shutil.copy(p, os.path.join(slides_dir, os.path.basename(p)))
        print(f"[plots] wrote {len(outputs)} semi-Lagrangian figures to {figdir}"
              + (f" and {slides_dir}" if slides_dir else ""))
        return outputs

    # Steady (oscillation off) studies get a 'steady_' figure prefix so their
    # exports never clobber the reversed-benchmark figure family in doc/slides.
    steady = _study_oscillation(study_dir) == "off"
    prefix = "steady_" if steady else ""

    inds = sorted({r.get("phaseIndicator", "") for r in rows if r.get("phaseIndicator")})
    if len(inds) > 1:
        # Study varies the phase indicator -> overlay comparison.
        outputs = [p for p in [indicator_comparison(rows, figdir)] if p]
    else:
        # Study varies the velocity-extension model -> one triptych per model.
        models = sorted({r.get("velocityExtension", "none") for r in rows})
        outputs = [p for p in (convergence_triptych(rows, m, figdir, prefix)
                               for m in models) if p]

    montage = interface_montage(study_dir, figdir, "none", prefix, steady)
    if montage:
        outputs.append(montage)

    # Maximal-deformation (t = T/2) convergence: the reversal-free reading of
    # the same data (needs the half.* columns from aggregate.py).
    md = maxdef_convergence(rows, figdir, prefix)
    if md:
        outputs.append(md)

    # Automated field atlas: alpha + |grad psi|-1 at t = 0, T/2, T, per (model, T).
    try:
        outputs += alpha_evolution(study_dir, figdir, prefix, steady)
        outputs += gradient_error_field(study_dir, figdir, prefix, steady)
        vh = volume_history_T8(study_dir, figdir)
        if vh:
            outputs.append(vh)
    except Exception as exc:  # noqa: BLE001
        print(f"[plots] field atlas skipped: {type(exc).__name__}: {exc}")

    if steady:
        td = timeseries_dashboard(study_dir, figdir, prefix)
        if td:
            outputs.append(td)
        # Shape error vs the marker-traced quasi-exact reference (there is no
        # analytic final interface without reversal).
        try:
            import marker_ref
            outputs += marker_ref.shape_report(study_dir, figdir, prefix)
        except Exception as exc:  # noqa: BLE001
            print(f"[plots] marker shape metric skipped: {type(exc).__name__}: {exc}")

    if slides_dir:
        for p in outputs:
            shutil.copy(p, os.path.join(slides_dir, os.path.basename(p)))
    dest = slides_dir if slides_dir else "(study only)"
    print(f"[plots] wrote {len(outputs)} figures to {figdir} and {dest}")
    return outputs


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else sys.argv[1])
