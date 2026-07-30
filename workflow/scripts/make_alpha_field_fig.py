#!/usr/bin/env python3
"""Volume-fraction FIELD montage for ONE 2D kinematic semi-Lagrangian study.

Renders the volume-fraction field alpha (the phase indicator the level set
carries) of the 2D reversed-vortex benchmark at t = 0, T/2 (maximum
deformation) and T (return to the initial circle), across the resolution
ladder (coarsest / median / finest by default):

    columns = resolution (coarse -> fine),  rows = t=0 / T/2 / T

Each panel shows the alpha field itself (graded fill, so the width of the
indicator transition region -- the numerical smearing of the method -- is
visible), the psi=0 level-set contour (the true interface location), and, on
the t=T row, the exact initial circle (dashed) so the reversal error can be
judged by eye. Written to the study's theme data/figures as

    alpha_field_<case>_<mesh>.png    (e.g. alpha_field_2Dvortex_hex.png)

Fields are read with VTK's vtkOpenFOAMReader (cell-centred data + cell centres
pulled to numpy, same robust reader as make_droplet_fields_fig.py). The 2D
ladders sweep CFL; one CFL is selected for the montage (default 0.5, matching
the 3D studies). Cases stuck at t=0 (failed solves) are skipped.

Usage:
    python3 make_alpha_field_fig.py <study_dir> [--theme T] [--cfl 0.5] [--nres 3]
"""
import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import paths

# Geometry of cases/2Dvortex (classic single-vortex benchmark): unit box,
# initial circle centre (0.5, 0.75), radius 0.15.
CIRCLE_CENTRE = (0.5, 0.75)
CIRCLE_R = 0.15
TIMES_FRAC = (0.0, 0.5, 1.0)    # t = 0, T/2, T (fractions of the end time)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _complete(case_dir):
    """True if the case advanced past t=0 (>=2 numeric time directories)."""
    n = 0
    for e in os.listdir(case_dir):
        if os.path.isdir(os.path.join(case_dir, e)):
            try:
                float(e); n += 1
            except ValueError:
                pass
    return n >= 2


def study_cases(study_dir, cfl=None):
    """[(h, label, case_dir, case, mesh, tokens)] coarse -> fine, one CFL only."""
    rows = []
    for cp in sorted(glob.glob(os.path.join(study_dir, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp))
        tok = d.get("tokens", {})
        case, mesh = d.get("case"), d.get("mesh")
        if cfl is not None:
            c = _f(tok.get("CFL"))
            if c is not None and abs(c - cfl) > 1e-12:
                continue
        if mesh == "poly":
            h = _f(tok.get("MAX_CELL_SIZE")); label = f"$h={h:g}$" if h else "?"
        else:
            n = _f(tok.get("N_CELLS")); h = (1.0 / n) if n else None
            label = f"$N={int(n)}$" if n else "?"
        case_dir = os.path.dirname(cp)
        if not h:
            continue
        if not _complete(case_dir):
            print(f"[alpha] skip incomplete case (only t=0): {case_dir}")
            continue
        rows.append((h, label, case_dir, case, mesh, tok))
    rows.sort(key=lambda r: r[0], reverse=True)     # coarse (large h) -> fine
    return rows


def read_case(case_dir, target_times):
    """{t: {x, y, alpha, psi}} of cell-centred data at the nearest written times."""
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(foam)
    reader.CreateCellToPointOff()          # keep native cell-centred data
    reader.SetSkipZeroTime(0)
    reader.Update()
    reader.EnableAllCellArrays()
    reader.Update()
    tv = reader.GetTimeValues()
    tvals = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())] if tv else []
    if not tvals:
        return {}

    def _internal(mb):
        it = mb.NewIterator(); it.InitTraversal()
        while not it.IsDoneWithTraversal():
            blk = it.GetCurrentDataObject()
            if blk is not None and blk.IsA("vtkUnstructuredGrid") \
               and blk.GetCellData().GetArray("alpha") is not None:
                return blk
            it.GoToNextItem()
        return None

    out = {}
    for t in target_times:
        tsel = min(tvals, key=lambda x: abs(x - t))
        reader.UpdateTimeStep(tsel)
        internal = _internal(reader.GetOutput())
        if internal is None:
            continue
        cc = vtk.vtkCellCenters(); cc.SetInputData(internal); cc.Update()
        pts = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())
        cd = internal.GetCellData()

        def arr(name):
            a = cd.GetArray(name)
            return vtk_to_numpy(a) if a is not None else None

        out[t] = {"x": pts[:, 0], "y": pts[:, 1],
                  "alpha": arr("alpha"), "psi": arr("psi"), "time": tsel}
    return out


def make(study_dir, theme="semi-lagrangian-level-set", cfl=0.5, nres=3):
    cases = study_cases(study_dir, cfl=cfl)
    if not cases:
        # ladder without a CFL axis (or a different sweep): take everything
        cases = study_cases(study_dir, cfl=None)
    if not cases:
        print(f"[alpha] no cases in {study_dir}; skip")
        return None
    case, mesh = cases[0][3], cases[0][4]
    # coarsest / median / finest of the ladder (nres columns)
    if len(cases) > nres:
        idx = np.unique(np.round(np.linspace(0, len(cases) - 1, nres)).astype(int))
        cases = [cases[i] for i in idx]
    end = _f(cases[0][5].get("END_TIME")) or 2.0
    targets = [frac * end for frac in TIMES_FRAC]
    recon = cases[0][5].get("SL_RECONSTRUCTION", "")

    ncol, nrow = len(cases), len(targets)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.6 * nrow),
                             squeeze=False)
    # Levels STRADDLE 0 and 1: the bulk field is exactly 0/1 (flat), and
    # tricontourf drops flat patches lying exactly ON a contour level (the top
    # interval boundary), leaving the drop interior unfilled. Half-interval
    # offsets put the flat values strictly inside intervals.
    levels = np.linspace(-0.025, 1.025, 22)
    mappable = None
    for j, (h, label, cdir, _c, _m, _tok) in enumerate(cases):
        fields = read_case(cdir, targets)
        for i, t in enumerate(targets):
            ax = axes[i][j]
            ax.set_xticks([]); ax.set_yticks([])
            f = fields.get(t)
            if f is None or f["alpha"] is None:
                ax.text(0.5, 0.5, "n/a", ha="center", va="center")
                continue
            # graded alpha fill: the transition-band width shows the smearing
            mappable = ax.tricontourf(f["x"], f["y"], np.clip(f["alpha"], 0, 1),
                                      levels=levels, cmap="Blues")
            # psi=0 = the level-set interface (smoother than the alpha=1/2 step)
            if f.get("psi") is not None:
                ax.tricontour(f["x"], f["y"], f["psi"], levels=[0.0],
                              colors="k", linewidths=1.0)
            # exact initial circle on the reversal row: the shape error by eye
            if i == nrow - 1:
                th = np.linspace(0, 2 * np.pi, 361)
                ax.plot(CIRCLE_CENTRE[0] + CIRCLE_R * np.cos(th),
                        CIRCLE_CENTRE[1] + CIRCLE_R * np.sin(th),
                        "r--", lw=1.0, label="exact")
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
            if i == 0:
                ax.set_title(label, fontsize=13)
            if j == 0:
                lab = {0: r"$t=0$", 1: r"$t=T/2$", 2: r"$t=T$"}.get(i, f"t={t:g}")
                ax.set_ylabel(lab, fontsize=13, rotation=90, labelpad=8)
    if mappable is not None:
        cb = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
        cb.set_label(r"volume fraction $\alpha$")
        cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    sub = f" — {recon}" if recon else ""
    fig.suptitle(rf"Volume fraction $\alpha$ — {case} ({mesh}, CFL {cfl:g}){sub}",
                 fontsize=15)
    out = os.path.join(paths.figs_dir(theme), f"alpha_field_{case}_{mesh}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[alpha] wrote {out}  ({ncol} resolutions x {nrow} times)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="semi-lagrangian-level-set",
                    help="docs theme receiving the figure (see paths.py)")
    ap.add_argument("--cfl", type=float, default=0.5,
                    help="CFL slice of the ladder to render (default 0.5)")
    ap.add_argument("--nres", type=int, default=3,
                    help="number of resolutions (columns) in the montage")
    args = ap.parse_args()
    raise SystemExit(0 if make(args.study_dir, theme=args.theme,
                               cfl=args.cfl, nres=args.nres) else 1)
