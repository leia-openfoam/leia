#!/usr/bin/env python3
"""Paper-quality psi=0 iso-surface montage for ONE 3D semi-Lagrangian study.

A study is one case x one mesh (e.g. uncachedConv3Dshear = 3D shear, hex). This
renders the zero level set of psi -- the fluid interface -- at t=0, T/2 and T on
the study's finest THREE resolutions, in a fixed canonical orientation, and tiles
them into a single montage:

    columns = resolution (coarse -> fine),  rows = t=0 / T/2 / T

written to the semi-Lagrangian theme's data/figures as

    isosurface_<case>_<mesh>.png     (e.g. isosurface_3Dshear_hex.png)

so both the reveal deck and the article can show the interface evolution + its
convergence with resolution. foamlib is used to locate the case; VTK's
vtkOpenFOAMReader (robust on polyhedra) reads psi and vtkContourFilter extracts the
iso-surface, rendered offscreen with a shared camera so all panels are comparable.

Usage:  python3 make_isosurface_fig.py <study_dir>
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import vtk

import paths

PANEL_PX = 900          # per-panel render resolution
TIMES_FRAC = (0.0, 0.5, 1.0)   # t = 0, T/2, T (fractions of the end time)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def study_cases(study_dir):
    """Return [(h, label, case_dir, case, mesh)] for a study, finest last."""
    rows = []
    for cp in sorted(glob.glob(os.path.join(study_dir, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp))
        tok = d.get("tokens", {})
        case, mesh = d.get("case"), d.get("mesh")
        if mesh == "poly":
            h = _f(tok.get("MAX_CELL_SIZE")); label = f"$h={h:g}$" if h else "?"
        else:
            n = _f(tok.get("N_CELLS")); h = (1.0 / n) if n else None
            label = f"$N={int(n)}$" if n else "?"
        if h:
            rows.append((h, label, os.path.dirname(cp), case, mesh))
    rows.sort(key=lambda r: r[0], reverse=True)     # coarse (large h) -> fine
    return rows


def render_iso(case_dir, target_times, camera, out_prefix):
    """Render the psi=0 iso-surface at each target time -> PNG files.
    Returns [(t, png|None)]. Reuses one offscreen window + a shared camera."""
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()                         # marker for the reader

    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(foam)
    reader.CreateCellToPointOn()                    # contour needs point data
    reader.SetSkipZeroTime(0)                       # keep t=0
    reader.Update()
    tv = reader.GetTimeValues()
    tvals = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())] if tv else []
    if not tvals:
        return [(t, None) for t in target_times]

    contour = vtk.vtkContourFilter()
    contour.SetValue(0, 0.0)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(contour.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.30, 0.55, 0.90)  # steel blue interface
    actor.GetProperty().SetInterpolationToPhong()

    ren = vtk.vtkRenderer()
    ren.SetBackground(1, 1, 1)
    ren.AddActor(actor)
    ren.AddLight(_light(0.3))
    win = vtk.vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.AddRenderer(ren)
    win.SetSize(PANEL_PX, PANEL_PX)

    results = []
    for t in target_times:
        tv = min(tvals, key=lambda x: abs(x - t))
        reader.UpdateTimeStep(tv)
        internal = _internal_mesh(reader.GetOutput())
        if internal is None:
            results.append((t, None)); continue
        internal.GetPointData().SetActiveScalars("psi")
        contour.SetInputData(internal)
        contour.Update()
        # add a domain outline once (first time) for spatial reference
        if not ren.GetActors2D().GetNumberOfItems():
            outline = vtk.vtkOutlineFilter(); outline.SetInputData(internal)
            om = vtk.vtkPolyDataMapper(); om.SetInputConnection(outline.GetOutputPort())
            oa = vtk.vtkActor(); oa.SetMapper(om)
            oa.GetProperty().SetColor(0.6, 0.6, 0.6); oa.GetProperty().SetLineWidth(1)
            ren.AddActor(oa)
        ren.ResetCamera()
        _apply_camera(ren.GetActiveCamera(), camera)
        win.Render()
        w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(win); w2i.Update()
        png = f"{out_prefix}_t{t:g}.png"
        wr = vtk.vtkPNGWriter(); wr.SetFileName(png)
        wr.SetInputConnection(w2i.GetOutputPort()); wr.Write()
        results.append((t, png if os.path.isfile(png) else None))
    win.Finalize()
    return results


def _internal_mesh(mb):
    """The internalMesh unstructured grid carrying the psi point array."""
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        blk = it.GetCurrentDataObject()
        if blk is not None and blk.IsA("vtkUnstructuredGrid") \
           and blk.GetPointData().GetArray("psi") is not None:
            return blk
        it.GoToNextItem()
    return None


def _light(intensity):
    lt = vtk.vtkLight(); lt.SetLightTypeToCameraLight()
    lt.SetPosition(1, 1, 1); lt.SetIntensity(intensity)
    return lt


def _apply_camera(cam, camera):
    az, el, zoom = camera
    cam.Azimuth(az); cam.Elevation(el)
    cam.OrthogonalizeViewUp()
    cam.Zoom(zoom)


def make(study_dir, camera=(30.0, 20.0, 1.25)):
    cases = study_cases(study_dir)
    if not cases:
        print(f"[iso] no cases in {study_dir}; skip"); return None
    case, mesh = cases[0][3], cases[0][4]
    finest3 = cases[-3:] if len(cases) >= 3 else cases   # coarse->fine of the finest 3
    end = None
    # end time from any case's controlDict token (END_TIME) via case_params
    d0 = json.load(open(os.path.join(finest3[0][2], "case_params.json")))
    end = _f(d0.get("tokens", {}).get("END_TIME")) or 3.0
    targets = [frac * end for frac in TIMES_FRAC]

    ncol, nrow = len(finest3), len(targets)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * nrow), squeeze=False)
    import tempfile
    tmp = tempfile.mkdtemp(prefix="iso_")
    for j, (h, label, cdir, _c, _m) in enumerate(finest3):
        panels = dict(render_iso(cdir, targets, camera, os.path.join(tmp, f"c{j}")))
        for i, t in enumerate(targets):
            ax = axes[i][j]; ax.set_xticks([]); ax.set_yticks([])
            png = panels.get(t)
            if png:
                ax.imshow(mpimg.imread(png))
            else:
                ax.text(0.5, 0.5, "n/a", ha="center", va="center")
            for s in ax.spines.values():
                s.set_visible(False)
            if i == 0:
                ax.set_title(label, fontsize=13)
            if j == 0:
                lab = {0: r"$t=0$", 1: r"$t=T/2$", 2: r"$t=T$"}.get(i, f"t={t:g}")
                ax.set_ylabel(lab, fontsize=13, rotation=90, labelpad=8)
    fig.suptitle(rf"$\psi=0$ interface — {case} ({mesh})", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(paths.figs_dir("semi-lagrangian-level-set"),
                       f"isosurface_{case}_{mesh}.png")
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"[iso] wrote {out}  ({ncol} resolutions x {nrow} times, camera={camera})")
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: make_isosurface_fig.py <study_dir>"); sys.exit(2)
    sys.exit(0 if make(sys.argv[1]) else 1)
