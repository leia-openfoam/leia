#!/usr/bin/env python3
"""Paper-quality psi=0 iso-surface montage for ONE 3D advection study.

A study is one case x one mesh (e.g. uncachedConv3Dshear = 3D shear, hex), but
it may hold several METHOD ARMS -- an SDPLS study is six (noSource/R/beta x two
linearizations). One montage is rendered PER ARM: the zero level set of psi --
the fluid interface -- at t=0, T/2 and T on that arm's finest THREE resolutions,
in a fixed canonical orientation, tiled as

    columns = resolution (coarse -> fine),  rows = t=0 / T/2 / T

written to the study's theme data/figures (default: the quadratic semi-Lagrangian
theme; other lines pass their own theme) as

    isosurface_<case>_<mesh>.png            single-arm study (unchanged name,
                                            referenced by the SL/LSL articles)
    isosurface_<case>_<mesh>_<arm>.png      multi-arm study, e.g.
                                            ..._euler_SDPLS_R_simpleImp.png

so the article can show the interface evolution and its convergence with resolution.
foamlib is used to locate the case; VTK's
vtkOpenFOAMReader (robust on polyhedra) reads psi and vtkContourFilter extracts the
iso-surface, rendered offscreen with a shared camera so all panels are comparable.

Usage:  python3 make_isosurface_fig.py <study_dir> [--theme <docs theme>]
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

import method_label
import paths

PANEL_PX = 900          # per-panel render resolution
TIMES_FRAC = (0.0, 0.5, 1.0)   # t = 0, T/2, T (fractions of the end time)


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _complete(case_dir):
    """True if the case advanced past t=0, i.e. has >=2 reconstructed numeric time
    directories. Guards the montage against a case that failed mid-solve (only t=0),
    which would otherwise be tiled as a broken column."""
    n = 0
    for e in os.listdir(case_dir):
        if os.path.isdir(os.path.join(case_dir, e)):
            try:
                float(e); n += 1
            except ValueError:
                pass
    return n >= 2


def study_cases(study_dir):
    """Return [(h, label, case_dir, case, mesh, arm)] for a study, finest last.
    Cases that did not advance past t=0 (failed solves) are skipped.

    `arm` is the method slug. A semi-Lagrangian study is one arm and behaves
    exactly as before; an SDPLS study is six (noSource/R/beta x two
    linearizations), and without this the "finest three cases" below would have
    mixed methods into a single montage and labelled it with none of them.
    """
    rows = []
    for cp in sorted(glob.glob(os.path.join(study_dir, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp))
        tok = d.get("tokens", {})
        case, mesh = d.get("case"), d.get("mesh")
        arm = method_label.method_slug(tok)
        if mesh == "poly":
            h = _f(tok.get("MAX_CELL_SIZE")); label = f"$h={h:g}$" if h else "?"
        else:
            n = _f(tok.get("N_CELLS")); h = (1.0 / n) if n else None
            label = f"$N={int(n)}$" if n else "?"
        case_dir = os.path.dirname(cp)
        if not h:
            continue
        if not _complete(case_dir):
            print(f"[iso] skip incomplete case (only t=0): {case_dir}")
            continue
        rows.append((h, label, case_dir, case, mesh, arm))
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
    actor.GetProperty().SetColor(0.66, 0.81, 0.97)  # light blue interface
    actor.GetProperty().SetAmbient(0.55)            # flatter, lighter shading
    actor.GetProperty().SetDiffuse(0.55)
    actor.GetProperty().SetInterpolationToPhong()

    ren = vtk.vtkRenderer()
    ren.SetBackground(1, 1, 1)
    ren.AddActor(actor)
    ren.AddLight(_light(0.55))
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
        # No domain bounding box: frame the interface itself so features are visible.
        ren.ResetCamera()
        _apply_camera(ren.GetActiveCamera(), camera)
        ren.ResetCameraClippingRange()
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
    """Position the camera. With an explicit view-up (camera[3]) the camera is placed on
    an orbit around the focal point at (azimuth, elevation) with that up-vector, so the
    chosen axis stays vertical in the image (SetViewUp alone degenerates when it is
    parallel to the default view direction). Without a view-up, azimuth/elevation are
    applied as increments from the default view. The caller resets the clipping range."""
    import math
    az, el, zoom = camera[0], camera[1], camera[2]
    viewup = camera[3] if len(camera) > 3 else None
    if viewup is not None:
        fp, pos0 = cam.GetFocalPoint(), cam.GetPosition()
        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(pos0, fp))) or 1.0
        a, e = math.radians(az), math.radians(el)
        d = (math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e))
        cam.SetPosition(fp[0] + dist * d[0], fp[1] + dist * d[1], fp[2] + dist * d[2])
        cam.SetViewUp(*viewup)
    else:
        cam.Azimuth(az); cam.Elevation(el)
        cam.OrthogonalizeViewUp()
    cam.Zoom(zoom)


def _camera_for(case):
    """Canonical camera per case: (azimuth, elevation, zoom, viewUp).

    The 3D shear domain is tall in z ([0,1]^2 x [0,2]); it is viewed from the side with
    z vertical so the sheared interface stretches vertically downward, as reported in the
    literature. The unit-cube deformation keeps a three-quarter view. z is the view-up in
    both, so the interface is upright."""
    c = (case or "").lower()
    if "shear" in c:
        return (-72.0, 8.0, 0.9, (0.0, 0.0, 1.0))   # zoom<1 fits the tall vertical profile
    # deformation: az=90 opens the double-scroll so both curled "ears" are visible at T/2
    return (90.0, 25.0, 1.0, (0.0, 0.0, 1.0))


def make(study_dir, camera=None, theme="semi-lagrangian-level-set"):
    """One montage per ARM. Returns the list of PNGs written.

    Single-arm studies keep the historical `isosurface_<case>_<mesh>.png` name
    byte-identically -- the SL and LSL articles reference it directly.
    """
    cases = study_cases(study_dir)
    if not cases:
        print(f"[iso] no cases in {study_dir}; skip"); return []
    arms = {}
    for row in cases:
        arms.setdefault(row[5], []).append(row)
    outs = []
    for arm, arm_cases in sorted(arms.items()):
        out = _make_one(arm_cases, arm, len(arms) > 1, camera, theme)
        if out:
            outs.append(out)
    return outs


def _make_one(cases, arm, tag_arm, camera, theme):
    case, mesh = cases[0][3], cases[0][4]
    if camera is None:
        camera = _camera_for(case)
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
    for j, (h, label, cdir, _c, _m, _a) in enumerate(finest3):
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
    title = rf"$\psi=0$ interface — {case} ({mesh})"
    if tag_arm:
        title += f" — {arm}"
    fig.suptitle(title, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    suffix = f"_{arm}" if tag_arm else ""
    out = os.path.join(paths.figs_dir(theme),
                       f"isosurface_{case}_{mesh}{suffix}.png")
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"[iso] wrote {out}  ({ncol} resolutions x {nrow} times, camera={camera})")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("study_dir")
    ap.add_argument("--theme", default="semi-lagrangian-level-set",
                    help="docs theme receiving the figure (see paths.py)")
    args = ap.parse_args()
    sys.exit(0 if make(args.study_dir, theme=args.theme) else 1)
# make() returns the LIST of montages (one per arm); an empty list is the
# failure case, which is what the exit status above reports.
