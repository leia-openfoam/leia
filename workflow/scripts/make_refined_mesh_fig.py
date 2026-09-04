#!/usr/bin/env python3
"""Mid-plane slice of the statically refined meshes: cells coloured by refinement level,
edges drawn, the psi = 0 interface in red. One panel per arm, one PNG.

Level is measured from the built mesh -- round(log2(h_coarsest / V_c^(1/3))) -- so the
same colouring works for hexRef8 levels and for a cfMesh graded polyhedral mesh; the
cellLevel file is not needed (and vtkOpenFOAMReader does not read it). Polyhedra are
kept whole (DecomposePolyhedraOff), otherwise the transition cells would be split into
pyramids with fractional volumes and mis-levelled.

Usage:  python3 workflow/scripts/make_refined_mesh_fig.py [--root .]
            [--studies stationaryDroplet3Drefined,stationaryDroplet3DrefinedL2,...]
            [--out docs/.../data/figures/refined_mesh_slice.png] [--time 0]
"""
import argparse
import csv
import glob
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.image as mpimg         # noqa: E402
import vtk                               # noqa: E402
from vtk.util import numpy_support as ns  # noqa: E402

PANEL_PX = 900
DEFAULT_STUDIES = ("stationaryDroplet3Drefined", "stationaryDroplet3DrefinedL2",
                   "stationaryDroplet3DrefinedBall", "polyDroplet3Drefined_r13p8")
# level 0 (coarsest) light, finer levels progressively darker blue-greys
LEVEL_RGB = [(0.97, 0.97, 0.95), (0.80, 0.86, 0.93), (0.58, 0.70, 0.86), (0.36, 0.52, 0.76)]


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def internal_mesh(mb):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        blk = it.GetCurrentDataObject()
        if blk is not None and blk.IsA("vtkUnstructuredGrid") and blk.GetNumberOfCells() > 0:
            return blk
        it.GoToNextItem()
    return None


def render_slice(case_dir, z, out_png, t=0.0, half_width=None, centre=None):
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(foam)
    reader.CreateCellToPointOn()
    reader.SetSkipZeroTime(0)
    if hasattr(reader, "DecomposePolyhedraOff"):
        reader.DecomposePolyhedraOff()
    reader.Update()
    tv = reader.GetTimeValues()
    tvals = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())] if tv else [0.0]
    reader.UpdateTimeStep(min(tvals, key=lambda x: abs(x - t)))
    grid = internal_mesh(reader.GetOutput())
    if grid is None:
        return None

    size = vtk.vtkCellSizeFilter()
    size.SetInputData(grid)
    size.ComputeVolumeOn(); size.ComputeAreaOff(); size.ComputeLengthOff(); size.ComputeVertexCountOff()
    size.Update()
    g = size.GetOutput()
    vol = ns.vtk_to_numpy(g.GetCellData().GetArray("Volume"))
    h = np.cbrt(np.abs(vol))
    level = np.rint(np.log2(h.max() / h)).clip(0, len(LEVEL_RGB) - 1)
    arr = ns.numpy_to_vtk(level.astype(np.float64), deep=1); arr.SetName("level")
    g.GetCellData().AddArray(arr)
    g.GetCellData().SetActiveScalars("level")

    # On a blockMesh box the mid-plane IS a face plane, so the honest picture is the set
    # of cell FACES lying on it: keep the cells whose centre is below the plane, take the
    # exposed boundary of that half, and select the faces sitting on the plane. Those are
    # the real faces (quads for hexahedra, the split quads of the hanging-node transition
    # cells), with the parent cell's level as cell data. A mesh with cell centres ON the
    # plane (a polyhedral cfMesh mesh) has no such face plane and falls back to a cutter.
    sl, z_cut = plane_faces(g, z, h.min())
    if sl is None:
        z_cut = z + 0.5 * h.min()
        plane = vtk.vtkPlane(); plane.SetOrigin(0.0, 0.0, z_cut); plane.SetNormal(0.0, 0.0, 1.0)
        cutter = vtk.vtkCutter(); cutter.SetCutFunction(plane); cutter.SetInputData(g); cutter.Update()
        sl = cutter.GetOutput()

    lut = vtk.vtkLookupTable(); lut.SetNumberOfTableValues(len(LEVEL_RGB))
    for i, c in enumerate(LEVEL_RGB):
        lut.SetTableValue(i, c[0], c[1], c[2], 1.0)
    lut.SetRange(0, len(LEVEL_RGB) - 1); lut.Build()
    fill_m = vtk.vtkPolyDataMapper(); fill_m.SetInputData(sl)
    fill_m.SetScalarModeToUseCellFieldData(); fill_m.SelectColorArray("level")
    fill_m.SetLookupTable(lut); fill_m.SetScalarRange(0, len(LEVEL_RGB) - 1)
    fill = vtk.vtkActor(); fill.SetMapper(fill_m); fill.GetProperty().SetAmbient(1.0); fill.GetProperty().SetDiffuse(0.0)

    edges = vtk.vtkExtractEdges(); edges.SetInputData(sl); edges.Update()
    edge_m = vtk.vtkPolyDataMapper(); edge_m.SetInputConnection(edges.GetOutputPort()); edge_m.ScalarVisibilityOff()
    edge = vtk.vtkActor(); edge.SetMapper(edge_m)
    edge.GetProperty().SetColor(0.25, 0.25, 0.25); edge.GetProperty().SetLineWidth(1.0)

    iso = None
    if sl.GetPointData().GetArray("psi") is not None:
        sl.GetPointData().SetActiveScalars("psi")
        cont = vtk.vtkContourFilter(); cont.SetInputData(sl); cont.SetValue(0, 0.0); cont.Update()
        iso_m = vtk.vtkPolyDataMapper(); iso_m.SetInputConnection(cont.GetOutputPort()); iso_m.ScalarVisibilityOff()
        iso = vtk.vtkActor(); iso.SetMapper(iso_m)
        iso.GetProperty().SetColor(0.80, 0.10, 0.10); iso.GetProperty().SetLineWidth(3.0)

    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
    ren.AddActor(fill); ren.AddActor(edge)
    if iso is not None:
        ren.AddActor(iso)
    win = vtk.vtkRenderWindow(); win.SetOffScreenRendering(1); win.AddRenderer(ren); win.SetSize(PANEL_PX, PANEL_PX)
    cam = ren.GetActiveCamera(); cam.ParallelProjectionOn()
    b = sl.GetBounds()
    cx, cy = (centre if centre else (0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])))
    hw = half_width if half_width else 0.5 * max(b[1] - b[0], b[3] - b[2])
    cam.SetFocalPoint(cx, cy, z_cut); cam.SetPosition(cx, cy, z_cut + 1.0); cam.SetViewUp(0, 1, 0)
    cam.SetParallelScale(hw)
    ren.ResetCameraClippingRange()
    win.Render()
    w2i = vtk.vtkWindowToImageFilter(); w2i.SetInput(win); w2i.Update()
    wr = vtk.vtkPNGWriter(); wr.SetFileName(out_png); wr.SetInputConnection(w2i.GetOutputPort()); wr.Write()
    win.Finalize()
    return out_png if os.path.isfile(out_png) else None


def _threshold_between(ds, name, lo, hi):
    thr = vtk.vtkThreshold(); thr.SetInputData(ds)
    thr.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, name)
    thr.SetLowerThreshold(lo); thr.SetUpperThreshold(hi)
    thr.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
    thr.Update()
    return thr.GetOutput()


def plane_faces(g, z, h_min):
    """Faces of the mesh lying on the plane z (None if z is not a face plane)."""
    cc = vtk.vtkCellCenters(); cc.SetInputData(g); cc.Update()
    zc = ns.vtk_to_numpy(cc.GetOutput().GetPoints().GetData())[:, 2]
    if np.any(np.abs(zc - z) < 1e-3 * h_min):
        return None, None
    below = (zc < z).astype(np.float64)
    arr = ns.numpy_to_vtk(below, deep=1); arr.SetName("below"); g.GetCellData().AddArray(arr)
    half = _threshold_between(g, "below", 0.5, 1.5)
    surf = vtk.vtkDataSetSurfaceFilter(); surf.SetInputData(half); surf.Update()
    pd = surf.GetOutput()
    fc = vtk.vtkCellCenters(); fc.SetInputData(pd); fc.Update()
    zf = ns.vtk_to_numpy(fc.GetOutput().GetPoints().GetData())[:, 2]
    on = (np.abs(zf - z) < 1e-3 * h_min).astype(np.float64)
    arr2 = ns.numpy_to_vtk(on, deep=1); arr2.SetName("onplane"); pd.GetCellData().AddArray(arr2)
    faces = _threshold_between(pd, "onplane", 0.5, 1.5)
    geo = vtk.vtkGeometryFilter(); geo.SetInputData(faces); geo.Update()
    return geo.GetOutput(), z


def arms(root, study):
    out = []
    for cp in sorted(glob.glob(os.path.join(root, "studies", study, "*_[0-9]*", "case_params.json"))):
        d = json.load(open(cp)); tok = d.get("tokens", {}); case_dir = os.path.dirname(cp)
        if not os.path.isfile(os.path.join(case_dir, "constant", "polyMesh", "owner")):
            continue
        band = {}
        bp = os.path.join(case_dir, "refinedBand.csv")
        if os.path.isfile(bp):
            rows = list(csv.DictReader(open(bp))); band = rows[0] if rows else {}
        out.append({"dir": case_dir, "tok": tok, "mesh": d.get("mesh"), "band": band, "study": study})
    out.sort(key=lambda a: (_f(a["tok"].get("N_CELLS")) or 0, a["tok"].get("REFINE_SOURCE", "")))
    return out


def label_of(a):
    tok, band = a["tok"], a["band"]
    n = _f(tok.get("N_CELLS")); lev = tok.get("REFINE_LEVELS", "?"); src = tok.get("REFINE_SOURCE", "")
    kind = "poly" if (a["mesh"] or "").startswith("poly") else "hex"
    s = f"{kind}, $N_{{\\mathrm{{fine}}}}={int(n) if n else '?'}$, {lev} level{'s' if lev != '1' else ''}"
    if src and src != "interface":
        s += f", {src} control"
    if band.get("nCells"):
        s += f"\n{int(float(band['nCells'])):,} cells, {float(band.get('fineLayersWorst', 'nan')):.1f} fine layers (worst)"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--studies", default=",".join(DEFAULT_STUDIES))
    ap.add_argument("--out", default="docs/semi-lagrangian-level-set/sl-level-set-article/"
                                     "data/figures/refined_mesh_slice.png")
    ap.add_argument("--time", type=float, default=0.0)
    ap.add_argument("--zoom", type=float, default=None,
                    help="half-width of the view in metres (default: the whole slice)")
    a = ap.parse_args()
    panels = []
    scratch = os.path.join(a.root, ".refined_mesh_fig")
    os.makedirs(scratch, exist_ok=True)
    for study in [s for s in a.studies.split(",") if s]:
        for arm in arms(a.root, study):
            tok = arm["tok"]
            L = _f(tok.get("DOMAIN_LENGTH")) or 0.006
            centre = (L / 2, L / 2)
            png = os.path.join(scratch, f"{study}_{os.path.basename(arm['dir'])}.png")
            got = render_slice(arm["dir"], L / 2, png, a.time, a.zoom, centre if a.zoom else None)
            if got:
                panels.append((label_of(arm), got))
                print(f"[refined_mesh_fig] {arm['dir']} -> {got}")
    if not panels:
        raise SystemExit("[refined_mesh_fig] no arm with a mesh found")
    ncol = min(4, len(panels)); nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.9 * nrow), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, (lab, png) in zip(axes.flat, panels):
        ax.imshow(mpimg.imread(png)); ax.set_title(lab, fontsize=8)
    fig.tight_layout()
    out = os.path.join(a.root, a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=170); plt.close(fig)
    print(f"[refined_mesh_fig] wrote {out} ({len(panels)} panels)")


if __name__ == "__main__":
    main()
