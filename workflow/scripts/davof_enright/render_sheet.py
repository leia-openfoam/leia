#!/usr/bin/env pvpython
"""The alpha = 0.5 surface of the 3D deformation runs (2026-11-06): one PNG per
(case, time) with ParaView's OpenFOAM reader on the decomposed or serial case,
cell-to-point interpolated alpha.water contoured at 0.5, a fixed camera.
    pvpython render_sheet.py <case dir> <time> <out.png> [<title>]
The case: cases/davof/sphereDeformation3D (LeVeque 1996, doi:10.1137/0733033;
Enright et al. 2002, doi:10.1006/jcph.2002.7166).
"""
import os, sys
from paraview.simple import *

case, t, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
title = sys.argv[4] if len(sys.argv) > 4 else ""
foam = os.path.join(case, "case.foam")
open(foam, "a").close()
r = OpenFOAMReader(FileName=foam)
r.CaseType = "Decomposed Case" if os.path.isdir(os.path.join(case, "processor0")) else "Reconstructed Case"
r.MeshRegions = ["internalMesh"]
r.CellArrays = ["alpha.water"]
r.CreateCellToPoint = 1
r.UpdatePipelineInformation()
times = list(r.TimestepValues) if hasattr(r.TimestepValues, "__len__") else [r.TimestepValues]
tt = min(times, key=lambda x: abs(x - t))
if abs(tt - t) > 1e-6:
    sys.stderr.write("render_sheet: no time %g in %s (have %s)\n" % (t, case, times)); sys.exit(2)
c = Contour(Input=r, ContourBy=["POINTS", "alpha.water"], Isosurfaces=[0.5])
c.ComputeNormals = 1
view = CreateView("RenderView")
view.ViewSize = [1200, 1000]
view.Background = [1, 1, 1]
try: view.UseColorPaletteForBackground = 0
except Exception: pass
view.OrientationAxesVisibility = 1
d = Show(c, view)
d.DiffuseColor = [0.25, 0.45, 0.85]
d.AmbientColor = [0.25, 0.45, 0.85]
ColorBy(d, None)
box = Outline(Input=r)
db = Show(box, view); db.DiffuseColor = [0.3, 0.3, 0.3]; db.AmbientColor = [0.3, 0.3, 0.3]
view.ViewTime = tt
view.CameraPosition = [2.6, 1.9, 2.4]
view.CameraFocalPoint = [0.5, 0.5, 0.5]
view.CameraViewUp = [0, 0, 1]
view.CameraViewAngle = 28
if title:
    tx = Text(Text=title); dt = Show(tx, view); dt.Color = [0, 0, 0]; dt.FontSize = 20; dt.WindowLocation = "Upper Left Corner"
Render(view)
SaveScreenshot(out, view, ImageResolution=[1200, 1000])
print("wrote", out, "time", tt)
