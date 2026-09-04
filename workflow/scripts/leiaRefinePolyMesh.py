#!/usr/bin/env python3
"""leiaRefinePolyMesh -- static local refinement of a cfMesh polyhedral mesh around
the interface, by RE-MESHING around the level set's zero iso-surface.

cfMesh has no in-place refiner and OpenFOAM-v2512's only one (hexRef8) is
hex-only, so the polyhedral analogue of leiaRefineHexMesh.py regenerates the mesh
with the criterion expressed as geometry: the psi = 0 iso-surface of the level set
initialised on the current mesh, plus a refinement thickness in fine cells. The
loop, per pass i = 1..N, with h_i = MAX_CELL_SIZE / 2^i:

    pMesh                                                  pass 0: uniform maxCellSize
    for pass i = 1..N:
        0/ := 0.org ; leiaSetFields        (1) psi, alpha on the CURRENT mesh
        postProcess -func isoInterface -time 0
                                           (2) the psi = 0 iso-surface as STL
                                               -> constant/triSurface/interface_pass<i>.stl
        system/meshDict += surfaceMeshRefinement { interface { surfaceFile;
            additionalRefinementLevels i; refinementThickness bandCells*h_i } }
            (--size-mode cellSize asks for h_i instead; MEASURED: cfMesh rounded a
            cellSize of half the base to TWO octree levels, so levels is the default)
        pMesh                              (3) the graded mesh, regenerated
    0/ := 0.org ; leiaSetFields            the FINAL fields: cfMesh has no fields, and
                                           nothing is ever mapped
    checkMesh ; postProcess -func writeCellVolumes ; band check
    -> refinement.csv, refinedBand.csv, exit 2 on FAIL

The iso-surface extracted from mesh i-1 is faceted to within h_{i-1} = 2 h_i; the
thickness bandCells*h_i (6 h_i by default) covers that error, and pass i+1
re-extracts it from the finer psi. REFINE_LEVELS therefore means the same as for
hex: the near-interface cell size is halved N times. cfMesh grades beyond the
requested thickness (measured: 11.8 complete fine layers for 6 requested at one
level, 1.78x the coarse cell count) -- a cost the band check reports as
fineExtentOverHfine, not a defect.

N_CELLS IS THE CAPILLARY-dt HANDLE and cfMesh's realised cell size is 0.63-0.76 x
the requested one, so the band check measures the fine spacing on the BUILT mesh
and prints the N_CELLS to pin (round(DOMAIN_LENGTH/hFine) from the median
interface-cell size, the same average the uniform polyhedral rungs were pinned
with). The mesh does not depend on N_CELLS, so pin it in the config and re-render:
the check FAILS on a pin more than 5 % off unless --allow-pin-mismatch is given
(the pinning pass itself).

Tokens (case_params.json, overridable): REFINE_LEVELS, REFINE_BAND_CELLS,
MAX_CELL_SIZE, DOMAIN_LENGTH, N_CELLS.
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr  # noqa: E402

FUNC = "isoInterface"
# cfMesh on glibc >= 2.39 aborts with `free(): invalid pointer` unless jemalloc is
# preloaded, but a STUDY-GLOBAL LD_PRELOAD segfaults the MPI solver at startup (empty log,
# rc 139 -- measured 2026-09-04 on this very smoke). The workaround is therefore scoped to
# the one tool that needs it: set LEIA_PMESH_PREFIX="env LD_PRELOAD=/path/libjemalloc.so.2"
# in the study's env_preamble and only the pMesh command is prefixed.
PMESH = (os.environ.get("LEIA_PMESH_PREFIX", "").strip() + " pMesh").strip()
MARK_BEGIN = "// >>> leiaRefinePolyMesh: written by workflow/scripts/leiaRefinePolyMesh.py, do not edit"
MARK_END = "// <<< leiaRefinePolyMesh"

# A `postProcess -func isoInterface` configuration (system/isoInterface, read like
# the etc/caseDicts/postProcessing entries -- no FoamFile header).
ISO_DICT = """// Written by workflow/scripts/leiaRefinePolyMesh.py -- the psi = 0 iso-surface of
// the level set on the current mesh, as an STL for cfMesh's surfaceMeshRefinement.
type            surfaces;
libs            (sampling);
writeControl    writeTime;
surfaceFormat   stl;
interpolationScheme cellPoint;
fields          (psi);
surfaces
{
    interface
    {
        type        isoSurface;
        isoField    psi;
        isoValue    0;
        interpolate true;
    }
}
"""

REFINEMENT_BLOCK = """{begin}
// pass {pass_}: {size_comment}, refinementThickness = {band} x (MAX_CELL_SIZE / 2^{pass_})
surfaceMeshRefinement
{{
    interface
    {{
        surfaceFile         "{stl}";
        {size_entry}
        refinementThickness {thickness:.10g};
    }}
}}
{end}
"""


def strip_marked(text):
    while MARK_BEGIN in text:
        b = text.index(MARK_BEGIN)
        e = text.index(MARK_END, b) + len(MARK_END)
        if e < len(text) and text[e] == "\n":
            e += 1
        text = text[:b] + text[e:]
    return text


def set_mesh_dict_refinement(case, pass_, stl_rel, cell_size, thickness, band, size_mode):
    path = os.path.join(case, "system", "meshDict")
    with open(path) as fh:
        text = strip_marked(fh.read())
    if size_mode == "levels":
        # cfMesh refines in whole octree levels; asking for the integer directly is what
        # REFINE_LEVELS means and avoids the mesher's own rounding of a cellSize.
        size_entry = f"additionalRefinementLevels {pass_};"
        size_comment = f"additionalRefinementLevels = {pass_} (one octree level per pass)"
    else:
        size_entry = f"cellSize            {cell_size:.10g};"
        size_comment = f"cellSize = MAX_CELL_SIZE / 2^{pass_} = {cell_size:.6g}"
    block = REFINEMENT_BLOCK.format(begin=MARK_BEGIN, end=MARK_END, pass_=pass_, band=band,
                                    stl=stl_rel, size_entry=size_entry, size_comment=size_comment,
                                    thickness=thickness)
    tail_marker = "// *****"
    i = text.rfind(tail_marker)
    text = (text[:i] + block + "\n" + text[i:]) if i >= 0 else (text.rstrip("\n") + "\n\n" + block)
    with open(path, "w") as fh:
        fh.write(text)


def clear_mesh_dict_refinement(case):
    path = os.path.join(case, "system", "meshDict")
    with open(path) as fh:
        text = fh.read()
    with open(path, "w") as fh:
        fh.write(strip_marked(text))


def extract_interface(case, pass_, dry):
    with open(os.path.join(case, "system", FUNC), "w") as fh:
        fh.write(ISO_DICT)
    shutil.rmtree(os.path.join(case, "postProcessing", FUNC), ignore_errors=True)
    lr.run(f"postProcess -func {FUNC} -time 0", f"log.{FUNC}.pass{pass_}", case, dry)
    dst_rel = os.path.join("constant", "triSurface", f"interface_pass{pass_}.stl")
    if dry:
        return dst_rel
    stls = lr.find_stl(case, FUNC)
    if not stls:
        # A writer that could not emit STL directly: convert whatever it wrote.
        import glob
        others = sorted(glob.glob(os.path.join(case, "postProcessing", FUNC, "**", "*.vt*"),
                                  recursive=True))
        if not others:
            raise SystemExit(f"{lr.TAG} postProcess -func {FUNC} wrote no surface; see "
                             f"log.{FUNC}.pass{pass_}")
        conv = os.path.join(case, "postProcessing", FUNC, "interface.stl")
        lr.run(f"surfaceMeshConvert {os.path.relpath(others[0], case)} {os.path.relpath(conv, case)}",
               f"log.surfaceMeshConvert.pass{pass_}", case, dry)
        stls = [conv]
    src = stls[0]
    if os.path.getsize(src) < 100:
        raise SystemExit(f"{lr.TAG} {src} is empty: the psi = 0 iso-surface was not found")
    os.makedirs(os.path.join(case, "constant", "triSurface"), exist_ok=True)
    shutil.copyfile(src, os.path.join(case, dst_rel))
    return dst_rel


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("case", nargs="?", default=".")
    ap.add_argument("--levels", type=int, help="halvings of the near-interface size (REFINE_LEVELS)")
    ap.add_argument("--band-cells", type=int,
                    help="refinement thickness in fine cells each side (REFINE_BAND_CELLS)")
    ap.add_argument("--max-cell-size", type=float, help="pass-0 cell size (MAX_CELL_SIZE)")
    ap.add_argument("--size-mode", choices=("levels", "cellSize"), default="levels",
                    help="how the near-interface size is requested from cfMesh: whole octree "
                         "levels (default; what REFINE_LEVELS means) or a cellSize the mesher rounds")
    ap.add_argument("--setfields", default=None,
                    help='field initialiser command, e.g. "leiaSetFields -alphaName alpha.water"')
    ap.add_argument("--alpha-name", default=None)
    ap.add_argument("--no-base-mesh", action="store_true",
                    help="reuse the existing mesh as pass 0 (do not run pMesh first)")
    ap.add_argument("--allow-pin-mismatch", action="store_true",
                    help="report, but do not fail on, an N_CELLS that is not the fine spacing")
    ap.add_argument("--keep-postprocessing", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    case, dry = a.case, a.dry_run
    toks = lr.tokens(case)
    levels = a.levels if a.levels is not None else lr.tok_int(toks, "REFINE_LEVELS", 1)
    band_cells = a.band_cells if a.band_cells is not None else lr.tok_int(toks, "REFINE_BAND_CELLS", 6)
    mcs = a.max_cell_size if a.max_cell_size is not None else lr.tok_float(toks, "MAX_CELL_SIZE")
    setfields = a.setfields or "leiaSetFields -alphaName alpha.water"
    alpha_name = a.alpha_name or lr.alpha_name_from(setfields)
    if levels < 1:
        raise SystemExit(f"{lr.TAG} --levels / REFINE_LEVELS must be >= 1 (got {levels})")
    if mcs is None or mcs <= 0:
        raise SystemExit(f"{lr.TAG} --max-cell-size / MAX_CELL_SIZE must be > 0 (got {mcs})")
    lr.say(f"poly refinement: levels={levels} bandCells={band_cells} maxCellSize={mcs:g} "
           f"-> fine {mcs / 2 ** levels:g}, alpha={alpha_name}, setfields={setfields!r}")

    if not a.no_base_mesh:
        clear_mesh_dict_refinement(case)
        if not dry:
            shutil.rmtree(os.path.join(case, "constant", "polyMesh"), ignore_errors=True)
        lr.run(PMESH, "log.pMesh.pass0", case, dry)

    rows = []
    for i in range(1, levels + 1):
        lr.say(f"---- pass {i}/{levels} ----")
        before = lr.n_cells(case) if not dry else -1
        lr.reset_zero(case, dry)
        lr.run(setfields, f"log.leiaSetFields.pass{i}", case, dry)
        stl_rel = extract_interface(case, i, dry)
        cell_size = mcs / 2 ** i
        thickness = band_cells * cell_size
        set_mesh_dict_refinement(case, i, stl_rel, cell_size, thickness, band_cells, a.size_mode)
        lr.run(PMESH, f"log.pMesh.pass{i}", case, dry)
        after = lr.n_cells(case) if not dry else -1
        lr.say(f"pass {i}: {before} -> {after} cells, cellSize {cell_size:g}, "
               f"thickness {thickness:g}, surface {stl_rel}")
        rows.append({"pass": i, "nCellsBefore": before, "nCellsAfter": after,
                     "cellSize": cell_size, "refinementThickness": thickness,
                     "surfaceFile": stl_rel})

    lr.say("---- final fields on the final mesh ----")
    lr.reset_zero(case, dry)
    lr.run(setfields, "log.leiaSetFields.final", case, dry)
    lr.run("checkMesh", "log.checkMesh", case, dry)
    lr.run("postProcess -func writeCellVolumes -time 0", "log.writeCellVolumes", case, dry)
    if dry:
        return 0
    res = lr.band_check(case, "poly", band_cells, alpha_name, toks,
                        allow_pin_mismatch=a.allow_pin_mismatch
                        or lr.tok_bool(toks, "REFINE_ALLOW_PIN_MISMATCH"))
    res.update(lr.checkmesh_stats(os.path.join(case, "log.checkMesh")))
    res["levels"] = levels
    res["sizeMode"] = a.size_mode
    res["maxCellSize"] = mcs
    res["cellSizeRequested"] = mcs / 2 ** levels
    lr.write_csv(os.path.join(case, "refinement.csv"), rows)
    lr.write_csv(os.path.join(case, "refinedBand.csv"), [res])
    lr.report(res)

    os.remove(os.path.join(case, "0", "V"))
    if not a.keep_postprocessing:
        shutil.rmtree(os.path.join(case, "postProcessing", FUNC), ignore_errors=True)
    return 0 if res["status"] != "FAIL" else 2


if __name__ == "__main__":
    sys.exit(main())
