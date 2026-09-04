#!/usr/bin/env python3
"""leiaRefineHexMesh -- static local refinement of a hex mesh around the interface.

Turns a level-set shape (fvSolution.levelSet.implicitSurface, read by
leiaSetFields) into a statically refined hexahedral mesh WITH re-initialised
initial fields on it. No solver change: leiaSemiLagrangianLevelSetTwoPhaseFoam
reads the result like any other mesh.

    blockMesh                                             (unless --no-base-mesh)
    for pass i = 1..N:
        0/ := 0.org ; leiaSetFields        (1) psi, alpha on the CURRENT mesh
        topoSet -dict system/topoSetDict.refine
                                           (2) refineCells = cells with 0 < alpha < 1
                                               (their faces are exactly the faces where
                                               snGrad(alpha) != 0 -- the capillary force
                                               support), dilated k = ceil(bandCells/2)
                                               times through faces
        refineHexMesh refineCells -overwrite
                                           (3) hexRef8 2x2x2 split; cellLevel/pointLevel
                                               are PERSISTED, so the 2:1 constraint holds
                                               across passes (OpenFOAM's plain refineMesh
                                               does not persist them -- never use it here)
    0/ := 0.org ; leiaSetFields            the FINAL fields: nothing mapped survives
    checkMesh ; postProcess -func writeCellVolumes ; band check
    -> refinement.csv (one row per pass), refinedBand.csv (the check), exit 2 on FAIL

How many dilations: the selection is made on the mesh BEFORE the split, so each
selected layer becomes two fine layers -- but face dilation grows a Manhattan
diamond, and along a sphere's diagonals one dilation buys ~1.15 fine cells rather
than 2 (measured at N = 60: 3/4/5 dilations -> 4.6/5.8/7.0 complete fine layers at
the worst point). The driver therefore starts at ceil(bandCells/2) dilations and,
using the psi it already has on the current mesh, adds dilations until every
UNSELECTED cell centre satisfies 2|psi|/h - 1 >= bandCells, i.e. until the split
is guaranteed to leave bandCells complete fine layers beyond the interface at the
worst point. psi only verifies the width; the criterion stays the snGrad(alpha)
seed. The stencil minimum is 4 layers (CPC curvature fit, one semi-Lagrangian
step, the explicit viscous transpose term -- the layer table in the article); the
final band check FAILS below it and WARNS below bandCells.

--source ball is the CONTROL arm: a sphereToCell ball of radius R + (1+k) h that
refines the whole droplet interior, to test whether the inner transition matters.

Tokens (case_params.json, overridable by flags): REFINE_LEVELS, REFINE_BAND_CELLS,
REFINE_SOURCE, DROPLET_RADIUS, DOMAIN_LENGTH, DOMAIN_HALF_LENGTH, N_CELLS_BASE,
N_CELLS (the capillary-dt handle; on the refined mesh it must be the FINE count and
the band check fails when it is not).

Standalone example, in a rendered case directory:
    python3 leiaRefineHexMesh.py --levels 1 --band-cells 6 \
        --setfields "leiaSetFields -alphaName alpha.water"
"""
import argparse
import math
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leia_refine as lr  # noqa: E402

TOPO_HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict.refine;
}}
// Written by workflow/scripts/leiaRefineHexMesh.py -- pass {pass_}, do not edit.
// Criterion: {criterion}
"""

SEED_INTERFACE = """    {{
        name    refineCells;
        type    cellSet;
        action  new;
        source  fieldToCell;
        field   {alpha};
        min     {lo:.1e};
        max     {hi:.9f};
    }}
"""

SEED_BALL = """    {{
        name    refineCells;
        type    cellSet;
        action  new;
        source  sphereToCell;
        origin  ({cx:.10g} {cy:.10g} {cz:.10g});
        radius  {radius:.10g};
    }}
"""

DILATION = """    {
        name    refineFaces;
        type    faceSet;
        action  new;
        source  cellToFace;
        set     refineCells;
        option  all;
    }
    {
        name    refineCells;
        type    cellSet;
        action  add;
        source  faceToCell;
        set     refineFaces;
        option  any;
    }
"""


def write_topo_set_dict(case, pass_, source, alpha_name, k, ball=None):
    if source == "interface":
        criterion = (f"seed 0 < {alpha_name} < 1 (the snGrad(alpha) support), "
                     f"then {k} face dilations")
        body = SEED_INTERFACE.format(alpha=alpha_name, lo=lr.ALPHA_LO, hi=lr.ALPHA_HI)
        body += DILATION * k
    elif source == "ball":
        cx, cy, cz, radius = ball
        criterion = f"CONTROL: ball of radius {radius:.6g} about ({cx:.6g} {cy:.6g} {cz:.6g})"
        body = SEED_BALL.format(cx=cx, cy=cy, cz=cz, radius=radius)
    else:
        raise ValueError(f"source must be interface or ball, got {source!r}")
    text = TOPO_HEADER.format(pass_=pass_, criterion=criterion) + "actions\n(\n" + body + ");\n"
    path = os.path.join(case, "system", "topoSetDict.refine")
    with open(path, "w") as fh:
        fh.write(text)
    return criterion


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("case", nargs="?", default=".")
    ap.add_argument("--levels", type=int, help="refinement passes N (token REFINE_LEVELS)")
    ap.add_argument("--band-cells", type=int,
                    help="target fine layers each side of the interface (REFINE_BAND_CELLS)")
    ap.add_argument("--source", choices=("interface", "ball"),
                    help="criterion: interface band (default) or the ball control (REFINE_SOURCE)")
    ap.add_argument("--setfields", default=None,
                    help='field initialiser command, e.g. "leiaSetFields -alphaName alpha.water"')
    ap.add_argument("--alpha-name", default=None, help="volume fraction field (from --setfields)")
    ap.add_argument("--no-base-mesh", action="store_true", help="do not run blockMesh first")
    ap.add_argument("--keep-sets", action="store_true", help="keep constant/polyMesh/sets")
    ap.add_argument("--allow-pin-mismatch", action="store_true",
                    help="report, but do not fail on, an N_CELLS that is not the fine count")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    case, dry = a.case, a.dry_run
    toks = lr.tokens(case)
    levels = a.levels if a.levels is not None else lr.tok_int(toks, "REFINE_LEVELS", 1)
    band_cells = a.band_cells if a.band_cells is not None else lr.tok_int(toks, "REFINE_BAND_CELLS", 6)
    source = a.source or toks.get("REFINE_SOURCE", "interface")
    setfields = a.setfields or "leiaSetFields -alphaName alpha.water"
    alpha_name = a.alpha_name or lr.alpha_name_from(setfields)
    if levels < 1:
        raise SystemExit(f"{lr.TAG} --levels / REFINE_LEVELS must be >= 1 (got {levels})")
    if band_cells < 1:
        raise SystemExit(f"{lr.TAG} --band-cells / REFINE_BAND_CELLS must be >= 1 (got {band_cells})")
    k0 = int(math.ceil(band_cells / 2.0))
    lr.say(f"hex refinement: levels={levels} bandCells={band_cells} (starting at {k0} "
           f"dilations/pass, added until psi proves the width), source={source}, "
           f"alpha={alpha_name}, setfields={setfields!r}")

    ball_inputs = None
    if source == "ball":
        R = lr.tok_float(toks, "DROPLET_RADIUS")
        L = lr.tok_float(toks, "DOMAIN_LENGTH")
        c = lr.tok_float(toks, "DOMAIN_HALF_LENGTH")
        nb = lr.tok_float(toks, "N_CELLS_BASE")
        if None in (R, L, c, nb) or nb <= 0:
            raise SystemExit(f"{lr.TAG} --source ball needs the tokens DROPLET_RADIUS, "
                             "DOMAIN_LENGTH, DOMAIN_HALF_LENGTH and N_CELLS_BASE in case_params.json")
        ball_inputs = (R, L / nb, c)

    if not a.no_base_mesh:
        # blockMesh rewrites points/faces/owner/neighbour/boundary but leaves a previous
        # run's cellLevel/pointLevel/refinementHistory in place, and hexRef8 READS them
        # (READ_IF_PRESENT): start from an empty polyMesh so a re-run cannot inherit
        # levels that belong to another mesh.
        if not dry:
            shutil.rmtree(os.path.join(case, "constant", "polyMesh"), ignore_errors=True)
        lr.run("blockMesh", "log.blockMesh", case, dry)
    if not dry:
        types = lr.boundary_types(case)
        if "empty" in types.values():
            raise SystemExit(f"{lr.TAG} the mesh has an empty patch (2D): refineHexMesh splits "
                             "2x2x2 and static interface refinement is 3D only")

    rows = []
    k_max = band_cells + 2
    for i in range(1, levels + 1):
        lr.say(f"---- pass {i}/{levels} ----")
        lr.reset_zero(case, dry)
        lr.run(setfields, f"log.leiaSetFields.pass{i}", case, dry)
        lr.run("postProcess -func writeCellVolumes -time 0", f"log.writeCellVolumes.pass{i}",
               case, dry)
        before = lr.n_cells(case) if not dry else -1
        k = k0
        layers = float("nan")
        while True:
            ball = None
            if ball_inputs:
                R, h_base, c = ball_inputs
                h_cur = h_base / 2 ** (i - 1)      # spacing near the interface BEFORE this pass
                ball = (c, c, c, R + (1 + k) * h_cur)
            criterion = write_topo_set_dict(case, i, source, alpha_name, k, ball)
            lr.run("topoSet -dict system/topoSetDict.refine", f"log.topoSet.pass{i}", case, dry)
            if dry:
                break
            selected = lr.read_set(case, "refineCells")
            layers = lr.worst_case_layers_after_split(case, selected, before)
            lr.say(f"pass {i}: {k} dilations -> {len(selected)} cells selected, worst-case "
                   f"complete fine layers after the split {layers:.2f} (requested {band_cells})")
            if layers >= band_cells or k >= k_max:
                break
            k += 1
        if not dry:
            os.remove(os.path.join(case, "0", "V"))   # not a field refineHexMesh should map
        n_sel = lr.set_size(case, "refineCells") if not dry else -1
        lr.run("refineHexMesh refineCells -overwrite", f"log.refineHexMesh.pass{i}", case, dry)
        after = lr.n_cells(case) if not dry else -1
        n_ref = -1
        if not dry:
            if (after - before) % 7:
                raise SystemExit(f"{lr.TAG} pass {i}: cell count grew by {after - before}, not a "
                                 "multiple of 7 -- a refined cell was not a hex")
            n_ref = (after - before) // 7
            lr.say(f"pass {i}: {before} cells, {n_sel} selected, {n_ref} refined "
                   f"(+{n_ref - n_sel} by consistentRefinement 2:1), {after} cells")
        rows.append({"pass": i, "nCellsBefore": before, "nSelected": n_sel, "nRefined": n_ref,
                     "growth2to1": (n_ref - n_sel) if not dry else "", "nCellsAfter": after,
                     "dilations": k, "fineLayersPredicted": layers, "source": source,
                     "criterion": criterion})

    lr.say("---- final fields on the final mesh ----")
    lr.reset_zero(case, dry)
    lr.run(setfields, "log.leiaSetFields.final", case, dry)
    lr.run("checkMesh", "log.checkMesh", case, dry)
    lr.run("postProcess -func writeCellVolumes -time 0", "log.writeCellVolumes", case, dry)
    if dry:
        return 0
    res = lr.band_check(case, "hex", band_cells, alpha_name, toks,
                        allow_pin_mismatch=a.allow_pin_mismatch)
    res.update(lr.checkmesh_stats(os.path.join(case, "log.checkMesh")))
    res["levels"] = levels
    lr.write_csv(os.path.join(case, "refinement.csv"), rows)
    lr.write_csv(os.path.join(case, "refinedBand.csv"), [res])
    lr.report(res)

    # Leave 0/ with exactly the initial fields, and no stale set labels (decomposePar
    # would otherwise try to decompose sets that no longer match the mesh).
    os.remove(os.path.join(case, "0", "V"))
    if not a.keep_sets:
        shutil.rmtree(os.path.join(case, "constant", "polyMesh", "sets"), ignore_errors=True)
    return 0 if res["status"] != "FAIL" else 2


if __name__ == "__main__":
    sys.exit(main())
