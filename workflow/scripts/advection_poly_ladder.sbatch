#!/bin/bash
#SBATCH -A special00004
#SBATCH -J leia-advpoly
#SBATCH -n 48
#SBATCH --mem-per-cpu=8000
#SBATCH -t 08:00:00
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#
# ONE rung of the POLYHEDRAL advection ladder: build ONE pMesh, then run every
# value-bound arm on THAT mesh.
#
# WHEN TO USE THIS RATHER THAN A SNAKEMAKE STUDY. Bit-identical meshes are NOT
# required for a polyhedral sweep, and config/advConv3DshearPoly.yaml is the normal
# route: let each arm build its own mesh and read the arms against a MEASURED
# mesh-noise floor. Use this script only when the effect under test is expected to sit
# AT that floor, so the mesh contribution has to be removed rather than bounded.
#
# The mesh really does vary. The workflow's `mesh` rule runs once per ARM, and
# MEASURED 2026-09-10 on the WSL laptop: two runs of pMesh from the IDENTICAL rendered
# case gave different meshes -- digests 7e9ee670c679 and e38b4e8e919b, and even
# different file sizes (8277379 against 8277373 bytes). cfMesh is therefore not
# reproducible run to run, not merely across a change of scale. A polyhedral parameter
# sweep driven by the workflow gives every arm its own mesh, and a cross-arm metric
# comparison carries that difference too. This script builds the mesh once and copies
# it into every arm, and prints the digest so the claim is checkable rather than
# asserted.
#
# --mem-per-cpu is mandatory on Lichtenberg. It is set high because pMesh is SERIAL
# and holds the whole mesh in ONE process, so it does not get cheaper as -n rises.
#
# Usage (submit one job per rung, from the repo root on the cluster):
#   SRC=studies/<meshStudy>/3Dshear_00000 NP=48 \
#     sbatch workflow/scripts/advection_poly_ladder.sbatch
#
#   SRC     a RENDERED case (snakemake --until generate_case), whose system/meshDict
#           already carries this rung's maxCellSize
#   NP      MPI ranks for the solver arms (default 48)
#   BOUNDS  space-separated valueBound list (default "none lipschitzCone stencilBounds")
#   OUTDIR  where the arms are written (default <SRC>/../../advPolyArms)

set +u
: "${SRC:?SRC must be a rendered case directory}"
NP=${NP:-48}
BOUNDS=${BOUNDS:-"none lipschitzCone stencilBounds"}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SRC=$(cd "$SRC" && pwd)
OUTDIR=${OUTDIR:-$(dirname "$(dirname "$SRC")")/advPolyArms}
mkdir -p "$OUTDIR"

# bashrc FIRST, then .leia_env: the bashrc RESETS WM_PROJECT_USER_DIR, so an export
# before it lands in the SHARED account prefix and overwrites another session's
# binaries. That happened on 2026-09-09.
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"
echo "prefix   = $WM_PROJECT_USER_DIR"
echo "src      = $SRC"
echo "maxCellSize = $(foamDictionary -entry maxCellSize -value "$SRC/system/meshDict")"

# ---- ONE mesh, built once -------------------------------------------------------
cd "$SRC" || exit 1
if [ ! -f constant/polyMesh/owner ]; then
    rm -rf 0 && cp -r 0.org 0 && find 0 -name '*.template' -delete
    ${LEIA_PMESH_PRELOAD:+env LD_PRELOAD=$LEIA_PMESH_PRELOAD} pMesh > log.pMesh 2>&1
    rc=$?
    if [ ! -f constant/polyMesh/owner ]; then
        echo "ABORT: pMesh produced no mesh (rc=$rc)"; tail -20 log.pMesh; exit 2
    fi
fi
DIG=$(md5sum constant/polyMesh/points | cut -c1-12)
NC=$(checkMesh -time 0 2>/dev/null | grep -m1 'cells:' | tr -dc '0-9')
echo "mesh     = $NC cells, digest $DIG"

# ---- every arm on THAT mesh -----------------------------------------------------
for B in $BOUNDS; do
    OUT="$OUTDIR/${B}"
    echo "=== arm $B -> $OUT ==="
    bash "$REPO/workflow/scripts/advect_bound_arm.sh" "$SRC" "$OUT" "$B" \
         unity cellOnly "$NP" leiaSemiLagrangeLevelSetFoam
    # Verify the ARTEFACT and the OUTCOME, never the return code.
    bash "$REPO/workflow/scripts/foam_log_state.sh" "$OUT/log.solve"
done

echo "=== census (every arm read at the LAST COMMON step) ==="
python3 "$REPO/workflow/scripts/value_bound_advection_census.py" \
    $(for B in $BOUNDS; do echo "$OUTDIR/$B"; done) \
    --label "poly mcs=$(foamDictionary -entry maxCellSize -value "$SRC/system/meshDict") cells=$NC" \
    --csv "$OUTDIR/census.csv"
