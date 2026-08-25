#!/bin/bash
# Classical reversed single-vortex (bulk) test — standalone run.
#   circle R=0.15 at (0.5,0.75) in the unit box, advected by the reversing
#   single-vortex field; returns to its initial shape at t = endTime.
# Edit system/fvSolution -> velocityExtension.type to try none |
# anisotropicDiffusion | pseudoTime.
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

rm -rf 0 && cp -r 0.org 0
find 0 -name '*.template' -delete 2>/dev/null || true

runApplication blockMesh
runApplication leiaSetFields
runApplication leiaLevelSetFoam
