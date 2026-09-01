#!/bin/bash
# Uniform-translation kinematic gate — standalone run.
#   circle R=0.15 at (0.25,0.5) in the unit box, carried by U = (1,0,0).
#   Exact solution: psi(x,t) = psi0(x - U t). Every error is transport error.
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

rm -rf 0 && cp -r 0.org 0
find 0 -name '*.template' -delete 2>/dev/null || true

runApplication blockMesh
runApplication leiaSetFields
runApplication leiaSemiLagrangeLevelSetFoam
