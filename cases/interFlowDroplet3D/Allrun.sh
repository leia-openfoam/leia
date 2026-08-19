#!/bin/bash
cd "${0%/*}" || exit 1
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

rm -rf 0 && cp -r 0.org 0

runApplication blockMesh
# The SHARED initialiser: identical initial alpha.water to the level-set and
# interFoam cases, so the three solvers differ only in their interface method.
runApplication leiaSetFields -alphaName alpha.water
runApplication interFlow
