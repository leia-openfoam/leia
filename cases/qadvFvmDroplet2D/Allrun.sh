#!/bin/bash
cd "${0%/*}" || exit 1
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

rm -rf 0 && cp -r 0.org 0

runApplication blockMesh
runApplication leiaSetFields -alphaName alpha.water
runApplication "$(getApplication)"
