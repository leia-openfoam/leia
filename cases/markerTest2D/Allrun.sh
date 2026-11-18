#!/bin/bash
# The marker test of step B (PDE plan section 6): graded circle psi0 = g(z)(|z| - R), q0 from 0.70 to
# 1.30 on the circle, moved by a linear flow; interfaceMarkers.csv reports q and psi at markers moved
# by RK4. Rendered by the workflow (cases/markerTest2D.parameter); standalone after rendering only.
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
rm -rf 0 && cp -r 0.org 0
runApplication blockMesh
runApplication leiaSetFields
runApplication leiaLevelSetFoam
