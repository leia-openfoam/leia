#!/usr/bin/env bash
# Polyhedral (Voronoi-type) mesh for a box case, WITHOUT cfMesh.
#
# WHY THIS EXISTS. The laptop's cfMesh build aborts with glibc heap corruption
# ("free(): invalid pointer" after "Finished smoothing mesh surface before
# mapping") on EVERY input tried -- including the transport suite's own unit-box
# configuration that works on Lichtenberg -- so `mesh: poly` studies cannot mesh
# locally. This script produces an equivalent polyhedral mesh from stock tools:
#   gmsh (tetrahedra, one physical surface "walls")  ->  gmshToFoam
#   ->  polyDualMesh <featureAngle> -overwrite       ->  checkMesh
# The dual of a tet mesh is a Voronoi-type polyhedral mesh: face orientations
# are scattered with no preferred angle, which is exactly the property the
# face-orientation study needs. Measured on the 6 mm box at lc = 1e-4:
# 166777 cells, one closed `walls` patch, max non-orthogonality 51 deg
# (average 12), max skewness 1.6 -- checkMesh clean.
#
# Usage:  make_poly_dual_mesh.sh <caseDir> <boxLength> <cellSize> [featureAngle]
# Writes constant/polyMesh in place; the case's 0.org must carry a `walls` patch.
set -euo pipefail
CASE=${1:?caseDir}; L=${2:?boxLength}; LC=${3:?cellSize}; ANGLE=${4:-45}
cd "$CASE"
cat > box.geo <<GEO
SetFactory("OpenCASCADE");
L = $L;
Box(1) = {0, 0, 0, L, L, L};
Characteristic Length{ PointsOf{ Volume{1}; } } = $LC;
Physical Surface("walls") = {1, 2, 3, 4, 5, 6};
Physical Volume("interior") = {1};
GEO
gmsh -3 box.geo -o box.msh -format msh2 -v 2 > log.gmsh 2>&1
gmshToFoam box.msh > log.gmshToFoam 2>&1
polyDualMesh "$ANGLE" -overwrite > log.polyDualMesh 2>&1
checkMesh > log.checkMesh 2>&1 || true
grep -E "cells:|non-orthogonality Max|Max skewness|Boundary openness" log.checkMesh
