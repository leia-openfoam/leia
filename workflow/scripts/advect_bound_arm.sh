#!/usr/bin/env bash
# Run ONE value-bound arm of an ADVECTION case on a GIVEN, SHARED mesh.
#
# WHY THIS EXISTS RATHER THAN A STUDY. The workflow's `mesh` rule runs once per ARM, and
# cfMesh is not bit-reproducible, so a polyhedral sweep gives every arm its own mesh -- and
# then a cross-arm metric comparison measures the MESHER, not the token. That is the same
# trap the clip gate table guards with a mesh digest, and there it can only report NOT
# TESTABLE. This script removes the trap instead of reporting it: the mesh is built ONCE by
# the workflow (`--until preprocess` on a single-arm config) and copied into every arm, so
# the only difference between arms is the dictionary entry. The digest is printed so the
# claim is checkable.
#
# WHY AN ADVECTION CASE. Kinematics before dynamics. `leiaSemiLagrangeLevelSetFoam` runs the
# transport with a PRESCRIBED velocity and no force in the loop, so a failure here is a
# transport failure and cannot be blamed on the force balance. MEASURED 2026-09-10: this rung
# falsified the distance-cone bound after one coupled case had made it look good.
#
# Pick the case for the displacement it produces IN THE CELLS THE BOUND CAN AFFECT. A shear
# or deformation gate has u = 0 on the walls, so it never tests the boundary stencils;
# `2Dtranslation` (|U| = 1 everywhere) is the O(1)-displacement gate.
#
# Usage:
#   workflow/scripts/advect_bound_arm.sh <meshed-src-case> <out-dir> <valueBound> \
#       [lipschitzMode] [onInadmissible] [np] [solver] [setfields-args]
#
#     <meshed-src-case>  a case that already has constant/polyMesh, e.g.
#                        studies/<meshStudy>/<case>_00000 after `--until preprocess`
#     <valueBound>       none | stencilBounds | lipschitzCone
#
# Prints the md5 of constant/polyMesh/points: arms are comparable only when it MATCHES.
# Read the outcome with workflow/scripts/foam_log_state.sh <out>/log.solve, never from rc.
#
# NO `set -u` and NO `set -e`: the OpenFOAM bashrc reads unbound variables and would abort
# the script there, leaving an EMPTY log that reads exactly like a solver that never ran.
SRC=${1:?meshed source case}; OUT=${2:?out dir}; BOUND=${3:?valueBound}
LMODE=${4:-unity}; INADM=${5:-cellOnly}; NP=${6:-4}; SOLVER=${7:-leiaSemiLagrangeLevelSetFoam}
SFARGS=${8:-}
SRC=$(cd "$SRC" && pwd)
rm -rf "$OUT"; mkdir -p "$OUT"
for f in 0.org constant system; do cp -r "$SRC/$f" "$OUT/"; done
for f in "$SRC"/*.fms "$SRC"/*.stl; do [ -e "$f" ] && cp "$f" "$OUT/"; done
cp "$SRC/case_params.json" "$OUT/" 2>/dev/null
cd "$OUT" || exit 1
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"
foamDictionary -entry levelSet/semiLagrangian/valueBound     -set "$BOUND" system/fvSolution >/dev/null 2>&1
foamDictionary -entry levelSet/semiLagrangian/lipschitzMode  -set "$LMODE" system/fvSolution >/dev/null 2>&1
foamDictionary -entry levelSet/semiLagrangian/onInadmissible -set "$INADM" system/fvSolution >/dev/null 2>&1
foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict >/dev/null 2>&1
# 0/ from 0.org plus the pre-processing, never from a finished 0/ (CLAUDE.md).
rm -rf 0 processor* && cp -r 0.org 0 && find 0 -name '*.template' -delete
leiaSetFields $SFARGS > log.setFields 2>&1 || { echo "ABORT setFields"; tail -5 log.setFields; exit 3; }
if [ "$NP" -gt 1 ]; then
  decomposePar -force > log.decomposePar 2>&1 || { echo "ABORT decomposePar"; exit 4; }
  mpirun -np "$NP" "$SOLVER" -parallel > log.solve 2>&1
else
  "$SOLVER" > log.solve 2>&1
fi
# md5 of the mesh points: the ONLY thing that makes a cross-arm comparison valid.
md5sum constant/polyMesh/points | cut -c1-12 | tr -d '\n'; echo "  <- mesh digest"
