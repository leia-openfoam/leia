#!/usr/bin/env bash
# Frozen-velocity power iteration of the semi-Lagrangian transport operator, one mesh.
#
# WHY. Lambda_c bounds ONE step in ONE cell. It does not decide repeated stability, and the
# anchored quadratic at x_i - C h on a uniform 1D grid is Lax-Wendroff: Lambda = 1 + C - C^2 > 1
# while |G|^2 = 1 - 4C^2(1-C^2)sin^4(theta/2) <= 1, i.e. L2-stable with negative weights. So the
# correlation between Lambda_max = 1.2608 on pMesh (which fails) and 1.05 on hexahedra (which do
# not) is correlation, not a mechanism. leiaTestTransportSpectrum measures the mechanism: with
# the velocity FROZEN and the clip OFF the transport stage is exactly linear, so its spectral
# radius and its finite-time power norms are well defined and measurable.
#
# This script builds ONE coarse mesh, regenerates 0/, decomposes and runs the power iteration.
# It is deliberately small: it is meant for the laptop at <= 8 ranks.
#
# Usage:
#   workflow/scripts/transport_spectrum_probe.sh <src-case> <mesher> <out-dir> <coarsen> <dt> [np] [seed]
#
#     <src-case>  a RENDERED case (studies/<study>/<case>_00000)
#     <mesher>    blockMesh | pMesh | cartesianMesh | snappyHexMesh | reuse
#     <coarsen>   blockMesh: the value for n_cells (n_cells_x is scaled by the same factor)
#                 cfMesh:    the value for maxCellSize
#                 reuse:     ignored
#     <dt>        deltaT, which fixes the departure displacement |d| = U dt. Choose it so that
#                 |d|/h MATCHES across the meshes being compared, or the comparison measures
#                 the Courant number instead of the stencil geometry.
#     [np]        MPI ranks, default 8. The mask and the stencil are decomposition independent
#                 by construction, so np is a cost choice -- but the seed is keyed on the cell
#                 CENTRE precisely so that it is identical at any np, which makes np a control.
#     [seed]      random (default) | checkerboard
#
# NO `set -u` and NO `set -e`: the OpenFOAM bashrc reads unbound variables and would abort the
# script there with an empty log. Errors are checked explicitly.

SRC=${1:?source rendered case}
MESHER=${2:?mesher}
OUT=${3:?output directory}
COARSEN=${4:?coarsen value (n_cells for blockMesh, maxCellSize for cfMesh)}
DT=${5:?deltaT}
NP=${6:-8}
SEED=${7:-random}

SRC=$(cd "$SRC" && pwd)
rm -rf "$OUT"; mkdir -p "$OUT"
for f in 0.org constant system; do cp -r "$SRC/$f" "$OUT/"; done
for f in "$SRC"/*.fms "$SRC"/*.stl; do [ -e "$f" ] && cp "$f" "$OUT/"; done
[ "$MESHER" = "reuse" ] || rm -rf "$OUT/constant/polyMesh"

cd "$OUT" || exit 1
LEIA_PMESH_PRELOAD=${LEIA_PMESH_PRELOAD:-}
# shellcheck disable=SC1090
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"

foamDictionary -entry deltaT -set "$DT" system/controlDict >/dev/null
foamDictionary -entry writeInterval -set "$DT" system/controlDict >/dev/null
foamDictionary -entry endTime -set "$DT" system/controlDict >/dev/null
# The power iteration needs a LINEAR operator; the clip is nonlinear. The app refuses to
# report a spectral radius otherwise, so force the entry off here and say so.
foamDictionary -entry levelSet/semiLagrangian/clipToStencilBounds -set false system/fvSolution >/dev/null
foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict >/dev/null 2>&1

case "$MESHER" in
  blockMesh)
      # n_cells_x keeps the box aspect ratio of the rendered dict.
      old=$(grep -oP '^n_cells\s+\K[0-9]+' system/blockMeshDict | head -1)
      oldx=$(grep -oP '^n_cells_x\s+\K[0-9]+' system/blockMeshDict | head -1)
      if [ -n "$old" ] && [ -n "$oldx" ]; then
          newx=$(( COARSEN * oldx / old ))
          sed -i -e "s/^n_cells  *[0-9]*;/n_cells   $COARSEN;/" \
                 -e "s/^n_cells_x  *[0-9]*;/n_cells_x $newx;/" system/blockMeshDict
          echo "  blockMesh: n_cells $old -> $COARSEN, n_cells_x $oldx -> $newx"
      fi
      blockMesh > log.mesh 2>&1 ;;
  pMesh|cartesianMesh)
      foamDictionary -entry maxCellSize -set "$COARSEN" system/meshDict >/dev/null
      echo "  cfMesh: maxCellSize $COARSEN, surfaceFile $(foamDictionary -entry surfaceFile -value system/meshDict)"
      if [ -n "$LEIA_PMESH_PRELOAD" ]; then
          env LD_PRELOAD="$LEIA_PMESH_PRELOAD" "$MESHER" > log.mesh 2>&1
      else
          "$MESHER" > log.mesh 2>&1
      fi ;;
  snappyHexMesh)
      blockMesh > log.blockMesh 2>&1
      mkdir -p constant/triSurface
      for f in *.stl; do [ -e "$f" ] && cp "$f" constant/triSurface/; done
      snappyHexMesh -overwrite > log.mesh 2>&1 ;;
  reuse) echo "  reusing the source mesh" ;;
  *) echo "ABORT: unknown mesher $MESHER"; exit 2 ;;
esac
rc=$?
[ $rc -eq 0 ] || { echo "ABORT: $MESHER rc=$rc"; tail -15 log.mesh; exit $rc; }

# 0/ from 0.org plus the pre-processing, never from a finished 0/ (CLAUDE.md).
rm -rf 0 && cp -r 0.org 0 && rm -f 0/*.template
leiaSetFields -alphaName alpha.water > log.setFields 2>&1 \
  || { echo "ABORT: leiaSetFields"; tail -15 log.setFields; exit 3; }

NCELLS=$(checkMesh -time 0 2>/dev/null | grep -m1 'cells:' | tr -dc '0-9')
echo "  cells = $NCELLS   deltaT = $DT   np = $NP   seed = $SEED"

if [ "$NP" -gt 1 ]; then
    decomposePar -force > log.decomposePar 2>&1 \
      || { echo "ABORT: decomposePar"; tail -15 log.decomposePar; exit 4; }
    mpirun -np "$NP" leiaTestTransportSpectrum -parallel -seed "$SEED" \
        "${@:8}" > log.spectrum 2>&1
else
    leiaTestTransportSpectrum -seed "$SEED" "${@:8}" > log.spectrum 2>&1
fi
src=$?

# Verify the ARTEFACT: the RESULT block must be present. Exit 1 from the app means
# rho > 1 + tol, which is a RESULT, not a failure.
grep -q 'RESULT' log.spectrum || { echo "ABORT: no RESULT block (rc=$src)"; tail -25 log.spectrum; exit 5; }
sed -n '/^leiaTestTransportSpectrum$/,/^  nIter/p' log.spectrum
grep -A7 '^RESULT' log.spectrum
grep -E '^(GROWS|BOUNDED)' log.spectrum
echo "probe done: $OUT (app rc=$src)"
