#!/usr/bin/env bash
# Build one mesh and measure the amplification bound of the semi-Lagrangian fit on it.
#
# WHY. Lambda_c = |1 - sum_j g_j| + sum_j |g_j| is the Lebesgue constant of the fit at the
# departure displacement d, and it depends on the stencil GEOMETRY alone -- no time steps, no
# flow, no field. So the question "which meshers produce cells that can amplify the level
# set?" is answered by ONE mesh pass per mesher, not by a coupled run. Lambda = 1 means the
# update is a convex combination of the stencil values and cannot create a new extremum.
# See workflow/README.md, "The amplification bound of the fit".
#
# WHAT THIS ANSWERS. Every polyhedral conclusion in this campaign rests on ONE mesher: the
# workflow's `mesh: poly` is pMesh, and there is no cartesianMesh or snappyHexMesh path. This
# script removes that limitation for the diagnostic, so the defect can be attributed to
# polyhedra, to cfMesh's boundary layer, or to small one-sided cells in general.
#
# cartesianMesh is the sharpest single arm: it is cfMesh's HEXAHEDRAL mesher reading the SAME
# meshDict and the SAME surface file as pMesh, so it isolates hex-versus-polyhedral WITHIN one
# mesher and one boundary-layer treatment.
#
# The run is serial and one step. It writes slFitAmplification, slFitOrder and slFitPivot into
# 0/, which workflow/scripts/sl_fit_amplification_census.py then tabulates.
#
# Usage:
#   workflow/scripts/fit_amplification_probe.sh <src-case> <mesher> <out-dir> [displacement]
#
#     <src-case>     a RENDERED case (studies/<study>/<case>_00000), never a case template
#     <mesher>       blockMesh | pMesh | cartesianMesh | snappyHexMesh
#     <out-dir>      created fresh; the probe never touches the source study
#     [displacement] the probe vector, default (-6.39947740050255e-07 0 0) = -U dt of the SI
#                    Popinet case at N = 64, so the numbers are directly comparable to the
#                    recorded hex 1.0527 and pMesh 1.2608
#
# maxCellSize and surfaceFile must already be set in the source case's system/meshDict for the
# two cfMesh meshers. A hex study renders them as -1.0 and none, so pass a POLY study as the
# source when probing cartesianMesh or pMesh, or set them by hand first.
# NO `set -u` and NO `set -e`: the OpenFOAM bashrc reads unbound variables
# (WM_PROJECT_DIR at etc/bashrc:184) and would abort the script there, leaving an EMPTY log
# that reads exactly like a mesher that produced nothing. Errors are checked explicitly below.

SRC=${1:?source rendered case}
MESHER=${2:?mesher: blockMesh | pMesh | cartesianMesh | snappyHexMesh}
OUT=${3:?output directory}
DISP=${4:-"(-6.39947740050255e-07 0 0)"}

SRC=$(cd "$SRC" && pwd)
rm -rf "$OUT"; mkdir -p "$OUT"
for f in 0.org constant system; do cp -r "$SRC/$f" "$OUT/"; done
# cfMesh reads the surface from the case directory; snappyHexMesh from constant/triSurface.
for f in "$SRC"/*.fms "$SRC"/*.stl; do [ -e "$f" ] && cp "$f" "$OUT/"; done
rm -rf "$OUT/constant/polyMesh"

cd "$OUT" || exit 1
# NOTE: sourcing the bashrc AFTER any WM_PROJECT_USER_DIR export would reset it; on a shared
# account source $HOME/.leia_env after the bashrc (CLUSTER.md).
# shellcheck disable=SC1090
LEIA_PMESH_PRELOAD=${LEIA_PMESH_PRELOAD:-}
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"

DT=$(foamDictionary -entry deltaT -value system/controlDict)

echo "probe: mesher=$MESHER  d=$DISP  dt=$DT"

case "$MESHER" in
  blockMesh)
      blockMesh > log.mesh 2>&1 ;;
  pMesh|cartesianMesh)
      mcs=$(foamDictionary -entry maxCellSize -value system/meshDict)
      sf=$(foamDictionary -entry surfaceFile -value system/meshDict)
      echo "  cfMesh: maxCellSize=$mcs surfaceFile=$sf"
      case "$mcs" in -*|"") echo "ABORT: maxCellSize is $mcs -- pass a POLY study as the source"; exit 2 ;; esac
      # The mesher's own environment workaround stays scoped to the mesher: a study-global
      # LD_PRELOAD once segfaulted the MPI solver (CLAUDE.md, provenance section).
      if [ -n "$LEIA_PMESH_PRELOAD" ]; then
          env LD_PRELOAD="$LEIA_PMESH_PRELOAD" "$MESHER" > log.mesh 2>&1
      else
          "$MESHER" > log.mesh 2>&1
      fi ;;
  snappyHexMesh)
      blockMesh > log.blockMesh 2>&1 || { echo "ABORT: background blockMesh failed"; exit 2; }
      mkdir -p constant/triSurface
      for f in *.stl; do [ -e "$f" ] && cp "$f" constant/triSurface/; done
      snappyHexMesh -overwrite > log.mesh 2>&1 ;;
  *)  echo "ABORT: unknown mesher $MESHER"; exit 2 ;;
esac
rc=$?
[ $rc -eq 0 ] || { echo "ABORT: $MESHER rc=$rc"; tail -20 log.mesh; exit $rc; }

# 0/ is regenerated from 0.org plus the pre-processing, NEVER copied from a finished 0/:
# the two-phase solver writes its projected pressure and recomputed phase indicator back into
# time 0, so a finished 0/ is not the initial state (CLAUDE.md).
rm -rf 0 && cp -r 0.org 0 && rm -f 0/*.template
leiaSetFields -alphaName alpha.water > log.setFields 2>&1 \
  || { echo "ABORT: leiaSetFields failed"; tail -20 log.setFields; exit 3; }

foamDictionary -entry levelSet/semiLagrangian/writeFitOrder -set true system/fvSolution >/dev/null
foamDictionary -entry levelSet/semiLagrangian/fitProbeDisplacement -set "$DISP" system/fvSolution >/dev/null
foamDictionary -entry endTime -set "$DT" system/controlDict >/dev/null
foamDictionary -entry writeInterval -set "$DT" system/controlDict >/dev/null

leiaSemiLagrangianLevelSetTwoPhaseFoam > log.probe 2>&1
src=$?
# Verify the ARTEFACT, never the exit code: the diagnostic must have written its field.
grep -m1 'amplification bound' log.probe || { echo "ABORT: no amplification line (rc=$src)"; tail -20 log.probe; exit 4; }
grep -m1 'ill-conditioned quadratic stencils' log.probe
postProcess -func writeCellVolumes -time 0 > log.writeCellVolumes 2>&1
ls 0/ | grep -c slFit | sed 's/^/  slFit fields written: /'
echo "probe done: $OUT"
