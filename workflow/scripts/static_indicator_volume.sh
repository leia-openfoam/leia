#!/bin/bash
# Static (t=0) phase-indicator volume convergence to the EXACT circle area.
#
# Initialises the bulk shear-vortex circle (implicitSphere r=0.15, centre (0.5,0.75),
# pseudo-2D domain [0,1]^2 x [-0.05,0.05], z-thickness 0.1) with each phase indicator
# (geometric, detrixheAslam) at several resolutions and records the initialised phase
# volume VOL_ALPHA_0 = sum(alpha_c V_c) written by leiaSetFields. No advection.
#
# The convergence order is then measured (make_indicator_static_fig.py) against the
# EXACT cylinder volume pi r^2 * 0.1 (equivalently the circle area pi r^2 with the
# pseudo-2D height divided out).
#
# Usage: bash workflow/scripts/static_indicator_volume.sh [out.csv]
set +eu; source $HOME/OpenFOAM/OpenFOAM-v2506/etc/bashrc; set -e

REPO=/home/tmaric/OpenFOAM/repos/leia
SRC=$REPO/studies/phaseIndicatorConvergence/2Dvortex_00000   # 2Dvortex / shear2D base
OUT=${1:-$REPO/workflow/scripts/static_indicator_volume.csv}
WORK=/home/tmaric/qa_deck/staticvol
rm -rf "$WORK"; mkdir -p "$WORK"
echo "indicator,N,h,vol0" > "$OUT"

run() { # indicator N
  local IND=$1 N=$2
  local C=$WORK/${IND}_$N
  rm -rf "$C"; cp -r "$SRC" "$C"; cd "$C"
  rm -rf processor* postProcessing *.csv
  for d in *; do [ -d "$d" ] && [[ "$d" =~ ^[0-9]+(\.[0-9]+)?$ ]] && rm -rf "$d"; done
  cp -r 0.org 0
  # mesh resolution: rewrite the hex block cell counts to (N N 1)
  sed -i -E "s/\([0-9]+ [0-9]+ 1\) simpleGrading/($N $N 1) simpleGrading/" system/blockMeshDict
  # phase indicator selection
  foamDictionary system/fvSolution -entry "levelSet/phaseIndicator/type" -set "$IND" >/dev/null
  blockMesh > log.blockMesh 2>&1
  leiaSetFields > log.setFields 2>&1
  local H=$(tail -1 leiaSetFields.csv | cut -d, -f1)
  local V=$(tail -1 leiaSetFields.csv | cut -d, -f2)
  echo "$IND,$N,$H,$V" >> "$OUT"
  printf "  %-14s N=%-4s h=%-11s VOL_ALPHA_0=%s\n" "$IND" "$N" "$H" "$V"
}

for IND in geometric detrixheAslam; do
  for N in 32 64 128 256; do run "$IND" "$N"; done
done
echo "STATIC-VOL-DONE -> $OUT"
