#!/usr/bin/env bash
#
# Unit gate for the semi-Lagrangian gradient-control source
# (docs/plan-halo-limited-gradient-control.md, D3). Serial, then 4 ranks with
# the band across a processor seam. Exits non-zero on any failure:
#
#     bash cases/slSourceUnit/Allrun.sh
#
# Requires a sourced OpenFOAM environment plus etc/leia-env.sh, and a built
# leiaTestSlSource (./Allwmake at the repo root).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

rm -rf processor* constant/polyMesh
blockMesh > log.blockMesh 2>&1
leiaTestSlSource > log.leiaTestSlSource.serial 2>&1 || { tail -30 log.leiaTestSlSource.serial; exit 1; }
tail -1 log.leiaTestSlSource.serial
decomposePar -force > log.decomposePar 2>&1
mpirun -np 4 leiaTestSlSource -parallel > log.leiaTestSlSource.np4 2>&1 || { tail -30 log.leiaTestSlSource.np4; exit 1; }
tail -1 log.leiaTestSlSource.np4
