#!/usr/bin/env bash
#
# Unit gate for the halo-limited velocity extension
# (docs/plan-halo-limited-gradient-control.md, D4). Serial, then 4 ranks with
# the interface next to a processor boundary. Exits non-zero on any failure:
#
#     bash cases/haloLimitedUnit/Allrun.sh
#
# Requires a sourced OpenFOAM environment plus etc/leia-env.sh, and a built
# leiaTestHaloLimited (./Allwmake at the repo root).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

rm -rf processor* constant/polyMesh
blockMesh > log.blockMesh 2>&1
leiaTestHaloLimited > log.leiaTestHaloLimited.serial 2>&1 || { tail -30 log.leiaTestHaloLimited.serial; exit 1; }
grep 'passed' log.leiaTestHaloLimited.serial
decomposePar -force > log.decomposePar 2>&1
mpirun -np 4 leiaTestHaloLimited -parallel > log.leiaTestHaloLimited.np4 2>&1 || { tail -30 log.leiaTestHaloLimited.np4; exit 1; }
grep 'passed' log.leiaTestHaloLimited.np4
