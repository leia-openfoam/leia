#!/usr/bin/env bash
#
# Analytic regression test for the SDPLS source term. Exits non-zero on any
# mismatch, so it is usable directly as a repo gate:
#
#     bash cases/sdplsSourceUnit/Allrun.sh
#
# Requires a sourced OpenFOAM environment and a built libleiaLevelSet +
# leiaTestSdplsSource (./Allwmake at the repo root).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

blockMesh > log.blockMesh 2>&1
leiaTestSdplsSource 2>&1 | tee log.leiaTestSdplsSource

# tee masks the exit status; recover it from the pipe.
exit "${PIPESTATUS[0]}"
