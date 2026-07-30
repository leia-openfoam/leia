#!/usr/bin/env bash
# Run the 1D plane redistancing test for every redistancer model.
# The initial psi is a signed distance with |grad psi| = 2 (implicitPlane
# 'gradient 2'); a correct redistancer recovers the EXACT signed distance
# (planeFootWave: globally machine-exact; see leiaTestRedistance.csv).

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

set -o errexit

for REDIST in PDE planeFootWave anchoredEikonal; do
    echo ">>> redistancer: $REDIST"
    sed -i "s/^\( *type *\)\(PDE\|planeFootWave\|anchoredEikonal\|noRedistancing\);/\1$REDIST;/" \
        system/fvSolution
    # PDE 'write true' leaves ~500 pseudo-time dirs that poison the next
    # variant (controlDict picks up the latest time) -- start pristine.
    ./Allclean >/dev/null 2>&1 || true
    rm -rf 0 && cp -r 0.org 0
    touch "$(basename ${PWD}).foam"
    mkdir -p results.variants
    blockMesh >"results.variants/log.blockMesh.$REDIST" 2>&1
    leiaSetFields >"results.variants/log.leiaSetFields.$REDIST" 2>&1
    leiaTestRedistance >"results.variants/log.leiaTestRedistance.$REDIST" 2>&1
    [ -f leiaTestRedistance.csv ] && cp leiaTestRedistance.csv "results.variants/leiaTestRedistance.$REDIST.csv"
    echo "    $(tail -1 leiaTestRedistance.csv)"
    sed -n '21,45p' 1/psi > "results.variants/psi.redistanced.$REDIST"
done

# restore the committed default
sed -i "s/^\( *type *\)\(planeFootWave\|anchoredEikonal\|noRedistancing\);/\1PDE;/" system/fvSolution
echo "Done. Outputs in results.variants/: leiaTestRedistance.<model>.csv, psi.redistanced.<model>, logs."
