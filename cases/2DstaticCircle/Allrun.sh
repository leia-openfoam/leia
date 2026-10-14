#!/bin/bash
# The static source test (STATUS 11.27), standalone, in a RENDERED case directory: the workflow
# renders the *.template files (workflow/scripts/materialize.py); this script then runs the steps
# of its mesh, preprocess and solve rules. Needs the OpenFOAM environment and etc/leia-env.sh
# (it exports LEIA_ROOT, the clone root).
cd "${0%/*}" || exit 1
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
: "${LEIA_ROOT:?source etc/leia-env.sh of the clone first}"

rm -rf 0 && cp -r 0.org 0
find 0 -name '*.template' -delete 2>/dev/null || true

runApplication blockMesh
# psi0: one variant of system/initFields (system/setExprFieldsDict).
runApplication setExprFields
# The reference psiEnd = d and alphaEnd = alpha_DA(d) (endReference translate at zero velocity).
python3 "$LEIA_ROOT/workflow/scripts/write_end_reference.py" > log.endReference.driver 2>&1 || exit 1
runApplication leiaLevelSetFoam
