#!/bin/bash
# The psi boundary-condition probe of the one-way 2Dtranslation (STATUS 11.19, 2026-09-29).
#
# Copies a rendered 2Dtranslation arm (system, constant, 0.org), sets the psi boundary
# condition of the four patches to one variant, rebuilds 0/ exactly as the workflow does
# (0.org, leiaSetFields, write_end_reference.py) and runs the SL kinematic solver in SERIAL
# with a write every 0.05 (writeControl runTime: the time-step sequence does not change).
#
# Variants (left = inflow, right = outflow, top/bottom = no flow):
#   exactAll         the exact signed distance on all four patches
#   exactNoOutflow   exact on left, top, bottom; zeroGradient on right
#   exactInflowOnly  exact on left; zeroGradient on right, top, bottom
#   zeroGradAll      zeroGradient on all four patches (the committed case)
#
# Usage: translation_bc_probe.sh <rendered arm dir> <probe dir> <variant>
# Then:  translation_bc_probe_scan.py <probe dir parent> <N> <probe dir name> ...
# FOAM_BASHRC overrides the OpenFOAM environment (default: v2606 in $HOME/OpenFOAM; the
# measurements of STATUS 11.19 ran on v2512).
src=$1; dst=$2; var=$3
[ -n "$src" ] && [ -n "$dst" ] && [ -n "$var" ] || { sed -n '2,20p' "$0"; exit 1; }
REPO=$(cd "$(dirname "$0")/../.." && pwd)
source "${FOAM_BASHRC:-$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc}"
. "$REPO/etc/leia-env.sh"
src=$(cd "$src" && pwd)
rm -rf "${dst:?}"; mkdir -p "$dst"
cp -r "$src/system" "$src/constant" "$src/0.org" "$dst/"
cd "$dst" || exit 1
python3 - "$var" <<'PY'
import sys
var = sys.argv[1]
p = "0.org/psi"; txt = open(p).read()
head = txt.split("boundaryField", 1)[0]
ex = """{
        type            uniformFixedValue;
        value           uniform 0;
        uniformValue
        {
            type        expression;
            expression  #{ sqrt(sqr(pos().x() - 0.25 - arg()) + sqr(pos().y() - 0.5)) - 0.15 #};
        }
    }"""
zg = "{ type zeroGradient; }"
m = {"exactAll": dict(left=ex, right=ex, top=ex, bottom=ex),
     "zeroGradAll": dict(left=zg, right=zg, top=zg, bottom=zg),
     "exactInflowOnly": dict(left=ex, right=zg, top=zg, bottom=zg),
     "exactNoOutflow": dict(left=ex, right=zg, top=ex, bottom=ex)}[var]
open(p, "w").write(head + "boundaryField\n{\n" + "".join(f"    {k} {v}\n" for k, v in m.items()) + "}\n")
PY
rm -rf 0 && cp -r 0.org 0
leiaSetFields > log.leiaSetFields 2>&1
python3 "$REPO/workflow/scripts/write_end_reference.py" > log.endReference.driver 2>&1 || { echo "endReference failed"; exit 3; }
foamDictionary -entry writeInterval -set 0.05 system/controlDict > /dev/null 2>&1
leiaSemiLagrangeLevelSetFoam > log.leiaSemiLagrangeLevelSetFoam 2>&1
"$REPO/workflow/scripts/foam_log_state.sh" log.leiaSemiLagrangeLevelSetFoam
