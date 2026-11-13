#!/usr/bin/env bash
# Integration identity check: compare the studies of a merged clone with those of a base
# clone, case by case, at tolerance 0.
#
#   integration_identity_check.sh <merged clone> <base clone> [<ghost reference study dir>]
#
# The two clones have run the same configs (workflow/scripts/integration/README.md gives the
# sequence). The script compares, with integration_identity_check.py:
#   1. studies/advConv2Dtranslation, studies/advConv3DshearHex,
#      studies/pdeLevelSet/s0Translation2D, studies/pdeLevelSet/a0Identity2Dvortex:
#      the cases are paired by their tokens;
#   2. studies/intPoly: the polyhedral arms on shared meshes, paired by directory name;
#   3. studies/ghostCellRefinedBoxN64 of the merged clone against the cases 00002 and 00003
#      of the reference study (default: the feature-branch clone leia-ghost).
# It writes one table per study into <merged clone>/studies/integrationIdentity/ and prints
# them. Exit code: 0 when every CSV pair is identical, 2 for a difference, 1 for a
# structural problem. No `set -e`: every study is compared, also after a failure.
M=${1:?merged clone}; B=${2:?base clone}
G=${3:-/work/scratch/tm83tomy/leia-ghost/studies/ghostCellRefinedBox3D}
PY="python3 $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/integration_identity_check.py"
OUT=$M/studies/integrationIdentity; mkdir -p "$OUT"
worst=0
run() {  # run <table name> <args of the python script ...>
    local name=$1; shift
    $PY "$@" --out "$OUT/$name.csv" > "$OUT/$name.md" 2>&1
    local rc=$?
    cat "$OUT/$name.md"; echo "[identity] $name rc=$rc"; echo
    [ $rc -gt $worst ] && worst=$rc   # 2 (a difference) outranks 1 (a structural problem)
}
for s in advConv2Dtranslation advConv3DshearHex pdeLevelSet/s0Translation2D pdeLevelSet/a0Identity2Dvortex; do
    run "$(basename $s)" "$M/studies/$s" "$B/studies/$s" --label "$(basename $s)"
done
run advConv3DshearPoly "$M/studies/intPoly" "$B/studies/intPoly" --match name --label advConv3DshearPoly
run ghostCellRefinedBoxN64 "$M/studies/ghostCellRefinedBoxN64" "$G" --only-b 00002,00003 \
    --show SL_RECONSTRUCTION,N_CELLS --label ghostCellRefinedBoxN64
echo "[identity] overall rc=$worst (0 identical, 2 difference, 1 structural problem)"
exit $worst
