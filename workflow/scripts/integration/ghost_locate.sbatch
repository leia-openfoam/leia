#!/bin/bash
#SBATCH -J leia-int-locate
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=6000
#SBATCH -t 03:00:00
#SBATCH -o int.locate.%j.out
#SBATCH -e int.locate.%j.err
# Locate the decomposition dependence of the ghost-cell reconstruction: one prepared case
# (mesh and 0/ from the workflow, no solver run), copied twice, solved in serial and on 4
# ranks with psi written at every step, then compared cell by cell. The pre-registered
# read-out is in config/integration/ghostDecompLocateN64.yaml.
#   sbatch --export=ALL,INT_SCRIPTS=<dir>,SRC=<prepared case>,OUT=<output dir>,DECOMP=<decomposeParDict> ghost_locate.sbatch
# DECOMP is the decomposeParDict of the 4-rank case. The workflow renders the prepared case in
# serial mode, and its decomposeParDict has one subdomain. Give the file of a 4-rank study of
# the same case (for example studies/ghostDecompN64np4/<case>/system/decomposeParDict).
# The last action is int.locate.done; the tables are in <OUT>/locate.md.
cd "$SLURM_SUBMIT_DIR" || exit 1
. "${INT_SCRIPTS:?}/env.sh"
rm -f int.locate.done
SOLVER=leiaSemiLagrangeLevelSetFoam
SRC=$(cd "${SRC:?}" && pwd); mkdir -p "${OUT:?}"; OUT=$(cd "$OUT" && pwd)
{
echo "src=$SRC commit=$(git rev-parse HEAD) date=$(date -Is) host=$(hostname)"
[ -f "$SRC/log.$SOLVER" ] && { echo "ABORT: a solver ran in the source case"; exit 2; }
for a in serial np4; do
    rm -rf "$OUT/$a"; mkdir -p "$OUT/$a"
    for f in 0 0.org constant system; do cp -r "$SRC/$f" "$OUT/$a/"; done
    cp "$SRC/case_params.json" "$OUT/$a/"
    ( cd "$OUT/$a" || exit 1
      rm -rf processor* [1-9]* 0.[0-9]*
      # psi at every step, 17 digits, time names that separate the steps
      foamDictionary -entry writeControl   -set timeStep system/controlDict > /dev/null
      foamDictionary -entry writeInterval  -set 1        system/controlDict > /dev/null
      foamDictionary -entry writeFormat    -set ascii    system/controlDict > /dev/null
      foamDictionary -entry writePrecision -set 17       system/controlDict > /dev/null
      foamDictionary -entry timePrecision  -set 14       system/controlDict > /dev/null
      # Optional: SLSET="<keyword> <value>" sets one entry of levelSet/semiLagrangian in
      # fvSolution, for an option that the case template does not render.
      if [ -n "$SLSET" ]; then
          set -- $SLSET
          foamDictionary -entry "levelSet/semiLagrangian/$1" -set "$2" system/fvSolution > /dev/null
          echo "$a: $1 = $(foamDictionary -entry "levelSet/semiLagrangian/$1" -value system/fvSolution)" >&2
      fi )
done
echo "initial state equal: $(diff -rq "$OUT/serial/0" "$OUT/np4/0" > /dev/null && diff -rq "$OUT/serial/constant" "$OUT/np4/constant" > /dev/null && diff -rq "$OUT/serial/system" "$OUT/np4/system" > /dev/null && echo yes || echo NO)"
echo "mesh digest: $(md5sum "$OUT/serial/constant/polyMesh/points" | cut -c1-12) $(md5sum "$OUT/np4/constant/polyMesh/points" | cut -c1-12)"
echo "psi digest:  $(md5sum "$OUT/serial/0/psi" | cut -c1-12) $(md5sum "$OUT/np4/0/psi" | cut -c1-12)"
echo "---- constant/polyMesh/boundary"; cat "$OUT/serial/constant/polyMesh/boundary"
cd "$OUT/serial" || exit 1
postProcess -func writeCellCentres -time 0 > log.writeCellCentres 2>&1; echo "writeCellCentres rc=$?"
$SOLVER > log.$SOLVER 2>&1; echo "serial solver rc=$?"
echo "serial state=$(bash "$SLURM_SUBMIT_DIR/workflow/scripts/foam_log_state.sh" log.$SOLVER)"
cd "$OUT/np4" || exit 1
cp "${DECOMP:?}" system/decomposeParDict; grep -v "^ *//" system/decomposeParDict | grep "numberOfSubdomains\|method\|^ *n "
decomposePar -force -cellDist > log.decomposePar 2>&1; echo "decomposePar rc=$?"
mpirun -np 4 $SOLVER -parallel > log.$SOLVER 2>&1; echo "np4 solver rc=$?"
echo "np4 state=$(bash "$SLURM_SUBMIT_DIR/workflow/scripts/foam_log_state.sh" log.$SOLVER)"
reconstructPar -fields '(psi)' > log.reconstructPar 2>&1; echo "reconstructPar rc=$?"
echo "---- boundary-face and bound lines"
grep -h "stencilBoundaryFaces\|lipschitzCone" "$OUT/serial/log.$SOLVER" "$OUT/np4/log.$SOLVER" | sort | uniq -c | cut -c1-200
echo "---- ghost data lines (serial, then np 4)"
grep -h "ghostCellQuadraticWeightedLeastSquares:" "$OUT/serial/log.$SOLVER" | sort | uniq -c
grep -h "ghostCellQuadraticWeightedLeastSquares:" "$OUT/np4/log.$SOLVER" | sort | uniq -c
PY=python3; $PY -c "import numpy" 2>/dev/null || PY=$(head -1 "$HOME/.local/bin/snakemake" | sed 's/^#!//')
echo "python: $PY"
$PY "$INT_SCRIPTS/ghost_locate_compare.py" "$OUT/serial" "$OUT/np4" --out "$OUT/locate" > "$OUT/locate.md" 2>&1; echo "compare rc=$?"
} > "$OUT/locate.job.log" 2>&1
# The two solved copies hold psi at every step (about 12,500 files). Pack them into one
# tar and remove them; locate.md, locate_steps.csv and locate.job.log stay next to it.
( cd "$OUT" && tar -cf fields.tar serial np4 && [ "$(tar -tf fields.tar | wc -l)" = "$(find serial np4 | wc -l)" ] \
  && rm -rf serial np4 && echo "packed serial and np4 into fields.tar" >> locate.job.log )
echo LOCATE_FINISHED > "$SLURM_SUBMIT_DIR/int.locate.done"
