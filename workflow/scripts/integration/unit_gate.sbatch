#!/bin/bash
#SBATCH -J leia-int-unit
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=3000
#SBATCH -t 01:00:00
#SBATCH -o int.unit.%j.out
#SBATCH -e int.unit.%j.err
# The exact-solution gate leiaTestGhostCellReconstruction in serial and on 4 ranks, in every
# case directory of CASES (prepared, decomposed to 4 subdomains, no solver run).
#   sbatch --export=ALL,INT_SCRIPTS=<dir>,CASES="<case dir> ..." unit_gate.sbatch
# The last action is int.unit.done; the numbers are in int.unit.log.
cd "$SLURM_SUBMIT_DIR" || exit 1
. "${INT_SCRIPTS:?}/env.sh"
rm -f int.unit.done
for d in $CASES; do
  ( cd $d || exit 1
    echo "== $d"; grep -h "reconstruction " system/fvSolution | head -2
    echo "subdomains: $(ls -d processor* | wc -l)"
    leiaTestGhostCellReconstruction > log.ghostGate.serial 2>&1; echo "serial rc=$?"
    tail -25 log.ghostGate.serial
    mpirun -np 4 leiaTestGhostCellReconstruction -parallel > log.ghostGate.np4 2>&1; echo "np4 rc=$?"
    grep -m1 "nProcs" log.ghostGate.np4
    tail -25 log.ghostGate.np4 )
done > int.unit.log 2>&1
echo UNIT_FINISHED > int.unit.done
