#!/bin/bash
#SBATCH -J leia-int-build
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=3000
#SBATCH -t 02:00:00
#SBATCH -o int.build.%j.out
#SBATCH -e int.build.%j.err
# Build ONE clone. Submit from the clone root:
#   sbatch <merged clone>/workflow/scripts/integration/build.sbatch
# The job writes int.build.log and, as its last action, int.build.done.
cd "$SLURM_SUBMIT_DIR" || exit 1
rm -f int.build.done
. "${INT_SCRIPTS:?}/env.sh"
export WM_NCOMPPROCS=24
{
  echo "clone=$PWD commit=$(git rev-parse HEAD) date=$(date -Is) host=$(hostname)"
  ./Allwmake
  echo "Allwmake rc=$?"
  python3 etc/leia-check-deps.py
  echo "leia-check-deps rc=$?"
  ls -la --time-style=full-iso platforms/*/lib/ platforms/*/bin/
  echo "ghostCell strings count: $(strings platforms/*/lib/libleiaSemiLagrangian.so | grep -c ghostCellQuadraticWeightedLeastSquares)"
  echo "footAnchored strings count: $(strings platforms/*/lib/libleiaSemiLagrangian.so | grep -c footAnchoredQuadraticWeightedLeastSquares)"
  which pMesh leiaSemiLagrangeLevelSetFoam leiaLevelSetFoam leiaSetFields
} > int.build.log 2>&1
echo BUILD_FINISHED > int.build.done
