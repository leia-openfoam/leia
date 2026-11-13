#!/bin/bash
#SBATCH -J leia-int-poly
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem-per-cpu=4000
#SBATCH -t 08:00:00
#SBATCH -o int.poly.%j.out
#SBATCH -e int.poly.%j.err
# One polyhedral arm on a GIVEN mesh with THIS clone's binaries (advect_bound_arm.sh).
# Submit from the clone root:
#   sbatch --export=ALL,INT_SCRIPTS=<dir>,SRC=<meshed case>,OUT=<out dir>,BOUND=<valueBound> poly_arm.sbatch
# The last action is <OUT>.done. Read the outcome with foam_log_state.sh <OUT>/log.solve.
cd "$SLURM_SUBMIT_DIR" || exit 1
. "${INT_SCRIPTS:?}/env.sh"
mkdir -p "$(dirname "$OUT")"; rm -f "$OUT.done"
bash workflow/scripts/advect_bound_arm.sh "$SRC" "$OUT" "$BOUND" unity cellOnly 48 \
    leiaSemiLagrangeLevelSetFoam > "$OUT.arm.log" 2>&1
echo "arm rc=$?" >> "$OUT.arm.log"
echo "state=$(bash workflow/scripts/foam_log_state.sh $OUT/log.solve)" >> "$OUT.arm.log"
echo ARM_FINISHED > "$OUT.done"
