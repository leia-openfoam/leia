#!/bin/bash
# The scheme matrix on every mesh, one array task per mesh, all at the same time.
# Submit from the clone root:  sbatch docs/analysis-amplification-sl/sympy/scheme_matrix.sbatch
# Needs a Python with numpy in $LEIA_AMP_PY (default: python3 on the PATH). No OpenFOAM.
#SBATCH --job-name=leia-schemes
#SBATCH --account=special00005
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --output=studies/schemeMatrix/%x-%A_%a.out
#SBATCH --error=studies/schemeMatrix/%x-%A_%a.err
MESHES=(square quadPert quad21 quad21patch cube hexPert hex21 hex21patch honeycomb honeycombPert)
m=${MESHES[$SLURM_ARRAY_TASK_ID]}
out=studies/schemeMatrix
mkdir -p $out
${LEIA_AMP_PY:-python3} -u docs/analysis-amplification-sl/sympy/scheme_matrix.py --mesh $m \
    --jobs ${SLURM_CPUS_PER_TASK:-4} --csv $out/scheme_matrix_$m.csv > $out/scheme_matrix_$m.txt
