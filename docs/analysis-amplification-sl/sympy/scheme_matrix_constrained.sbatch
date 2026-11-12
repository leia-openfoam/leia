#!/bin/bash
# The +c variant of the scheme matrix (hanging points as ghost values) on the four meshes
# with a 2:1 transition, one array task per mesh, all at the same time.
# Submit from the clone root:
#   sbatch docs/analysis-amplification-sl/sympy/scheme_matrix_constrained.sbatch
# Needs a Python with numpy in $LEIA_AMP_PY (default: python3 on the PATH). No OpenFOAM.
#SBATCH --job-name=leia-schemes-c
#SBATCH --account=special00004
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --output=studies/schemeMatrix/%x-%A_%a.out
#SBATCH --error=studies/schemeMatrix/%x-%A_%a.err
MESHES=(quad21 quad21patch hex21 hex21patch)
m=${MESHES[$SLURM_ARRAY_TASK_ID]}
out=studies/schemeMatrix
mkdir -p $out
${LEIA_AMP_PY:-python3} -u docs/analysis-amplification-sl/sympy/scheme_matrix.py --mesh $m \
    --variant constrained --jobs ${SLURM_CPUS_PER_TASK:-4} \
    --csv $out/scheme_matrix_${m}_constrained.csv > $out/scheme_matrix_${m}_constrained.txt
