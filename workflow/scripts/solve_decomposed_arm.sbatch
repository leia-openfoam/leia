#!/bin/bash
#SBATCH -A special00004
#SBATCH -J leia-advsolve
#SBATCH --mem-per-cpu=4000
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#
# Run the SOLVER on ONE already-meshed and already-decomposed case.
#
# WHY THIS EXISTS. The snakemake SLURM executor completed mesh, preprocess and
# decompose for the 3D advection studies of 2026-09-10 and then ran the `solve` rule
# LOCALLY instead of submitting it, so it died on the login node with no solver log at
# all -- while the 2D studies of the same batch got their rule_solve SLURM jobs
# normally. Rather than leave prepared cases unused, this submits the remaining step
# directly. The case is untouched apart from the run itself.
#
# --mem-per-cpu is mandatory on Lichtenberg. -n comes from the submit line, and it
# MUST equal numberOfSubdomains in the case's decomposeParDict.
#
# Usage:
#   sbatch -n <np> -t <time> --export=ALL,CASE=<case-dir>,SOLVER=<app> \
#          workflow/scripts/solve_decomposed_arm.sbatch
#
# Read the outcome with workflow/scripts/foam_log_state.sh <case>/log.<SOLVER>,
# never from the return code: a deadlock leaves the job alive with a truncated log.

set +u
: "${CASE:?CASE must be an already-decomposed case directory}"
SOLVER=${SOLVER:-leiaSemiLagrangeLevelSetFoam}

source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"
cd "$CASE" || exit 1

NP=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict 2>/dev/null)
echo "case   = $CASE"
echo "prefix = $WM_PROJECT_USER_DIR"
echo "np     = $SLURM_NTASKS (decomposeParDict says $NP)"
if [ -n "$NP" ] && [ "$NP" != "$SLURM_NTASKS" ]; then
    echo "ABORT: -n $SLURM_NTASKS does not match numberOfSubdomains $NP"; exit 2
fi
ls -d processor0 >/dev/null 2>&1 || { echo "ABORT: case is not decomposed"; exit 3; }

# --cpu-bind=none is NOT cosmetic: without it every rank is confined to one CPU and
# OpenMPI's busy-wait engine collapses (measured 11000x slower, no error message).
# The tell is ExecutionTime/ClockTime: ~1 means the ranks are on separate cores.
srun --ntasks="$SLURM_NTASKS" --overlap --cpu-bind=none "$SOLVER" -parallel \
     > "log.$SOLVER" 2>&1
rc=$?
echo "solver rc=$rc"
bash "$(dirname "$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")")/workflow/scripts/foam_log_state.sh" "log.$SOLVER" 2>/dev/null \
  || tail -5 "log.$SOLVER"
