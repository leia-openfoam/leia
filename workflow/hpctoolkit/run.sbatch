#!/bin/bash
#SBATCH -J leia-hpc-run
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH --mem-per-cpu=12000
#SBATCH -t 04:00:00
#SBATCH -o hpcprofile/%x.%j.out
#SBATCH -e hpcprofile/%x.%j.err
#
# One solver run of the HPCToolkit profile, in a COPY of a prepared and never solved case.
#
#   sbatch -J leia-hpc-<run> -n <np> [-c <cpus>] \
#       --export=ALL,SRC=<prepared case>,RUN=<run name>,NP=<np>,HPC=<0|1> \
#       workflow/hpctoolkit/run.sbatch
#
# SRC  the prepared case, relative to the clone root (studies/ghostCellProfileNp1/...)
# RUN  the name of the run; the copy is studies/ghostCellProfileRuns/<RUN>
# NP   1: serial, the solver command of the workflow's solve rule for mode serial
#      >1: the launcher of profiles/slurm and the option -parallel
# HPC  1: the solver runs under hpcrun (-e REALTIME@5000, one sample per 5 ms of wall clock
#         per thread), then hpcstruct and hpcprof run in this job; 0: the plain solver
#
# The marker <copy>/.hpc_run_done is the last file that the job writes.
set -o pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from the clone root}" || exit 1
REPO=$PWD
: "${SRC:?}" "${RUN:?}" "${NP:?}" "${HPC:?}"
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. "$PWD/etc/leia-env.sh" || exit 1
export PATH="$HOME/.local/bin:$PATH"
HPCTK=$HOME/hpctoolkit/bin
SOLVER=leiaSemiLagrangeLevelSetFoam
OUT=studies/ghostCellProfileRuns/$RUN

[ -e "$OUT" ] && { echo "refuse: $OUT exists" >&2; exit 1; }
[ -f "$SRC/.leia_preprocessed" ] || { echo "refuse: $SRC is not prepared" >&2; exit 1; }
if [ -e "$SRC/$SOLVER.csv" ] || [ -e "$SRC/log.$SOLVER" ]; then
    echo "refuse: $SRC was solved; its 0/ is not the initial state" >&2; exit 1
fi
mkdir -p studies/ghostCellProfileRuns
cp -r "$SRC" "$OUT" || exit 1
cd "$OUT" || exit 1

PAR=""; LAUNCH=""; END=0/psiEnd
if [ "$NP" -gt 1 ]; then
    PAR="-parallel"
    LAUNCH="srun --ntasks=$NP --overlap --cpu-bind=none --kill-on-bad-exit=1"
    END=processor0/0/psiEnd
fi
[ -f "$END" ] || { echo "refuse: $END is missing" >&2; exit 1; }
PRE=""
[ "$HPC" = 1 ] && PRE="$HPCTK/hpcrun -e REALTIME@5000 -o $PWD/hpctoolkit-measurements"

{
echo "solver=$SOLVER"; echo "np=$NP"; echo "hpc=$HPC"; echo "slurmJobId=$SLURM_JOB_ID"
echo "host=$(hostname)"; echo "cpus=$(nproc)"; echo "commit=$(git -C "$REPO" rev-parse HEAD)"
echo "command=$LAUNCH $PRE $SOLVER $PAR"
echo "start=$(date +%s.%N)"
} > .leia_launch
if [ "$NP" -gt 1 ]; then
    $LAUNCH $PRE $SOLVER $PAR > log.$SOLVER 2>&1
else
    /usr/bin/time -v -o time.$SOLVER $PRE $SOLVER > log.$SOLVER 2>&1
fi
rc=$?
{
echo "end=$(date +%s.%N)"; echo "rc=$rc"
echo "state=$(bash "$REPO/workflow/scripts/foam_log_state.sh" log.$SOLVER 2>/dev/null)"
} >> .leia_launch

if [ "$HPC" = 1 ]; then
    J=$(nproc)
    echo "hpcstructStart=$(date +%s)" >> .leia_launch
    $HPCTK/hpcstruct -j "$J" hpctoolkit-measurements > log.hpcstruct 2>&1
    echo "hpcstructRc=$?" >> .leia_launch
    $HPCTK/hpcprof -j "$J" -o hpctoolkit-database hpctoolkit-measurements > log.hpcprof 2>&1
    echo "hpcprofRc=$?" >> .leia_launch
    echo "hpcprofEnd=$(date +%s)" >> .leia_launch
fi
sleep 20   # sacct needs the step record
sacct -j "$SLURM_JOB_ID" --format=JobID%24,JobName%32,Elapsed,MaxRSS,AveRSS,NCPUS,State -P > sacct.txt 2>&1
echo "HPC_RUN_DONE $(date)" > .hpc_run_done
