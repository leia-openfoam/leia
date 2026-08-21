#!/bin/bash
#SBATCH -J leia-curv
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem-per-cpu=4000
#SBATCH -t 12:00:00
#SBATCH -o leia-sink.%j.out
#SBATCH -e leia-sink.%j.err
#
# Sink-measurement driver: restart existing stationary-droplet arms with sigma = 0
# and measure the per-step decay of max|U|. Rationale and method are documented in
# workflow/scripts/sink_decay_restart.sh; the short version is that every
# stability exponent we have quoted is a local slope on a NET gain
# A = -sink + amplifier whose sink has never been measured, and A changes sign
# under refinement so those slopes are inflated by cancellation.
#
# This is deliberately NOT a snakemake study. A study is one (case, mesh, mode)
# built from templates at t = 0; this is a restart from another study's written
# field, which the DAG has no way to express. The driver runs the solver directly
# with mpirun inside its own allocation, sequentially over the arms.
#
#   sbatch -n 8  -t 04:00:00 --export=ALL,SRC=studies/upwindConvection2D,\
#          DST=studies/sinkDecay2D,NSTEPS=3000,SIGMA=0,GLOB='stationaryDroplet2D_0000[0-3]' \
#          run-sink-decay.sbatch
#
#   sbatch -n 32 -t 12:00:00 --export=ALL,SRC=studies/filterOffAmplifier3D,\
#          DST=studies/sinkDecay3D,NSTEPS=3000,SIGMA=0,GLOB='stationaryDroplet3D_0000[0-3]' \
#          run-sink-decay.sbatch
#
# Re-submitting is safe and resumes: the driver skips any destination arm whose
# log already contains OpenFOAM's terminating `End`, so a job killed at its wall
# limit continues at the next arm rather than redoing finished ones. Completeness
# is judged on the LOG, never on a post-processing output that only appears after
# a run ends -- that mistake destroyed a live 16 h interFoam arm.
set -o pipefail   # NOTE: no `set -u` -- OpenFOAM etc/bashrc reads WM_PROJECT_DIR unset

REPO="${REPO:-/work/scratch/tm83tomy/leia-curvature}"
cd "$REPO" || exit 1

module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/curvature-v2512/etc/bashrc"
export PATH="$HOME/.local/bin:$PATH"

echo "===== BEGIN sink-decay  $(date) on $(hostname) ====="
echo "  REPO=$REPO  SRC=$SRC  DST=$DST  NSTEPS=$NSTEPS  SIGMA=$SIGMA  GLOB=$GLOB"
echo "  ntasks=$SLURM_NTASKS  solver=$(which leiaSemiLagrangianLevelSetTwoPhaseFoam)"

bash workflow/scripts/sink_decay_restart.sh \
     "$SRC" "$DST" "$NSTEPS" "$SIGMA" "$GLOB" "${NP:-$SLURM_NTASKS}"
rc=$?

echo "===== END sink-decay rc=$rc  $(date) ====="
echo "SINK-DECAY-DONE"
