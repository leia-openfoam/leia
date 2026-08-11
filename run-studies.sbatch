#!/bin/bash
#SBATCH -J leia-studies
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=4000
#SBATCH -t 7-00:00:00
#SBATCH -o leia-studies.%j.out
#SBATCH -e leia-studies.%j.err
#
# Orchestrator for a study GROUP on Lichtenberg.
#
# Why a batch job and not just `make ... &` on a login node: a login-node run of
# `make studies-euler` was SIGTERMed after ~5.5 h, mid-sweep (login nodes reap
# long-lived user processes). This job only orchestrates -- one core, little
# memory -- while snakemake submits the real work as separate jobs.
#
#   sbatch run-studies.sbatch                              # TARGET=studies-euler
#   sbatch --export=ALL,TARGET=studies-sdpls run-studies.sbatch
#   sbatch --export=ALL,STUDY=sdplsConv3Dshear run-studies.sbatch   # one study
#
# Resuming is free: snakemake skips cases whose outputs already exist, so
# re-submitting after any interruption continues where it stopped.
set -o pipefail   # NOTE: no `set -u` -- OpenFOAM etc/bashrc reads WM_PROJECT_DIR unset

cd /work/scratch/tm83tomy/leia || exit 1
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc"
export PATH="$HOME/.local/bin:$PATH"

# STUDY runs exactly one config; TARGET runs a whole group.
if [ -n "$STUDY" ]; then
    echo "===== BEGIN study $STUDY  $(date) on $(hostname) ====="
    make studies-one STUDY="$STUDY" PROFILE=profiles/slurm
    echo "===== END study $STUDY rc=$?  $(date) ====="
else
    TARGET="${TARGET:-studies-euler}"
    echo "===== BEGIN $TARGET  $(date) on $(hostname) ====="
    make "$TARGET" PROFILE=profiles/slurm
    echo "===== END $TARGET rc=$?  $(date) ====="
fi
echo "ALL-STUDIES-DONE"
