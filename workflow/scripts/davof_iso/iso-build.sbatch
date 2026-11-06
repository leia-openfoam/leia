#!/bin/bash
#SBATCH -J davof-iso-build
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=3600
#SBATCH -t 00:40:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The fork worktree iso/$WT (default TwoPhaseFlow; created from the clone if missing) at
# ISO_SHA (davof-isoAdvector) rebuilt (Allwmake, default flags), then the preparation
# studies of the verification submitted from iso/leia (PREP="study1 study2 ..."; empty:
# none) with LEIA_TPF_DIR = that worktree. Claude, 2026-11-06/07.
W=/work/scratch/tm83tomy/davof/agents/swept
I=$W/iso
: ${ISO_SHA:?}
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
echo "node $(hostname), job $SLURM_JOB_ID, ISO_SHA $ISO_SHA, $(date)"
WT=${WT:-TwoPhaseFlow}
git -C $W/TwoPhaseFlow fetch -q origin || echo "fetch failed"
if [ ! -d $I/$WT ]; then git -C $W/TwoPhaseFlow worktree add -q --detach $I/$WT $ISO_SHA || { echo "worktree $I/$WT failed"; exit 1; }; echo "worktree $I/$WT created"; fi
cd $I/$WT && git checkout -q -- . && git checkout -q --detach $ISO_SHA || { echo "checkout $ISO_SHA failed"; exit 1; }
export WM_NCOMPPROCS=24; unset FOAM_EXTRA_CXXFLAGS
find . -name 'Allwmake*' -o -name 'Allwclean*' | xargs chmod +x
t0=$(date +%s); ./Allwmake > build.log 2>&1; rc=$?
echo "build at $(git log --oneline -1 | cut -c1-80) rc=$rc in $(( $(date +%s) - t0 )) s, errors $(grep -c -E 'error:' build.log)"
[ $rc -eq 0 ] || { grep -B2 -A8 -E 'error:' build.log | grep -v '^g++' | head -40; exit 1; }
. etc/davof-env.sh; ls -la platforms/*/bin/advectorVoF platforms/*/bin/initShearFlow3D platforms/*/bin/initVortex2D 2>/dev/null | cut -c1-110
cd $I/leia && git pull -q --ff-only origin davof-swept && echo "leia $(git log --oneline -1 | cut -c1-80)"
for v in $(env | grep -o '^SLURM_[^=]*'); do unset $v; done   # the child jobs' srun must not inherit this job's TRES
for s in $PREP; do
  J=$(sbatch --parsable -J leia-davof-$s --export=ALL,STUDY=davof/$s,LEIA_TPF_DIR=$I/$WT run-studies.sbatch 2>/dev/null | tail -1)
  echo "$J leia-davof-$s (the preparation study $s of the isoAdvector verification)" >> $W/.my_jobs
  echo "submitted $J for $s"
done
echo "ISO-BUILD-DONE $(date)"
