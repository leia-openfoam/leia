#!/bin/bash
#SBATCH -J d2-hp-build
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:30:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
# The profiling build of the pieces path (diag2, 2026-11-16): a detached worktree of this clone at HEAD in prof/leia,
# -O3 (the default WM_COMPILE_OPTION) plus FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer", its own platforms/.
# The production build in the clone is not touched.   cd diag2 && sbatch workflow-copy of this file
W=${SLURM_SUBMIT_DIR:?}
P=$W/prof
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
mkdir -p $P
if [ ! -d $P/leia ]; then git -C $W/leia worktree add --detach $P/leia $(git -C $W/leia rev-parse HEAD) || exit 1; fi
cd $P/leia || exit 1
echo "worktree $(git log --oneline -1) on $(hostname), $(date)"
export FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer"
export WM_NCOMPPROCS=${SLURM_CPUS_PER_TASK:-24}
echo "WM_COMPILE_OPTION=$WM_COMPILE_OPTION WM_OPTIONS=$WM_OPTIONS FOAM_EXTRA_CXXFLAGS='$FOAM_EXTRA_CXXFLAGS' gcc $(g++ --version | head -1)"
./Allwmake > build.log 2>&1; rc=$?
echo "--- one compile line (davofSweptAdvectionFlux.C):"
grep -m1 'davofSweptAdvectionFlux.C' build.log | tr ' ' '\n' | grep -E '^-(O|g|f|std|m|D)' | tr '\n' ' '; echo
echo "--- errors: $(grep -c -i ' error' build.log)"
. etc/leia-env.sh
echo "davofAdvectorFoam: $(command -v davofAdvectorFoam)"; ls -la platforms/*/lib/libleiaDavof.so platforms/*/bin/davofAdvectorFoam
echo "debug_line sections in libleiaDavof.so: $(readelf -S platforms/*/lib/libleiaDavof.so | grep -c debug_line)"
mkdir -p $P/struct-cache
[ -d /work/scratch/tm83tomy/davof/agents/swept/prof/struct-cache ] && cp -rn /work/scratch/tm83tomy/davof/agents/swept/prof/struct-cache/. $P/struct-cache/ 2>/dev/null
echo "HP-BUILD-DONE rc=$rc $(date)"
exit $rc
