#!/bin/bash
#SBATCH -J davof-hp-build
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The profiling build (Claude, 2026-11-04): a detached git worktree of the fork at
# davof-swept HEAD (720cfc6) in prof/TwoPhaseFlow, built with the default
# WM_COMPILE_OPTION=Opt (c++OPT = -O3) plus FOAM_EXTRA_CXXFLAGS="-g
# -fno-omit-frame-pointer" (wmake/rules/linux64Gcc/c++ appends FOAM_EXTRA_CXXFLAGS
# to c++FLAGS), installed in the worktree's own platforms/ (etc/davof-env.sh).
# The default build in ../TwoPhaseFlow is not touched.
W=/work/scratch/tm83tomy/davof/agents/swept
P=$W/prof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
mkdir -p $P
if [ ! -d $P/TwoPhaseFlow ]; then git -C $W/TwoPhaseFlow worktree add --detach $P/TwoPhaseFlow 720cfc6 || exit 1; fi
cd $P/TwoPhaseFlow || exit 1
echo "worktree $(git log --oneline -1) on $(hostname), $(date)"
export FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer"
export WM_NCOMPPROCS=${SLURM_CPUS_PER_TASK:-24}
echo "WM_COMPILER=$WM_COMPILER WM_COMPILE_OPTION=$WM_COMPILE_OPTION WM_OPTIONS=$WM_OPTIONS WM_COMPILE_CONTROL='$WM_COMPILE_CONTROL' FOAM_EXTRA_CXXFLAGS='$FOAM_EXTRA_CXXFLAGS' gcc $(g++ --version | head -1)"
wmakeLnInclude -u src/VoF > /dev/null 2>&1
./Allwmake > build.log 2>&1; rc=$?
echo "--- one compile line (isoAdvectionSwept.C):"
grep -m1 'isoAdvectionSwept.C' build.log | tr ' ' '\n' | grep -E '^-(O|g|f|std|m|D)' | tr '\n' ' '; echo
echo "--- errors: $(grep -c -i ' error' build.log)"
. ./etc/davof-env.sh
echo "advectorVoF: $(command -v advectorVoF)"; ls -la platforms/*/lib/libVoF.so platforms/*/bin/advectorVoF
echo "debug_line sections in libVoF.so: $(readelf -S platforms/*/lib/libVoF.so | grep -c debug_line)"
strings platforms/*/lib/libVoF.so | grep -q sweptFaceCells || { echo "libVoF lacks sweptFaceCells"; rc=1; }
echo "HP-BUILD-DONE rc=$rc $(date)"
exit $rc
