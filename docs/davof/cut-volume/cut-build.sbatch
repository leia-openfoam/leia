#!/bin/bash
#SBATCH -J cut-build
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 01:30:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
# the minimal build of the cut workspace (branch davof-ccs-cut): levelSetImplicitSurfaces, leiaCore (for
# leiaPerturbMesh), leiaPerturbMesh, davofTestCutVolume; into this clone's own platforms/
cd "${SLURM_SUBMIT_DIR:?}/leia" || exit 1
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export WM_NCOMPPROCS=${SLURM_CPUS_PER_TASK:-8}
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
ROOT=$PWD
. etc/leia-env.sh || exit 1
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN" && echo "$WM_PROJECT_VERSION" > "$ROOT/platforms/$WM_OPTIONS/.openfoam-version"
echo "===== cut build $(date): $(git log --oneline -1) on $(git branch --show-current) ====="
wmake $ROOT/src/levelSetImplicitSurfaces || exit 1
wmakeLnInclude -u $ROOT/src/leiaLevelSet > /dev/null
etc/leia-stamp.sh src/leiaLevelSet leiaCore || exit 1
wmake $ROOT/src/leiaLevelSet || exit 1
wmake $ROOT/applications/utilities/preProcessing/leiaPerturbMesh || exit 1
wmake $ROOT/applications/test/davofTestCutVolume || exit 1
ls -la $FOAM_USER_APPBIN
echo "CUT-BUILD-DONE rc=$? $(date)"