# Environment of every integration-identity job. Source it from the clone root.
# No `set -e` and no `set -u`: the OpenFOAM bashrc reads unbound variables.
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. "$PWD/etc/leia-env.sh" || exit 1
export PATH="$HOME/.local/bin:$PATH"
