#!/bin/bash
#SBATCH -J leia-hpc-prep
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=24000
#SBATCH -t 02:00:00
#SBATCH -o hpcprofile/%x.%j.out
#SBATCH -e hpcprofile/%x.%j.err
#
# Prepare the cases of the HPCToolkit profile (config/gates/ghostCellProfileNp{1,4}.yaml).
# The workflow runs the mesh rule, the preprocess rule and, for the parallel study, the
# decompose rule. It runs with the local profile inside this one job, so this job is the
# only job to record in .my_jobs. The solve rule does not run.
#
#   mkdir -p hpcprofile && sbatch workflow/hpctoolkit/prep.sbatch      (from the clone root)
set -o pipefail
cd "${SLURM_SUBMIT_DIR:?submit with sbatch from the clone root}" || exit 1
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. "$PWD/etc/leia-env.sh" || exit 1
export PATH="$HOME/.local/bin:$PATH"
rm -f hpcprofile/prep.done
rc=0
for pair in ghostCellProfileNp1:preprocess ghostCellProfileNp4:decompose; do
    cfg=${pair%%:*}; rule=${pair##*:}
    snakemake --workflow-profile profiles/local --nolock --keep-going --cores 1 --jobs 1 \
        --rerun-triggers mtime --configfile config/gates/$cfg.yaml --until $rule \
        > hpcprofile/prep.$cfg.log 2>&1
    r=$?; echo "snakemake rc=$r" >> hpcprofile/prep.$cfg.log; [ $r -ne 0 ] && rc=$r
done
echo "PREP_FINISHED rc=$rc $(date)" > hpcprofile/prep.done
