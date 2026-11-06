#!/bin/bash
#SBATCH -J davof-iso-pvc
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=3600
#SBATCH -t 00:40:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# planeVolumeConsistency on the written states of the DAVOF rows of the timed verification
# (F-3, the residual of vofCutCell's planes against ccs's on the same states and normals;
# the plicRDF rows write no alphaf, so the app (MUST_READ alphaf) runs on the DAVOF states).
# Normals: areaNormal (the DAVOF identity) and field (the ring normal the run used);
# positioning vofCutCell (positionTol 1e-8, plicRDF's) and ccs (1e-14 of V_c); the latest written
# time (the initial time holds no interfaceNormal). Claude, 2026-11-06.
W=/work/scratch/tm83tomy/davof/agents/swept
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $W/iso/TwoPhaseFlow/etc/davof-env.sh
export OMPI_MCA_btl=^openib
echo "node $(hostname), job $SLURM_JOB_ID, $(which planeVolumeConsistency), $(date)"
CELLS=${CELLS:-"Rotation_N32_co0/RING_r1 Rotation_N64_co0/RING_r1 Rotation_N32_co0/RINGvcc_r1 Deformation_N32_co0/RING_r1"}
for c in $CELLS; do
  d=$W/iso/time/$c; [ -d $d ] || { echo "MISSING $d"; continue; }
  par=""; np=1; [ -d $d/processor0 ] && { np=$(ls -d $d/processor* | wc -l); par="-parallel"; }
  for pos in vofCutCell ccs; do
    out=$d/postProcessing/planeVolumeConsistency/pvc_$pos.dat
    (cd $d && srun --ntasks=$np --cpu-bind=cores --kill-on-bad-exit=1 planeVolumeConsistency $par -latestTime -normals '(areaNormal field)' -positioning $pos -out $out > log.pvc_$pos 2>&1); rc=$?
    echo "== $c $pos (np $np) rc=$rc"; grep -v '^#' $out 2>/dev/null | cut -c1-400; grep -E 'HELD|FALSIFIED|E1|E2|E3|C1|C2' log.pvc_$pos 2>/dev/null | head -8 | cut -c1-200
  done
done
echo "ISO-PVC-DONE $(date)"
