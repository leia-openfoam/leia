#!/bin/bash
#SBATCH -J leia-int-drv
#SBATCH -A special00005
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=12000
#SBATCH -t 24:00:00
#SBATCH -o int.drv.%j.out
#SBATCH -e int.drv.%j.err
# One-core driver: runs the studies of CFGS (a space-separated list of config files) in this
# clone. Submit from the clone root:
#   sbatch --export=ALL,INT_SCRIPTS=<dir>,TAG=<tag>,CFGS="config/a.yaml config/b.yaml",PROFILE=profiles/slurm[,EXTRA="--until preprocess"] driver.sbatch
# The studies start 150 s apart (two snakemake starts in one worktree in the same second
# fail). Each study writes int.<tag>.<study>.log. The last action is int.<tag>.done.
cd "$SLURM_SUBMIT_DIR" || exit 1
. "${INT_SCRIPTS:?}/env.sh"
TAG=${TAG:?}; PROFILE=${PROFILE:-profiles/slurm}
rm -f int.$TAG.done
for c in $CFGS; do
    s=$(basename $c .yaml)
    ( make studies-one-file CFG=$c PROFILE=$PROFILE TASKS=10000 SNAKE_EXTRA="$EXTRA" > int.$TAG.$s.log 2>&1
      echo "make rc=$?" >> int.$TAG.$s.log ) &
    sleep 150
done
wait
echo DRIVER_FINISHED > int.$TAG.done
