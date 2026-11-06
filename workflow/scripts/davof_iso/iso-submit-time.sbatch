#!/bin/bash
#SBATCH -J davof-iso-submit
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=2000
#SBATCH -t 00:20:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Submits the timed verification jobs (iso-time.sbatch) for the cells whose preparation
# study is complete (15 cases with advectorVoF.csv); run after the preparation orchestrators
# (--dependency=afterany). FLOWS, NS, COS select the cells; N 128 cells are split per
# repetition (REP0=REP1=r) plus one job for the profiled runs. Claude, 2026-11-06.
# When this script runs as a job, --export=ALL would carry this job's SLURM_TRES_PER_TASK
# (cpu=4 on the submit allocation) into the timed jobs, whose own allocation sets
# SLURM_CPUS_PER_TASK=1: srun then refuses the step ("cpus-per-task set by two different
# environment variables"; the shear batch of 2026-11-06 12:14 lost every row that way).
# The timed jobs get a clean environment.
for v in $(env | grep -o '^SLURM_[^=]*' ; env | grep -o '^SBATCH_[^=]*'); do unset $v; done
W=/work/scratch/tm83tomy/davof/agents/swept
I=$W/iso; S=$I/leia/studies/davof; SB=$I/leia/workflow/scripts/davof_iso/iso-time.sbatch
FLOWS=${FLOWS:-"Rotation Shear Deformation"}; NS=${NS:-"32 64 128"}; COS=${COS:-"0 1 2"}
# MESHCLASSES (default hex): the mesh classes; a class other than hex has N 32 and 64 (OF Co
# 0.2, 0.5 and 1 like hex: the standing protocol of 2026-11-07), and its job names carry the class.
MESHCLASSES=${MESHCLASSES:-hex}   # tet: dropped by decision (Tomislav, 2026-11-08), do not submit it
declare -A NP=([32]=1 [64]=8 [128]=64)
# the 2D flows (FLOW *2D, the study iso<FLOW>_...): N 64 / 128 / 256 on 1 / 1 / 8 ranks
declare -A NP2D=([64]=1 [128]=1 [256]=8)
for MC in $MESHCLASSES; do for F in $FLOWS; do for N in $NS; do
  [ $MC != hex ] && [ $N -eq 128 ] && continue
  STB=iso${F}3D; R=${NP[$N]}; case $F in *2D) STB=iso$F; R=${NP2D[$N]};; esac
  st=$S/${STB}_${MC}N$N; n=$(ls $st/*_0*/advectorVoF.csv 2>/dev/null | wc -l); nc=$(ls -d $st/*_0* 2>/dev/null | wc -l)
  if [ "$n" -lt 1 ] || [ "$n" -lt "$nc" ]; then echo "SKIP $F $MC N$N: $n of $nc prepared"; continue; fi
  for CO in $COS; do
    name=davof-iso-$F-N$N-co$CO; [ $MC != hex ] && name=davof-iso-$F-$MC-N$N-co$CO
    if [ $N -eq 128 ] && [ "$R" -ge 64 ]; then
      for r in 1 2 3; do
        J=$(sbatch --parsable -J $name-r$r -n $R -t 03:00:00 --export=ALL,FLOW=$F,MESHCLASS=$MC,N=$N,CO=$CO,REP0=$r,REP1=$r,NOPROF=1 $SB 2>/dev/null | tail -1)
        echo "$J $name-r$r (the timed verification, repetition $r)" >> $W/.my_jobs; echo "submitted $J $name-r$r"
      done
      J=$(sbatch --parsable -J $name-prof -n $R -t 03:00:00 --export=ALL,FLOW=$F,MESHCLASS=$MC,N=$N,CO=$CO,REP0=prof,REP1=prof $SB 2>/dev/null | tail -1)
      echo "$J $name-prof (the timed verification, the profiled runs)" >> $W/.my_jobs; echo "submitted $J $name-prof"
    else
      lim=01:30:00; [ $N -eq 32 ] && lim=00:30:00
      J=$(sbatch --parsable -J $name -n $R -t $lim --export=ALL,FLOW=$F,MESHCLASS=$MC,N=$N,CO=$CO $SB 2>/dev/null | tail -1)
      echo "$J $name (the timed verification: 5 rows x 3 reps + profiled runs)" >> $W/.my_jobs; echo "submitted $J $name"
    fi
  done
done; done; done
echo "ISO-SUBMIT-DONE $(date)"
