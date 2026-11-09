#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
# The cost of the face-fraction corrector (leia davof-iter, 2026-11-09; REPORT-swept.md
# "The iterative face-fraction corrector", P4): the prepared cases of a study
# (rendered, meshed, initialised, decomposed by the workflow) copied per corrector setting
# and run on ONE exclusive i01 node at the study's rank count, END steps of the fixed dt,
# the function objects at the end only, REPS repetitions with the order rotated; then one
# profiled run (OpenFOAM's profiling) of none and mean1, not timed. The CSV: the wall
# clock per run and the median of the per-step ExecutionTime differences.
#   sbatch -J iter-time-N64 -n 8 -t 06:00:00 --export=ALL,STUDY=iter_rot_N64,CO=0.0462996198086,END=0.25 workflow/scripts/iter_time.sbatch
W=/work/scratch/tm83tomy/davof/agents/iter
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/leia && . etc/leia-env.sh || exit 1
: ${STUDY:?} ${CO:?}
REPS=${REPS:-2}
END=${END:-0.25}
NP=$SLURM_NTASKS
SETTINGS=${SETTINGS:-"none mean1 mean2 mean4"}
echo "leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP"
T=$W/time/${STUDY}_co${CO}
rm -rf $T; mkdir -p $T
pick() { /usr/bin/python3 - "$S/$STUDY" "$1" "$CO" <<'EOF'
import json, glob, os, sys
st, corr, co = sys.argv[1], sys.argv[2], sys.argv[3]
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c) or not os.path.isfile(c + "/case_params.json"): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_FACE_FRACTION_CORRECTOR"] == corr and t["TARGET_ADVECTIVE_CO"] == co:
        print(c); break
EOF
}
declare -A SRC
for s in $SETTINGS; do SRC[$s]=$(pick $s); echo "source $s ${SRC[$s]}"; [ -n "${SRC[$s]}" ] || exit 1; done
prep() { # setting rep [prof] -> case dir
  local s=$1 r=$2 prof=$3 src=${SRC[$1]} c=$T/${1}_r$2${3:+_prof}
  mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
  else
    cp -r $src/0 $c/
  fi
  local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
  local n=$(python3 -c "print(max(1, round($END/$dt)))")
  foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
  foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null
  for fo in areaFractionError volumeFractionError; do
    foamDictionary -entry functions/$fo/writeControl -set onEnd $c/system/controlDict > /dev/null 2>&1
    foamDictionary -entry functions/$fo/executeControl -set onEnd $c/system/controlDict > /dev/null 2>&1
  done
  # the corrector's per-step line stays (it is the diagnostic); the profiling for the
  # profiled run only
  if [ -n "$prof" ]; then foamDictionary -entry profiling/active -set true $c/system/controlDict > /dev/null 2>&1; fi
  echo $c
}
csv=$T/timing.csv
echo "study,co,np,setting,rep,profiled,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
run() { # setting rep prof
  local s=$1 r=$2 prof=$3
  local c=$(prep $s $r $prof)
  local par=""; [ $NP -gt 1 ] && par="-parallel"
  local t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${s}_r$r davofAdvectorFoam $par > log.davofAdvectorFoam 2>&1); local rc=$?
  local t1=$(date +%s.%N)
  /usr/bin/python3 - $c $s $r "${prof:-0}" $rc $t0 $t1 >> $csv <<'EOF'
import sys, re, os, statistics
c, s, r, prof, rc, t0, t1 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
log = open(c + "/log.davofAdvectorFoam").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log))
st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
e = os.environ
print(",".join(str(x) for x in (e["STUDY"], e["CO"], e["SLURM_NTASKS"], s, r, 1 if prof != "0" else 0, steps, dt, et,
      ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0), "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
  echo "$s r$r ${prof:+prof} rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s"
}
arr=($SETTINGS); n=${#arr[@]}
for r in $(seq 1 $REPS); do
  for k in $(seq 0 $((n-1))); do run ${arr[$(( (k + r - 1) % n ))]} $r; done
done
for s in none mean1; do [[ " $SETTINGS " == *" $s "* ]] && run $s 1 prof; done
# the profiling files of the profiled runs
for s in none mean1; do
  c=$T/${s}_r1_prof; [ -d $c ] || continue
  f=$(ls -d $c/processor0/*/uniform/profiling $c/*/uniform/profiling 2>/dev/null | tail -1)
  [ -n "$f" ] && cp $f $T/profiling_$s.txt && echo "profiling of $s: $T/profiling_$s.txt"
done
sleep 5
sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps.txt
echo ITER-TIME-DONE
