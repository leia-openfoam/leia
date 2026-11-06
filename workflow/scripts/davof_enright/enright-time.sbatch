#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The timed 3D deformation (Enright / LeVeque) comparison of the swept-face advection with
# plicRDF + isoAdvector (2026-11-06; the pre-registration E-1 .. E-5 in
# config/davof/sphereDeformation3D_hexSweptN32.yaml). One exclusive i01 node per (N, CO):
# REPS repetitions of the two methods (ring = the swept-face advection with the ring normal
# in the study setting, rdf = plicRDF + isoAdvector), the order alternated per repetition,
# each solve one srun step on NP ranks (--cpu-bind=cores); the case copies come from the
# preparation studies sphereDeformation3D_hex{Swept,RDF}N<N> (rendered, meshed, initialised,
# decomposed), endTime set to T = 3 (3 N 2.5 / Co steps of the fixed dt, a multiple of four),
# the fields and the diagnostics written at T/4, T/2, 3T/4 and T.
#   sbatch -J davof-enright-N<N>-co<CO> -n <NP> -t <limit> --export=ALL,N=..,CO=0|1 enright-time.sbatch
# REP0 / REP1 (optional): only the repetitions REP0..REP1 (one job per repetition for the
# long cells); their files carry _r<REP0>. METHODS (optional): only these methods.
W=/work/scratch/tm83tomy/davof/agents/enright
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/leia && . etc/leia-env.sh || exit 1
: ${N:?} ${CO:?}
REPS=${REPS:-3}
NP=$SLURM_NTASKS
END=3; WI=0.75
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-90); leia $(git -C $W/leia log --oneline -1 | cut -c1-90)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP, N $N, CO index $CO"
echo "advectorVoF $(command -v advectorVoF)"
T=$W/time/N${N}_co${CO}
if [ -z "$REP0" ]; then rm -rf $T
elif [ -n "$METHODS" ]; then for m in $METHODS; do rm -rf $T/${m}_r$REP0; done
else rm -rf $T/*_r$REP0; fi
mkdir -p $T
SUF=${REP0:+_r$REP0}${METHODS:+_${METHODS// /-}}
pick() { /usr/bin/python3 - "$S/sphereDeformation3D_hex$1N$N" "$2" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/sphereDeformation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort()
print(cs[co][1])
EOF
}
declare -A SRC
SRC[ring]=$(pick Swept davof-planeSurf-ringArea)
SRC[rdf]=$(pick RDF plicRDF)
for m in ring rdf; do echo "source $m ${SRC[$m]}"; done
prep() { # method rep -> case dir
  local m=$1 r=$2 src=${SRC[$1]} c=$T/${1}_r$2
  mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
  else
    cp -r $src/0 $c/
  fi
  foamDictionary -entry endTime -set $END $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  echo "$src" > $c/.source
  echo $c
}
csv=$T/timing$SUF.csv
echo "N,np,co_index,method,rep,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
orders=("ring rdf" "rdf ring" "ring rdf")
for r in $(seq ${REP0:-1} ${REP1:-$REPS}); do
  for m in ${orders[$(( (r-1) % 3 ))]}; do
    if [ -n "$METHODS" ] && [[ " $METHODS " != *" $m "* ]]; then continue; fi
    c=$(prep $m $r)
    [ $r -eq ${REP0:-1} ] && echo "$m: $(grep -E '^(deltaT|endTime|writeInterval)' $c/system/controlDict | tr -s ' ' | tr '\n' ' ')"
    par=""; [ $NP -gt 1 ] && par="-parallel"
    t0=$(date +%s.%N)
    (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${m}_r$r advectorVoF $par > log.advectorVoF 2>&1); rc=$?
    t1=$(date +%s.%N)
    (cd $c && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1)
    /usr/bin/python3 - $c $m $r $rc $t0 $t1 >> $csv <<'EOF'
import sys, re, os, statistics
c, m, r, rc, t0, t1 = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log))
st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
e = os.environ
print(",".join(str(x) for x in (e["N"], e["SLURM_NTASKS"], e["CO"], m, r, steps, dt, et,
      ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0), "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
    echo "$m r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s, steps $(grep -c '^Time = ' $c/log.advectorVoF), $(grep -c 'ABORT' $c/log.advectorVoF) abort lines, interface cells at the end: $(grep -o 'davof: [0-9]* interface cells\|plicRDF.*interface' $c/log.advectorVoF | tail -1 | cut -c1-40)"
  done
done
sleep 5
sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps$SUF.txt
cat $T/sacct_steps$SUF.txt
echo ENRIGHT-TIME-DONE
