#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The rerun of the reference setting (Claude, 2026-11-04; fork davof-swept 407f969, the main
# clone's default build): the rendered study cases copied WITHOUT the keys of the reference
# setting in alpha.water (strip), so that the fork's defaults act; scored with the converter.
#   MODE=plane : planeTranslation3D_hexSweptWispM9Ccs2, N 20..160, Co 0.2 / 0.95, both arms and rows (32 cases, PAR at a time)
#   MODE=cmp   : FLOW=Rot|LineX|LineDiag N=.. CO=0|1 (the comparison's cells): ring REPS repetitions, exact one (the control), -n NP
#   MODE=spin  : sphereRotation3D_hexSpinCoPyr (N 32 / 64, Co 0.2 / 0.5), PAR at a time
#   MODE=par   : sphereRotation3D_hexSweptParN32np{1,2,4,8,16}, 20 steps, both rows, the decomposition comparison against np 1
#   sbatch -J davof-ref-<tag> -n NP -t LIMIT --export=ALL,MODE=..,FLOW=..,N=..,CO=.. ref-run.sbatch
W=/work/scratch/tm83tomy/davof/agents/swept
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/TwoPhaseFlow && . etc/davof-env.sh || exit 1
: ${MODE:?}
NP=$SLURM_NTASKS; REPS=${REPS:-3}; PAR=${PAR:-8}
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-70); leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), job $SLURM_JOB_ID, MODE $MODE np $NP $(date)"
echo "advectorVoF $(command -v advectorVoF)"
KEYS='sweptWispSurface|sweptWispMatchTol|sweptWispTol|sweptConsistent|sweptVolumeCorrection|sweptScaleToPhi|positioning|positioningTol|positioningMaxIter|anchorSkipConverged|sweptPreTest|sweptStateShortcut|sweptTraceCache'
strip() { sed -i -E "/^\s*($KEYS)\s/d" $1/system/fvSolution; local n=$(grep -c -E "^\s*($KEYS)\s" $1/system/fvSolution); [ $n -eq 0 ] || echo "STRIP FAILED in $1: $n keys left"; }
copy() { # src dst (the fields of 0 or processor*/0, constant, system; no keys)
  mkdir -p $2; cp -r $1/system $1/constant $2/
  if [ -d $1/processor0 ]; then for p in $1/processor*; do mkdir -p $2/$(basename $p); cp -r $p/0 $p/constant $2/$(basename $p)/; done; else cp -r $1/0 $2/; fi
  strip $2; echo "$1" > $2/.source
}
metrics() { (cd $1 && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1); }
tok() { /usr/bin/python3 -c "import json,sys; t=json.load(open(sys.argv[1]+'/case_params.json'))['tokens']; print(' '.join(str(t.get(k,'')) for k in sys.argv[2:]))" "$@"; }
case $MODE in
plane)
  # PTAG: a suffix of the output directory (a split submission); NSET: only these N (default all four)
  T=$W/ref/plane${PTAG:-}; rm -rf $T; mkdir -p $T; NSET=${NSET:-"20 40 80 160"}
  run1() { local src=$1 c=$T/$(basename $1); copy $src $c; (cd $c && advectorVoF > log.advectorVoF 2>&1; echo "$(basename $c) rc=$? steps $(grep -c '^Time = ' log.advectorVoF) $(grep -o 'ClockTime = [^ ]* s' log.advectorVoF | tail -1)"); metrics $c; }
  export -f run1 copy strip metrics; export T W KEYS
  for c in $S/planeTranslation3D_hexSweptWispM9Ccs2/planeTranslation3D_0*; do
    [ -f $c/case_params.json ] || continue
    read N co row arm <<< "$(tok $c N_CELLS TARGET_ADVECTIVE_CO TPF_ROW TPF_VELOCITY_ARM)"
    N=${N%.*}; [[ " $NSET " == *" $N "* ]] || continue
    case $co in 0.2|0.95) echo $c;; esac
  done | xargs -P $PAR -I{} bash -c 'run1 {}'
  echo "plane: $(ls -d $T/planeTranslation3D_0* | wc -l) cases"
  ;;
spin)
  T=$W/ref/spin; rm -rf $T; mkdir -p $T
  run1() { local src=$1 c=$T/$(basename $1); copy $src $c; (cd $c && advectorVoF > log.advectorVoF 2>&1; echo "$(basename $c) rc=$? steps $(grep -c '^Time = ' log.advectorVoF) $(grep -o 'timeStepTime = [^ ]* s' log.advectorVoF | tail -1)"); metrics $c; }
  export -f run1 copy strip metrics; export T W KEYS
  ls -d $S/sphereRotation3D_hexSpinCoPyr/sphereRotation3D_0* | xargs -P $PAR -I{} bash -c 'run1 {}'
  ;;
cmp)
  : ${FLOW:?} ${N:?} ${CO:?}
  case $FLOW in Rot) END=1;; *) END=0.45;; esac
  T=$W/ref/cmp/${FLOW}_N${N}_co${CO}; rm -rf $T; mkdir -p $T
  pick() { /usr/bin/python3 - "$S/sphereRotation3D_hexCmp${FLOW}SweptN$N" "$1" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/sphereRotation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[co][1])
EOF
  }
  declare -A SRC; SRC[ring]=$(pick davof-planeSurf-ringArea); SRC[exact]=$(pick davof-planeSurf-exact)
  for m in ring exact; do echo "source $m ${SRC[$m]}"; done
  prep() { local m=$1 r=$2 c=$T/${1}_r$2; copy ${SRC[$m]} $c
    local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict); local n=$(python3 -c "print(max(1, round($END/$dt)))")
    foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
    foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
    foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null; echo $c; }
  csv=$T/timing.csv; echo "flow,N,np,co_index,method,rep,corr,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
  par=""; [ $NP -gt 1 ] && par="-parallel"
  for r in $(seq 1 $REPS); do for m in ring exact; do
    [ $m = exact ] && [ $r -gt 1 ] && continue
    c=$(prep $m $r); t0=$(date +%s.%N)
    (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${m}_r$r advectorVoF $par > log.advectorVoF 2>&1); rc=$?
    t1=$(date +%s.%N); metrics $c
    /usr/bin/python3 - $c $m $r $rc $t0 $t1 "$FLOW" "$N" "$NP" "$CO" >> $csv <<'EOF'
import sys, re, os, statistics
c, m, r, rc, t0, t1, flow, N, np_, co = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10]
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log)); st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
print(",".join(str(x) for x in (flow, N, np_, co, m, r, "pyramid", steps, dt, et, ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0),
      "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
    echo "$m r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s; $(grep -o 'faces by path: [^;]*' $c/log.advectorVoF | tail -1)"
  done; done
  sleep 5; sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps.txt; cat $T/sacct_steps.txt
  ;;
par)
  T=$W/ref/par; rm -rf $T; mkdir -p $T
  for k in 1 2 4 8 16; do for i in 0 1; do
    src=$S/sphereRotation3D_hexSweptParN32np$k/sphereRotation3D_0000$i; c=$T/np${k}_$i; copy $src $c
    par=""; [ $k -gt 1 ] && par="-parallel"
    (cd $c && srun --ntasks=$k --cpu-bind=cores --kill-on-bad-exit=1 advectorVoF $par > log.advectorVoF 2>&1; echo "np$k row $i rc=$? steps $(grep -c '^Time = ' log.advectorVoF) rank fallbacks $(grep -o 'rank fallbacks [0-9]*' log.advectorVoF | awk '{s+=$3} END {print s+0}')")
    [ $k -gt 1 ] && (cd $c && reconstructPar -withZero > log.reconstructPar 2>&1)
  done; done
  for i in 0 1; do for k in 2 4 8 16; do
    /usr/bin/python3 $W/leia/workflow/scripts/davof_decomposition_compare.py --ref $T/np1_$i --case $T/np${k}_$i --out $T/cmp_np${k}_$i.csv > $T/cmp_np${k}_$i.log 2>&1 || echo "compare np$k row $i failed"
    /usr/bin/python3 - $T/cmp_np${k}_$i.csv $k $i <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
g = lambda r, k: float(r[k]) if r.get(k) not in (None, "") else float("nan")
print("np%s row %s: steps %d; step 1 max|dalpha| %.2e max|dalphaf| %.2e; max over steps |dalpha| %.2e |dalphaf| %.2e; nMissing %d nDup %d" % (
    sys.argv[2], sys.argv[3], len(rows), g(rows[0], "maxAlphaDiff"), g(rows[0], "maxAlphafDiff"), max(g(r, "maxAlphaDiff") for r in rows), max(g(r, "maxAlphafDiff") for r in rows),
    sum(int(float(r.get("nMissing", 0) or 0)) for r in rows), sum(int(float(r.get("nDup", 0) or 0)) for r in rows)))
EOF
  done; done
  ;;
*) echo "unknown MODE $MODE"; exit 1;;
esac
echo "REF-RUN-DONE $MODE $(date)"
