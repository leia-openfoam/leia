#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The HPCToolkit profiles of the reference setting (Claude, 2026-11-04): hp-run.sbatch for the
# profiling build prof/TwoPhaseFlow at 407f969, the cases of the first profile (the comparison's
# rotation, OpenFOAM Co 1) WITHOUT the keys of the reference setting in alpha.water (the fork's
# defaults act), the same hpcrun events, hpcstruct and hpcprof; the databases under prof/runs2/.
#   sbatch -J davof-hp2-<MODE> -n NP -t LIMIT --export=ALL,MODE=serial|par8|par64 hp-run2.sbatch
W=/work/scratch/tm83tomy/davof/agents/swept
P=$W/prof
S=$W/leia/studies/davof
HP=/home/tm83tomy/spack/opt/spack/linux-x86_64/hpctoolkit-2025.0.1-ky3i4lswv3ynbyvk6ygnunxtn6tsfa67/bin
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $P/TwoPhaseFlow && . etc/davof-env.sh || exit 1
: ${MODE:?}
NP=$SLURM_NTASKS
# RUNDIR (optional): the runs directory under prof (default runs2); ADDKEYS (optional): key=value+key=value
# inserted into the alpha.water dictionary after reconstructionScheme (the opt-in keys of a round)
T=$P/${RUNDIR:-runs2}/$MODE
if [ -n "$ANALYSE" ]; then
  # the analysis alone (the measurement job ran with MEASURE_ONLY)
  for c in $T/*/; do
    M=$c/hpctoolkit-measurements; [ -d $M ] || continue
    ls $M | grep -c '\.hpcrun$' | xargs echo "$(basename $c): hpcrun files"
    $HP/hpcstruct -j $(nproc) -c $P/struct-cache $M > $c/log.hpcstruct 2>&1 || echo "hpcstruct failed for $c"
    $HP/hpcprof -j $(nproc) -o $c/hpctoolkit-database $M > $c/log.hpcprof 2>&1 || echo "hpcprof failed for $c"
    ls $c/hpctoolkit-database 2>/dev/null | tr '\n' ' '; echo
  done
  echo "HP-RUN-DONE $MODE (analysis) $(date)"; exit 0
fi
rm -rf $T; mkdir -p $T
echo "fork worktree $(git -C $P/TwoPhaseFlow log --oneline -1 | cut -c1-70); leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP"
echo "advectorVoF $(command -v advectorVoF)"; ldd $(command -v advectorVoF) | grep libVoF
echo "$($HP/hpcrun --version | head -1); perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid)"
case $MODE in
  serial) RUNS="ring32:full:nohp ring32:full ring64:20:nohp ring64:20 ring64:full";;
  par8)   RUNS="ring64:20:nohp ring64:full";;
  par64)  RUNS="ring128:20:nohp ring128:full";;
  *) echo "unknown MODE $MODE"; exit 1;;
esac
EV="-e CPUTIME@1000"
if $HP/hpcrun -e ix86arch::LLC_MISSES@100000 -o $T/.evtest /bin/true > $T/.evtest.log 2>&1; then EV="$EV -e ix86arch::LLC_MISSES@100000"; else echo "LLC_MISSES not available: $(tail -1 $T/.evtest.log)"; fi
echo "hpcrun events: $EV"
KEYS='sweptWispSurface|sweptWispMatchTol|sweptWispTol|sweptConsistent|sweptVolumeCorrection|sweptScaleToPhi|positioning|positioningTol|positioningMaxIter|anchorSkipConverged|sweptPreTest|sweptStateShortcut|sweptTraceCache'
pick() { /usr/bin/python3 - "$S/sphereRotation3D_hexCmpRot$1N$2" "$3" <<'EOF'
import json, glob, os, sys
st, row = sys.argv[1], sys.argv[2]
cs = []
for c in sorted(glob.glob(st + "/sphereRotation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[1][1])
EOF
}
prep() { # method N steps|full tag -> case dir
  local m=$1 N=$2 steps=$3 tag=$4 src c=$T/$4
  src=$(pick Swept $N davof-planeSurf-ringArea)
  mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
    [ $(ls -d $src/processor* | wc -l) -eq $NP ] || { echo "decomposition of $src is not $NP"; exit 1; }
  else
    cp -r $src/0 $c/
  fi
  local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
  local n; if [ $steps = full ]; then n=$(python3 -c "print(max(1, round(1.0/$dt)))"); else n=$steps; fi
  foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
  foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null
  sed -i -E "/^\s*($KEYS)\s/d" $c/system/fvSolution
  if [ -n "$ADDKEYS" ]; then for kv in ${ADDKEYS//+/ }; do
    sed -i "s/^\(\s*\)reconstructionScheme davof;/&\n\1${kv%%=*} ${kv#*=};/" $c/system/fvSolution
  done; fi
  echo "$src" > $c/.source; echo $c
}
csv=$T/timing.csv
echo "mode,np,method,N,steps,hpcrun,wall_s,exec_s,clock_s,n_steps,mean_step_s,median_step_s,first20_mean_s,node" > $csv
par=""; [ $NP -gt 1 ] && par="-parallel"
for spec in $RUNS; do
  IFS=: read mN steps hp <<< "$spec"
  m=${mN%%[0-9]*}; N=${mN##*[a-z]}
  tag=${m}${N}_${steps}${hp:+_$hp}
  c=$(prep $m $N $steps $tag)
  echo "=== $tag from $(cat $c/.source) $(date); keys left: $(grep -c -E "^\s*($KEYS)\s" $c/system/fvSolution)"
  t0=$(date +%s.%N)
  if [ "$hp" = nohp ]; then
    (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=$tag advectorVoF $par > log.advectorVoF 2>&1); rc=$?
  else
    (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=$tag $HP/hpcrun $EV -o $c/hpctoolkit-measurements advectorVoF $par > log.advectorVoF 2>&1); rc=$?
  fi
  t1=$(date +%s.%N)
  /usr/bin/python3 - "$c" "$m" "$N" "$steps" "${hp:-hp}" "$t0" "$t1" "$MODE" "$NP" >> $csv <<'EOF'
import sys, re, os, statistics
c, m, N, steps, hp, t0, t1, mode, np_ = sys.argv[1:10]
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
st = [float(x[2]) for x in ct if x[2]]
n = len(re.findall(r"\nTime = ", log))
print(",".join(str(x) for x in (mode, np_, m, N, steps, "0" if hp == "nohp" else "1", "%.2f" % (float(t1) - float(t0)),
      ct[-1][0] if ct else "", ct[-1][1] if ct else "", n, "%.4f" % statistics.mean(st) if st else "", "%.4f" % statistics.median(st) if st else "",
      "%.4f" % statistics.mean(st[:20]) if st else "", os.uname().nodename)))
EOF
  echo "$tag rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s; last: $(grep -o 'timeStepTime = [^ ]* s' $c/log.advectorVoF | tail -1); $(grep -m1 -o 'sweptPreTest [a-z]*, sweptStateShortcut [0-9], sweptTraceCache [0-9]' $c/log.advectorVoF)"
  [ $rc -eq 0 ] || { echo "FAILED $tag"; tail -20 $c/log.advectorVoF; }
done
echo "=== analysis $(date)"
# ANALYSE=only: no runs, hpcstruct + hpcprof of the measurements already in $T (a second 30-minute job)
[ -n "$MEASURE_ONLY" ] && { cat $csv; echo "HP-RUN-MEASURED $MODE $(date)"; exit 0; }
for c in $T/*/; do
  M=$c/hpctoolkit-measurements; [ -d $M ] || continue
  ls $M | grep -c '\.hpcrun$' | xargs echo "$(basename $c): hpcrun files"
  $HP/hpcstruct -j $(nproc) -c $P/struct-cache $M > $c/log.hpcstruct 2>&1 || echo "hpcstruct failed for $c"
  $HP/hpcprof -j $(nproc) -o $c/hpctoolkit-database $M > $c/log.hpcprof 2>&1 || echo "hpcprof failed for $c"
  ls $c/hpctoolkit-database 2>/dev/null | tr '\n' ' '; echo
done
cat $csv
echo "HP-RUN-DONE $MODE $(date)"
