#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# HPCToolkit 2025.0.1 profiles of the rotation, OpenFOAM Co 1 (Claude, 2026-11-04):
# the swept-face advection (ring normal, scale correction, the final wisp rule)
# and plicRDF + isoAdvector, prepared from the comparison's study cases exactly
# as davof-cmp-time.sbatch does, run with the profiling build prof/TwoPhaseFlow.
#   sbatch -J davof-hp-<MODE> -n NP [-c C --mem-per-cpu=M] -t LIMIT --export=ALL,MODE=serial|par8|par64 hp-run.sbatch
# A run spec is <method><N>:<steps|full>[:nohp]; nohp = no hpcrun (the baseline).
# hpcrun samples CPUTIME every 1000 us and LLC misses every 100000 on every rank;
# hpcstruct (cached) and hpcprof follow in the job.
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
T=$P/runs/$MODE
rm -rf $T; mkdir -p $T
echo "fork worktree $(git -C $P/TwoPhaseFlow log --oneline -1); leia $(git -C $W/leia log --oneline -1)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP, cpus/task ${SLURM_CPUS_PER_TASK:-1}"
echo "advectorVoF $(command -v advectorVoF)"; ldd $(command -v advectorVoF) | grep libVoF
echo "$($HP/hpcrun --version | head -1); perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid)"
case $MODE in
  serial) RUNS="ring32:full:nohp ring32:full rdf32:full:nohp rdf32:full ring64:20:nohp ring64:20 ring64:full rdf64:20:nohp rdf64:full ring128:20 rdf128:20";;
  par8)   RUNS="ring64:20:nohp ring64:full rdf64:20:nohp rdf64:full";;
  par64)  RUNS="ring128:20:nohp ring128:full rdf128:20:nohp rdf128:full";;
  *) echo "unknown MODE $MODE"; exit 1;;
esac
# perf-event availability on this node
EV="-e CPUTIME@1000"
if $HP/hpcrun -e ix86arch::LLC_MISSES@100000 -o $T/.evtest /bin/true > $T/.evtest.log 2>&1; then EV="$EV -e ix86arch::LLC_MISSES@100000"; else echo "LLC_MISSES not available: $(tail -1 $T/.evtest.log)"; fi
echo "hpcrun events: $EV"
pick() { /usr/bin/python3 - "$S/sphereRotation3D_hexCmpRot$1N$2" "$3" <<'EOF'
import json, glob, os, sys
st, row = sys.argv[1], sys.argv[2]
cs = []
for c in sorted(glob.glob(st + "/sphereRotation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort()
print(cs[1][1])
EOF
}
prep() { # method N steps|full tag -> case dir
  local m=$1 N=$2 steps=$3 tag=$4 src c=$T/$4
  case $m in ring) src=$(pick Swept $N davof-planeSurf-ringArea);; rdf) src=$(pick RDF $N plicRDF);; esac
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
  if [ $m = ring ]; then sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" $c/system/fvSolution; fi
  echo "$src" > $c/.source
  echo $c
}
csv=$T/timing.csv
echo "mode,np,method,N,steps,hpcrun,wall_s,exec_s,clock_s,n_steps,mean_step_s,median_step_s,first20_mean_s,node" > $csv
par=""; [ $NP -gt 1 ] && par="-parallel"
for spec in $RUNS; do
  IFS=: read mN steps hp <<< "$spec"
  m=${mN%%[0-9]*}; N=${mN##*[a-z]}
  tag=${m}${N}_${steps}${hp:+_$hp}
  c=$(prep $m $N $steps $tag)
  echo "=== $tag from $(cat $c/.source) $(date)"
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
  echo "$tag rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s; last: $(grep -o 'timeStepTime = [^ ]* s' $c/log.advectorVoF | tail -1)"
  [ $rc -eq 0 ] || { echo "FAILED $tag"; tail -20 $c/log.advectorVoF; }
done
echo "=== analysis $(date)"
for c in $T/*/; do
  M=$c/hpctoolkit-measurements; [ -d $M ] || continue
  ls $M | grep -c '\.hpcrun$' | xargs echo "$(basename $c): hpcrun files"
  $HP/hpcstruct -j $(nproc) -c $P/struct-cache $M > $c/log.hpcstruct 2>&1 || echo "hpcstruct failed for $c"
  $HP/hpcprof -j $(nproc) -o $c/hpctoolkit-database $M > $c/log.hpcprof 2>&1 || echo "hpcprof failed for $c"
  ls $c/hpctoolkit-database 2>/dev/null | tr '\n' ' '; echo
done
cat $csv
echo "HP-RUN-DONE $MODE $(date)"
