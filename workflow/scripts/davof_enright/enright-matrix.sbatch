#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The matrix runner of the decisions 1-3 tests (2026-11-07): every line of SPEC is one run
#   name|source case|np|steps (or end:<t>)|reps|key=value;key=value;...
# copied from the prepared source (0/ or processor*/), endTime = steps * deltaT (or end:<t>),
# one field write at the end (timeStep control), writePrecision 17, the keys set in the
# alpha.water dictionary, run with srun on np ranks (--cpu-bind=cores), timed, the converter
# run after it; PARALLEL=1 runs the lines concurrently (identity tests), else sequentially
# (timing). Output under $W/matrix/$TAG/<name>[_r<rep>]/ and $W/matrix/$TAG/timing.csv.
#   sbatch -J davof-mx-<tag> -n <max np> [--exclusive --constraint=i01] -t .. --export=ALL,SPEC=..,TAG=..,PARALLEL=0|1 enright-matrix.sbatch
W=${ENRIGHT_W:-/work/scratch/tm83tomy/davof/agents/enright}
APP=${APP:-advectorVoF}
: ${SPEC:?} ${TAG:?}
PARALLEL=${PARALLEL:-0}
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"; cd $W/leia && . etc/leia-env.sh || exit 1
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-90); leia $(git -C $W/leia log --oneline -1 | cut -c1-90)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), job $SLURM_JOB_ID, ntasks $SLURM_NTASKS, SPEC $SPEC, TAG $TAG, PARALLEL $PARALLEL"
D=$W/matrix/$TAG; mkdir -p $D
csv=$D/timing.csv; [ -f $csv ] || echo "name,rep,np,steps,dt,wall_s,clock_s,exec_s,median_step_s,rc,abort,node" > $csv
setkey() { # case key value
  if grep -qE "^\s*$2\s" $1/system/fvSolution; then sed -i -E "s/^(\s*)$2\s+[^;]*;/\1$2 $3;/" $1/system/fvSolution
  else sed -i -E "s/^(\s*)advectionScheme\s+isoAdvection;/\1advectionScheme isoAdvection;\n\1$2 $3;/" $1/system/fvSolution; fi
}
runone() { # name src np steps rep keys
  local name=$1 src=$2 np=$3 steps=$4 rep=$5 keys=$6
  local c=$D/$name; [ "$rep" != "1" ] && c=$D/${name}_r$rep
  rm -rf $c; mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $np -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
    [ -d $c/processor0 ] || { echo "$name: no decomposed source for np $np"; return; }
  else
    if [ -d $src/0 ]; then cp -r $src/0 $c/; else
      # a decomposed source run in serial: the initial state reconstructed
      for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
      (cd $c && reconstructPar -withZero -time 0 > log.reconstructPar0 2>&1) && rm -rf $c/processor*
      [ -d $c/0 ] || { echo "$name: reconstructPar of the initial state failed"; return; }
    fi
  fi
  local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
  local endT; local nst
  if [[ "$steps" == end:* ]]; then endT=${steps#end:}; nst=$(python3 -c "print(int(round($endT/$dt)))"); else nst=$steps; endT=$(python3 -c "print(repr($steps*$dt))"); fi
  foamDictionary -entry endTime -set $endT $c/system/controlDict > /dev/null
  foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $nst $c/system/controlDict > /dev/null
  foamDictionary -entry writePrecision -set 17 $c/system/controlDict > /dev/null
  for fo in areaFractionError volumeFractionError; do
    foamDictionary -entry functions/$fo/writeControl -set timeStep $c/system/controlDict > /dev/null 2>&1
    foamDictionary -entry functions/$fo/writeInterval -set $nst $c/system/controlDict > /dev/null 2>&1
    foamDictionary -entry functions/$fo/writePrecision -set 17 $c/system/controlDict > /dev/null 2>&1
  done
  IFS=';' read -ra KV <<< "$keys"
  for kv in "${KV[@]}"; do [ -z "$kv" ] && continue; setkey $c ${kv%%=*} ${kv#*=}; done
  echo "$name r$rep np $np steps $nst dt $dt keys: $keys"
  local par=""; [ $np -gt 1 ] && par="-parallel"
  local t0=$(date +%s.%N)
  # stdin closed: srun must not read the spec lines of the enclosing loop
  (cd $c && srun --ntasks=$np --cpu-bind=cores --exclusive --kill-on-bad-exit=1 --job-name=$name $APP $par < /dev/null > log.advectorVoF 2>&1); local rc=$?
  local t1=$(date +%s.%N)
  (cd $c && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1)
  /usr/bin/python3 - $c $name $rep $np $dt $rc $t0 $t1 >> $csv <<'EOF'
import sys, re, os, statistics
c, name, rep, np_, dt, rc, t0, t1 = sys.argv[1:9]
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log)); st = [float(x[2]) for x in ct if x[2]]
ab = re.search(r"ABORT at step (\d+)", log)
print(",".join(str(x) for x in (name, rep, np_, steps, dt, "%.2f" % (float(t1) - float(t0)), ct[-1][1] if ct else "", ct[-1][0] if ct else "",
      "%.4f" % statistics.median(st) if st else "", rc, ab.group(1) if ab else "", os.uname().nodename)))
EOF
  echo "$name r$rep rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s $(grep -c '^Time = ' $c/log.advectorVoF) steps, abort $(grep -c ABORT $c/log.advectorVoF)"
}
export -f runone setkey; export D W APP csv
exec 3< $SPEC
while IFS='|' read -r -u 3 name src np steps reps keys; do
  [ -z "$name" ] && continue; [[ "$name" == \#* ]] && continue
  reps=${reps:-1}
  for r in $(seq 1 $reps); do
    if [ "$PARALLEL" = "1" ]; then runone "$name" "$src" "$np" "$steps" "$r" "$keys" & else runone "$name" "$src" "$np" "$steps" "$r" "$keys"; fi
  done
done
exec 3<&-
wait
echo MATRIX-DONE
