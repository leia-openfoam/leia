#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The port's timed comparison (docs/davof/PORT-DESIGN.md 4.4, 2026-11-06), the fork's
# davof-cmp-time.sbatch (2026-11-03) with BOTH routes on one exclusive i01 node per
# (FLOW, N, CO): the prepared leia-route cases of the studies port_cmp<FLOW>{Swept,RDF}N<N>_leia
# (rendered, meshed, initialised, decomposed; END_TIME 0.02), copied per method and route;
# the fork-route twin of a case differs in the three tokenised lines only (application,
# advectionScheme, the function objects' libs), rewritten here by sed: the acceptance's own
# premise. Methods: ring, exact (the DAVOF arm, davofSwept against the fork's isoAdvection at
# c5c00c5) and rdf (plicRDF + isoAdvector: OpenFOAM v2606's against TwoPhaseFlow's), each on
# both routes; REPS repetitions of ring (the timed cell), one of exact and rdf; the order of
# the methods rotated per repetition; endTime n dt, n = round(END/dt); the fields and the
# end-only diagnostics at the last step.
#   sbatch -J port-cmp-<FLOW>-N<N>-co<CO> -n <NP> -t <limit> --export=ALL,FLOW=Rot|LineX|LineDiag,N=32|64|128,CO=0|1 port-cmp-time.sbatch
W=/work/scratch/tm83tomy/davof/agents/port
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/TwoPhaseFlow && . etc/davof-env.sh || exit 1
cd $W/leia && . etc/leia-env.sh || exit 1
: ${FLOW:?} ${N:?} ${CO:?}
REPS=${REPS:-3}
NP=$SLURM_NTASKS
case $FLOW in Rot) END=1;; *) END=0.45;; esac
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-70); leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP"
echo "advectorVoF $(command -v advectorVoF); davofAdvectorFoam $(command -v davofAdvectorFoam)"
T=$W/cmp-time/${FLOW}_N${N}_co${CO}
rm -rf $T; mkdir -p $T
pick() { /usr/bin/python3 - "$S/port_cmp${FLOW}$1N${N}_leia" "$2" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/sphereRotation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort()
print(cs[co][1])
EOF
}
declare -A SRC
SRC[ring]=$(pick Swept davof-planeSurf-ringArea)
SRC[exact]=$(pick Swept davof-planeSurf-exact)
SRC[rdf]=$(pick RDF plicRDF)
for m in ring exact rdf; do echo "source $m ${SRC[$m]}"; done
# the fork-route twin: the three tokenised lines (materialize.py DAVOF_* tokens)
toFork() {
  # foamDictionary reformats the dictionaries it edits (spaces inside the parentheses),
  # so the patterns tolerate white space; warnings go to stderr (prep's stdout is the path)
  sed -i 's/^application\s\+davofAdvectorFoam;/application     advectorVoF;/' $1/system/controlDict
  sed -i 's/(\s*"libleiaDavofFunctionObjects.so"\s*)/("libpostProcess.so")/' $1/system/controlDict
  sed -i 's/^\(\s*\)advectionScheme\s\+\(davofSwept\|isoAdvector\);/\1advectionScheme      isoAdvection;/' $1/system/fvSolution
  local n=$(cat $1/system/controlDict $1/system/fvSolution | grep -c 'application\s\+advectorVoF;\|"libpostProcess.so"\|advectionScheme\s\+isoAdvection;')
  [ "$n" -ge 3 ] || echo "TOFORK: $n of 3 lines rewritten in $1" >&2
  true
}
prep() { # method route rep -> case dir
  local m=$1 route=$2 r=$3 src=${SRC[$1]} c=$T/${1}_${2}_r$3
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
  foamDictionary -entry writePrecision -set 17 $c/system/controlDict > /dev/null
  [ $route = fork ] && toFork $c
  echo $c
}
csv=$T/timing.csv
echo "flow,N,np,co_index,method,route,rep,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
orders=("ring exact rdf" "exact rdf ring" "rdf ring exact")
for r in $(seq 1 $REPS); do
  for m in ${orders[$(( (r-1) % 3 ))]}; do
    [ $r -gt 1 ] && [ $m != ring ] && continue
    for route in leia fork; do
      c=$(prep $m $route $r)
      par=""; [ $NP -gt 1 ] && par="-parallel"
      if [ $route = leia ]; then app=davofAdvectorFoam; log=log.davofAdvectorFoam; else app=advectorVoF; log=log.advectorVoF; fi
      t0=$(date +%s.%N)
      (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${m}_${route}_r$r $app $par > $log 2>&1); rc=$?
      t1=$(date +%s.%N)
      (cd $c && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1)
      /usr/bin/python3 - $c $m $route $r $rc $t0 $t1 $log >> $csv <<'EOF'
import sys, re, os, statistics
c, m, route, r, rc, t0, t1, logn = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
log = open(c + "/" + logn).read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log))
st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
e = os.environ
print(",".join(str(x) for x in (e["FLOW"], e["N"], e["SLURM_NTASKS"], e["CO"], m, route, r, steps, dt, et,
      ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0), "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
      echo "$m $route r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s"
    done
  done
done
sleep 5
sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps.txt
cat $T/sacct_steps.txt
echo CMP-TIME-DONE
