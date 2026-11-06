#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Phase 3 of the port (docs/davof/PORT-DESIGN.md, 2026-11-06): the Enright agent's timed 3D
# deformation comparison (enright-time.sbatch, 2026-11-06) run on BOTH routes from the
# leia-route preparations port_enright{Swept,RDF}N<N>_leia: ring (the swept-face advection
# with the ring normal) and rdf (plicRDF + isoAdvector) on the leia route and on the fork
# route, the fork twin of a prepared case by the three tokenised lines (as
# davof_cmp/port-cmp-time.sbatch); prism (the plane-clipped loop, the Enright agent's third
# row) on the fork route only from its own preparation sphereDeformation3D_hexPrismN<N>
# (that loop is not ported). REPS repetitions of ring, one of the others; endTime T = 3,
# the fields and the diagnostics at T/4, T/2, 3T/4, T.
#   sbatch -J port-enright-N<N>-co<CO> -n <NP> -t <limit> --export=ALL,N=..,CO=0|1 port-enright-time.sbatch
# Optional: ROUTES="leia fork" (a third route, leiaGF, is the leia binary with `grazingCut fork`,
# the opt-in davofCutCellPLIC, PORT-DESIGN 1.2 item 1); METHODS="ring rdf prism"; REP0/REP1
# the repetitions to run (then only those runs of those methods and routes are removed, the
# others of the cell stay, and the timing file carries the suffix _r<REP0>_<methods>_<routes>).
W=/work/scratch/tm83tomy/davof/agents/port
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/TwoPhaseFlow && . etc/davof-env.sh || exit 1
cd $W/leia && . etc/leia-env.sh || exit 1
: ${N:?} ${CO:?}
REPS=${REPS:-3}
NP=$SLURM_NTASKS
END=3; WI=0.75
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-70); leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP, N $N, CO index $CO"
echo "advectorVoF $(command -v advectorVoF); davofAdvectorFoam $(command -v davofAdvectorFoam)"
ROUTES=${ROUTES:-leia fork}; METHODS=${METHODS:-ring rdf prism}
T=$W/enright-time/N${N}_co${CO}
if [ -z "$REP0" ]; then rm -rf $T
else for m in $METHODS; do for ro in $ROUTES; do for r in $(seq $REP0 ${REP1:-$REP0}); do rm -rf $T/${m}_${ro}_r$r; done; done; done; fi
mkdir -p $T
SUF=${REP0:+_r${REP0}_${METHODS// /-}_${ROUTES// /-}}
pick() { /usr/bin/python3 - "$1" "$2" "$CO" <<'EOF'
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
SRC[ring]=$(pick $S/port_enrightSweptN${N}_leia davof-planeSurf-ringArea)
SRC[rdf]=$(pick $S/port_enrightRDFN${N}_leia plicRDF)
SRC[prism]=$(pick $S/sphereDeformation3D_hexPrismN${N} davof-ringArea 2>/dev/null)
for m in ring rdf prism; do echo "source $m ${SRC[$m]:-none}"; done
toFork() {
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
  foamDictionary -entry endTime -set $END $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry writePrecision -set 17 $c/system/controlDict > /dev/null
  # the prism row's preparation is already a fork-route case; the others become one by the three lines
  [ $route = fork ] && [ $m != prism ] && toFork $c
  # the leiaGF route: the leia binary with the fork's grazing-cut threshold (the opt-in)
  if [ $route = leiaGF ]; then
    sed -i 's/^\(\s*\)grazingCut\s\+openfoam;/\1grazingCut fork;/' $c/system/fvSolution
    grep -q 'grazingCut\s\+fork;' $c/system/fvSolution || echo "LEIAGF: grazingCut not set in $c" >&2
  fi
  echo "$src" > $c/.source
  echo $c
}
csv=$T/timing$SUF.csv
echo "N,np,co_index,method,route,rep,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
orders=("ring rdf prism" "rdf prism ring" "prism ring rdf")
for r in $(seq ${REP0:-1} ${REP1:-$REPS}); do
  for m in ${orders[$(( (r-1) % 3 ))]}; do
    [ $r -gt 1 ] && [ $m != ring ] && continue
    [[ " $METHODS " == *" $m "* ]] || continue
    routes="$ROUTES"; [ $m = prism ] && routes="fork"
    [ $m = prism ] && [ -z "${SRC[prism]}" ] && continue
    [ $m = prism ] && [[ " $ROUTES " != *" fork "* ]] && continue
    for route in $routes; do
      c=$(prep $m $route $r)
      par=""; [ $NP -gt 1 ] && par="-parallel"
      if [ $route != fork ]; then app=davofAdvectorFoam; log=log.davofAdvectorFoam; else app=advectorVoF; log=log.advectorVoF; fi
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
print(",".join(str(x) for x in (e["N"], e["SLURM_NTASKS"], e["CO"], m, route, r, steps, dt, et,
      ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0), "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
      echo "$m $route r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s, steps $(grep -c '^Time = ' $c/$log), $(grep -c 'ABORT' $c/$log) abort lines"
    done
  done
done
sleep 5
sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps$SUF.txt
cat $T/sacct_steps$SUF.txt
echo ENRIGHT-TIME-DONE
