#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Item 5 b (Claude, 2026-11-04): the interface-weighted decomposition renewed every K steps with
# redistributePar (OpenFOAM v2606, weightField in decomposeParDict), the run continued from the
# written state (the fork restarts bit for bit, job 55931978):
#   sbatch -J davof-redist-<TAG> -n NP -t LIMIT --export=ALL,TAG=..,N=64|128,W=..,LAYERS=..,K=100,STEPS=full|n,KEYS=.. opt-redist.sbatch
# One turn: the initial weighted decomposition (opt-par's), then cycles of K steps; the wall clock
# of the solves and of the redistributions apart; the per-step history over the turn.
W_=/work/scratch/tm83tomy/davof/agents/swept
P=$W_/prof
S=$W_/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $P/TwoPhaseFlow && . etc/davof-env.sh || exit 1
: ${TAG:?} ${N:?}
NP=$SLURM_NTASKS; K=${K:-100}; STEPS=${STEPS:-full}; W=${W:-100}; LAYERS=${LAYERS:-2}
echo "node $(hostname), job $SLURM_JOB_ID, TAG $TAG N $N np $NP W $W LAYERS $LAYERS K $K STEPS $STEPS KEYS ${KEYS:-none}; fork $(git log --oneline -1 | cut -c1-70)"
T=$P/par/$TAG; rm -rf $T; mkdir -p $T
src=$(/usr/bin/python3 - "$S/sphereRotation3D_hexCmpRotSweptN$N" <<'EOF'
import json, glob, os, sys
cs = []
for c in sorted(glob.glob(sys.argv[1] + "/*_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-ringArea": cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[1][1])
EOF
)
c=$T/case; mkdir -p $c; cp -r $src/0 $src/constant $src/system $c/
dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
if [ $STEPS = full ]; then n=$(python3 -c "print(max(1, round(1.0/$dt)))"); else n=$STEPS; fi
foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
foamDictionary -entry writeInterval -set $K $c/system/controlDict > /dev/null
foamDictionary -entry writePrecision -set 17 $c/system/controlDict > /dev/null
foamDictionary -entry startFrom -set latestTime $c/system/controlDict > /dev/null
foamDictionary -entry purgeWrite -set 2 $c/system/controlDict > /dev/null
sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" $c/system/fvSolution
if [ -n "$KEYS" ]; then for kv in ${KEYS//+/ }; do
  sed -i "s/^\(\s*\)reconstructionScheme davof;/&\n\1${kv%%=*} ${kv#*=};/" $c/system/fvSolution
done; fi
foamDictionary -entry numberOfSubdomains -set $NP $c/system/decomposeParDict > /dev/null
foamDictionary -entry method -set scotch $c/system/decomposeParDict > /dev/null
/usr/bin/python3 -I $P/cellweight.py $c 0 $W $LAYERS
foamDictionary -entry weightField -set cellWeight $c/system/decomposeParDict > /dev/null
(cd $c && decomposePar -force > log.decomposePar 2>&1; echo "decomposePar rc=$?")
/usr/bin/python3 -I $P/rankcells.py $c 0
tsolve=0; tredist=0; step=0; cyc=0; allsteps=""
while [ $step -lt $n ]; do
  m=$(( n - step )); [ $m -gt $K ] && m=$K
  e=$(python3 -c "print(repr(($step+$m)*$dt))")
  foamDictionary -entry endTime -set "$e" $c/system/controlDict > /dev/null
  t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 advectorVoF -parallel > log.solve.$cyc 2>&1); rc=$?
  t1=$(date +%s.%N); tsolve=$(python3 -c "print($tsolve + $t1 - $t0)")
  st=$(grep -o 'timeStepTime = [^ ]*' $c/log.solve.$cyc | awk '{print $3}' | tr '\n' ' '); allsteps="$allsteps $st"
  ns=$(grep -c '^Time = ' $c/log.solve.$cyc); step=$(( step + ns )); cyc=$(( cyc + 1 ))
  [ $rc -eq 0 ] && [ $ns -gt 0 ] || { echo "solve failed (rc $rc, $ns steps) in cycle $cyc"; tail -5 $c/log.solve.$((cyc-1)); break; }
  if [ $step -lt $n ]; then
    latest=$(ls -d $c/processor0/[0-9]* | grep -v '/0$' | sort -t/ -k9 -g | tail -1 | xargs basename)
    /usr/bin/python3 -I $P/cellweight_par.py $c $latest $W $LAYERS
    t2=$(date +%s.%N)
    (cd $c && srun --ntasks=$NP --cpu-bind=cores redistributePar -parallel -overwrite -latestTime > log.redist.$cyc 2>&1); rr=$?
    t3=$(date +%s.%N); tredist=$(python3 -c "print($tredist + $t3 - $t2)")
    [ $rr -eq 0 ] || { echo "redistributePar failed rc $rr"; tail -8 $c/log.redist.$cyc; break; }
    /usr/bin/python3 -I $P/rankcells.py $c $latest | head -1
  fi
done
/usr/bin/python3 - "$allsteps" $tsolve $tredist $cyc $TAG <<'EOF'
import sys, statistics as st
t = [float(x) for x in sys.argv[1].split()]
b = [st.mean(t[i:i+50]) for i in range(0, len(t), 50)]
print("%s: %d steps in %d cycles; solve wall %.1f s, redistributions %.1f s (%.1f s each); mean step %.4f median %.4f; bins50 %s" % (sys.argv[5], len(t), int(sys.argv[4]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[3])/max(1, int(sys.argv[4]) - 1), st.mean(t), st.median(t), " ".join("%.2f" % x for x in b)))
EOF
echo "REDIST-DONE $TAG $(date)"
