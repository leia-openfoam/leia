#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The decomposition of the swept-face advection (Claude, 2026-11-04, item 5 a):
#   sbatch -J davof-par-<TAG> -n NP -t LIMIT --export=ALL,TAG=..,N=64|128,DECOMP=scotch|weighted,W=..,LAYERS=..,REPS=..,STEPS=full|n opt-par.sbatch
# The rotation (ring row, OF Co 1, the final wisp rule, the optimisation defaults of the prof build),
# the undecomposed study case decomposed here: scotch over the whole mesh (the comparison's) or scotch
# with the cell weights cellWeight (cellweight.py: W for the cells within LAYERS face layers of a
# partial cell at t = 0, 1 elsewhere); the partial cells per rank, REPS runs of STEPS steps, the
# per-step time, the 50-step-bin history, the in-code timers.
W_=/work/scratch/tm83tomy/davof/agents/swept
P=$W_/prof
S=$W_/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $P/TwoPhaseFlow && . etc/davof-env.sh || exit 1
: ${TAG:?} ${N:?} ${DECOMP:?}
NP=$SLURM_NTASKS; REPS=${REPS:-1}; STEPS=${STEPS:-full}; W=${W:-100}; LAYERS=${LAYERS:-2}
echo "node $(hostname), job $SLURM_JOB_ID, TAG $TAG N $N np $NP DECOMP $DECOMP W $W LAYERS $LAYERS REPS $REPS STEPS $STEPS; fork $(git log --oneline -1 | cut -c1-70)"
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
foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null
sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" $c/system/fvSolution
# KEYS=key=value+key=value into the alpha.water.* dictionary (the optimisation keys)
if [ -n "$KEYS" ]; then for kv in ${KEYS//+/ }; do
  sed -i "s/^\(\s*\)reconstructionScheme davof;/&\n\1${kv%%=*} ${kv#*=};/" $c/system/fvSolution
done; fi
echo "KEYS: ${KEYS:-none}"
foamDictionary -entry numberOfSubdomains -set $NP $c/system/decomposeParDict > /dev/null
foamDictionary -entry method -set scotch $c/system/decomposeParDict > /dev/null
if [ $DECOMP = weighted ]; then
  /usr/bin/python3 -I $P/cellweight.py $c 0 $W $LAYERS
  foamDictionary -entry weightField -set cellWeight $c/system/decomposeParDict > /dev/null
fi
(cd $c && decomposePar -force > log.decomposePar 2>&1; echo "decomposePar rc=$?"); grep -E 'Processor [0-9]+$|Number of cells' $c/log.decomposePar | head -0
/usr/bin/python3 -I $P/rankcells.py $c 0
csv=$T/timing.csv; echo "tag,N,np,decomp,W,layers,rep,steps,wall_s,mean_step_s,median_step_s,max_step_s,bins50" > $csv
for r in $(seq 1 $REPS); do
  rm -rf $c/processor*/[1-9]* $c/processor*/0.[0-9]*
  t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 advectorVoF -parallel > log.advectorVoF.r$r 2>&1); rc=$?
  t1=$(date +%s.%N)
  /usr/bin/python3 - $c/log.advectorVoF.r$r $t0 $t1 $r $TAG $N $NP $DECOMP $W $LAYERS >> $csv <<'EOF'
import sys, re, statistics as st
log = open(sys.argv[1]).read()
t = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)]
b = [st.mean(t[i:i+50]) for i in range(0, len(t), 50)]
rows = re.findall(r"wall time max over (\d+) ranks: advect (\S+) s, curvedFlux (\S+) s, exchange (\S+) s \(wait (\S+)", log)
print(",".join(str(x) for x in (sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[4], len(t), "%.1f" % (float(sys.argv[3]) - float(sys.argv[2])), "%.4f" % st.mean(t), "%.4f" % st.median(t), "%.4f" % max(t[:-1]) if len(t) > 1 else "", " ".join("%.3f" % x for x in b))))
print("rep %s: steps %d, wall %.1f s, mean step %.4f, median %.4f, bins50 %s; in-code advect %.4f exchange %.4f wait %.4f (means)" % (sys.argv[4], len(t), float(sys.argv[3]) - float(sys.argv[2]), st.mean(t), st.median(t), " ".join("%.2f" % x for x in b),
      st.mean([float(x[1]) for x in rows]) if rows else 0, st.mean([float(x[3]) for x in rows]) if rows else 0, st.mean([float(x[4]) for x in rows]) if rows else 0), file=sys.stderr)
EOF
  echo "rep $r rc=$rc; $(grep -o 'faces by path: [^;]*' $c/log.advectorVoF.r$r | tail -1)"
done
cat $csv | cut -d, -f7-12
echo "PAR-DONE $TAG $(date)"
