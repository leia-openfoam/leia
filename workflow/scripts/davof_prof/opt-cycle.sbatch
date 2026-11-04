#!/bin/bash
#SBATCH -A special00004
#SBATCH -p deflt_short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3800
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# One optimisation cycle of the swept-face advection (Claude, 2026-11-04):
#   sbatch -J davof-opt-<TAG> --export=ALL,TAG=..,SHA=<commit>,REF=<reference tag> opt-cycle.sbatch
# 1. the profiling build of prof/TwoPhaseFlow at SHA (-O3 -g -fno-omit-frame-pointer; skipped when SHA is empty),
# 2. the identity checks: 20-step serial copies (writePrecision 17) of the rotation N 32 and 64 (ring row, OF Co 1),
#    the plane (normal arm, Co 0.95, N 40, exact row of hexSweptWispM9) and the diagonal line N 32 (ring row), the
#    written fields compared bit for bit (cmp) and by max |difference| with the REF tag, the log counters diffed,
# 3. the serial timing: the rotation N 64, STEPS steps, REPS repetitions (median and spread of the per-step time),
# 4. a PSTEPS-step hpcrun profile of the same case, analysed with hpan.py (rules_swept.py).
W=/work/scratch/tm83tomy/davof/agents/swept
P=$W/prof
S=$W/leia/studies/davof
HP=/home/tm83tomy/spack/opt/spack/linux-x86_64/hpctoolkit-2025.0.1-ky3i4lswv3ynbyvk6ygnunxtn6tsfa67/bin
: ${TAG:?}
REF=${REF:-before}; STEPS=${STEPS:-200}; REPS=${REPS:-3}; PSTEPS=${PSTEPS:-100}
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
echo "node $(hostname), job $SLURM_JOB_ID, TAG $TAG SHA ${SHA:-none} REF $REF, $(date)"
if [ -n "$SHA" ]; then
  git -C $W/TwoPhaseFlow fetch -q origin davof-swept || echo "fetch failed"
  cd $P/TwoPhaseFlow && git checkout -q -- . && git checkout -q --detach $SHA || { echo "checkout $SHA failed"; exit 1; }
  export FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer" WM_NCOMPPROCS=${SLURM_CPUS_PER_TASK:-24}
  find . -name 'Allwmake*' -o -name 'Allwclean*' | xargs chmod +x
  wmakeLnInclude -u src/VoF > /dev/null 2>&1
  t0=$(date +%s); ./Allwmake > build.log 2>&1; rc=$?
  echo "build at $(git log --oneline -1) rc=$rc in $(( $(date +%s) - t0 )) s, errors $(grep -c -i ' error' build.log)"
  [ $rc -eq 0 ] || { grep -B2 -A8 -i ' error' build.log | head -40; exit 1; }
fi
cd $P/TwoPhaseFlow && . etc/davof-env.sh || exit 1
echo "binary $(command -v advectorVoF), fork $(git -C $P/TwoPhaseFlow log --oneline -1)"
T=$P/opt/$TAG; rm -rf $T; mkdir -p $T/verify $T/timing
pick() { /usr/bin/python3 - "$1" "$2" "$3" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[co][1])
EOF
}
pickplane() { /usr/bin/python3 - "$S/planeTranslation3D_hexSweptWispM9" <<'EOF'
import json, glob, os, sys
for c in sorted(glob.glob(sys.argv[1] + "/planeTranslation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-exact" and t["TPF_VELOCITY_ARM"] == "normal" and float(t["TARGET_ADVECTIVE_CO"]) == 0.95 and int(float(t["N_CELLS"])) == 40:
        print(c); break
EOF
}
prep() { # src dst steps [write interval]
  mkdir -p "$2"; cp -r "$1/0" "$1/constant" "$1/system" "$2/"
  local dt=$(foamDictionary -entry deltaT -value "$2/system/controlDict")
  foamDictionary -entry endTime -set "$(python3 -c "print(repr($3*$dt))")" "$2/system/controlDict" > /dev/null
  foamDictionary -entry writeControl -set timeStep "$2/system/controlDict" > /dev/null
  foamDictionary -entry writeInterval -set ${4:-$3} "$2/system/controlDict" > /dev/null
  foamDictionary -entry writePrecision -set 17 "$2/system/controlDict" > /dev/null
  sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" "$2/system/fvSolution"
  # KEYS=key=value+key=value ('+' separated: sbatch --export splits commas): inserted into the
  # alpha.water.* dictionary after reconstructionScheme
  if [ -n "$KEYS" ]; then for kv in ${KEYS//+/ }; do
    sed -i "s/^\(\s*\)reconstructionScheme davof;/&\n\1${kv%%=*} ${kv#*=};/" "$2/system/fvSolution"
  done; fi
  echo "$1" > "$2/.source"
}
echo "KEYS: ${KEYS:-none}"
# ---- 2. identity checks: the rotation N 32 and 64, the plane N 40 normal (exact row, M9), the
# diagonal N 32, and the plane gate (hexSweptPyr: N 20 and 40, both arms, ring row, Co 0.5, the
# final wisp rule matched + consistent; sweptVolumeCorrection scale as everywhere here)
pickpyr() { /usr/bin/python3 - "$S/planeTranslation3D_hexSweptPyr" "$1" "$2" <<'EOF'
import json, glob, os, sys
for c in sorted(glob.glob(sys.argv[1] + "/planeTranslation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-ringArea" and t["TPF_VELOCITY_ARM"] == sys.argv[3] and float(t["TARGET_ADVECTIVE_CO"]) == 0.5 and int(float(t["N_CELLS"])) == int(sys.argv[2]):
        print(c); break
EOF
}
VSET=${VSET:-"R32 R64 P D P20n P20o P40n P40o"}
prep $(pick $S/sphereRotation3D_hexCmpRotSweptN32 davof-planeSurf-ringArea 1) $T/verify/R32 20
prep $(pick $S/sphereRotation3D_hexCmpRotSweptN64 davof-planeSurf-ringArea 1) $T/verify/R64 20
prep $(pickplane) $T/verify/P 20
prep $(pick $S/sphereRotation3D_hexCmpLineDiagSweptN32 davof-planeSurf-ringArea 1) $T/verify/D 20
prep $(pickpyr 20 normal) $T/verify/P20n 20
prep $(pickpyr 20 oblique) $T/verify/P20o 20
prep $(pickpyr 40 normal) $T/verify/P40n 20
prep $(pickpyr 40 oblique) $T/verify/P40o 20
for c in $VSET; do
  (cd $T/verify/$c && advectorVoF > log.advectorVoF 2>&1; echo "verify $c rc=$? $(grep -c '^Time = ' log.advectorVoF) steps; $(grep -o 'pre-test check: [^;]*' log.advectorVoF | tail -1)")
done
if [ -d $P/opt/$REF/verify ] && [ "$REF" != "$TAG" ]; then
  echo "=== identity against $REF"
  for c in $VSET; do
    a=$T/verify/$c; b=$P/opt/$REF/verify/$c
    ta=$(ls -d $a/[0-9]* | grep -v '/0$' | sort -t/ -k2 -g | tail -1); tb=$b/$(basename $ta)
    [ -d "$tb" ] || { echo "$c: no $tb in $REF"; continue; }
    same=0; diff=""
    for f in $(cd $ta && ls | grep -v uniform); do if cmp -s $ta/$f $tb/$f; then same=$((same+1)); else diff="$diff $f"; fi; done
    echo "$c: $same identical field files; differing:${diff:- none}"
    /usr/bin/python3 - $ta $tb <<'EOF'
import sys, re
def field(p):
    s = open(p).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\((.*?)\n\)", s, re.S)
    if not m: return None
    return [float(x) for x in m.group(2).split()]
for f in ("alpha.water", "alphaf.water"):
    try:
        a, b = field(sys.argv[1] + "/" + f), field(sys.argv[2] + "/" + f)
    except FileNotFoundError:
        continue
    if a is None or b is None: continue
    d = [abs(x - y) for x, y in zip(a, b)]
    n = sum(1 for x in d if x > 0)
    print("    %-13s n %d, differing %d, max|diff| %.3e, max|diff| relative to max|a| %.3e" % (f, len(a), n, max(d), max(d)/max(1e-300, max(abs(x) for x in a))))
EOF
    # the per-step log lines without the timings
    diff <(grep -v -E 'ExecutionTime|ClockTime|timeStepTime|wall time|^Date|^Host|^PID|^Case|^nProcs|^Build|^Exec|^Time = |^Create|^Reading|Starting|Writing|^End|^$' $a/log.advectorVoF) \
         <(grep -v -E 'ExecutionTime|ClockTime|timeStepTime|wall time|^Date|^Host|^PID|^Case|^nProcs|^Build|^Exec|^Time = |^Create|^Reading|Starting|Writing|^End|^$' $b/log.advectorVoF) > $a/log.diff
    echo "    log lines differing (ignoring timings): $(grep -c '^[<>]' $a/log.diff); first: $(grep -m1 '^[<>]' $a/log.diff | cut -c1-200)"
  done
fi
[ -n "$VERIFY_ONLY" ] && { echo "OPT-CYCLE-DONE $TAG (verify only) $(date)"; exit 0; }
# ---- 3. timing
src=$(pick $S/sphereRotation3D_hexCmpRotSweptN64 davof-planeSurf-ringArea 1)
for r in $(seq 1 $REPS); do
  prep $src $T/timing/r$r $STEPS
  t0=$(date +%s.%N); (cd $T/timing/r$r && advectorVoF > log.advectorVoF 2>&1); rc=$?; t1=$(date +%s.%N)
  /usr/bin/python3 - $T/timing/r$r $t0 $t1 $r <<'EOF' | tee -a $T/timing/timing.txt
import sys, re, statistics as st
log = open(sys.argv[1] + "/log.advectorVoF").read()
t = [float(x) for x in re.findall(r"timeStepTime = (\S+) s", log)][1:-1]
sh = re.findall(r"faces by path: full (\d+), empty (\d+), pieces (\d+)", log)
pre = re.findall(r"pre-test: ring cells skipped (\d+), clipped (\d+)", log)
print("rep %s: steps %d mean %.4f median %.4f min %.4f max %.4f s/step; wall %.1f s%s%s" % (sys.argv[4], len(t), st.mean(t), st.median(t), min(t), max(t), float(sys.argv[3]) - float(sys.argv[2]),
      ("; last step faces full/empty/pieces %s/%s/%s" % sh[-1]) if sh else "", ("; ring cells skipped/clipped %s/%s" % pre[-1]) if pre else ""))
EOF
done
/usr/bin/python3 - $T/timing/timing.txt <<'EOF'
import sys, re, statistics as st
m = [float(x) for x in re.findall(r"mean (\S+) median", open(sys.argv[1]).read())]
print("TIMING %s: mean s/step over reps: median %.4f, min %.4f, max %.4f (spread %.1f%%)" % (sys.argv[1].split("/")[-3], st.median(m), min(m), max(m), 100*(max(m)-min(m))/st.median(m)))
EOF
# ---- 4. profile
prep $src $T/profile $PSTEPS
(cd $T/profile && $HP/hpcrun -e CPUTIME@1000 -e ix86arch::LLC_MISSES@100000 -o hpctoolkit-measurements advectorVoF > log.advectorVoF 2>&1; echo "profile rc=$?")
$HP/hpcstruct -j 24 -c $P/struct-cache $T/profile/hpctoolkit-measurements > $T/profile/log.hpcstruct 2>&1
$HP/hpcprof -j 24 -o $T/profile/hpctoolkit-database $T/profile/hpctoolkit-measurements > $T/profile/log.hpcprof 2>&1
n=$(grep -c '^Time = ' $T/profile/log.advectorVoF)
/usr/bin/python3 -I $P/hpan.py --db $T/profile/hpctoolkit-database --out $T/profile/an --steps $n --rules $P/rules_swept.py --src 'TwoPhaseFlow|OpenFOAM' --tree-min 0.4 > $T/profile/an.log 2>&1
echo "=== profile categories (s/step, %), $n steps"; tail -3 $T/profile/an.log
/usr/bin/python3 - $T/profile/an.cat.csv <<'EOF'
import csv, sys
try: rows = list(csv.DictReader(open(sys.argv[1])))
except Exception as e: print("no categories:", e); sys.exit(0)
for r in sorted(rows, key=lambda r: -float(r["cputime_s"])):
    if r["category"].startswith("(") or float(r["pct"]) < 0.3: continue
    print("  %-78s %8s s/step %6s %%" % (r["category"][:78], r["s_per_step_per_rank"], r["pct"]))
EOF
grep -o 'faces by path: [^;]*\|pre-test: [^;]*\|trace cache: [^;]*\|cells: pieces [0-9]*' $T/profile/log.advectorVoF | tail -4
cat $T/timing/timing.txt
echo "OPT-CYCLE-DONE $TAG $(date)"
