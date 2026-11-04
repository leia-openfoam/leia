#!/bin/bash
#SBATCH -J davof-opt-restart
#SBATCH -A special00004
#SBATCH -p deflt_short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=8000
#SBATCH -t 00:25:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The restart identity of the swept-face advection (Claude, 2026-11-04): the rotation N 32 (ring row,
# OF Co 1) run 40 steps in one go (A) against 20 steps, a stop, and 20 more from the written state (B),
# writePrecision 17; the fields at step 40 compared bit for bit. A precondition of redistributing a
# parallel run between steps (item 5 b). Uses the default build of $W/TwoPhaseFlow.
W=/work/scratch/tm83tomy/davof/agents/swept
S=$W/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/prof/TwoPhaseFlow && . etc/davof-env.sh || exit 1
echo "binary $(command -v advectorVoF), fork $(git log --oneline -1 | cut -c1-80)"
T=$W/prof/opt/restart; rm -rf $T; mkdir -p $T
src=$(/usr/bin/python3 - "$S/sphereRotation3D_hexCmpRotSweptN32" <<'EOF'
import json, glob, os, sys
cs = []
for c in sorted(glob.glob(sys.argv[1] + "/*_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-ringArea": cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[1][1])
EOF
)
for c in A B; do
  mkdir -p $T/$c; cp -r $src/0 $src/constant $src/system $T/$c/
  dt=$(foamDictionary -entry deltaT -value $T/$c/system/controlDict)
  foamDictionary -entry writeControl -set timeStep $T/$c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set 20 $T/$c/system/controlDict > /dev/null
  foamDictionary -entry writePrecision -set 17 $T/$c/system/controlDict > /dev/null
  foamDictionary -entry startFrom -set latestTime $T/$c/system/controlDict > /dev/null
  sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" $T/$c/system/fvSolution
done
e40=$(python3 -c "print(repr(40*$dt))"); e20=$(python3 -c "print(repr(20*$dt))")
foamDictionary -entry endTime -set $e40 $T/A/system/controlDict > /dev/null
(cd $T/A && advectorVoF > log.advectorVoF 2>&1; echo "A rc=$? $(grep -c '^Time = ' log.advectorVoF) steps")
foamDictionary -entry endTime -set $e20 $T/B/system/controlDict > /dev/null
(cd $T/B && advectorVoF > log.1 2>&1; echo "B part 1 rc=$? $(grep -c '^Time = ' log.1) steps; written: $(ls -d [0-9]* | tr '\n' ' ')")
foamDictionary -entry endTime -set $e40 $T/B/system/controlDict > /dev/null
(cd $T/B && advectorVoF > log.2 2>&1; echo "B part 2 rc=$? $(grep -c '^Time = ' log.2) steps from $(grep -m1 '^Time = ' log.2)")
ta=$(ls -d $T/A/[0-9]* | sort -t/ -k8 -g | tail -1); tb=$T/B/$(basename $ta)
echo "compare $ta vs $tb"
same=0; diff=""
for f in $(cd $ta && ls | grep -v uniform); do if cmp -s $ta/$f $tb/$f; then same=$((same+1)); else diff="$diff $f"; fi; done
echo "step 40: $same identical field files; differing:${diff:- none}"
/usr/bin/python3 - $ta $tb <<'EOF'
import sys, re
def field(p):
    s = open(p).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\((.*?)\n\)", s, re.S)
    return [float(x) for x in m.group(2).split()] if m else None
for f in ("alpha.water", "alphaf.water"):
    a, b = field(sys.argv[1] + "/" + f), field(sys.argv[2] + "/" + f)
    d = [abs(x - y) for x, y in zip(a, b)]
    print("  %-13s differing %d of %d, max|diff| %.3e" % (f, sum(1 for x in d if x > 0), len(a), max(d)))
EOF
echo "A step-40 line: $(grep '^t = ' $T/A/log.advectorVoF | tail -1 | cut -c1-120)"
echo "B step-40 line: $(grep '^t = ' $T/B/log.2 | tail -1 | cut -c1-120)"
echo RESTART-DONE
