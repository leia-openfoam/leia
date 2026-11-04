#!/bin/bash
#SBATCH -J davof-opt-dbg
#SBATCH -A special00004
#SBATCH -p deflt_short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=3000
#SBATCH -t 00:25:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# A debug build in its own worktree (prof/TwoPhaseFlow-dbg at SHA) and the two N 20 plane cases with
# the pre-test check, so that the running cycle chain's worktree is not touched (Claude, 2026-11-04).
W=/work/scratch/tm83tomy/davof/agents/swept
P=$W/prof
S=$W/leia/studies/davof
: ${SHA:?}
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
git -C $W/TwoPhaseFlow fetch -q origin
D=$P/TwoPhaseFlow-dbg
if [ ! -d $D ]; then git -C $W/TwoPhaseFlow worktree add -q --detach $D $SHA || exit 1; else (cd $D && git checkout -q -- . && git checkout -q --detach $SHA) || exit 1; fi
cd $D && find . -name 'Allwmake*' -o -name 'Allwclean*' | xargs chmod +x
export FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer" WM_NCOMPPROCS=24
wmakeLnInclude -u src/VoF > /dev/null 2>&1
./Allwmake > build.log 2>&1; rc=$?; echo "dbg build at $(git log --oneline -1 | cut -c1-70) rc=$rc errors $(grep -c -i ' error' build.log)"
[ $rc -eq 0 ] && [ $(grep -c -i ' error' build.log) -eq 0 ] || { grep -B2 -A6 -i ' error' build.log | head -30; exit 1; }
. etc/davof-env.sh
T=$P/opt/dbg; rm -rf $T; mkdir -p $T
pickpyr() { /usr/bin/python3 - "$S/planeTranslation3D_hexSweptPyr" "$1" "$2" <<'EOF'
import json, glob, os, sys
for c in sorted(glob.glob(sys.argv[1] + "/planeTranslation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-ringArea" and t["TPF_VELOCITY_ARM"] == sys.argv[3] and float(t["TARGET_ADVECTIVE_CO"]) == 0.5 and int(float(t["N_CELLS"])) == int(sys.argv[2]):
        print(c); break
EOF
}
pickrot() { /usr/bin/python3 - "$S/sphereRotation3D_hexCmpRotSweptN32" <<'EOF'
import json, glob, os, sys
cs = []
for c in sorted(glob.glob(sys.argv[1] + "/*_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == "davof-planeSurf-ringArea": cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[1][1])
EOF
}
for arm in ${CASES:-normal oblique rot32}; do
  if [ $arm = rot32 ]; then src=$(pickrot); else src=$(pickpyr 20 $arm); fi
  c=$T/P20$arm; mkdir -p $c; cp -r $src/0 $src/constant $src/system $c/
  dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
  foamDictionary -entry endTime -set "$(python3 -c "print(repr(3*$dt))")" $c/system/controlDict > /dev/null
  sed -i "s/^\(\s*\)sweptVolumeCorrection\s.*/\1sweptVolumeCorrection scale;/" $c/system/fvSolution
  for kv in ${KEYS//+/ }; do sed -i "s/^\(\s*\)reconstructionScheme davof;/&\n\1${kv%%=*} ${kv#*=};/" $c/system/fvSolution; done
  (cd $c && advectorVoF > log.advectorVoF 2>&1; echo "== P20 $arm rc=$? $(grep -c '^Time = ' log.advectorVoF) steps")
  grep 'sweptPreTestCheck:' $c/log.advectorVoF | grep -v shortcut | head -8
  grep 'sweptPreTestCheck:' $c/log.advectorVoF | grep shortcut | head -6
  grep -o 'pre-test[^;]*\|cells: pieces [^;]*\|beyond layer [^,]*\|degenerate [0-9]*' $c/log.advectorVoF | tail -4
  echo "-- boundary patches:"; grep -A3 -E '^\s+[a-zA-Z]+\s*$' $c/constant/polyMesh/boundary | grep -E 'type|nFaces|startFace|^\s+[a-zA-Z]+\s*$' | paste - - - - | head -8
  echo "-- U at the boundary:"; grep -B1 -A4 'inlet\|outlet\|left\|right' $c/0/U | head -20
done
echo DBG-DONE
