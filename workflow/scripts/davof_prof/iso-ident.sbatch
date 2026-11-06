#!/bin/bash
#SBATCH -J davof-iso-ident
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3800
#SBATCH -t 01:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The identity check of the davof-isoAdvector package (Claude, 2026-11-06): the fork branch
# davof-isoAdvector (ISO_SHA) and davof-swept c5c00c5 (SWEPT_SHA) built side by side with the
# default flags (-O3) in $W/iso/{TwoPhaseFlow,TwoPhaseFlow-swept} (worktrees of the main clone);
# 20-step copies (writePrecision 17) of the plane-clipped loop's rows run with both binaries,
# the written fields compared bit for bit (cmp) and by max |difference|, the per-step log
# lines diffed. Rows: R32 / R32x (rotation N 32, OF Co 1, ring / exact normal), D32 (the
# diagonal line N 32, ring), P20n / P20o (the plane gate N 20, Co 0.5, both arms, ring),
# RDF32 (plicRDF, rotation N 32), ISO32 (isoAlpha, rotation N 32, hexComparators), R32p8
# (R32 at 8 ranks). The swept cases get clippingSurface none / fluxIntegration prism /
# fitAnchor faces (the plane-clipped loop, as ref/cmpC); the plicRDF and isoAlpha cases lose
# the keys positioning* (the new binary's plicRDF reads them, the old one ignored them).
# RDF32ccs: plicRDF with positioning ccs on the new binary alone (a smoke test of (b)).
W=/work/scratch/tm83tomy/davof/agents/swept
S=$W/leia/studies/davof
I=$W/iso
: ${ISO_SHA:?} ${SWEPT_SHA:=c5c00c5}
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
echo "node $(hostname), job $SLURM_JOB_ID, ISO_SHA $ISO_SHA SWEPT_SHA $SWEPT_SHA, $(date)"
mkdir -p $I
git -C $W/TwoPhaseFlow fetch -q origin || echo "fetch failed"
build() { # dir sha
  if [ ! -d $1 ]; then git -C $W/TwoPhaseFlow worktree add --detach $1 $2 || { echo "worktree $1 failed"; exit 1; }; fi
  cd $1 && git checkout -q -- . && git checkout -q --detach $2 || { echo "checkout $2 in $1 failed"; exit 1; }
  export WM_NCOMPPROCS=24; unset FOAM_EXTRA_CXXFLAGS
  find . -name 'Allwmake*' -o -name 'Allwclean*' | xargs chmod +x
  rm -rf src/VoF/lnInclude; find src solver apps -path '*Make/linux64*' -name '*.dep' -delete 2>/dev/null
  t0=$(date +%s); ./Allwmake > build.log 2>&1; rc=$?
  echo "build $1 at $(git log --oneline -1 | cut -c1-70) rc=$rc in $(( $(date +%s) - t0 )) s, errors $(grep -c -E 'error:' build.log)"
  [ $rc -eq 0 ] || { grep -B2 -A8 -E 'error:' build.log | grep -v '^g++' | head -40; exit 1; }
  ls -la platforms/*/bin/advectorVoF | cut -c1-120
  cd $W
}
if [ -z "$NO_BUILD" ]; then build $I/TwoPhaseFlow-swept $SWEPT_SHA; build $I/TwoPhaseFlow $ISO_SHA; fi
pick() { /usr/bin/python3 - "$1" "$2" "$3" "${4:-}" "${5:-}" <<'EOF'
import json, glob, os, sys
st, row, co, N, arm = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
cs = []
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c) or not os.path.exists(c + "/case_params.json"): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] != row: continue
    if N and int(float(t["N_CELLS"])) != int(N): continue
    if arm and t.get("TPF_VELOCITY_ARM") != arm: continue
    cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort()
# co: an index into the sorted Co values, or a Co value (with a decimal point): the closest
if "." in co: k = min(range(len(cs)), key=lambda i: abs(cs[i][0] - float(co)))
else: k = int(co)
sys.stderr.write("pick %s %s: Co %s -> %s (%s)\n" % (os.path.basename(st), row, [x[0] for x in cs], cs[k][0], os.path.basename(cs[k][1])))
print(cs[k][1])
EOF
}
T=$I/ident; rm -rf $T; mkdir -p $T
prep() { # src dst steps [plane|rdf]
  mkdir -p "$2"; cp -r "$1/0" "$1/constant" "$1/system" "$2/"
  local dt=$(foamDictionary -entry deltaT -value "$2/system/controlDict")
  foamDictionary -entry endTime -set "$(python3 -c "print(repr($3*$dt))")" "$2/system/controlDict" > /dev/null
  foamDictionary -entry writeControl -set timeStep "$2/system/controlDict" > /dev/null
  foamDictionary -entry writeInterval -set $3 "$2/system/controlDict" > /dev/null
  foamDictionary -entry writePrecision -set 17 "$2/system/controlDict" > /dev/null
  if [ "$4" = plane ]; then
    sed -i -E 's/^(\s*)clippingSurface\s+\w+;/\1clippingSurface none;/; s/^(\s*)fitAnchor\s+\w+;/\1fitAnchor faces;/; s/^(\s*)fluxIntegration\s+\w+;/\1fluxIntegration prism;/' "$2/system/fvSolution"
  elif [ "$4" = rdf ]; then
    sed -i -E '/^\s*(positioning|positioningTol|positioningMaxIter)\s/d' "$2/system/fvSolution"
  fi
  echo "$1" > "$2/.source"
}
prep $(pick $S/sphereRotation3D_hexCmpRotSweptN32 davof-planeSurf-ringArea 1) $T/src/R32 20 plane
prep $(pick $S/sphereRotation3D_hexCmpRotSweptN32 davof-planeSurf-exact 1) $T/src/R32x 20 plane
prep $(pick $S/sphereRotation3D_hexCmpLineDiagSweptN32 davof-planeSurf-ringArea 1) $T/src/D32 20 plane
prep $(pick $S/planeTranslation3D_hexSweptPyr davof-planeSurf-ringArea 0.5 20 normal) $T/src/P20n 20 plane
prep $(pick $S/planeTranslation3D_hexSweptPyr davof-planeSurf-ringArea 0.5 20 oblique) $T/src/P20o 20 plane
prep $(pick $S/sphereRotation3D_hexCmpRotRDFN32 plicRDF 1) $T/src/RDF32 20 rdf
prep $(pick $S/sphereRotation3D_hexComparators isoAlpha 0 32) $T/src/ISO32 20 rdf
cp -r $T/src/R32 $T/src/R32p8
cp -r $T/src/RDF32 $T/src/RDF32ccs; sed -i "s/^\(\s*\)reconstructionScheme plicRDF;/&\n\1positioning ccs;/" $T/src/RDF32ccs/system/fvSolution
for c in R32 R32x D32 P20n P20o RDF32 ISO32; do echo "$c keys: $(grep -E '^\s*(reconstructionScheme|normalModel|clippingSurface|fluxIntegration|fitAnchor|positioning)\s' $T/src/$c/system/fvSolution | tr -s ' ' | tr '\n' ' ')"; done
echo "RDF32ccs keys: $(grep -E '^\s*(reconstructionScheme|positioning)\s' $T/src/RDF32ccs/system/fvSolution | tr -s ' ' | tr '\n' ' ')"
runset() { # A|B
  local tag=$1 clone=$2
  (cd $clone && . etc/davof-env.sh) || exit 1
  cd $clone && . etc/davof-env.sh && cd $W
  echo "=== $tag: $(command -v advectorVoF) ($(git -C $clone log --oneline -1 | cut -c1-60))"
  for c in R32 R32x D32 P20n P20o RDF32 ISO32; do
    cp -r $T/src/$c $T/$tag/$c 2>/dev/null || { mkdir -p $T/$tag; cp -r $T/src/$c $T/$tag/$c; }
    (cd $T/$tag/$c && advectorVoF > log.advectorVoF 2>&1; echo "$tag $c rc=$? $(grep -c '^Time = ' log.advectorVoF) steps, $(grep -o 'After conservative bounding: min(alpha) = [^,]*' log.advectorVoF | tail -1)")
  done
  c=R32p8; mkdir -p $T/$tag; cp -r $T/src/$c $T/$tag/$c
  cat > $T/$tag/$c/system/decomposeParDict <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 8;
method simple;
coeffs { n (2 2 2); }
EOF
  (cd $T/$tag/$c && decomposePar -force > log.decomposePar 2>&1 && srun --ntasks=8 --cpu-bind=cores --kill-on-bad-exit=1 advectorVoF -parallel > log.advectorVoF 2>&1; echo "$tag $c rc=$? $(grep -c '^Time = ' log.advectorVoF) steps"; reconstructPar -latestTime > log.reconstructPar 2>&1)
  if [ $tag = B ]; then c=RDF32ccs; cp -r $T/src/$c $T/$tag/$c; (cd $T/$tag/$c && advectorVoF > log.advectorVoF 2>&1; echo "$tag $c rc=$? $(grep -c '^Time = ' log.advectorVoF) steps; $(grep -m1 'positioning ccs:' log.advectorVoF | cut -c1-160)"); fi
}
runset A $I/TwoPhaseFlow-swept
runset B $I/TwoPhaseFlow
echo "=== identity B (davof-isoAdvector) against A (davof-swept $SWEPT_SHA): the last written time"
cmpcase() { # a b label
  local a=$1 b=$2
  local ta=$(ls -d $a/[0-9]* | grep -v '/0$' | sort -t/ -k2 -g | tail -1); local tb=$b/$(basename $ta)
  [ -d "$tb" ] || { echo "$3: no $tb"; return; }
  local same=0 diff=""
  for f in $(cd $ta && ls | grep -v uniform); do if cmp -s $ta/$f $tb/$f; then same=$((same+1)); else diff="$diff $f"; fi; done
  echo "$3: $same identical field files; differing:${diff:- none}"
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
    print("    %-13s n %d, differing %d, max|diff| %.3e" % (f, len(a), n, max(d)))
EOF
  diff <(grep -E 'After conservative bounding|^davof: [0-9]+ interface|alphaf set by a plane|bounding extended' $a/log.advectorVoF) \
       <(grep -E 'After conservative bounding|^davof: [0-9]+ interface|alphaf set by a plane|bounding extended' $b/log.advectorVoF) > $b/log.diff
  echo "    per-step log lines differing: $(grep -c '^[<>]' $b/log.diff); first: $(grep -m1 '^[<>]' $b/log.diff | cut -c1-160)"
}
for c in R32 R32x D32 P20n P20o RDF32 ISO32 R32p8; do cmpcase $T/A/$c $T/B/$c $c; done
echo "=== RDF32ccs (B) against RDF32 (B): the ccs positioning of plicRDF against vofCutCell"
cmpcase $T/B/RDF32 $T/B/RDF32ccs RDF32ccs
echo "ISO-IDENT-DONE $(date)"
