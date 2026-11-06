#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -c 1
#SBATCH --mem-per-cpu=4000
#SBATCH -t 04:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The diagnosis of the fragmentation on the 3D deformation (Tomislav, 2026-11-06, the
# predictions D-1 .. D-6 in REPORT-swept "The 3D Enright deformation"): N 32, Co 0.5, serial,
# the full T = 3, every run a copy of the prepared production case with the keys below set in
# the alpha solver dictionary (fork branch davof-swept-enright):
#   frag        the swept row + sweptFragmentLog (the per-cell history, (1); the transport unchanged)
#   exactCaps   + sweptExactCaps (control a: the exact piece volume where the cap system is not solvable)
#   reinit      + sweptReinitFaceFractions (control c: the face fractions from the previous planes every step)
#   reinitPrism the prism row + sweptReinitFaceFractions (control d)
#   wispTol     wispTol 1e-2 in the reconstruction (control e)
#   exactReinit + sweptExactCaps + sweptReinitFaceFractions (the two together; not a tuning, the
#               cross-check of D-2 and D-3)
# N (default 32) and CO (0 = Co 0.5, 1 = Co 1) select the prepared case; NP > 1 runs the
# decomposed copy in parallel (N 64 only where a control needs it).
#   sbatch -J davof-enright-diag3 --export=ALL,N=32,CO=0 enright-diag.sbatch
W=${ENRIGHT_W:-/work/scratch/tm83tomy/davof/agents/enright}
S=${STUDIES:-$W/leia/studies/davof}
APP=${APP:-advectorVoF}
N=${N:-32}; CO=${CO:-0}; NP=${NP:-1}; TAG=${TAG:-diag3}
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"; cd $W/leia && . etc/leia-env.sh || exit 1
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-90); leia $(git -C $W/leia log --oneline -1 | cut -c1-90)"
pick() { /usr/bin/python3 - "$S/sphereDeformation3D_hex$1N$N" "$2" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/sphereDeformation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[co][1])
EOF
}
SRC_SWEPT=$(pick Swept davof-planeSurf-ringArea); SRC_PRISM=$(pick Prism davof-ringArea)
echo "sources: swept $SRC_SWEPT; prism $SRC_PRISM"
D=$W/$TAG/N${N}_co$CO; mkdir -p $D
run() { # name src keys...
  local name=$1 src=$2; shift 2
  local c=$D/$name; rm -rf $c; mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done; else cp -r $src/0 $c/; fi
  local endT=3
  if [ -n "$STEPS" ]; then endT=$(python3 -c "print(repr($STEPS*$(foamDictionary -entry deltaT -value $c/system/controlDict)))"); fi
  foamDictionary -entry endTime -set $endT $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set 0.75 $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set 0.75 $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set 0.75 $c/system/controlDict > /dev/null
  for kv in "$@"; do
    local k=${kv%%=*} v=${kv#*=}
    if grep -qE "^\s*$k\s" $c/system/fvSolution; then sed -i -E "s/^(\s*)$k\s+[^;]*;/\1$k $v;/" $c/system/fvSolution
    else sed -i -E "s/^(\s*)advectionScheme\s+isoAdvection;/\1advectionScheme isoAdvection;\n\1$k $v;/" $c/system/fvSolution; fi
  done
  echo "$name: $(grep -E '^\s*(fluxIntegration|wispTol|sweptFragmentLog|sweptExactCaps|sweptReinitFaceFractions)\s' $c/system/fvSolution | tr -s ' ' | tr '\n' ' ')"
  local par=""; [ $NP -gt 1 ] && par="-parallel"
  ( cd $c && t0=$(date +%s.%N); srun --ntasks=$NP --cpu-bind=cores --exclusive --job-name=$name $APP $par > log.advectorVoF 2>&1; rc=$?; t1=$(date +%s.%N)
    echo "$name rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s steps $(grep -c '^Time = ' log.advectorVoF) abort $(grep -c ABORT log.advectorVoF) $(grep '^t = ' log.advectorVoF | tail -1 | cut -c1-120)"
    python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1 ) &
}
# RUNS (optional) selects the runs; the two *Classify runs add sweptClassify (the first
# unbounded cell with its faces' fluxes) and the validation line of sweptExactCaps
RUNS=${RUNS:-"frag exactCaps reinit reinitPrism wispTol exactReinit"}
want() { [[ " $RUNS " == *" $1 "* ]]; }
want frag && run frag $SRC_SWEPT sweptFragmentLog=true
want exactCaps && run exactCaps $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true
want reinit && run reinit $SRC_SWEPT sweptFragmentLog=true sweptReinitFaceFractions=true
want reinitPrism && run reinitPrism $SRC_PRISM sweptFragmentLog=true sweptReinitFaceFractions=true
want wispTol && run wispTol $SRC_SWEPT sweptFragmentLog=true wispTol=1e-2
want exactReinit && run exactReinit $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true sweptReinitFaceFractions=true
want exactCapsClassify && run exactCapsClassify $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true sweptClassify=true
want exactReinitClassify && run exactReinitClassify $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true sweptReinitFaceFractions=true sweptClassify=true
# the coverage check (which path sums the ring's pieces to the swept volume): the production
# row and the exact-caps row with sweptCoverageCheck
want covFrag && run covFrag $SRC_SWEPT sweptFragmentLog=true sweptCoverageCheck=true
want covExact && run covExact $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true sweptCoverageCheck=true
# the piece dump (the closure/exact disagreements above 1e-6 V_c) with the fragment log's
# bounding class; and the same on the rotation N 32 case of the comparison (ROT_SRC, STEPS)
want exactCapsDump && run exactCapsDump $SRC_SWEPT sweptFragmentLog=true sweptExactCaps=true sweptPieceDumpTol=1e-6
want rotExactDump && [ -n "$ROT_SRC" ] && run rotExactDump $ROT_SRC sweptFragmentLog=true sweptExactCaps=true sweptPieceDumpTol=1e-6
wait
echo DIAG-DONE
