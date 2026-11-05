#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --constraint=i01
#SBATCH --mem-per-cpu=3000
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The clipping rules of the swept-face flux against the pieces (2026-11-05, fork
# davof-swept-joint df125bb1 on the reference setting 407f969; the note
# leia/docs/davof/davof-swept-joint). The rendered study cases of the swept agent
# (read only) are copied WITHOUT the reference-setting keys (the fork's defaults
# act, as the reference rerun did) and WITH the two keys of a variant:
#   pieces     sweptClip pieces, sweptHexClip false   (the reference)
#   joint      sweptClip joint
#   hybrid     sweptClip hybrid
#   piecesHex  sweptClip pieces, sweptHexClip true    (bit for bit the pieces, expected)
#   hybridHex  sweptClip hybrid, sweptHexClip true
# The fields are written at precision 17 (the bit-for-bit check of the hex clip).
#   MODE=plane : planeTranslation3D_hexSweptWispM9Ccs2, N in NSET (20 40), Co 0.2 / 0.95, both arms and rows,
#                every variant (PAR cases at a time, serial)
#   MODE=cmp   : FLOW=Rot|LineDiag|LineX N=.. CO=0|1: the ring row REPS repetitions of pieces / joint / hybrid in a
#                rotated order, one of piecesHex / hybridHex, the exact row one repetition of every variant; -n NP
#   sbatch -J joint-<tag> -n NP -t LIMIT --export=ALL,MODE=..,FLOW=..,N=..,CO=.. joint-run.sbatch
W=/work/scratch/tm83tomy/davof/agents/joint
S=/work/scratch/tm83tomy/davof/agents/swept/leia/studies/davof
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/TwoPhaseFlow && . etc/davof-env.sh || exit 1
: ${MODE:?}
NP=${SLURM_NTASKS:-1}; REPS=${REPS:-3}; PAR=${PAR:-8}
VARIANTS=${VARIANTS:-"pieces joint hybrid piecesHex hybridHex"}
echo "fork $(git -C $W/TwoPhaseFlow log --oneline -1 | cut -c1-90); leia $(git -C $W/leia log --oneline -1 | cut -c1-70)"
echo "node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), job $SLURM_JOB_ID, MODE $MODE np $NP $(date)"
echo "advectorVoF $(command -v advectorVoF)"
KEYS='sweptWispSurface|sweptWispMatchTol|sweptWispTol|sweptConsistent|sweptVolumeCorrection|sweptScaleToPhi|positioning|positioningTol|positioningMaxIter|anchorSkipConverged|sweptPreTest|sweptStateShortcut|sweptTraceCache|sweptClip|sweptHexClip'
strip() { sed -i -E "/^\s*($KEYS)\s/d" $1/system/fvSolution; local n=$(grep -c -E "^\s*($KEYS)\s" $1/system/fvSolution); [ $n -eq 0 ] || echo "STRIP FAILED in $1: $n keys left"; }
clipOf() { case $1 in pieces|piecesHex) echo pieces;; joint) echo joint;; hybrid|hybridHex) echo hybrid;; esac; }
hexOf() { case $1 in piecesHex|hybridHex) echo true;; *) echo false;; esac; }
addkeys() { # case variant: the two keys after the fluxIntegration line, the write precision
  sed -i "s/^\(\s*\)fluxIntegration\s.*/&\n\1sweptClip $(clipOf $2);\n\1sweptHexClip $(hexOf $2);/" $1/system/fvSolution
  foamDictionary -entry writePrecision -set 17 $1/system/controlDict > /dev/null
  grep -q "sweptClip $(clipOf $2);" $1/system/fvSolution || echo "ADDKEYS FAILED in $1"
}
copy() { # src dst variant
  mkdir -p $2; cp -r $1/system $1/constant $2/
  if [ -d $1/processor0 ]; then for p in $1/processor*; do mkdir -p $2/$(basename $p); cp -r $p/0 $p/constant $2/$(basename $p)/; done; else cp -r $1/0 $2/; fi
  strip $2; addkeys $2 $3; echo "$1" > $2/.source; echo "$3" > $2/.variant
}
metrics() { (cd $1 && python3 $W/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1); }
tok() { /usr/bin/python3 -c "import json,sys; t=json.load(open(sys.argv[1]+'/case_params.json'))['tokens']; print(' '.join(str(t.get(k,'')) for k in sys.argv[2:]))" "$@"; }
case $MODE in
plane)
  T=$W/runs/plane${PTAG:-}; rm -rf $T; mkdir -p $T; NSET=${NSET:-"20 40"}
  # (one assignment per local: a single "local a=$1 b=$a" expands every word before it assigns)
  run1() { local src=$1; local v=$2; local c=$T/$v/$(basename $src); copy $src $c $v; (cd $c && advectorVoF > log.advectorVoF 2>&1; echo "$v $(basename $c) rc=$? steps $(grep -c '^Time = ' log.advectorVoF) $(grep -o 'ClockTime = [^ ]* s' log.advectorVoF | tail -1)"); metrics $c; }
  export -f run1 copy strip addkeys clipOf hexOf metrics; export T W KEYS
  for c in $S/planeTranslation3D_hexSweptWispM9Ccs2/planeTranslation3D_0*; do
    [ -f $c/case_params.json ] || continue
    read N co row arm <<< "$(tok $c N_CELLS TARGET_ADVECTIVE_CO TPF_ROW TPF_VELOCITY_ARM)"
    N=${N%.*}; [[ " $NSET " == *" $N "* ]] || continue
    case $co in 0.2|0.95) for v in $VARIANTS; do echo "$c $v"; done;; esac
  done | xargs -P $PAR -L 1 bash -c 'run1 $0 $1'
  echo "plane: $(ls -d $T/*/planeTranslation3D_0* | wc -l) runs"
  ;;
cmp)
  : ${FLOW:?} ${N:?} ${CO:?}
  case $FLOW in Rot) END=1;; *) END=0.45;; esac
  T=$W/runs/cmp/${FLOW}_N${N}_co${CO}; [ -z "$REP0" ] && rm -rf $T; mkdir -p $T
  pick() { /usr/bin/python3 - "$S/sphereRotation3D_hexCmp${FLOW}SweptN$N" "$1" "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/sphereRotation3D_0*")):
    if "." in os.path.basename(c): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[co][1])
EOF
  }
  declare -A SRC; SRC[ring]=$(pick davof-planeSurf-ringArea); SRC[exact]=$(pick davof-planeSurf-exact)
  for m in ring exact; do echo "source $m ${SRC[$m]}"; done
  prep() { local row=$1; local v=$2; local r=$3; local c=$T/${row}_${v}_r$r; copy ${SRC[$row]} $c $v
    local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict); local n=$(python3 -c "print(max(1, round($END/$dt)))")
    foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
    foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
    foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null; echo $c; }
  SUF=${REP0:+_r$REP0}
  csv=$T/timing$SUF.csv; echo "flow,N,np,co_index,row,variant,rep,steps,dt,end_time,clock_s,exec_s,wall_s,median_step_s,node" > $csv
  par=""; [ $NP -gt 1 ] && par="-parallel"
  orders=("pieces joint hybrid" "joint hybrid pieces" "hybrid pieces joint")
  runone() { local row=$1 v=$2 r=$3; local c=$(prep $row $v $r); local t0=$(date +%s.%N)
    (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=${row}_${v}_r$r advectorVoF $par > log.advectorVoF 2>&1); local rc=$?
    local t1=$(date +%s.%N); metrics $c
    /usr/bin/python3 - $c $row $v $r $rc $t0 $t1 "$FLOW" "$N" "$NP" "$CO" >> $csv <<'EOF'
import sys, re, os, statistics
c, row, v, r, rc, t0, t1, flow, N, np_, co = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11]
log = open(c + "/log.advectorVoF").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
steps = len(re.findall(r"\nTime = ", log)); st = [float(x[2]) for x in ct if x[2]]
cd = open(c + "/system/controlDict").read()
dt = re.search(r"\ndeltaT\s+([^;]+);", cd).group(1); et = re.search(r"\nendTime\s+([^;]+);", cd).group(1)
print(",".join(str(x) for x in (flow, N, np_, co, row, v, r, steps, dt, et, ct[-1][1] if ct else "", ct[-1][0] if ct else "", "%.2f" % (t1 - t0),
      "%.4f" % statistics.median(st) if st else "", os.uname().nodename)) + ("" if rc == 0 else ",rc=%d" % rc))
EOF
    echo "$row $v r$r rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s; $(grep -o 'half-space clips [0-9]*' $c/log.advectorVoF | tail -1); $(grep -o -E '(joint|hybrid) clip: [^;]*' $c/log.advectorVoF | tail -1)"; }
  DOHEX=${DOHEX:-1}; DOEXACT=${DOEXACT:-1}   # the hex variants and the exact row (off for the 64-rank cells)
  for r in $(seq ${REP0:-1} ${REP1:-$REPS}); do
    for v in ${orders[$(( (r-1) % 3 ))]}; do runone ring $v $r; done
    if [ $r -eq 1 ]; then
      [ $DOHEX -eq 1 ] && for v in piecesHex hybridHex; do runone ring $v $r; done
      if [ $DOEXACT -eq 1 ]; then
        for v in pieces joint hybrid; do runone exact $v $r; done
        [ $DOHEX -eq 1 ] && for v in piecesHex hybridHex; do runone exact $v $r; done
      fi
    fi
  done
  # the bit-for-bit check of the hex clip: the written fields of pieces against piecesHex, hybrid against hybridHex
  for row in ring exact; do for pair in "pieces piecesHex" "hybrid hybridHex"; do set -- $pair
    a=$T/${row}_$1_r1; b=$T/${row}_$2_r1; last=$(ls -d $a/[0-9]* 2>/dev/null | sort -g | tail -1); [ -n "$last" ] || continue; tdir=$(basename $last)
    if [ -d $a/processor0 ]; then n=0; for p in $a/processor*; do pb=$b/$(basename $p); cmp -s $p/$tdir/alpha.water $pb/$tdir/alpha.water || n=$((n+1)); done; echo "hex check $row $1 vs $2 at $tdir: $n of $(ls -d $a/processor* | wc -l) ranks differ in alpha.water"
    else cmp -s $a/$tdir/alpha.water $b/$tdir/alpha.water && echo "hex check $row $1 vs $2 at $tdir: alpha.water IDENTICAL" || echo "hex check $row $1 vs $2 at $tdir: alpha.water DIFFERS"; fi
  done; done
  sleep 5; sacct -j $SLURM_JOB_ID -n -P -o JobID,JobName,Elapsed,ElapsedRaw,NNodes,NTasks,State > $T/sacct_steps$SUF.txt; cat $T/sacct_steps$SUF.txt
  ;;
*) echo "unknown MODE $MODE"; exit 1;;
esac
echo "JOINT-RUN-DONE $MODE $(date)"
