#!/bin/bash
#SBATCH -J davof-iso-dec
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --mem-per-cpu=3600
#SBATCH -t 03:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Decisions 2 and 3 in the plane-clipped loop (Tomislav's item 3, 2026-11-07): one cell
# (FLOW Shear | Deformation, N, CO index) of the hexahedral standing matrix, the DAVOF ring +
# ccs row, the nine combinations sweptBoundedFaceFraction {none, scale, recompute} x
# interfaceCellRule {normalOnly, normalAndAlpha, normalAndAlphaReset}, scored as the standing
# matrix (tpf_advection_metrics.py). The keys are appended to the rendered fvSolution by sed
# (a foamDictionary rewrite perturbs a run at round-off, 2026-11-07). Binary: TPF_DIR
# (iso/TwoPhaseFlow-e). Output iso/dec/<FLOW>_N<N>_co<CO>/<bff>-<rule>/. Claude, 2026-11-07.
W=/work/scratch/tm83tomy/davof/agents/swept
I=$W/iso; S=$I/leia/studies/davof; TPF=${TPF_DIR:-$I/TwoPhaseFlow-e}
: ${FLOW:?} ${N:?} ${CO:?}
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
export LEIA_TPF_DIR=$TPF; cd $I/leia && . etc/leia-env.sh || exit 1
export OMPI_MCA_btl=^openib
NP=$SLURM_NTASKS
case $FLOW in Shear|Deformation) END=3; WI=0.75;; *) echo "FLOW $FLOW?"; exit 1;; esac
echo "node $(hostname), job $SLURM_JOB_ID, np $NP, FLOW $FLOW N $N CO $CO, fork $(git -C $TPF log --oneline -1 | cut -c1-70), $(date)"
src=$(/usr/bin/python3 - "$S/iso${FLOW}3D_hexN$N" davof-ringArea "$CO" <<'EOF'
import json, glob, os, sys
st, row, co = sys.argv[1], sys.argv[2], int(sys.argv[3])
cs = []
for c in sorted(glob.glob(st + "/*_0*")):
    if "." in os.path.basename(c) or not os.path.exists(c + "/case_params.json"): continue
    t = json.load(open(c + "/case_params.json"))["tokens"]
    if t["TPF_ROW"] == row: cs.append((float(t["TARGET_ADVECTIVE_CO"]), c))
cs.sort(); print(cs[co][1])
EOF
)
echo "source $src"
T=$I/dec/${FLOW}_N${N}_co${CO}; mkdir -p $T
par=""; [ $NP -gt 1 ] && par="-parallel"
for bff in none scale recompute; do for rule in normalOnly normalAndAlpha normalAndAlphaReset; do
  c=$T/$bff-$rule; rm -rf $c; mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done; else cp -r $src/0 $c/; fi
  foamDictionary -entry endTime -set $END $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/areaFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  foamDictionary -entry functions/volumeFractionError/writeInterval -set $WI $c/system/controlDict > /dev/null
  sed -i "/^\s*positioning\s\+ccs;/a\        sweptBoundedFaceFraction $bff;\n        interfaceCellRule $rule;" $c/system/fvSolution
  grep -c "sweptBoundedFaceFraction $bff;" $c/system/fvSolution > /dev/null || { echo "key not set in $c"; exit 1; }
  echo "$src" > $c/.source; printf "fork_commit=%s\nfork_dir=%s\nbff=%s\nrule=%s\n" "$(git -C $TPF rev-parse --short HEAD)" "$TPF" "$bff" "$rule" > $c/.stamp
  t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 advectorVoF $par > log.advectorVoF 2>&1); rc=$?
  t1=$(date +%s.%N)
  (cd $c && python3 $W/iso/leia/workflow/scripts/tpf_advection_metrics.py --case . --out advectorVoF.csv > log.metrics 2>&1)
  echo "$bff-$rule rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s, $(grep -c '^Time = ' $c/log.advectorVoF) steps; $(grep -o 'ABORT at step [0-9]*, t = [^:]*: [a-zA-Z() ]*' $c/log.advectorVoF | head -1); $(grep 'sum(alpha\*V)' $c/log.advectorVoF | tail -1 | grep -o 'dev = [^%]*%\|min(alpha1) = [^ ]*' | tr '\n' ' ')"
  echo "   $(grep -m1 'bounded face fractions' $c/log.advectorVoF | cut -c1-160)"
  echo "   $(grep 'interface-cell rule' $c/log.advectorVoF | tail -1 | grep -o 'interface-cell rule.*' | cut -c1-200)"
done; done
echo "ISO-DEC-DONE $FLOW N$N co$CO $(date)"
