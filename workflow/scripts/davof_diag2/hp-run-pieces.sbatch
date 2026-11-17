#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
# HPCToolkit 2025.0.1 profiles of the pieces path (diag2, 2026-11-16), deformation3D / shear2D, OF Co 0.5, the cases of the
# studies diag2_c_* (pieces), the stored standing plicRDF and ccsB_*_iso (Stage B), run with the profiling build prof/leia.
#   sbatch -J d2-hp-<MODE> -n NP -t LIMIT --export=ALL,MODE=serial|par8|par64|s2d  workflow/scripts/davof_diag2/hp-run-pieces.sbatch
# Run spec <method>:<family>:<N>:<steps>; the first 20 steps of the pieces serial run are also run with the production build.
W=/work/scratch/tm83tomy/davof/agents/diag2
P=$W/prof
S=$W/leia/studies/davof
SB=/work/scratch/tm83tomy/davof/agents/ccs/leia/studies/davof
K=/work/scratch/tm83tomy/davof/agents/enright/kiss/leia/studies/davof
HP=/home/tm83tomy/spack/opt/spack/linux-x86_64/hpctoolkit-2025.0.1-ky3i4lswv3ynbyvk6ygnunxtn6tsfa67/bin
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
: ${MODE:?}
NP=$SLURM_NTASKS
T=$P/runs/$MODE
rm -rf $T; mkdir -p $T
echo "profiled build $(git -C $P/leia log --oneline -1); node $(hostname) $(lscpu | grep -m1 'Model name' | sed 's/  */ /g'), cores $(nproc --all), job $SLURM_JOB_ID, np $NP"
echo "$($HP/hpcrun --version | head -1); perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid)"
case $MODE in
  serial) RUNS="pieces:deformation3D:64:100 rdf:deformation3D:64:100 iso:deformation3D:64:100"; IDENT=1;;
  par8)   RUNS="pieces:deformation3D:64:200 rdf:deformation3D:64:200"; IDENT=0;;
  par64)  RUNS="pieces:deformation3D:128:200 rdf:deformation3D:128:200"; IDENT=0;;
  s2d)    RUNS="pieces:shear2D:128:500"; IDENT=0;;
  *) echo "unknown MODE"; exit 1;;
esac
# the timers: CPUTIME and REALTIME together if accepted
EV="-e CPUTIME@1000"
if $HP/hpcrun -e CPUTIME@1000 -e REALTIME@1000 -o $T/.evtest /bin/true > $T/.evtest.log 2>&1; then EV="$EV -e REALTIME@1000"; else echo "REALTIME with CPUTIME not accepted: $(tail -1 $T/.evtest.log)"; fi
echo "hpcrun events: $EV"
pick() { /usr/bin/python3 - "$1" "$2" "$3" "$S" "$SB" "$K" <<'EOF'
import csv, glob, sys
m, fam, N, S, SB, K = sys.argv[1:7]
pre = "isoDeformation3D" if fam == "deformation3D" else "isoShear2D"
if m == "pieces":
    print(glob.glob("%s/diag2_c_%s_N%s/%s_00001" % (S, fam, N, pre))[0]); sys.exit()
d = "%s/standing_%s_N%s" % (K, fam, N) if m == "rdf" else "%s/ccsB_%s_N%s_iso" % (SB, fam, N)
db = glob.glob(d + "/*_database.csv")[0]
for r in csv.DictReader(open(db)):
    if r["TARGET_ADVECTIVE_CO"] != "0.5": continue
    if m == "rdf" and (r["TPF_ROW"] != "plicRDF" or r["TPF_REPETITION"] != "1"): continue
    if m == "iso" and r["TPF_ROW"] != "isoAlpha": continue
    print(glob.glob("%s/*_%s" % (d, r["index"]))[0]); break
EOF
}
prep() { # method family N steps tag -> case dir
  local m=$1 fam=$2 N=$3 steps=$4 tag=$5 src c=$T/$5
  src=$(pick $m $fam $N); [ -d "$src" ] || { echo "no source for $m $fam $N" >&2; return 1; }
  mkdir -p $c; cp -r $src/system $src/constant $c/
  if [ $NP -gt 1 ]; then
    for p in $src/processor*; do mkdir -p $c/$(basename $p); cp -r $p/0 $p/constant $c/$(basename $p)/; done
    [ $(ls -d $src/processor* | wc -l) -eq $NP ] || { echo "decomposition of $src is not $NP" >&2; return 1; }
  else
    cp -r $src/0 $c/
  fi
  local dt=$(foamDictionary -entry deltaT -value $c/system/controlDict)
  local n=$steps
  foamDictionary -entry endTime -set "$(python3 -c "print(repr($n*$dt))")" $c/system/controlDict > /dev/null
  foamDictionary -entry writeControl -set timeStep $c/system/controlDict > /dev/null
  foamDictionary -entry writeInterval -set $n $c/system/controlDict > /dev/null
  echo "$src" > $c/.source
  echo $c
}
csv=$T/timing.csv
echo "mode,np,method,family,N,steps,hpcrun,wall_s,n_steps,mean_step_s,median_step_s,first20_mean_s,node" > $csv
par=""; [ $NP -gt 1 ] && par="-parallel"
logstat() { /usr/bin/python3 - "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" >> $csv <<'EOF'
import sys, re, os, statistics
c, m, fam, N, steps, hp, wall, mode, np_ = sys.argv[1:10]
log = open(c + "/log.solver").read()
ct = re.findall(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s(?:\s+timeStepTime = (\S+) s)?", log)
st = [float(x[2]) for x in ct if x[2]]
n = len(re.findall(r"\nTime = ", log))
print(",".join(str(x) for x in (mode, np_, m, fam, N, steps, hp, wall, n, "%.4f" % statistics.mean(st) if st else "", "%.4f" % statistics.median(st) if st else "", "%.4f" % statistics.mean(st[:20]) if st else "", os.uname().nodename)))
EOF
}
# ---- the identity check: the production build, the first 20 steps of the pieces serial case, and the profiled build without hpcrun
if [ $IDENT = 1 ]; then
  c=$(prep pieces deformation3D 64 20 ident_prod)
  ( . $W/leia/etc/leia-env.sh && cd $c && srun --ntasks=1 --cpu-bind=cores davofAdvectorFoam > log.solver 2>&1 ); echo "ident_prod rc=$? $(command -v davofAdvectorFoam)"
  logstat $c pieces deformation3D 64 20 prodbuild 0 $MODE 1
  c=$(prep pieces deformation3D 64 20 ident_prof)
  ( . $P/leia/etc/leia-env.sh && cd $c && srun --ntasks=1 --cpu-bind=cores davofAdvectorFoam > log.solver 2>&1 ); echo "ident_prof rc=$?"
  logstat $c pieces deformation3D 64 20 profbuild 0 $MODE 1
fi
. $P/leia/etc/leia-env.sh || exit 1
echo "davofAdvectorFoam $(command -v davofAdvectorFoam)"
for spec in $RUNS; do
  IFS=: read m fam N steps <<< "$spec"
  tag=${m}_${fam}_N${N}_${steps}
  c=$(prep $m $fam $N $steps $tag) || { echo "FAILED prep $tag"; continue; }
  echo "=== $tag from $(cat $c/.source) $(date)"
  t0=$(date +%s.%N)
  (cd $c && srun --ntasks=$NP --cpu-bind=cores --kill-on-bad-exit=1 --job-name=$tag $HP/hpcrun $EV -o $c/hpctoolkit-measurements davofAdvectorFoam $par > log.solver 2>&1); rc=$?
  t1=$(date +%s.%N)
  logstat $c $m $fam $N $steps hp "$(python3 -c "print('%.2f' % ($t1-$t0))")" $MODE $NP
  echo "$tag rc=$rc wall $(python3 -c "print('%.1f' % ($t1-$t0))") s; last: $(grep -o 'timeStepTime = [^ ]* s' $c/log.solver | tail -1)"
  [ $rc -eq 0 ] || { echo "FAILED $tag"; tail -20 $c/log.solver; }
  if [ $NP -gt 1 ]; then
    tl=$(ls $c/processor0 | grep -v constant | sort -g | tail -1)
    /usr/bin/python3 $W/leia/workflow/scripts/davof_prof/rankcells.py $c $tl > $c/rankcells.txt 2>&1; head -2 $c/rankcells.txt | cut -c1-300
  fi
done
echo "=== analysis $(date)"
for c in $T/*/; do
  M=$c/hpctoolkit-measurements; [ -d $M ] || continue
  ls $M | grep -c '\.hpcrun$' | xargs echo "$(basename $c): hpcrun files"
  $HP/hpcstruct -j $(nproc) -c $P/struct-cache $M > $c/log.hpcstruct 2>&1 || echo "hpcstruct failed for $c"
  $HP/hpcprof -j $(nproc) -o $c/hpctoolkit-database $M > $c/log.hpcprof 2>&1 || echo "hpcprof failed for $c"
  ls $c/hpctoolkit-database 2>/dev/null | tr '\n' ' '; echo
done
cat $csv
echo "HP-RUN-DONE $MODE $(date)"
