#!/bin/bash
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=8000
#SBATCH -t 02:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# The pictures of the 3D deformation runs (2026-11-06): for every case of the timed study
# (and the diagnostic controls) the written times T/4, T/2, 3T/4, T are reconstructed
# (decomposed cases: reconstructPar -time <t> -fields '(alpha.water)') and the alpha = 0.5
# surface is drawn by render_alpha_mc.py (marching cubes + matplotlib; the cluster's ParaView
# modules carry no pvpython) into $W/pics/<N>_co<CO>_<method>_t<t>.png.
#   sbatch -J davof-enright-render enright-render.sbatch
W=${ENRIGHT_W:-/work/scratch/tm83tomy/davof/agents/enright}
module purge
module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
cd $W/leia && . etc/leia-env.sh || exit 1
PY=/shared/spapps_2026_01/gcc-11.5.0/python-3.11.14-xvq4x57r4fxsq6nqrtdvjkriwjkwgce4/bin/python3
R=$W/leia/workflow/scripts/davof_enright/render_alpha_mc.py
mkdir -p $W/pics
render() { # case N label
  local c=$1 N=$2 lab=$3
  for t in 0.75 1.5 2.25 3; do
    if [ -d $c/processor0/$t ] && [ ! -f $c/$t/alpha.water ]; then
      (cd $c && reconstructPar -time $t -fields '(alpha.water)' > log.reconstructPar.$t 2>&1) || echo "reconstructPar $c $t failed"
    fi
    [ -f $c/$t/alpha.water ] || continue
    $PY $R $c $t $N $W/pics/${lab}_t$t.png "$lab, t = $t" 2>&1 | tail -1
  done
}
for d in $W/time/N*_co*; do
  N=$(basename $d | sed 's/N\([0-9]*\)_co.*/\1/'); co=$(basename $d | sed 's/.*_co//')
  for m in ring rdf; do [ -d $d/${m}_r1 ] && render $d/${m}_r1 $N N${N}_co${co}_$m; done
done
for d in $W/diag/*_co1 $W/diag/*_co05; do [ -d $d ] && render $d 32 diag_$(basename $d); done
ls -la $W/pics
echo RENDER-DONE
