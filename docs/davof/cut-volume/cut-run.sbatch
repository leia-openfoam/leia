#!/bin/bash
#SBATCH -J cut-vol
#SBATCH -A special00004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=4000
#SBATCH -t 04:00:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
# The cell cut on warped faces (REPORT-swept.md): the meshes of the CCS unit gate (hex1, pert, poly; built as in
# ccs_unit.sbatch) plus the unperturbed box, the test application davofTestCutVolume on each.
cd "${SLURM_SUBMIT_DIR:?}" || exit 1
ROOT=$PWD/leia
module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. $ROOT/etc/leia-env.sh
U=$PWD/cases-${SLURM_JOB_ID}
mkdir -p $U
EXTRA=$PWD/ccs_failures_20_per_mesh.txt
hdr() { printf 'FoamFile\n{\n version 2.0;\n format ascii;\n class dictionary;\n object %s;\n}\n' "$1"; }
mkcase() {  # name
  d=$U/$1; mkdir -p $d/system $d/constant
  { hdr controlDict; printf 'application x;\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\nendTime 0;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n'; } > $d/system/controlDict
  { hdr fvSchemes; printf 'ddtSchemes{default Euler;}\ngradSchemes{default Gauss linear;}\ndivSchemes{default none;}\nlaplacianSchemes{default Gauss linear corrected;}\ninterpolationSchemes{default linear;}\nsnGradSchemes{default corrected;}\n'; } > $d/system/fvSchemes
  { hdr fvSolution; printf 'solvers{}\n'; } > $d/system/fvSolution
}
block() {  # name nx ny nz lz
  mkcase $1
  { hdr blockMeshDict; printf 'scale 1;\nvertices ((0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 %s) (1 0 %s) (1 1 %s) (0 1 %s));\nblocks (hex (0 1 2 3 4 5 6 7) (%s %s %s) simpleGrading (1 1 1));\nedges ();\nboundary (walls {type patch; faces ((0 4 7 3) (2 6 5 1) (1 5 4 0) (3 7 6 2) (0 3 2 1) (4 5 6 7));});\nmergePatchPairs ();\n' $5 $5 $5 $5 $2 $3 $4; } > $U/$1/system/blockMeshDict
  (cd $U/$1 && blockMesh > log.blockMesh 2>&1; echo "blockMesh $1 rc=$?")
}
block hex1 1 1 1 1
block pert 16 16 32 2
block box 16 16 32 2
(cd $U/pert && leiaPerturbMesh > log.perturb 2>&1; echo "perturb rc=$?"; [ -d 0/polyMesh ] && { cp -r 0/polyMesh/ constant/; rm -rf 0/polyMesh; }; ls constant/polyMesh | head -3)
mkcase poly
cp $ROOT/cases/davof/isoShear3D/box1x1x2-one.stl $U/poly/box.stl
{ hdr meshDict; printf 'maxCellSize 0.0787;\nsurfaceFile "box.stl";\n'; } > $U/poly/system/meshDict
(cd $U/poly && pMesh > log.pMesh 2>&1; echo "pMesh rc=$?"; checkMesh -constant 2>&1 | grep -i "cells:\|polyhedra\|hexahedra\|prisms\|tetra\|Mesh OK\|Failed")
for t in hex1 pert poly box; do
  (cd $U/$t && davofTestCutVolume -tag $t -nCells 200 -nRandom 100 -extra $EXTRA > ../../cutvol-${SLURM_JOB_ID}-$t.log 2>&1; echo "cutvol $t rc=$?") &
done
wait
echo "CUT-RUN-DONE $(date)"