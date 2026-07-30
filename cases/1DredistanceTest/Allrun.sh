#!/usr/bin/env bash
# Standalone DEBUG runner (single redistancer, as configured in fvSolution).
# CANONICAL, reproducible path: the snakemake `utilities` rule,
#   snakemake --workflow-profile profiles/local utilities
# All-models comparison on this case: ./Allrun_variants.sh

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

set -o verbose

rm -rf 0 && cp -r 0.org 0  
touch "$(basename ${PWD}).foam" && \
set -o errexit
blockMesh 


leiaSetFields 
leiaTestRedistance

sed -n '21,45p' 1/psi
