#!/usr/bin/env bash

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

# if mpi_call is unset, initialise it
if [[ -z ${mpi_call+x} ]]; then
  mpi_call="mpirun -np 4"
fi

set -o verbose

rm -rf 0 && cp -r 0.org 0  
find 0/ -name "*.template" -delete && \
touch "$(basename ${PWD}).foam" && \
set -o errexit
blockMesh 


decomposePar -force
${mpi_call} leiaSetFields -parallel
${mpi_call} leiaLevelSetFoam -parallel

