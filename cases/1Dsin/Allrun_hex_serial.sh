#!/usr/bin/env bash
set -o verbose

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

restore0Dir
touch "$(basename ${PWD}).foam" && \
set -o errexit
blockMesh 


leiaSetFields 
leiaToyProblem

