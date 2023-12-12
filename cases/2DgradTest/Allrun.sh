#!/usr/bin/env bash

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

set -o verbose

rm -rf 0 && cp -r 0.org 0  
touch "$(basename ${PWD}).foam" && \
set -o errexit

blockMesh 
leiaSetFields 
leiaTestGradScheme

