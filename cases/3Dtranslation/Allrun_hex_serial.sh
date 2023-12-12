#!/usr/bin/env bash

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

set -o verbose

rm -rf 0 && cp -r 0.org 0  
touch "$(basename ${PWD}).foam" && \
set -o errexit
blockMesh 


leiaSetFields 

# Succesfull execution is optional. Purpose: error calculation
./system/init_End/init_End.serial || { echo "Failed initialising endTime state." && exit 1; }

leiaLevelSetFoam 

