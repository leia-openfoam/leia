#!/usr/bin/env bash
set -o verbose

rm -rf 0 && cp -r 0.org 0  
touch "$(basename ${PWD}).foam" && \
set -o errexit
blockMesh 

checkMesh | tee log.checkMesh 
# Fetch minCellVolume from checkMesh
minVol=$(awk '/Min volume/{print $4}' log.checkMesh | sed 's/\.$//')
# Calc deltaX with python, because #eval #calc are buggy with small float
deltaX=$(python3 -c "dX=${minVol}**(1/2); print(f'{dX:e}')")
echo $deltaX
# Set deltaX entry in controlDict
sed -i "/^deltaX\>/s/\<[^ \t]*;$/${deltaX};/" system/controlDict

leiaSetFields 
leiaLevelSetFoam 

