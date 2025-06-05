#!/bin/bash

#------------------------------------------------------------------------------

h=$(awk '/\<height\>/{print $2}' system/blockMeshDict | tr -d ';')
export h
export w=$(echo "scale=3 ; 40 / $h" | bc)

echo "h=$h,w=$w"

envsubst <blockMeshDict.template> system/blockMeshDict

#------------------------------------------------------------------------------
