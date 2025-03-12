#!/bin/bash
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

PATTERN=$1
SOLVER=$2

# Submits a parameter study to SLURM, cases fit PATTERN and use SOLVER 
for case in $PATTERN*; 
do
    cd $case
    ./Allrun # Make sure the runApplication line is commented out
    sbatch ../$SOLVER.sbatch
    sleep 1
    cd ..
done
