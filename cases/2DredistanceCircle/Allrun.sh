#!/usr/bin/env bash
# Standalone DEBUG runner for the static redistancing gate (single variant).
#
# CANONICAL, reproducible path (all variants, aggregation, figures, tables):
#   snakemake --workflow-profile profiles/local --configfile config/redistanceCircle2D.yaml
# (or `make studies-grl` from the repo root). This script only renders ONE
# variant of the templates with sed for quick interactive debugging.
#
# Usage: ./Allrun.sh [REDISTANCER] [N_CELLS]     (defaults: planeFootWave 64)

. ${WM_PROJECT_DIR}/bin/tools/RunFunctions

set -o errexit

REDISTANCER=${1:-planeFootWave}
N=${2:-64}

for t in system/*.template; do
    sed -e "s/@!N_CELLS!@/$N/g" \
        -e "s/@!END_TIME!@/1/g" \
        -e "s/@!PROFILE!@/tanh/g" \
        -e "s/@!TANH_EPS!@/0.5/g" \
        -e "s/@!TANH_LIMIT!@/0.1/g" \
        -e "s/@!PHASE_INDICATOR!@/detrixheAslam/g" \
        -e "s/@!REDISTANCER!@/$REDISTANCER/g" \
        -e "s/@!REDIST_TRIGGER!@/interval/g" \
        -e "s/@!REDIST_THRESHOLD!@/-1/g" \
        -e "s/@!REDIST_INTERVAL!@/1/g" \
        "$t" > "${t%.template}"
done

rm -rf 0 && cp -r 0.org 0
touch "$(basename ${PWD}).foam"
blockMesh
leiaSetFields
leiaTestRedistance

echo "--- leiaTestRedistance.csv ($REDISTANCER, N=$N):"
cat leiaTestRedistance.csv
