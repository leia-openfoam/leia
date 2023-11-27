#! /usr/bin/env bash

readonly SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
readonly DIR_PATH="$(dirname $SCRIPT_PATH)"


study_cases.py --from-json $1 | xargs -I{} bash -c "
    cd {} 
    python ${DIR_PATH}/plot_fieldsAtFixedTime.py
    python ${DIR_PATH}/plot_zeroLevelSetEvolution.py
    "
