#!/usr/bin/env bash

readonly SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
readonly DIR_PATH="$(dirname $SCRIPT_PATH)"

mkdir -p plots
study_cases.py . | xargs -I{} bash -c "cd {} && python ${DIR_PATH}/plot.py && cp plot.png ../plots/{}.png"
cp study_*.json plots/
