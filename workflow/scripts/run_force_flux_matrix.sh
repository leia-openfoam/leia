#!/usr/bin/env bash
# Run a materialized translating-droplet force-flux matrix while preserving a
# log and exit status for every model, including models that fail before endTime.
set -u

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "usage: $0 STUDY_DIR [JOBS]" >&2
    exit 2
fi

study_dir=$(realpath "$1")
jobs=${2:-4}
solver=leiaSemiLagrangianLevelSetTwoPhaseFoam

if [[ ! -d "$study_dir" ]]; then
    echo "study directory does not exist: $study_dir" >&2
    exit 2
fi

run_case()
{
    local case_dir=$1
    local status=0

    if [[ -e "$case_dir/log.solve" || -e "$case_dir/run.status" ]]; then
        echo "$(basename "$case_dir"): existing run artifacts; skipped" >&2
        return 0
    fi

    (
        cd "$case_dir" || exit 125
        leiaSetFields -alphaName alpha.water > log.setFields 2>&1 || exit $?
        # A violently unstable model can drive adjustTimeStep toward zero.
        # Bound wall-clock time without hiding the incomplete physical time.
        timeout --signal=TERM 300 "$solver" > log.solve 2>&1
    ) || status=$?

    printf '%s\n' "$status" > "$case_dir/run.status"
    echo "$(basename "$case_dir"): status=$status"
}

export -f run_case
export solver

find "$study_dir" -mindepth 1 -maxdepth 1 -type d \
    -name 'transISTDroplet2D_*' -print0 \
| sort -z \
| xargs -0 -r -n1 -P "$jobs" bash -c 'run_case "$1"' _

