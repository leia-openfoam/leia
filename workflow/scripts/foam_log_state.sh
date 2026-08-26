#!/usr/bin/env bash
# THE ONE WAY to classify an OpenFOAM solver log. Never grep a log ad hoc.
#
# WHY THIS FILE EXISTS. Every OpenFOAM solver prints, in its STARTUP banner:
#
#     trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
#
# so any grep for a bare 'Floating point exception' (or bare 'core dumped',
# which SLURM attaches to unrelated kill messages) classifies EVERY healthy run
# as diverged. This exact false positive has now bitten three separate times:
# the original Snakefile divergence classifier (fixed there, with the full
# post-mortem in workflow/Snakefile ~line 190), a dt-bisection probe, and a
# completion-waiter that declared three live studies finished -- the same class
# of mistake that once deleted a live 16-hour interFoam arm. Ad-hoc
# reimplementations keep resurrecting the bug; this helper is the single
# implementation every script and every interactive check must call instead.
#
# States (stdout, first word) and exit codes:
#   COMPLETED       0   OpenFOAM's terminating '^End$' is present
#   DIVERGED        3   solver died: sigFpe/sigSegv handler stack, a genuine
#                       '... (core dumped)' signal line, or FOAM FATAL -- the
#                       run is OVER and the outcome is a RESULT (a blow-up is
#                       data, not garbage)
#   LAUNCH_FAILURE  4   MPI/SLURM could not even start or place the ranks --
#                       the run never happened and must NOT be recorded
#   RUNNING         1   none of the above and the log grew recently
#   STALLED         2   none of the above and the log has not grown for
#                       --stall seconds (default 900) -- check the process
#                       (pgrep/squeue) before touching anything
#   MISSING         5   no log file
#
# Also printed: steps=<n Time-lines> age=<seconds since last write>.
#
# Usage:
#   foam_log_state.sh <log> [--stall SECONDS]
#   foam_log_state.sh <log> --wait [--poll SECONDS]   # block until not RUNNING
#
# The patterns are EXACTLY the Snakefile classifier's; if one changes, change
# both, and say so in the commit.

set -u

DIVERGED_RE='sigFpe::sigHandler|sigSegv::sigHandler|Floating point exception \(core dumped\)|Segmentation fault \(core dumped\)|FOAM FATAL'
LAUNCHFAIL_RE='Unable to create step|not enough slots|More processors requested than permitted|Unable to satisfy cpu bind|CPU binding outside|srun: error: Unable to|error: Task launch|mpirun.*(error|Error)'

LOG=${1:?usage: foam_log_state.sh <log> [--stall S] [--wait [--poll S]]}
shift
STALL=900; WAIT=0; POLL=30
while [ $# -gt 0 ]; do
    case "$1" in
        --stall) STALL=$2; shift 2;;
        --wait)  WAIT=1; shift;;
        --poll)  POLL=$2; shift 2;;
        *) echo "unknown option $1" >&2; exit 64;;
    esac
done

classify() {
    [ -f "$LOG" ] || { echo "MISSING steps=0 age=-1"; return 5; }
    local steps age
    steps=$(grep -c '^Time = ' "$LOG" 2>/dev/null || echo 0)
    age=$(( $(date +%s) - $(stat -c %Y "$LOG") ))
    if grep -qE '^End$' "$LOG"; then
        echo "COMPLETED steps=$steps age=$age"; return 0
    fi
    # LAUNCH_FAILURE before DIVERGED: an aborted launch also spews MPI noise
    if grep -qE "$LAUNCHFAIL_RE" "$LOG" && [ "$steps" -eq 0 ]; then
        echo "LAUNCH_FAILURE steps=$steps age=$age"; return 4
    fi
    if grep -qE "$DIVERGED_RE" "$LOG"; then
        echo "DIVERGED steps=$steps age=$age"; return 3
    fi
    if [ "$age" -gt "$STALL" ]; then
        echo "STALLED steps=$steps age=$age"; return 2
    fi
    echo "RUNNING steps=$steps age=$age"; return 1
}

if [ "$WAIT" -eq 1 ]; then
    while true; do
        classify; rc=$?
        [ $rc -ne 1 ] && exit $rc
        sleep "$POLL"
    done
fi
classify
