#!/usr/bin/env bash
#
# MEASURE THE SINK.
#
# Every stability exponent quoted for the stationary droplet is a local slope on
# the NET per-step gain A = ln(u_end/u_0)/nSteps, and A is a difference of two
# things we have never separated:
#
#     A(h)  =  -sink(h)  +  amplifier(h)
#
# A changes sign under refinement in BOTH dimensions (3D between R/h 12.7 and
# 15.8, 2D between N 64 and 128), so no power law describes A globally and the
# slope measured just above the root is inflated by the cancellation. Fitting the
# additive form to 4 points is unidentifiable: a constant-sink fit returns
# amplifier ~ h^-1.36 in 3D with a 12.6% residual and DEGENERATES in 2D
# (q -> 0.01, sink ~ source ~ 0.104, 33.6% residual). The way out is to measure
# the sink independently instead of fitting it.
#
# METHOD. Restart an existing arm from its written field at t ~ 0.02 -- so the
# leading eigenmode is already established, no artificial perturbation has to be
# seeded and no mode-selection assumption is made -- and set sigma = 0. With no
# capillary force the interface perturbation cannot feed back into momentum, so
# what remains is exactly the sink: viscous diffusion plus the removal of the
# non-solenoidal part by the pressure projection. The per-step decay rate of
# max|U| over the restart window IS sink(h).
#
# Then amplifier(h) = A(h) + sink(h) becomes a four-point fit rather than a
# two-point slope, and sink(h) is independently falsifiable against nu*k^2*dt:
#   * a geometry-set mode, k ~ 1/R, gives sink ~ dt ~ h^+1.5  (weakens on refinement)
#   * a mesh-set mode,     k ~ 1/h, gives sink ~ h^-2*h^1.5 = h^-0.5  (strengthens)
# The two damping arms already imply sink ~ h^+2.18, which is incompatible with
# the constant sink the growing arms want -- so at least one is wrong, and this
# run says which. Air is the relevant phase: nu_air = 1.53e-5 is 15x nu_water,
# and the currents are strongest on the light side.
#
# THE CONTROL IS NOT OPTIONAL. A restart that perturbs the state, loses the BDF2
# history (OpenFOAM's `backward` restarts from a single level, so the first step
# after a restart is Euler) or mis-reads the band would produce a decay that has
# nothing to do with the sink. Every ladder therefore gets sigma-UNCHANGED
# restart arms whose measured gain must reproduce the parent arm's local gain
# over the same window. If the control fails, the sigma = 0 numbers mean nothing.
#
# WHAT THE OUTPUT MEANS. The solver truncates and rewrites its metrics CSV, and
# createDropletMetricsFile.H takes the liquid volume AT STARTUP as its reference,
# so on a restart `phaseVolumeRelError` is measured against the RESTART state,
# not against t = 0. It is still the right quantity for "did the interface move
# while the velocity decayed". `zeroSetRadialL2` is referenced to the analytic
# sphere/circle and so remains absolute across the restart -- prefer it.
#
# Usage:
#   sink_decay_restart.sh <srcStudyDir> <dstStudyDir> <nSteps> <sigma> <armGlob> [np]
#
# e.g. sink_decay_restart.sh studies/filterOffAmplifier3D studies/sinkDecay3D \
#          3000 0 'stationaryDroplet3D_0000[0-3]' 32
#
# Creates the restart cases and runs them SEQUENTIALLY inside the current
# allocation. Safe to re-invoke: an existing destination arm with a finished log
# is skipped, so a job that hits its wall limit resumes at the next arm.

set -u
set -o pipefail

SRC=${1:?srcStudyDir}
DST=${2:?dstStudyDir}
NSTEPS=${3:?nSteps}
SIGMA=${4:?sigma}
GLOB=${5:?armGlob}
NP=${6:-}

SOLVER=leiaSemiLagrangianLevelSetTwoPhaseFoam

mkdir -p "$DST"
echo "=== sink_decay_restart: $SRC -> $DST  nSteps=$NSTEPS sigma=$SIGMA glob='$GLOB'"

shopt -s nullglob
for src in $SRC/$GLOB; do
    [ -d "$src" ] || continue
    arm=$(basename "$src")
    dst="$DST/$arm"

    # --- liveness/idempotence: never touch an arm that already finished. Test on
    # --- the LOG, not on a post-processing output that only appears at the end.
    if [ -f "$dst/log.$SOLVER" ] && grep -q "^End$" "$dst/log.$SOLVER" 2>/dev/null; then
        echo "--- $arm: already complete, skipping"
        continue
    fi

    np=$NP
    [ -n "$np" ] || np=$(ls -d "$src"/processor* 2>/dev/null | wc -l)
    if [ "$np" -lt 1 ]; then echo "!!! $arm: no processor dirs, skipping"; continue; fi

    # --- the restart instant: the latest written time in processor0, excluding
    # --- `constant` and excluding 0 (we deliberately do not copy t = 0, so that
    # --- `startFrom latestTime` resolves to the field we did copy and no time
    # --- string has to be formatted to matching precision).
    t0=$(ls -d "$src"/processor0/[0-9]* 2>/dev/null | xargs -n1 basename \
         | awk '$1+0 > 0' | sort -g | tail -1)
    if [ -z "$t0" ]; then echo "!!! $arm: no written time > 0, skipping"; continue; fi

    dt=$(python3 -c "import json,sys;print(json.load(open('$src/case_params.json'))['tokens']['MAX_DELTA_T'])")
    tEnd=$(python3 -c "print(repr(float('$t0') + $NSTEPS*float('$dt')))")

    echo "--- $arm: np=$np  restart at t=$t0  dt=$dt  ->  endTime=$tEnd  ($NSTEPS steps)"

    rm -rf "$dst"
    mkdir -p "$dst"
    cp -a "$src/constant" "$src/system" "$dst/"
    [ -f "$src/case_params.json" ] && cp -a "$src/case_params.json" "$dst/"
    for p in "$src"/processor*; do
        pn=$(basename "$p")
        mkdir -p "$dst/$pn"
        cp -a "$p/constant" "$dst/$pn/"
        cp -a "$p/$t0"      "$dst/$pn/"
    done

    # --- sigma. A plain (untemplated) entry in constant/transportProperties.
    sed -i -E "s/^([[:space:]]*sigma[[:space:]]+).*;/\1${SIGMA};/" "$dst/constant/transportProperties"
    grep -E "^[[:space:]]*sigma" "$dst/constant/transportProperties" | sed 's/^/       sigma now: /'

    # --- time control. deltaT is untouched: the parent runs with
    # --- `adjustTimeStep no` and deltaT == MAX_DELTA_T, so the restart advances on
    # --- exactly the parent's time grid and per-step rates are directly comparable.
    sed -i -E "s/^([[:space:]]*startFrom[[:space:]]+).*;/\1latestTime;/"    "$dst/system/controlDict"
    sed -i -E "s/^([[:space:]]*endTime[[:space:]]+).*;/\1${tEnd};/"         "$dst/system/controlDict"
    sed -i -E "s/^([[:space:]]*writeInterval[[:space:]]+).*;/\1${tEnd};/"   "$dst/system/controlDict"

    ( cd "$dst" && mpirun -np "$np" "$SOLVER" -parallel > "log.$SOLVER" 2>&1 )
    rc=$?
    nrow=$( [ -f "$dst/$SOLVER.csv" ] && wc -l < "$dst/$SOLVER.csv" || echo 0 )
    echo "--- $arm: exit $rc, $(grep -c '^Time = ' "$dst/log.$SOLVER" 2>/dev/null) steps, $nrow csv rows"
done
echo "=== sink_decay_restart done"
