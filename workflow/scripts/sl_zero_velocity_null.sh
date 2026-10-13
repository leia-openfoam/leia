#!/usr/bin/env bash
# Is the semi-Lagrangian step the IDENTITY at zero velocity? Research instrument.
#
# WHY THIS EXISTS. A weighted least-squares fit is not in general an interpolant, so it need
# not return psi_c when evaluated at the cell centre. If it does not, the scheme moves the
# level set even when nothing moves -- an error source with no mesh geometry in it, and a
# candidate mechanism for the amplification measured in METHOD.md 8.2.
#
# The CODE REFUTES that for the anchored value fits: the right-hand side is the difference
# psi_i - psi_c over a constant-free basis and evaluateRaw starts from psi_c
# (uncachedQuadraticWeightedLeastSquaresReconstruction.C:688,712,749-766). So this gate is
# run to find the stages where the identity BREAKS -- the reprojecting reconstruction, the
# renormalization, the value bounds -- not to confirm the fit. It also retires the claim at
# config/kinematicTranslation2D.yaml:66-68 that "a LSQ fit does not interpolate its nodal
# values", which is false on this code.
#
# WHAT THIS MEASURES. With velocityModel velocity (0 0 0) the departure foot is x_d = x_c
# exactly (pointValueScheme.C:430), so psi^{n+1} must equal psi^n to round-off. The read-out
# is L_INF_E_PSI, which the shared error writer restricts to the interface band at u = 0
# because it zeroes every cell with alpha in {0,1} and unchanged sign
# (advectionErrors.H:23-45). E_GEOM_ALPHA_REL, E_VOL_ALPHA_REL and E_BOUND_ALPHA are
# reported with it, never alone.
#
# VALIDITY GATE. Two traps, both already recorded:
#   - `type none` ABORTS: velocityModel::velocity is notImplemented (velocityModel.C:121-125).
#     Zero velocity is `type translation; velocity (0 0 0)`.
#   - `adjustTimeStep no` is MANDATORY: maxDeltaT is maxCo/(0.5 gMax(sumPhi/V)), i.e. 0/0 at
#     zero flux.
#
# Usage:
#   workflow/scripts/sl_zero_velocity_null.sh <rendered-case> [out-dir] [np] [nsteps]
#
# Render the case first, which is what config/gates/slZeroVelocityNull2D.yaml is for:
#   make studies-one-file CFG=config/gates/slZeroVelocityNull2D.yaml   # --until preprocess
#
# NO `set -u` and NO `set -e`: the OpenFOAM bashrc reads unbound variables and would abort
# the script there with an empty log. Errors are checked explicitly.

SRC=${1:?rendered case (studies/<study>/<case>_00000)}
OUT=${2:-/tmp/slZeroVelocityNull}
NP=${3:-4}
NSTEPS=${4:-50}

SRC=$(cd "$SRC" && pwd)
# shellcheck disable=SC1090
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
[ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"
REPO=$(cd "$(dirname "$0")/../.." && pwd)
# shellcheck disable=SC1091
. "$REPO/etc/leia-env.sh"

RECONS="uncachedQuadraticWeightedLeastSquares quadraticWeightedLeastSquares quadraticTaylor
        linearWeightedLeastSquares linearTaylor defectCorrectedIDW
        signedDistanceQuadraticWeightedLeastSquares"
BOUNDS="none stencilBounds lipschitzCone"
RENORMS="strain none"

rm -rf "$OUT"; mkdir -p "$OUT"
RESULT="$OUT/results.csv"
echo "ARM,RECONSTRUCTION,VALUE_BOUND,RENORMALIZATION,STATE,STEPS,L_INF_E_PSI,E_GEOM_ALPHA_REL,E_VOL_ALPHA_REL,E_BOUND_ALPHA" > "$RESULT"

run_arm () {
    local tag=$1 recon=$2 bound=$3 renorm=$4
    local d="$OUT/$tag"
    rm -rf "$d"; mkdir -p "$d"
    for f in 0.org constant system; do cp -r "$SRC/$f" "$d/"; done

    # 0/ from 0.org plus the pre-processing, NEVER from a finished 0/ (CLAUDE.md).
    ( cd "$d" || exit 1
      rm -rf 0 && cp -r 0.org 0 && rm -f 0/*.template

      foamDictionary -entry velocityModel/type     -set translation   system/fvSolution >/dev/null 2>&1
      foamDictionary -entry velocityModel/velocity -set "(0 0 0)"     system/fvSolution >/dev/null 2>&1
      foamDictionary -entry levelSet/semiLagrangian/reconstruction   -set "$recon"  system/fvSolution >/dev/null 2>&1
      foamDictionary -entry levelSet/semiLagrangian/valueBound       -set "$bound"  system/fvSolution >/dev/null 2>&1
      foamDictionary -entry levelSet/semiLagrangian/renormalization  -set "$renorm" system/fvSolution >/dev/null 2>&1

      # A fixed step: at zero flux the Courant-based dt is 0/0.
      foamDictionary -entry adjustTimeStep -set no      system/controlDict >/dev/null 2>&1
      foamDictionary -entry deltaT         -set 1e-3    system/controlDict >/dev/null 2>&1
      foamDictionary -entry endTime        -set "$(awk "BEGIN{print $NSTEPS*1e-3}")" system/controlDict >/dev/null 2>&1
      foamDictionary -entry writeInterval  -set "$(awk "BEGIN{print $NSTEPS*1e-3}")" system/controlDict >/dev/null 2>&1
      foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict >/dev/null 2>&1

      leiaSetFields > log.setFields 2>&1 || exit 3
      if [ "$NP" -gt 1 ]; then
          decomposePar -force > log.decomposePar 2>&1 || exit 4
          mpirun -np "$NP" leiaSemiLagrangeLevelSetFoam -parallel > log.solver 2>&1
      else
          leiaSemiLagrangeLevelSetFoam > log.solver 2>&1
      fi ) >/dev/null 2>&1

    local state steps csv
    state=$("$REPO/workflow/scripts/foam_log_state.sh" "$d/log.solver" 2>/dev/null | head -1)
    steps=$(grep -c "^Time = " "$d/log.solver" 2>/dev/null)
    csv="$d/leiaSemiLagrangeLevelSetFoam.csv"
    if [ -f "$csv" ]; then
        python3 - "$csv" "$tag" "$recon" "$bound" "$renorm" "$state" "$steps" >> "$RESULT" <<'PY'
import csv, sys
path, tag, recon, bound, renorm, state, steps = sys.argv[1:8]
rows = list(csv.DictReader(open(path)))
def worst(col):
    v = [abs(float(r[col])) for r in rows if r.get(col) not in (None, "")]
    return max(v) if v else float("nan")
print(",".join([tag, recon, bound, renorm, state, steps,
                "%.6e" % worst("L_INF_E_PSI"),
                "%.6e" % worst("E_GEOM_ALPHA_REL"),
                "%.6e" % worst("E_VOL_ALPHA_REL"),
                "%.6e" % worst("E_BOUND_ALPHA")]))
PY
    else
        echo "$tag,$recon,$bound,$renorm,${state:-NO_CSV},${steps:-0},nan,nan,nan,nan" >> "$RESULT"
    fi
    tail -1 "$RESULT"
}

echo "=== A. reconstruction (bound none, renormalization none) ==="
for r in $RECONS; do run_arm "recon_$r" "$r" none none; done
echo "=== B. value bound (production reconstruction) ==="
for b in $BOUNDS; do run_arm "bound_$b" uncachedQuadraticWeightedLeastSquares "$b" none; done
echo "=== C. renormalization (production reconstruction) ==="
for n in $RENORMS; do run_arm "renorm_$n" uncachedQuadraticWeightedLeastSquares none "$n"; done

echo
echo "results: $RESULT"
column -s, -t < "$RESULT"
