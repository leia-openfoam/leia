#!/bin/bash
# The zero-flow source gate E1.1 to E1.3 (config/gates/sourceZeroFlow2D.yaml): the SL source step
# alone, at u = 0. Serial, small (2D, N = 100): runs on the laptop.
#
#   workflow/scripts/source_zero_flow_gate.sh <studies_dir> [config] [mu_t_end] [scales]
#
# Defaults: config/gates/sourceZeroFlow2D.yaml (E1), mu t = 40, scales "1 1.5". Gate E2.1 runs
# config/gates/sourceZeroFlowUpwind2D.yaml to mu t = 80.
#
# 1. The workflow renders and initialises the six cases (GC_MU x SL_SOURCE_BAND_CELLS).
# 2. Every case is copied twice: psi0 = d and psi0 = 1.5 d (the internal field scaled; the zero set
#    and alpha do not change). In each copy: u = (0 0 0), adjustTimeStep no, deltaT = 1e-2/mu,
#    endTime = mu_t_end/mu (100 steps per unit of mu t), psi written every 50 steps.
# 3. The kinematic SL solver runs; the log classifier reports the state.
# No `set -e`/`-u`: OpenFOAM's etc/bashrc does not survive them.
set -o pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
dir=$1
config=${2:-config/gates/sourceZeroFlow2D.yaml}
mut=${3:-40}
scales=${4:-"1 1.5"}
[ -n "$dir" ] || { echo "usage: $0 <studies_dir> [config] [mu_t_end] [scales]" >&2; exit 2; }
study=$(sed -n 's/^study_name:[[:space:]]*//p' "$ROOT/$config")
mkdir -p "$dir"
set --
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. "$ROOT/etc/leia-env.sh" || exit 1
cd "$ROOT" || exit 1
snakemake --workflow-profile profiles/local --configfile "$config" \
    --config studies_dir="$dir" --until preprocess --cores 4 > "$dir/prepare.log" 2>&1 \
    || { echo "prepare failed, see $dir/prepare.log" >&2; exit 1; }
for c in "$dir/$study"/2Dtranslation_000??; do
    mu=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['tokens']['GC_MU'])" "$c/case_params.json")
    for scale in $scales; do
        d="${c}_s${scale}"
        rm -rf "${d:?}"; cp -a "$c" "$d"
        rm -f "$d/leiaSemiLagrangeLevelSetFoam.csv"
        foamDictionary -entry velocityModel/velocity -set "(0 0 0)" "$d/system/fvSolution" > /dev/null
        foamDictionary -entry adjustTimeStep -set no "$d/system/controlDict" > /dev/null
        foamDictionary -entry deltaT -set "$(python3 -c "print(repr(1e-2/$mu))")" "$d/system/controlDict" > /dev/null
        foamDictionary -entry endTime -set "$(python3 -c "print(repr($mut/$mu))")" "$d/system/controlDict" > /dev/null
        foamDictionary -entry writeControl -set timeStep "$d/system/controlDict" > /dev/null
        foamDictionary -entry writeInterval -set 50 "$d/system/controlDict" > /dev/null
        if [ "$scale" != 1 ]; then
            python3 - "$d/0/psi" "$scale" <<'PY'
import re, sys
p, s = sys.argv[1], float(sys.argv[2])
t = open(p).read()
m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", t)
a = m.end(); b = t.index(")", a)
vals = [float(x) * s for x in t[a:b].split()]
open(p, "w").write(t[:a] + "\n" + "\n".join(repr(v) for v in vals) + "\n" + t[b:])
PY
        fi
        ( cd "$d" && timeout 3600 leiaSemiLagrangeLevelSetFoam > log.leiaSemiLagrangeLevelSetFoam 2>&1 )
        echo "$(basename "$d") mu=$mu: $("$ROOT/workflow/scripts/foam_log_state.sh" "$d/log.leiaSemiLagrangeLevelSetFoam" | head -1)"
    done
done
