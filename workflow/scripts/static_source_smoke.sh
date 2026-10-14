#!/bin/bash
# The laptop smoke of the static source test (config/staticSource/, STATUS 11.27): the 4-rank
# check of CLAUDE.md before the cluster. Every study config of stage A and B runs at N = 64 to
# tau = 2 (256 steps of 4096 cells, so the RAM use is negligible) at np 4; Control, GcNone and
# GcBand also run serially, for the seam check (serial against np 4). A pace run measures the
# seconds per step at N = 128 and 256 for the solve_runtime of the cluster configs. The dead-zone
# run is linearQ in the band on the signed distance at N = 256 to tau = 1: with minIter 0 of the
# psi solver this case never changed psi (STATUS 11.27, the void data of commit 010b553a).
#
#   workflow/scripts/static_source_smoke.sh <smoke dir> [cores]
#
# The derived configs are written into <smoke dir>/configs (copies of the study configs with the
# smoke axes; no --config merging). Every snakemake call stops at aggregate and runs under
# timeout; the calls are sequential (one .snakemake per worktree). Then:
#   static_source_metrics.py selftest | check-init | iterations | null | seam | analyse | pace
# No `set -e`/`-u`: OpenFOAM's etc/bashrc does not survive them.
smoke=$1
cores=${2:-12}
[ -n "$smoke" ] || { echo "usage: $0 <smoke dir> [cores]" >&2; exit 2; }
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$smoke/configs" || exit 1
smoke=$(cd "$smoke" && pwd)
set --
source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
. "$ROOT/etc/leia-env.sh" || exit 1
cd "$ROOT" || exit 1
export PATH="$HOME/.local/bin:$PATH"

python3 - "$smoke" <<'PY' || exit 1
import os, sys, yaml
smoke = sys.argv[1]
def derive(name, suffix, **over):
    d = yaml.safe_load(open(f"config/staticSource/{name}.yaml"))
    d["study_name"] = f"{name}_{suffix}"
    d["studies_dir"] = os.path.join(smoke, "studies")
    d["solve_runtime"] = 30
    ax = d["axes_override"]
    ax["N_CELLS"] = over.pop("N_CELLS", [64])
    ax["END_TIME"] = over.pop("END_TIME", [2])
    ax["WRITE_EVERY_TAU"] = over.pop("WRITE_EVERY_TAU", [1])
    for k, v in over.pop("axes", {}).items():
        ax[k] = v
    d.update(over)
    yaml.safe_dump(d, open(os.path.join(smoke, "configs", f"{d['study_name']}.yaml"), "w"), sort_keys=False)
for name in ("staticSourceControl2D", "staticSourceNullR2D", "staticSourceNullStrain2D",
             "staticSourceGcNone2D", "staticSourceGcBand2D", "staticSourceBeta2D"):
    derive(name, "np4")
for name in ("staticSourceControl2D", "staticSourceGcNone2D", "staticSourceGcBand2D"):
    derive(name, "serial", mode="serial", np=1)
derive("staticSourceGcBand2D", "pace", N_CELLS=[128, 256], END_TIME=[1], WRITE_EVERY_TAU=[1],
       axes={"INIT_FIELD": ["alg"], "GC_LAW": ["linearQ", "softWall"]})
derive("staticSourceGcBand2D", "deadzone", N_CELLS=[256], END_TIME=[1], WRITE_EVERY_TAU=[1],
       axes={"INIT_FIELD": ["sdf"], "GC_LAW": ["linearQ"]})
PY

for cfg in "$smoke"/configs/*.yaml; do
    name=$(basename "$cfg" .yaml)
    echo "$(date +%T) start $name"
    timeout 3h snakemake --workflow-profile profiles/local --configfile "$cfg" --until aggregate \
        --cores "$cores" --resources tasks="$cores" --nolock > "$smoke/$name.log" 2>&1
    echo "$(date +%T) end $name rc=$?"
    sleep 3    # the concurrent-start race of snakemake in one worktree (process rules)
done
echo "SMOKE-RUNS-DONE"
