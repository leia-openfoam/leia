#!/bin/bash
# The flux-projection gate on Lichtenberg (STATUS 11.25): velocityModel.fluxProjection none against
# helmholtz on IDENTICAL copies of each prepared polyhedral case (same mesh, same 0/), so that the
# cfMesh noise of separate meshes cannot enter the comparison.
#
#   cd <clone> && workflow/scripts/flux_projection_gate.sh submit <config> <studies_dir> <hours>
#
# submit  one prepare job (the workflow up to decompose, profiles/local inside the job), then two
#         solve jobs per case that depend on it; every job id goes to <clone>/.my_jobs.
# prep    (inside the prepare job) the workflow up to decompose
# solve   (inside a solve job) <case dir> <none|helmholtz>: copy the case to <case>_<proj>, set the
#         entry, run leiaSemiLagrangeLevelSetFoam on the job's tasks
# No `set -e`/`-u`: OpenFOAM's etc/bashrc does not survive them (CLAUDE.md, shell traps).
set -o pipefail
cmd=$1; shift
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
[ -f "$ROOT/etc/leia-env.sh" ] || { echo "flux_projection_gate.sh: $ROOT is not a leia clone" >&2; exit 1; }

env_setup() {
    module purge
    module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv
    export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_mtl=^ofi,psm2
    set --
    source "$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc"
    [ -f "$HOME/.leia_env" ] && . "$HOME/.leia_env"
    . "$ROOT/etc/leia-env.sh" || exit 1
    export PATH="$HOME/.local/bin:$PATH"
}

case "$cmd" in
submit)
    cfg=$1; dir=$2; hours=$3
    study=$(python3 -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))['study_name'])" "$cfg")
    np=$(python3 -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))['np'])" "$cfg")
    ncase=$(python3 -c "import yaml,sys; print(len(yaml.safe_load(open(sys.argv[1]))['axes_override']['MAX_CELL_SIZE']))" "$cfg")
    logs="$dir/$study.logs"; mkdir -p "$logs"
    cd "$ROOT" || exit 1
    prep=$(sbatch --parsable -A special00004 -J fluxgate-prep -N 1 -n 1 -c 2 --mem-per-cpu=24000 \
           -t 06:00:00 -o "$logs/prep.%j.out" -e "$logs/prep.%j.err" \
           workflow/scripts/flux_projection_gate.sh prep "$cfg" "$dir" 2>/dev/null | tail -1)
    echo "$prep fluxgate-prep $study" >> "$ROOT/.my_jobs"; echo "prepare job $prep"
    for i in $(seq 0 $((ncase - 1))); do
        c="$dir/$study/3Dshear_$(printf '%05d' "$i")"
        for proj in none helmholtz; do
            j=$(sbatch --parsable -A special00004 -J fluxgate-solve -N 1 -n "$np" --mem-per-cpu=4000 \
                -t "$hours":00:00 --dependency=afterok:"$prep" \
                -o "$logs/solve_$(basename "$c")_$proj.%j.out" -e "$logs/solve_$(basename "$c")_$proj.%j.err" \
                workflow/scripts/flux_projection_gate.sh solve "$c" "$proj" 2>/dev/null | tail -1)
            echo "$j fluxgate-solve $study $(basename "$c") $proj" >> "$ROOT/.my_jobs"
            echo "solve job $j $(basename "$c") $proj"
        done
    done
    ;;
prep)
    cfg=$1; dir=$2
    env_setup
    cd "$ROOT" || exit 1
    echo "prepare $(date) on $(hostname)"
    snakemake --workflow-profile profiles/local --configfile "$cfg" --config studies_dir="$dir" \
        --until decompose --cores 2
    echo "prepare rc=$? $(date)"
    ;;
solve)
    c=$1; proj=$2
    env_setup
    [ -f "$c/processor0/constant/polyMesh/owner" ] || { echo "not prepared: $c" >&2; exit 1; }
    d="${c}_${proj}"
    rm -rf "${d:?}"; cp -a "$c" "$d"; rm -f "$d/leiaSemiLagrangeLevelSetFoam.csv"
    foamDictionary -entry velocityModel/fluxProjection -set "$proj" "$d/system/fvSolution" > /dev/null 2>&1
    cd "$d" || exit 1
    echo "solve $proj $(date) on $(hostname), $SLURM_NTASKS tasks"
    srun --ntasks="$SLURM_NTASKS" --cpu-bind=none leiaSemiLagrangeLevelSetFoam -parallel \
        > log.leiaSemiLagrangeLevelSetFoam 2>&1
    echo "solve rc=$? $(date)"
    "$ROOT/workflow/scripts/foam_log_state.sh" log.leiaSemiLagrangeLevelSetFoam
    ;;
*)
    echo "usage: $0 submit <config> <studies_dir> <hours> | prep <config> <dir> | solve <case> <proj>" >&2
    exit 2
    ;;
esac
