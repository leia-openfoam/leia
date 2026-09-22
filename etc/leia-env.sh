# etc/leia-env.sh -- source this file AFTER OpenFOAM's etc/bashrc, in every shell
# that builds or runs a leia binary:
#
#     source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && . ./etc/leia-env.sh
#
# It installs and finds THIS clone's binaries in <clone>/platforms/<WM_OPTIONS>/
# {bin,lib} (git-ignored). Two clones never share binaries, so a rebuild in one
# clone cannot change what another clone runs (MEASURED 2026-09-09: a build from
# one clone landed in the shared account default $HOME/OpenFOAM/<user>-v2512 and a
# second clone then ran a library about 200 commits ahead of its own source; see
# STATUS.md section 9 and docs/plan-library-split-and-build-policy.md).
#
# Sourced by: Allwmake, Allwclean, the Snakefile shell helper (every workflow job,
# local or SLURM, after the profile's or config's env_preamble), run-studies.sbatch.
# Idempotent: sourcing it twice gives the same PATH and LD_LIBRARY_PATH.
if [ -z "$WM_OPTIONS" ]; then
    echo "leia-env.sh: WM_OPTIONS is empty; source OpenFOAM's etc/bashrc first" >&2
    return 1 2>/dev/null || exit 1
fi
LEIA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LEIA_ROOT
export WM_PROJECT_USER_DIR="$LEIA_ROOT"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
# Drop every OpenFOAM USER dir that etc/bashrc (the account default) or an earlier
# source put on the search paths, then put this clone's dirs first. OpenFOAM's own
# installation ($WM_PROJECT_DIR) and ThirdParty ($WM_THIRD_PARTY_DIR) are kept:
# the linker resolves transitive library dependencies through LD_LIBRARY_PATH.
_leia_strip() {
    echo "$1" | tr ':' '\n' | awk -v of="$WM_PROJECT_DIR/" -v tp="${WM_THIRD_PARTY_DIR:-/nonexistent}/" \
        -v mine="$LEIA_ROOT/platforms/" '
        length($0) == 0            { next }
        index($0, of) == 1         { print; next }
        index($0, tp) == 1         { print; next }
        index($0, mine) == 1       { next }
        /\/OpenFOAM\/[^\/]+-v2512\/platforms\// { next }
        { print }' | paste -sd:
}
export PATH="$FOAM_USER_APPBIN:$(_leia_strip "$PATH")"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$(_leia_strip "$LD_LIBRARY_PATH")"
unset -f _leia_strip
