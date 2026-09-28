# etc/leia-env.sh -- source this file AFTER OpenFOAM's etc/bashrc, in every shell
# that builds or runs a leia binary:
#
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . ./etc/leia-env.sh
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
if [ -z "$WM_OPTIONS" ] || [ -z "$WM_PROJECT_VERSION" ]; then
    echo "leia-env.sh: WM_OPTIONS is empty; source OpenFOAM's etc/bashrc first" >&2
    return 1 2>/dev/null || exit 1
fi
LEIA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LEIA_ROOT
# The account default that etc/bashrc put first on the search paths, captured BEFORE
# it is overridden. When an earlier source (this file, or a TwoPhaseFlow clone's
# etc/davof-env.sh) already replaced it by a clone root, fall back to the bashrc's own
# formula (etc/bashrc: WM_PROJECT_USER_DIR="$HOME/$WM_PROJECT/${USER:-user}-$WM_PROJECT_VERSION").
_leia_acct="${WM_PROJECT_USER_DIR:-}"
if [ -z "$_leia_acct" ] || [ -f "$_leia_acct/etc/leia-env.sh" ] || [ -f "$_leia_acct/etc/davof-env.sh" ]; then
    _leia_acct="$HOME/${WM_PROJECT:-OpenFOAM}/${USER:-user}-$WM_PROJECT_VERSION"
fi
export WM_PROJECT_USER_DIR="$LEIA_ROOT"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
# v2512 and v2606 share WM_OPTIONS=linux64GccDPInt32Opt, so a platforms/ built against
# one is loaded under the other without any error from the loader. Allwmake writes the
# stamp and Allwclean removes the directory; refuse a mismatch here (the Snakefile's
# sh() sources this file with `|| exit 1`, so a job fails loudly instead of running
# binaries of the wrong OpenFOAM).
_leia_stamp="$LEIA_ROOT/platforms/$WM_OPTIONS/.openfoam-version"
if [ -f "$_leia_stamp" ] && [ "$(cat "$_leia_stamp")" != "$WM_PROJECT_VERSION" ]; then
    echo "leia-env.sh: $LEIA_ROOT/platforms/$WM_OPTIONS was built with OpenFOAM-$(cat "$_leia_stamp"), but OpenFOAM-$WM_PROJECT_VERSION is sourced; run ./Allwclean && ./Allwmake under the version you want" >&2
    unset _leia_stamp _leia_acct
    return 1 2>/dev/null || exit 1
fi
unset _leia_stamp
# Drop every OpenFOAM USER dir that etc/bashrc (the account default) or an earlier
# source put on the search paths, then put this clone's dirs first. OpenFOAM's own
# installation ($WM_PROJECT_DIR) and ThirdParty ($WM_THIRD_PARTY_DIR) are kept:
# the linker resolves transitive library dependencies through LD_LIBRARY_PATH.
# The account default is matched by its exact prefix (any OpenFOAM version); the
# regular expression additionally removes retired or renamed user dirs of the form
# $HOME/OpenFOAM/<name>-vNNNN/platforms (e.g. curvature-v2512 on the cluster).
_leia_strip() {
    echo "$1" | tr ':' '\n' | awk -v of="$WM_PROJECT_DIR/" -v tp="${WM_THIRD_PARTY_DIR:-/nonexistent}/" \
        -v mine="$LEIA_ROOT/platforms/" -v acct="$_leia_acct/platforms/" '
        length($0) == 0            { next }
        index($0, of) == 1         { print; next }
        index($0, tp) == 1         { print; next }
        index($0, mine) == 1       { next }
        index($0, acct) == 1       { next }
        /\/OpenFOAM\/[^\/]+-v[0-9][0-9][0-9][0-9][^\/]*\/platforms\// { next }
        { print }' | paste -sd:
}
export PATH="$FOAM_USER_APPBIN:$(_leia_strip "$PATH")"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$(_leia_strip "$LD_LIBRARY_PATH")"
unset -f _leia_strip
unset _leia_acct

# cfMesh (pMesh, cartesianMesh) is a third-party mesher built once per machine, like
# OpenFOAM itself, not per clone: $LEIA_CFMESH_DIR, in the layout of a
# WM_PROJECT_USER_DIR. Default: the per-version build $HOME/OpenFOAM/cfmesh-$WM_PROJECT_VERSION
# when it exists (cfmesh-v2606 on the laptop since 2026-09-28), else the shared
# $HOME/OpenFOAM/cfmesh (v2512). Appended AFTER the clone's own directories, so it can
# never shadow a leia binary. MEASURED 2026-09-23: stripping the shared user
# directories also removed the only pMesh on the cluster's PATH, and the polyhedral
# mesh rule failed with `pMesh: command not found` (WP3 gate, STATUS.md 10.3).
if [ -z "$LEIA_CFMESH_DIR" ]; then
    LEIA_CFMESH_DIR="$HOME/OpenFOAM/cfmesh"
    [ -d "$HOME/OpenFOAM/cfmesh-$WM_PROJECT_VERSION/platforms/$WM_OPTIONS/bin" ] && LEIA_CFMESH_DIR="$HOME/OpenFOAM/cfmesh-$WM_PROJECT_VERSION"
fi
export LEIA_CFMESH_DIR
_leia_cf="$LEIA_CFMESH_DIR/platforms/$WM_OPTIONS"
if [ -d "$_leia_cf/bin" ]; then
    # A build stamped for another OpenFOAM version would load the wrong libOpenFOAM;
    # warn (a build without a stamp is trusted, as before).
    if [ -f "$_leia_cf/.openfoam-version" ] && [ "$(cat "$_leia_cf/.openfoam-version")" != "$WM_PROJECT_VERSION" ]; then
        echo "leia-env.sh: WARNING: $LEIA_CFMESH_DIR is built for OpenFOAM-$(cat "$_leia_cf/.openfoam-version"), not $WM_PROJECT_VERSION; polyhedral (pMesh) cases will not run" >&2
    fi
    case ":$PATH:" in *":$_leia_cf/bin:"*) ;; *) export PATH="$PATH:$_leia_cf/bin" ;; esac
    case ":$LD_LIBRARY_PATH:" in *":$_leia_cf/lib:"*) ;; *) export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$_leia_cf/lib" ;; esac
fi
unset _leia_cf

# TwoPhaseFlow (interFlow, libVoF, libsurfaceForces): the geometric-VOF benchmark arm of
# the method-comparison studies (`solver: interFlow`, config/interFlowDroplet{2D,3D}.yaml),
# found on PATH only. Like cfMesh a sibling build, not part of this clone: $LEIA_TPF_DIR,
# default the TwoPhaseFlow clone next to this one (<workspace>/{leia,TwoPhaseFlow}), built
# by its own etc/davof-env.sh into its platforms/$WM_OPTIONS. Appended AFTER leia's own
# directories, so it can never shadow a leia binary; absent, or stamped for another
# OpenFOAM version: skipped.
LEIA_TPF_DIR="${LEIA_TPF_DIR:-$(dirname "$LEIA_ROOT")/TwoPhaseFlow}"; export LEIA_TPF_DIR
_leia_tpf="$LEIA_TPF_DIR/platforms/$WM_OPTIONS"
if [ -d "$_leia_tpf/bin" ]; then
    if [ -f "$_leia_tpf/.openfoam-version" ] && [ "$(cat "$_leia_tpf/.openfoam-version")" != "$WM_PROJECT_VERSION" ]; then
        echo "leia-env.sh: WARNING: $LEIA_TPF_DIR is built for OpenFOAM-$(cat "$_leia_tpf/.openfoam-version"), not $WM_PROJECT_VERSION; not added to PATH" >&2
    else
        case ":$PATH:" in *":$_leia_tpf/bin:"*) ;; *) export PATH="$PATH:$_leia_tpf/bin" ;; esac
        case ":$LD_LIBRARY_PATH:" in *":$_leia_tpf/lib:"*) ;; *) export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$_leia_tpf/lib" ;; esac
    fi
fi
unset _leia_tpf
