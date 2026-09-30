#!/bin/bash
# The leia CI, runnable anywhere: build with OpenFOAM, check every declared target, run the
# one-way 2Dtranslation case, check its numbers. The same script runs in the CI container
# (.github/workflows/build.yml) and on a workstation.
#
# Usage: .github/scripts/ci-build-and-smoke.sh [<OpenFOAM etc/bashrc>]
#        default: $LEIA_FOAM_BASHRC (set in the CI image), else $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc
#
# No `set -e`: OpenFOAM's etc/bashrc does not survive errexit (CLAUDE.md, shell traps), so every
# step checks its own exit code and the script stops at the first failure.
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# A failure prints one line and, under GitHub Actions, an error annotation, which the public
# check-runs API returns without login (the job logs need admin rights).
fail() {
    echo "CI: FAIL, $1"
    if [ -n "$GITHUB_ACTIONS" ]; then
        detail=""
        [ -n "$2" ] && [ -f "$2" ] && detail="%0A$(tail -n 30 "$2" | sed 's/%/%25/g' | awk '{printf "%s%%0A", $0}')"
        echo "::error title=leia CI::$1$detail"
    fi
    exit 1
}
BASHRC="${1:-${LEIA_FOAM_BASHRC:-$HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc}}"
[ -f "$BASHRC" ] || { echo "CI: no OpenFOAM bashrc at $BASHRC"; exit 2; }
# Clear the positional parameters first: OpenFOAM's etc/bashrc evaluates its arguments and
# sources any readable file among them, so a path in $1 made it source itself until the
# stack overflowed (MEASURED 2026-09-30: bash segmentation fault, exit 139).
set --
. "$BASHRC"
cd "$ROOT" || exit 2
echo "CI: OpenFOAM-$WM_PROJECT_VERSION ($WM_OPTIONS), leia $(git describe --always 2>/dev/null || echo unknown)"

# 1. Build. The full log is kept; its tail is printed.
export WM_NCOMPPROCS="${WM_NCOMPPROCS:-$(nproc)}"
./Allwmake > log.Allwmake 2>&1
rc=$?
tail -n 25 log.Allwmake
[ "$rc" -eq 0 ] || fail "./Allwmake rc=$rc" log.Allwmake

# 2. Every target of every Make/files exists (wmake does not always fail on one broken app).
. ./etc/leia-env.sh || fail "etc/leia-env.sh"
python3 etc/leia-check-build.py > log.checkBuild 2>&1; rc=$?; cat log.checkBuild
[ "$rc" -eq 0 ] || fail "a declared target is missing" log.checkBuild

# 3. Smoke run: the one-way translation, serial (cases/2Dtranslation/Allrun.sh, N = 128).
cd "$ROOT/cases/2Dtranslation" || exit 2
sh ./Allclean > /dev/null 2>&1
timeout 900 bash ./Allrun.sh   # the case scripts are not executable in git (mode 644)
rc=$?
state="$("$ROOT/workflow/scripts/foam_log_state.sh" log.leiaSemiLagrangeLevelSetFoam)"
echo "CI: Allrun.sh rc=$rc; solver log: $state"
case "$state" in
    COMPLETED*) ;;
    *) fail "the smoke run did not complete: $state" log.leiaSemiLagrangeLevelSetFoam ;;
esac

# 4. Its numbers. E_GEOM_ALPHA_REL = 2 at t = 0 shows that the exact end reference is in use (the
# start and end circles are disjoint); the end values must stay near the measured ones.
python3 - > log.ciNumbers 2>&1 <<'PY'
import csv, sys
rows = list(csv.DictReader(open("leiaSemiLagrangeLevelSetFoam.csv")))
g0, gT = float(rows[0]["E_GEOM_ALPHA_REL"]), float(rows[-1]["E_GEOM_ALPHA_REL"])
v0, vT = float(rows[0]["E_VOL_ALPHA_REL"]), float(rows[-1]["E_VOL_ALPHA_REL"])
t = float(rows[-1]["TIME"])
checks = [
    ("end time 0.5", abs(t - 0.5) < 1e-9),
    ("E_GEOM_ALPHA_REL(0) = 2 (the end reference is in use)", abs(g0 - 2.0) < 1e-6),
    ("E_VOL_ALPHA_REL(0) < 1e-12", v0 < 1e-12),
    ("E_GEOM_ALPHA_REL(T) < 3e-3", gT < 3e-3),
    ("E_VOL_ALPHA_REL(T) < 1e-3", vT < 1e-3),
    ("E_BOUND_ALPHA = 0", max(float(r["E_BOUND_ALPHA"]) for r in rows) < 1e-12),
]
print(f"CI: {len(rows) - 1} steps; E_GEOM_ALPHA_REL {g0:.6f} -> {gT:.4e}; E_VOL_ALPHA_REL {v0:.2e} -> {vT:.4e}")
for name, ok in checks:
    print(f"CI: {'PASS' if ok else 'FAIL'}  {name}")
sys.exit(0 if all(ok for _, ok in checks) else 1)
PY
rc=$?; cat log.ciNumbers
[ "$rc" -eq 0 ] || fail "the smoke numbers" log.ciNumbers
sh ./Allclean > /dev/null 2>&1
echo "CI: PASS"
