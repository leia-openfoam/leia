#!/bin/bash
#------------------------------------------------------------------------------
# droplet_mechanism_probe.sh
#
# Reproduces the raw data behind the stationary-droplet divergence-mechanism
# figure (sl_droplet_mechanism.png) and discussion in the article.
#
# It runs the stationary 2D water/air droplet with the two-phase SL solver at two
# resolutions -- N=64 (stable) and N=128 (diverging) -- plus a controlled
# confirmatory run: N=128 with the interface FROZEN (semi-Lagrangian advection
# skipped via SL_FREEZE_INTERFACE=1), which holds psi at the clean initial signed
# distance. The solver writes, per step, the parasitic-current amplitude, the
# band curvature error vs the exact 1/R, and min|grad psi| in the band
# (leiaSemiLagrangianLevelSetTwoPhaseFoam.csv; see writeDropletMetrics.H).
#
# The three raw CSVs are then downsampled to the compact, committed inputs of the
# figure:
#   docs/.../sl-level-set-article/data/mechanism/{stable_N64,diverging_N128,frozen_N128}.csv
#
# NOTE: SL_FREEZE_INTERFACE is a diagnostic-only gate in slAlphaEqn.H; the
# production solver never sets it. This probe is documentation of how the
# mechanism data were produced -- it is not part of the standard Snakemake sweep.
#
# Usage (in WSL, OpenFOAM sourced):  bash workflow/scripts/droplet_mechanism_probe.sh
#------------------------------------------------------------------------------
set -e
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
CASE="$REPO/cases/stationaryDroplet2D"
OUT="$REPO/docs/semi-lagrangian-level-set/sl-level-set-article/data/mechanism"
mkdir -p "$OUT"

wmake "$REPO/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam"

runcase () {  # $1=N $2=maxDt $3=endTime $4=tag  ($SL_FREEZE_INTERFACE honoured)
  local N=$1 DT=$2 ET=$3 TAG=$4 R="/tmp/mech_$4"
  rm -rf "$R"; mkdir -p "$R"
  cp -r "$CASE/0.org" "$CASE/constant" "$CASE/system" "$R/"
  cd "$R"; cp -r 0.org 0
  sed -i "s/^n_cells .*/n_cells $N;/;s/@!N_CELLS!@/$N/g" system/blockMeshDict
  sed -i "s/^maxDeltaT .*/maxDeltaT       $DT;/;s/@!MAX_DELTA_T!@/$DT/g" system/controlDict
  sed -i "s/^endTime .*/endTime         $ET;/;s/@!END_TIME!@/$ET/g;s/@!MAX_CO!@/0.1/g" system/controlDict
  blockMesh >/dev/null 2>&1
  leiaSetFields -alphaName alpha.water >/dev/null 2>&1
  leiaSemiLagrangianLevelSetTwoPhaseFoam > log.solver 2>&1 || true
  cd "$REPO"
}

runcase 64  2.12e-5 0.1  N64
runcase 128 7.5e-6  0.1  N128
SL_FREEZE_INTERFACE=1 runcase 128 7.5e-6 0.1 N128frozen

# Downsample to ~600 log-spaced rows each (keeps the committed CSVs small).
python3 - "$OUT" <<'PY'
import csv, sys, math, os
out = sys.argv[1]
srcs = {"stable_N64": "/tmp/mech_N64",
        "diverging_N128": "/tmp/mech_N128",
        "frozen_N128": "/tmp/mech_N128frozen"}
for name, d in srcs.items():
    p = os.path.join(d, "leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
    rows = list(csv.DictReader(open(p)))
    n = len(rows)
    keep = sorted(set([0, n-1] + [int(round(i*(n-1)/599)) for i in range(600)]))
    with open(os.path.join(out, name+".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        for i in keep:
            w.writerow(rows[i])
    print(f"[mechanism] {name}: {n} -> {len(keep)} rows")
PY
echo "mechanism CSVs written to $OUT"
