# Integration identity check

These files test one claim: a merged tree gives the same metric CSVs as the tree before the
merge, for existing studies with default tokens. The first use was the branch
`integration/sl-pde` against the hub commit 87a70b22 (2026-11-13, Lichtenberg).

All jobs run on Lichtenberg. Nothing runs on the laptop.

## Files

| file | function |
|---|---|
| `env.sh` | the environment of every job (modules, OpenFOAM-v2606, `etc/leia-env.sh` of the clone) |
| `build.sbatch` | builds one clone and records the artefacts in `int.build.log` |
| `driver.sbatch` | one-core driver; runs `make studies-one-file` for a list of configs |
| `poly_arm.sbatch` | one polyhedral arm on a given mesh with the binaries of the clone |
| `unit_gate.sbatch` | `leiaTestGhostCellReconstruction` in serial and on 4 ranks |
| `../integration_identity_check.py` | compares two study directories at tolerance 0 |
| `../integration_identity_check.sh` | compares all study pairs of two clones |
| `config/integration/*.yaml` | the mesh study of the polyhedral rung, the unit-gate mesh, the N = 64 ghost-cell rung |

## Sequence

`M` is the merged clone. `B` is the base clone (detached at the commit before the merge).
`S=$M/workflow/scripts/integration`. Append every job id to `.my_jobs` of the clone that
submits it. The child job ids of a driver are in its `int.<tag>.<study>.log`.

1. Build both clones. In each clone root:
   `sbatch --export=ALL,INT_SCRIPTS=$S $S/build.sbatch`.
   Read `int.build.log`: `Allwmake rc=0`, `leia-check-deps rc=0`, new timestamps, and the
   count of the string `ghostCellQuadraticWeightedLeastSquares` (more than 0 in `M`, 0 in `B`).
2. If the profile of `B` names another SLURM account, copy `profiles/slurm` of `B` to an
   untracked directory (`profiles/slurm-int5`), change the account there, and use that
   directory as `PROFILE` for `B`. Do not edit a tracked file of `B`.
3. Run the hexahedral studies and the Eulerian studies in each clone root:
   `sbatch --export=ALL,INT_SCRIPTS=$S,TAG=std,PROFILE=<profile>,CFGS="config/advConv2Dtranslation.yaml config/advConv3DshearHex.yaml config/pdeLevelSet/s0Translation2D.yaml config/pdeLevelSet/a0Identity2Dvortex.yaml" $S/driver.sbatch`
4. In `M` only, build the polyhedral meshes, the unit-gate mesh and the ghost-cell rung:
   - `TAG=polymesh`, `CFGS=config/integration/advConv3DshearPolyMesh.yaml`, `EXTRA="--until preprocess"`
   - `TAG=unitmesh`, `CFGS=config/integration/ghostCellUnitGateN32.yaml`, `EXTRA="--until decompose"`
   - `TAG=ghost`, `CFGS=config/integration/ghostCellRefinedBoxN64.yaml`
5. Polyhedral arms. For each mesh case `c` of `$M/studies/advConv3DshearPolyMesh` and each
   `valueBound` `b` (`none`, `lipschitzCone`, `stencilBounds`), in EACH clone root:
   `sbatch --export=ALL,INT_SCRIPTS=$S,SRC=$M/studies/advConv3DshearPolyMesh/$c,OUT=studies/intPoly/${c}_$b,BOUND=$b $S/poly_arm.sbatch`.
   Both clones read the same mesh case. No solver ran in that case, and each clone makes its
   own `0/` from `0.org` with its own `leiaSetFields`.
6. Unit gate, in `M`:
   `sbatch --export=ALL,INT_SCRIPTS=$S,CASES="<the two case directories of studies/ghostCellUnitGateN32>" $S/unit_gate.sbatch`.
7. Compare: `bash $M/workflow/scripts/integration_identity_check.sh $M $B`.

## Read-out

PASS: every CSV pair is identical in every column except `ELAPSED_CPU_TIME` and
`ELAPSED_CLOCK_TIME`. A case that diverged in both trees with equal CSVs is identical.
A launch failure is not a result: repair it and run the case again.
