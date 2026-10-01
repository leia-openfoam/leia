---
title: "leiaRedistancedLevelSetFoam retired: the redistancing line runs in leiaLevelSetFoam"
description: "Settled 2026-09-30: leiaLevelSetFoam with eulerian advection, velocityExtension none and a redistancer reproduces leiaRedistancedLevelSetFoam bit for bit in 10 arms on 4 ranks, so the second solver is deleted; the same test on the semi-Lagrangian pair fails, and leiaSemiLagrangeLevelSetFoam stays."
aliases: [retired redistanced solver, leiaRedistancedLevelSetFoam]
kind: decision
status: settled
part: advection
tags: [decision, part/advection]
date: 2026-09-30
date_settled: 2026-09-30
decided_by: ["author decision 2026-09-30", "retirement gate of STATUS 11.22 (scratch configs, axes in STATUS)"]
code: [applications/solvers/leiaLevelSetFoam/leiaLevelSetFoam.C, src/leiaLevelSet/advection/eulerianAdvection.C, Allwmake, workflow/scripts/aggregate.py, workflow/scripts/make_grl_fig.py, workflow/scripts/paths.py]
sources: ["STATUS 11.22 (L4859-L4939)", "commit 6bd5b7b7", "commit 2a53364b"]
---
# leiaRedistancedLevelSetFoam retired: the redistancing line runs in leiaLevelSetFoam

> Settled 2026-09-30 on the author's request. `leiaLevelSetFoam` with `eulerian` advection, `velocityExtension none` and a redistancer reproduces `leiaRedistancedLevelSetFoam` bit for bit: 10 arms on 4 ranks, every metric column identical at every step, the final psi and alpha byte-identical on every rank ([STATUS 11.22 items 1 to 2](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4866-L4891)). The second solver is deleted, and the nine redistancing studies select `leiaLevelSetFoam` with an explicit theme (commit 6bd5b7b7). The same test on the semi-Lagrangian pair fails: `leiaSemiLagrangeLevelSetFoam` stays ([STATUS 11.22 item 5](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4910-L4939), [[models/level-set-advection]]).

## The question

The author asked on 2026-09-30: `leiaLevelSetFoam` can run with no source term, no velocity extension and a redistancer, so why a second solver? The old solver assembled `ddt(psi) + div(phi psi) - psi div(phi) = S_SDPLS(psi)` inline, with the deferred-correction loop, and then called the gated redistancer ([leiaRedistancedLevelSetFoam.C L159-L184](https://github.com/leia-openfoam/leia/blob/06e15ecf/applications/solvers/leiaRedistancedLevelSetFoam/leiaRedistancedLevelSetFoam.C#L159-L184), the last commit that has it). `leiaLevelSetFoam` hands the same step to its `eulerian` model, which builds the same equation with the extension flux ([eulerianAdvection.C L90-L100](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/src/leiaLevelSet/advection/eulerianAdvection.C#L90-L100)); with `velocityExtension none` that flux is the prescribed flux. Then it calls the same redistancer and a volume correction whose default does nothing ([leiaLevelSetFoam.C L168-L226](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/applications/solvers/leiaLevelSetFoam/leiaLevelSetFoam.C#L168-L226)). Both solvers have the option `-fluxCorrection`, and both set the adaptive time step from the advecting flux, which is the prescribed flux with `velocityExtension none`. The commit that added the second solver says only that it "exercises the redistancer line on its own" (2a53364b, 2026-07-30): a reason of organisation, not of numerics.

## The measurement that decided it

Pre-registered: bit-identical in every metric column at every step, and the final fields byte-identical. The runs went through the workflow up to the solve rule, before (the old solver, the tree of 06e15ecf) and after (the changed tree, with the old solver deleted and its binary removed); laptop, OpenFOAM-v2606, np 4.

| arm | metric | value | where |
|---|---|---|---|
| 2Dvortex N = 32, T = 0.5; `noRedistancing`, `PDE`, `planeFootWave`, `anchoredEikonal` | solver CSV (14 columns) and `gradPsiError.csv`, `compare_metrics_csv.py --tol 0` | PASS in all 4 arms, 74 steps; final psi and alpha byte-identical, 8 of 8 files | MEASURED, [STATUS 11.22 item 2](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4877-L4891) |
| 2Dvortex N = 32; SDPLS `R` and `beta`, each `strictNegativeSpLinearImplicit` and `explicit` | the same | PASS in all 4 arms, 74 steps; fields byte-identical | MEASURED, same |
| 3Dshear hex N = 24, T = 0.75; `PDE`, `planeFootWave` | the same | PASS in both arms, 85 steps; fields byte-identical | MEASURED, same |
| control: `PDE` against `planeFootWave`, 2Dvortex | the same comparison | FAIL in 7 columns, so the redistancer acts and the identity is not trivial | MEASURED, same |
| after the change | build | 25 executables and 14 libraries, 0 missing; the CI script passes | MEASURED, [STATUS 11.22 item 3](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4892-L4902) |

## What it does not cover

1. The curated data tables, the plan documents and the older STATUS sections keep the name `leiaRedistancedLevelSetFoam`: they record what ran. `aggregate.py`, `make_grl_fig.py` and `paths.py` still read the old name for those archived studies.
2. The retirement gate is a transition gate: it needs the old solver, so it cannot run again after 6bd5b7b7. Its axes are in STATUS 11.22; its scratch studies are not kept.
3. Two traps, fixed in the same commit: a clone that built the old solver keeps its untracked object folder, and `wmake` stopped on the missing `Make/files` (`Allwmake` rc 2); `make_grl_fig.py` divided by zero on two rows at the same h, on the committed tables alone ([STATUS 11.22 item 4](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4903-L4909)).
4. The semi-Lagrangian pair is NOT equivalent; on 2026-10-01 the author settled that `leiaSemiLagrangeLevelSetFoam` stays as the kinematic solver of the SL two-phase solver ([[decisions/kinematic-solver-per-flow-solver]]), and the stale band is fixed (6f63418a). Three items were open: `semiLagrangianAdvection` has no trace-velocity option (no `projectedFlux`, the production default); the SL solver computes alpha with the narrow band of psi^n; the Eulerian `L_INF_E_PSI` has a sign-test defect ([STATUS 11.22 item 5](https://github.com/leia-openfoam/leia/blob/6bd5b7b7/STATUS.md#L4910-L4939), [[models/level-set-advection]]).
5. I was wrong in the chat answer of 2026-09-30: I wrote that both solvers came in one commit. `leiaLevelSetFoam` dates from 2021 (a0d38c2e); 2a53364b added the redistanced solver and changed `leiaLevelSetFoam`.

## Related

[[models/level-set-advection]] - [[models/redistancer]] - [[concepts/redistancing-geometric-grl]] - [[studies/grl-pre-print]] - [[concepts/eulerian-fv-transport]] - [[concepts/bit-identity-and-inertness-gates]] - [[hubs/advection]] - [[decision-log]]

## Log

### 2026-09-30
SETTLED: the retirement gate passed in 10 arms; the solver is deleted in 6bd5b7b7. Entered in [[decision-log#2026-09]].

### 2026-10-01
Item 4: the SL kinematic solver stays by the author's rule; the stale band is fixed ([[decisions/kinematic-solver-per-flow-solver]]).
