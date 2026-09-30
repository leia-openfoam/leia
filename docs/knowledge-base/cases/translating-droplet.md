---
title: "The translating droplet (2D and 3D)"
description: "The Galilean test of the coupled solver, a water droplet carried by a uniform air stream: a closed box until 2026-09-02, repaired and gated since; the method gate stops it at 0.05 s; beyond that the outlet triggers a fast growth and a slower interior growth remains open (2026-09-27)"
aliases: [translatingDroplet2D, translatingDroplet3D, translating arm, translating droplet]
kind: case
status: settled
part: mass-flux
tags: [case, part/mass-flux]
date: 2026-09-29
date_settled: 2026-09-27
decided_by: [config/translatingFreeStreamGate2D.yaml, config/gates/methodGate2D.yaml]
code: [cases/translatingDroplet2D, cases/translatingDroplet2D.parameter, cases/translatingDroplet3D, cases/translatingDroplet3D.parameter, cases/translatingDroplet3D_poly.parameter, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createDropletMetricsFile.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H, workflow/scripts/make_gate_summary.py]
sources: [STATUS 0, STATUS 11.13, STATUS 11.14, STATUS 11.15, SL article sec:translating, SL article sec:translating-late, CLAUDE wrong-setup rule, CLAUDE step 5 boundary paragraph, METHOD 8.1 row CURVATURE_EXTENSION, G2 translating arm]
---
# The translating droplet (2D and 3D)

> Verdict (2026-09-29). The translating droplet is the Galilean test of the coupled solver. A water droplet of R = 1 mm moves with a uniform air stream of U0 = 0.05 m/s. The exact solution is a rigid translation ([SL article, sec:translating](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2072-L2100), `sec:translating`). Until 2026-09-02 the 2D mesh had no inlet and no outlet, so every earlier result on the case is VOID ([STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L14-L46), [[retractions/closed-box-translating-droplet]]). Commit 440107f split the boundary into inlet, outlet and slip walls. The pre-registered free-stream gate then passed, with a first-step continuity error of 2.32e-20 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L34-L53)). The 3D case had no `dropletReferenceVelocity` until 2026-09-27, so its disturbance columns reported the translation itself ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3734-L3741)). The 2D method gate runs the case to 0.05 s with `cellCentreInverse`. The pre-registered horizon was 0.1 s; the cut came after the baseline diverged at every rung ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3614-L3640)). Beyond 0.05 s, the outlet triggers a fast growth in the 10 mm box. A slower interior growth of about 30 1/s remains in a 20 mm box ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3834-L3847), [[retractions/late-translating-instability-is-the-outlet]]). A 40 mm box completes 0.3 s. Its degradation is 2.3x, 5.4x and 2.1x smaller at N = 142 than at N = 100, on two rungs only ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3950-L3965)). The mechanism of the interior growth is open.

## What it is

The 2D case, `cases/translatingDroplet2D`:

1. The box is 10 mm by 10 mm and one cell thick: one `DOMAIN_LENGTH` of 0.01 m sets both x and y ([`blockMeshDict.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/system/blockMeshDict.template#L42-L60)).
2. The patches are `inlet` (left), `outlet` (right), `walls` (top and bottom, type `wall`) and `frontAndBack` (empty) ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/system/blockMeshDict.template#L66-L123)).
3. The conditions: `U` is `(U0 0 0)` at the inlet, zeroGradient at the outlet and slip on the walls ([`U.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/0.org/U.template#L36-L38)). `p_rgh` is fixed to 0 at the outlet only ([`p_rgh`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/0.org/p_rgh#L23-L25)). `alpha.water` is 0 at the inlet ([`alpha.water`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/0.org/alpha.water#L24-L26)). `psi` is zeroGradient on every patch.
4. The droplet is an `implicitSphere` of R = 1 mm at x = L/2 + `DROPLET_OFFSET_X` ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/system/fvSolution.template#L399-L410), [`default.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/default.parameter#L1009-L1028)). The default offset 0 puts it at the box centre; the gate sets -2.5 mm, so it starts at x = 2.5 mm.
5. The case layer: N = 32, 64, 128, 256, `END_TIME` 0.1 s, `dt = 0.010861 N^-1.5` ([`translatingDroplet2D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D.parameter#L10-L40)). The SL article uses N = 128 (R/h = 12.8), a density ratio of 838.8 and 0.2323 of the Brackbill limit. Every filter is off ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2098-L2100)).

The 3D case, `cases/translatingDroplet3D`, is a 6 mm box with the droplet at the centre. It runs 0.02 s, so the droplet moves 1 mm. Its downstream edge stops 1 mm (12.7 cells at N = 76) before the outlet ([`translatingDroplet3D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet3D.parameter#L1-L37)). It has six patches, `left` (inlet) and `right` (outlet) of type `patch` and four slip walls ([`blockMeshDict.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet3D/system/blockMeshDict.template#L54-L94)). The names are those of cfMesh, so one set of fields serves the hexahedral and the polyhedral mesh ([`translatingDroplet3D_poly.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet3D_poly.parameter#L1-L37)).

The metrics. The solver writes one row per step to `leiaSemiLagrangianLevelSetTwoPhaseFoam.csv` ([`createDropletMetricsFile.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createDropletMetricsFile.H#L1-L45)). The entry `dropletReferenceVelocity` (default `(0 0 0)`, [same](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createDropletMetricsFile.H#L53-L54)) sets the disturbance `U' = U - U_ref` and the reference centre `x_0 + U_ref t` ([`writeDropletMetrics.H`](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L27), [L56](https://github.com/leia-openfoam/leia/blob/d1e3414/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/writeDropletMetrics.H#L56)). The solution does not read it. The gate reads this vector at T, in L2 and L1 only ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L206-L238)):

| gate entry | definition |
|---|---|
| shape error | `zeroSetRadialL2 / R`: the zero-set crossings against the translated circle |
| volume error | `abs(phaseVolumeRelError)` at T and at T/2 |
| band gradient error | `gradPsiL2ErrorBand`: L2 of `abs(grad psi) - 1` in the band |
| spurious current | `l2MagUPrime` (L2) and `meanMagUPrime` (L1) of `abs(U - U_ref)` |
| pressure jump error | `abs(pLaplace - sigma/R) / (sigma/R)` |
| curvature error | `kErrL2Band * R` |
| travelled fraction error | `abs((x_c(T) - x_c(0)) / (U0 T) - 1)` |

`maxMagUPrime` is an L_inf norm; the gate does not score it ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L47-L54), [[concepts/error-vector-and-read-out-instants]]).

## Why it matters

In its own frame the translating droplet is the stationary droplet. A method that is right at rest and wrong in a moving frame fails here and nowhere else ([`translatingDroplet3D.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet3D.parameter#L3-L8)). The stationary droplet cannot see an error that is proportional to the bulk velocity ([SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2074-L2076), [[cases/stationary-droplet]]). On the repaired case, the curvature error is the source of the disturbance. The translation and the density ratio act on the amplifier ([[concepts/parasitic-current-mechanism]], [[concepts/density-ratio-amplifier]]). Two repository rules come from this case: "A wrong setup voids its data" ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L705-L749), [[concepts/wrong-setup-voids]]) and "Is the interface still inside the domain?" ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L528-L544)).

## Evidence

| claim | number | where |
|---|---|---|
| the closed box annihilated the stream | first-step continuity error 1.00e-05 against 2.32e-20 repaired; `mean\|U-U0\|` 5.011e-02 against 1.160e-05; whole-domain mean Ux -7.3e-07 | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L34-L46) |
| the free-stream gate passes (N = 128, 200 steps) | `mean\|U-U0\|` 1.2e-05 to 7.6e-05, about 50x below the artefact | MEASURED, [STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L48-L53), criteria in [`translatingFreeStreamGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/translatingFreeStreamGate2D.yaml#L13-L22) |
| the disturbance after one step | L1 2.32e-4 and L2 1.76e-3 of `U_ref` | MEASURED, [SL article](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2102-L2105) |
| the 3D reference velocity was missing | `meanMagUPrime = 0.0500 = U0` in the 3D gate smoke; four trace studies are VOID for the disturbance, zero-set and centroid columns | MEASURED, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3734-L3741), fix in [`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet3D/system/fvSolution.template#L196-L203) |
| the gate arm at 0.1 s | the baseline diverged at 0.0904, 0.0942 and 0.0772 s at N = 100, 142, 200 | MEASURED, [STATUS 11.13](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3614-L3626) |
| the gate arm at 0.05 s | the baseline completes every rung; 15 re-runs reproduce the first 0.05 s of their 0.1 s runs byte for byte | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3782-L3785) |
| the baseline vector at 0.05 s, N = 100, 142, 200 | shape/R 2.68e-2, 1.83e-2, 3.06e-3 (least-squares order 3.12, pairwise 1.08 and 5.23); `l2MagUPrime` 1.95e-3, 1.26e-3, 8.93e-4 m/s (1.13); volume 1.99e-3, 1.79e-3, 2.05e-4 (3.27); curvature 3.3, 2.4, 3.6 % (-0.10) | MEASURED, archive [`summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/summary.csv#L12-L14), [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/orders.csv#L23-L32) |
| the decomposition check of the arm | serial against np 4 at most 3.47e-7 (`meanMagUPrime`), tolerance 1e-5 | MEASURED, [`seam.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/seam.csv#L6), [[concepts/coupled-face-density-defect]] |
| the late instability is not a seam effect | serial 0.0775 s; np 4 0.0904 s before the fix and 0.0868 s after it | MEASURED, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3660-L3670) |
| one-change discriminators, N = 100, one resolution | `curvatureExtension none` 0.0695 s (20 % earlier); `footIntegrator rk2` 0.0842 s (no effect); midpoint force centring 0.0593 s (32 % earlier); density ratio 1 completes | MEASURED, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3743-L3761) |
| the outlet triggers the fast growth | 10 mm box: L2 `abs(U-U0)` 1.79e-3 at 0.06 s to 1.35e-1 at 0.08 s; the 20 mm box decays to 6.73e-4 | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3787-L3802) |
| the long-box ladder to 0.1 s | zero-set L2 1.38e-4, 7.08e-5, 3.99e-5 m (orders 1.90, 1.68); centroid orders 1.89, 1.69; curvature error 1.5 to 2.2 %, not convergent | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3818-L3832) |
| the interior growth, 20 mm box | N = 100 diverges at 0.2127 s, growth from 0.10 s at about 30 1/s with the droplet more than 8 mm from the outlet; N = 142 diverges at 0.2194 s | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3834-L3847) |
| the start-position test (x0 = 5 mm) | jump at 0.1372 s with the centroid at 12.15 mm, so neither prediction holds; the slow phase is identical to 0.12 s | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3866-L3881) |
| the 40 mm box, N = 100 | completes 0.3 s; L2 3.6e-4 to 8.3e-3; volume error 1.5e-1; lead 3.86 mm | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3950-L3959) |
| the 40 mm box shrinks with h | at 0.3 s 2.3x (current), 5.4x (volume), 2.1x (lead) smaller at N = 142; two rungs, no order | MEASURED, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3961-L3965) |

All the late-instability runs are laptop runs. Their libraries carry the stamp `...-183-g935fd4e-dirty`, their solvers come from b1798c3, and the raw output is preserved in the git-ignored `runs/gcls-laptop-20260927` ([archive README](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/README.md#L17), [`translating_runs.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/laptop/translating_runs.csv#L2-L17)).

## Why it failed, or why we think so

1. The closed box (VOIDED 2026-09-02). The blockMeshDict put all four sides into one `walls` patch. OpenFOAM ignored the field entries `inlet` and `outlet`, which matched no mesh patch. Slip gives `U.n = 0` on the x faces, so the projection removed the stream on step 1. Then `maxMagUPrime` measured the removed stream at about 2 U0 ([STATUS 0](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L23-L46)).
2. The missing 3D reference velocity (fixed 2026-09-27). The writer used `(0 0 0)`, so the disturbance and zero-set columns measured the translation and the displacement ([STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3734-L3741)).
3. The late instability (open). The reading that fits every run: a slow interior growth in time, and an amplification that increases as the droplet approaches the outlet. The jump occurs when the two together cross a threshold ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3873-L3881)). HYPOTHESIS for the interior part: a discretisation error, because it shrinks with h between two rungs ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3961-L3963)). It needs the density contrast. dec002f falsified the mass-momentum consistency as its dominant term ([[retractions/mass-momentum-consistency-dominant-term]]).

## Decisions

- Commit 440107f split the boundary, and the free-stream gate ran before any production arm ([`translatingFreeStreamGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/translatingFreeStreamGate2D.yaml#L1-L22)).
- The gate arm ends at 0.05 s (2.5 R of travel) with `cellCentreInverse`, `DROPLET_OFFSET_X -0.0025` and the `twoPhaseCoupling` block ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L126-L154)). The record states the cut as a change after seeing data. The `none` of 8a9b85a was reverted, because its evidence is VOID ([STATUS 11.13](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3631-L3640)).
- The 3D gate arm runs 0.02 s, 1 R of travel, at N = 60, 78, 102 ([`methodGate3D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate3D.yaml#L97-L111)). It has not run.
- `CURVATURE_EXTENSION` on this case is undecided: no valid measurement compares `none` with `cellCentreInverse` ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/d1e3414/METHOD.md#L393), [[decisions/curvature-extension-cell-centre-inverse]]).
- Case-dependent values live in the configuration files, not in the agent guide. The translating arm once inherited a global default that the guide called case-specific ([CLAUDE.md](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L170-L176), 2026-09-27).

## Open questions

1. The mechanism of the interior growth at the water/air density ratio. The third 40 mm rung (N = 200, 160 000 cells, 78 000 steps) belongs on the cluster; no order exists yet ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3963-L3965)).
2. A box-length token for the gates and ladders is an open author decision ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3800-L3802)). The template has one `DOMAIN_LENGTH` for x and y. The 20 mm and 40 mm boxes came from a hand edit of the rendered mesh file (a `domainX` entry in `runs/gcls-laptop-20260927/latecheck/box40/system/blockMeshDict`, read 2026-09-29).
3. Not tested: an outlet condition other than fixed `p_rgh` with zeroGradient U, and the semi-implicit capillary force ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3880-L3881), [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3760-L3761), [[models/semi-implicit-capillary-force]]).
4. Two stale texts. The 2D template still says "MEASURED 2026-09-01: ... the leading edge reaches the outlet at t = 0.08" ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D/system/fvSolution.template#L405-L407)). That is a closed-box measurement without a VOID marker. The header of `translatingDroplet2D.parameter` describes the stationary droplet ([L1-L4](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/translatingDroplet2D.parameter#L1-L4)).
5. A number conflict. The 40 mm box at N = 100 and t = 0.3 s has two volume values. One STATUS table and its sentence give 10 % ([L3941-L3945](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3941-L3945)). The next table gives 1.5e-1 ([L3956](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3956)), the SL article gives 15 % ([L2306-L2309](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/semi-lagrangian-level-set/sl-level-set-article/semiLagrangianLevelSet.tex#L2306-L2309)), and the archived history reads 0.147 ([`translating_histories.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/laptop/translating_histories.csv#L4599)). This note uses 1.5e-1.

## Related

- Hubs: [[hubs/mass-flux]], [[hubs/verification]].
- Retractions: [[retractions/closed-box-translating-droplet]], [[retractions/late-translating-instability-is-the-outlet]], [[retractions/mass-momentum-consistency-dominant-term]].
- Mechanism: [[concepts/parasitic-current-mechanism]], [[concepts/density-ratio-amplifier]], [[concepts/force-time-centring]], [[concepts/cell-centre-inverse-curvature]].
- Mass flux: [[models/mass-flux]], [[concepts/rholent-mass-flux]], [[concepts/coupled-face-density-defect]].
- Method: [[concepts/method-gates]], [[concepts/wrong-setup-voids]], [[concepts/error-vector-and-read-out-instants]].
- Cases: [[cases/popinet-translating-droplet]], [[cases/stationary-droplet]], [[cases/benchmark-cases]].
- Studies: [[studies/sl-quadratic-pre-print]], [[studies/method-gate-2d-campaign-2026-09]].

## Log

### 2026-09-29
Created from STATUS 0, 11.13 to 11.15, the SL article, the case files, the gate configs and the gcls data archive.
