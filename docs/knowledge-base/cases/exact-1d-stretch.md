---
title: "The exact 1D stretch (1Dstretch)"
description: "The only case with a closed form for the whole level-set field, u = alpha x with q = exp(-alpha t): the first arm of both method gates; the baseline converges at order 1.9, but for every candidate the script sets the closed form to None, so the arm cannot fail a candidate that completes (found 2026-09-28; the repair E0.1 is open)"
aliases: [1Dstretch, exact1D, exact 1D arm, uniaxial stretching gate]
kind: case
status: settled
part: verification
tags: [case, part/verification]
date: 2026-09-29
date_settled: 2026-09-28
decided_by: [config/gates/methodGate2D.yaml, config/gates/methodGate3D.yaml, config/sdpls1Dstretch.yaml]
code: [cases/1Dstretch, cases/1Dstretch.parameter, workflow/scripts/make_gate_summary.py, workflow/Snakefile.gate, config/gates/methodGate2D.yaml]
sources: [STATUS 11.5, STATUS 11.8, STATUS 11.12, STATUS 11.17, G2 exact1D arm, make_gate_summary.py L186-L198 and L479-L492, technical report sec:gate and sec:relocation and tab:experiments, sdpls1Dstretch header]
---
# The exact 1D stretch (1Dstretch)

> Verdict (2026-09-29). `cases/1Dstretch` is the only case with a closed form for the interface and for the whole level-set field ([`1Dstretch.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch.parameter#L1-L17)). Its error is therefore a true error. Both method gates run it first, at N = 32 to 256 and to T = 1 s. The other arms start only after it passes ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L89-L100)). The baseline band mean of `abs(grad psi)` reaches e^-1. Its relative error is 1.22e-4 at N = 32 and 2.53e-6 at N = 256, least-squares order 1.87 ([`summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/summary.csv#L2-L5), [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/orders.csv#L2)). For every candidate with an extension or a source, the script sets the closed form to None ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L186-L198)). The 1D check then passes a candidate on completion alone, so the arm cannot fail a candidate that completes ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L479-L492), [STATUS 11.17](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L4068-L4073)). The candidates' band means are recorded, 0.46 to 0.59 for HL0 and 0.998 to 0.697 for S1, but no closed form judges them ([`HL0/summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/HL0/summary.csv#L2-L5), [`S1/summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/S1/summary.csv#L2-L5)). The repair E0.1 integrates the closed form of each candidate per band cell; it has not run ([[concepts/gradient-control-next-experiments]]).

## What it is

The geometry and the flow:

1. The mesh is one-dimensional: N cells in x on [-2, 2] m, one cell in y and z, and `empty` transverse patches ([`blockMeshDict.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/blockMeshDict.template#L17-L35)).
2. Both x patches are outflow patches for this velocity, so zeroGradient is the correct condition; nothing enters ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/blockMeshDict.template#L60-L89)).
3. The velocity is `uniaxialStrain`, `u = alpha (x - x_0) e_x` with alpha = `STRAIN_RATE` = 1 1/s and `x_0 = 0` ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/fvSolution.template#L110-L132)). Its divergence is alpha: a 1D flow cannot be solenoidal. The solver adds `-Sp(divPhi, psi)`, so it solves the advective form.
4. The initial surface is the plane `psi_0 = x - x_i0` at `x_i0` = 0.25 m, an exact signed distance ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/fvSolution.template#L167-L194)). The interface reaches 0.25 e = 0.68 m at T = 1 s; an `END_TIME` above about 2.07 s needs a longer box ([`blockMeshDict.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/blockMeshDict.template#L23-L27)).

The closed forms ([`1Dstretch.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch.parameter#L8-L12)):

| arm | psi(x, t) | q = dpsi/dx | interface x_i(t) |
|---|---|---|---|
| no source (the gate baseline) | `x_0 + (x - x_0) exp(-alpha t) - x_i0` | `exp(-alpha t)` | `x_0 + (x_i0 - x_0) exp(alpha t)` |
| SDPLS source R | `x - x_i(t)` | 1 | the same |

The gate arm `exact1D` sets `END_TIME 1`, `CFL 0.5`, `implicitPlane`, `detrixheAslam`, `noRedistancing` and `N_LAYERS 3` at np 4 ([`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L89-L100)). The line tokens of the SL line (the uncached quadratic fit, the `projectedFlux` trace) replace the case layer's `cellCentred` trace ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L38-L46), [STATUS 11.8](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3403-L3405)). The 3D gate runs the same arm, so every 3D run carries its own exact check ([`methodGate3D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate3D.yaml#L9-L10)).

The read-out is `qBandMean`, the band mean of `abs(grad psi)` (`NARROW_MEAN_MAG_GRAD_PSI`) at T, and `qError = abs(qBandMean - exp(-alpha T)) / exp(-alpha T)` ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L186-L198)). The volume change and `abs(q - 1)` are physics in this case, not errors ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L187-L188)). The check `exact1d_check` has priority 100 and writes `exact1D.pass`; the other arms wait for it ([`Snakefile.gate`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/Snakefile.gate#L109-L121)). It passes when every rung completes and `qError` is None or at most 1e-3 at the finest rung ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L479-L492)).

## Why it matters

It is the only true error of the gate; every other kinematic case must use the surrogate `abs(abs(grad psi) - 1)` ([`1Dstretch.parameter`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch.parameter#L3-L6)). It also shows why the error vector has three entries. Without a source the zero contour is exact, while the distance property falls by a factor e. A shape-only or volume-only score then calls that arm perfect ([`sdpls1Dstretch.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/sdpls1Dstretch.yaml#L26-L30), [[concepts/error-vector-and-read-out-instants]]). And the strain field is an exact constant, so any deviation of a source or an extension is discretisation, never modelling error ([`fvSolution.template`](https://github.com/leia-openfoam/leia/blob/d1e3414/cases/1Dstretch/system/fvSolution.template#L112-L118)).

## Evidence

| claim | number | where |
|---|---|---|
| the baseline converges to the closed form | `qError` 1.22e-4, 3.98e-5, 1.01e-5, 2.53e-6 at N = 32, 64, 128, 256 (82 to 295 steps, np 4); pairwise orders 1.61, 1.98, 1.99; least squares 1.87 | MEASURED, [`summary.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/summary.csv#L2-L5), [`orders.csv`](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/orders.csv#L2) |
| the smoke and the cluster checks | 2D smoke: `qError` 1.6e-6, PASS; the first cluster campaign: six checks PASS, baseline 2.5e-6 | MEASURED, [STATUS 11.5](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3307-L3309), [STATUS 11.12](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3552-L3553) |
| the closed form is None for every candidate | `q_exact = exp(-alpha T)` only if `VELOCITY_EXTENSION none`, `SL_SOURCE none` and `SDPLS_SOURCE noSource` | MEASURED (a reading of the script), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L192-L198) |
| the 1D check cannot fail a completed candidate | PASS if no rung failed and `qError` is None or at most 1e-3; the log says "oracle pending: source or extension active" | MEASURED (a reading of the script), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L484-L491) |
| HL0 relocates the normal strain | band mean 0.589, 0.499, 0.510, 0.459 against 0.49, the band average of `exp(-K)` over d = 0.5 h, 1.5 h, 2.5 h | MEASURED and DERIVED, [technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L296-L298), [[concepts/extension-strain-relocation]] |
| the SDPLS arms at N = 64 | noSource 0.367872 (exact 0.367879), R 1.000000, Rdiv 0.367872; the interface error is +0.0002 h in all three | MEASURED, [`sdpls1Dstretch.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/sdpls1Dstretch.yaml#L12-L15) |
| Rdiv is inert in 1D | the two `fvm::div` terms cancel discretely and Rdiv assembles the zero matrix | DERIVED, [`sdpls1Dstretch.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/sdpls1Dstretch.yaml#L17-L24) |

The band means of the first campaign at T = 1 s. The fixed gate reproduced the kinematic arms byte for byte, so these values hold for both runs ([STATUS 11.15](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3897-L3899)). E0.1 gives the predictions ([technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L434-L442)):

| candidate | what it adds | N = 32 | N = 64 | N = 128 | N = 256 | E0.1 prediction | where |
|---|---|---|---|---|---|---|---|
| baseline | nothing | 0.3678 | 0.3679 | 0.3679 | 0.3679 | exp(-1) = 0.3679 | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/baseline/summary.csv#L2-L5) |
| HL0 | halo-limited extension, R = 1 h | 0.589 | 0.499 | 0.510 | 0.459 | about 0.49 | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/HL0/summary.csv#L2-L5) |
| HL1q | HL0 and the linearQ source | 0.799 | 0.764 | 0.651 | 0.526 | about 0.5 (mu = alpha) | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/HL1q/summary.csv#L2-L5) |
| HL1z | HL0 and the linearZ source | 0.765 | 0.729 | 0.629 | 0.516 | not recorded | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/HL1z/summary.csv#L2-L5) |
| HL2 | HL0 and the soft wall | 1.102 | 1.031 | 0.762 | 0.715 | not recorded | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/HL2/summary.csv#L2-L5) |
| S1 | the soft wall, no extension | 0.998 | 0.926 | 0.746 | 0.697 | 0.925 | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/S1/summary.csv#L2-L5) |
| FP0 | the closestPoint extension | 0.992 | 1.088 | 1.083 | 0.959 | not recorded | [csv](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-level-set-article/data/archive/shared-method-config-2026-09-01-192-g1150e68/gate/summaries/FP0/summary.csv#L2-L5) |

The candidate tokens are in `config/candidates/<name>.yaml` ([[concepts/method-gates]]). Every candidate completed every rung, so every one passed the 1D check. The record states no order for these band means. HL0 stays between 0.46 and 0.59, around its prediction of about 0.49. HL1q approaches its prediction as N grows. S1 falls below its prediction at N = 128 and 256.

## Why it failed, or why we think so

The gate had no closed form for a candidate. The script computes the exact value only for the plain method and leaves `qError` empty otherwise ([`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L192-L198)). `qError` is in the scored vector, but the regression test skips an empty value ([same](https://github.com/leia-openfoam/leia/blob/d1e3414/workflow/scripts/make_gate_summary.py#L404-L407)). So the arm reported a band mean for every candidate but could not fail anyone ([technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L154-L158)). The closed form of a candidate depends on its source rate F and on the strain K(d/R) that its extension transmits. Nobody derived it before the campaign ([STATUS 11.17](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L4068-L4073)).

## Decisions

- The exact 1D arm runs first, and no other arm starts before it passes (criterion 5, [`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/d1e3414/config/gates/methodGate2D.yaml#L28)).
- The SL line runs the production trace `projectedFlux` in this arm, set as a line token ([STATUS 11.8](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L3403-L3405)).

## Open questions

1. E0.1 (hours): integrate, per band cell, `dq/dt = q (F - alpha K(d/R))` and `dd/dt = alpha d (1 - c^2)`, with `K(t) = 1 - (1 - t^4)(1 + t^4)^(-3/2)`, `c = (1 + t^4)^(-1/4)` and `t = d/R` ([technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L284-L290), [STATUS 11.17](https://github.com/leia-openfoam/leia/blob/d1e3414/STATUS.md#L4071-L4073)). Then compare `qBandMean` with it. The predictions are in the table above; E0.1 falsifies nothing yet, and it makes the 1D gate real ([[concepts/gradient-control-next-experiments]]).
2. The Tier-I tests of the dossiers did not run: uniform translation, rigid rotation, planar strain (q_0 = 0.92, 1.00, 1.08) and a perturbed scaling. The zero-flow static test did not run either ([technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L148-L152)).
3. The 1D arm is the only exact rung; the kinematic arms carry no curvature column ([technical report](https://github.com/leia-openfoam/leia/blob/d1e3414/docs/gradient-controlled-level-set/gcls-technical-report/gclsTechnicalReport.tex#L164)). See the second blind spot in [[concepts/method-gates]].

## Related

- Hubs: [[hubs/verification]], [[hubs/gradient-control]].
- Method: [[concepts/method-gates]], [[concepts/error-vector-and-read-out-instants]], [[concepts/richardson-ladders-and-orders]].
- Gradient control: [[concepts/gradient-control-next-experiments]], [[concepts/extension-strain-relocation]], [[concepts/why-the-candidates-failed]], [[models/velocity-extension]], [[models/sl-source]], [[models/sdpls-source]].
- Campaign and studies: [[studies/method-gate-2d-campaign-2026-09]], [[studies/sdpls-pre-print]], [[studies/gcls-pre-print]].
- Cases: [[cases/benchmark-cases]], [[cases/kinematic-advection-cases]].

## Log

### 2026-09-29
Created from the case files, the gate configs, the scoring script, the gcls data archive, STATUS 11.5 to 11.17 and the technical report.
