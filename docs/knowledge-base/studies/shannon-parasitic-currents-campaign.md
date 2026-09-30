---
title: "The Shannon parasitic-currents campaign"
description: "The map of docs/plan-shannon-parasitic-currents.md (2026-08-19 to 08-25): the evidence base, the per-step gain as the order parameter, the two-factor law, the Popinet reframing with its retraction, the interFoam pairing audit, the clamp and freeze verdicts, and the Polya sections that set the execution order."
aliases: []
kind: study
status: settled
part: surface-tension
tags: [study, part/surface-tension]
date: 2026-09-28
code: [docs/plan-shannon-parasitic-currents.md, workflow/scripts/per_step_gain.py, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam]
sources: [PSH sections 0-6, STATUS 4 (08-19 to 08-31), kb-raw C4, kb-raw C7, kb-raw C8]
---
# The Shannon parasitic-currents campaign

> **Verdict (2026-09-28).** `docs/plan-shannon-parasitic-currents.md` (1100 lines, 2026-08-19 to 08-25, last commit bf2170c) is the second campaign document against the parasitic-current instability. It reduced the score to one number, the per-step gain g = ln(u_end/u_0)/nSteps, and to one non-absorbable quantity, the alpha-weighted face gradient of the curvature ([section 1, line 891](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L891)); it measured the amplifier bare and wrote the two-factor law max U(T) = u_0(h) exp G(h) ([section 0d, line 244](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L244)); it exonerated the force lag, the momentum time order and the convective scheme; and it closed the attribution with three experiments: clamp, freeze and operator swap, which leave the phase-coherent refill of the curvature by the fit as the pump ([section 0g, line 830](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L830)). Its first reading of the Popinet result ("the force is explicit at t^n") is retracted inside the document ([[retractions/force-at-n-not-n-plus-1]]), and its 2D convergence claim was withdrawn with the psi-filter seam bug ([[retractions/psi-filter-seam-bug]]). What it points to, the variational pairing, has no implementation ([[concepts/variational-capillary-force]]).

## What it is

Written "from measured data only, no new runs" on 2026-08-19, then extended with the results of the runs it ordered. Sections:

| section | content | line |
|---|---|---|
| 0 | The evidence base: what converges (delivered non-gradient content +2.04, band 2h amplitude +2.05, max U at t=0 +3.58), what does not (end-state currents), the constant per-step amplification of the 3D ladder (r dt = 4.0 to 4.9e-4), the exclusions (pressure-velocity coupling, normal corrugation, the K term, solenoidality, the gradient drift), the bugs found (psi filter decomposition-dependent, the 93.5 % inverse fill, the t=0 Young-Laplace solve) | [8](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L8) |
| 0b | The order parameter g; u_0 about h^3.6 in both dimensions; g converges in 2D (h^+0.77) and diverges in 3D (h^-1.90, three points) | [73](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L73) |
| 0c | First results post-fix (2026-08-19): the seam fix reverses the 2D convergence claim; the fixed-dt test; the gate `psiOuterCorrectorsGain3D` exonerates the force lag | [117](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L117) |
| 0d | The amplifier bare (filter off) and the two-factor law; the filter's damping flips sign; BDF2 and the convective scheme on a matched window; the 2D volume error at fourth order (2026-08-20) | [244](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L244) |
| 0e | Popinet (2009) reframes the object: the m=8 capillary wave; the retraction "our force is at n+1, the scheme is symplectic" (2026-08-23) | [398](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L398), [530](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530) |
| 0f | Why interFoam stays bounded: pairing, not accuracy; the operator-swap falsification `variationalForce2D` (2026-08-25) | [585](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L585) |
| 0g | SAAMPLE read against the campaign; the clamp verdict; the freeze-kappa verdict (2026-08-25) | [705](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L705), [807](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L807), [830](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L830) |
| 1 to 6 | The Polya sections: cut it down (two quantities matter), problems already solved (interFoam, interFlow), say it differently (is the discrete equilibrium a fixed point), break it up (sub-map gains), flip it (the specification g at most 4e-5), make it bigger (g as a solver column, a gain gate) | [891](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L891) to [1036](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L1036) |
| -- | Execution order by information per unit cost, phases 0 to 3 | [1059](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L1059) |

The curated outputs of its runs are in the method-comparison data folder (`filterOffAmplifier2D_errors.csv`, `kickOriginGate2D_errors.csv`, `mufGrid2D_*`, `parasitic_mode_maxloc_Rh25.csv`, the `oscTraceGate2D` and `oscOuterCorrectors2D` CSVs) and in `REPRODUCE.md` ([[studies/method-comparison]]).

## Why it matters

The campaign wrote the scoring rules that every later gate uses: compare arms at equal step counts, read the gain and not the blow-up time, report the whole vector ([[concepts/error-vector-and-read-out-instants]]); and it produced the mechanism statement of [[concepts/parasitic-current-mechanism]]: the source is the curvature error, the amplifier grows with the translation and the density ratio, the pump is the phase-coherent response of the refitted curvature to the interface displacement.

## Where in the code

`workflow/scripts/per_step_gain.py`, the columns `kappaMeanBand`, `kappaStdDevBand`, `kappaMeanActiveFaces`, `kappaStdDevActiveFaces` of the droplet metrics, the environment hooks `SL_FREEZE_KAPPA`, `SL_FREEZE_RHOPHI`, the `curvatureClamp` entry, the force models `divGradAlphaSnGradAlpha` (the operator swap) and `constantCurvatureSurfaceTension`, `config/filterOffAmplifier3D.yaml`, `config/psiOuterCorrectorsGain3D.yaml`, `config/upwindConvection2D.yaml`, `config/kappaClamp2D.yaml`, `config/freezeKappa2D.yaml`.

## Evidence

| claim | number | where |
|---|---|---|
| The seam fix reverses the 2D convergence claim | post-fix max U 2.24e-5, 2.81e-5, 8.83e-5 at N=64, 128, 256 (orders -0.32, -1.65); shape and volume still converge (+3.41 / +1.38, +5.00 / +4.09) | [section 0c, line 119](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L119), MEASURED; the maximum norm is what the record has for this claim |
| The gain changes sign with resolution in both dimensions | 2D gAvg -7.61e-4, -6.41e-5, +7.22e-5; 3D -1.91e-4, -8.35e-5, +2.70e-4 at R/h = 10.0, 12.7, 15.8 | [section 0c, line 172](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L172), MEASURED |
| The force lag is exonerated | `psiOuterCorrectors` changes the endpoint by +0.4 % and -0.1 % | [section 0c, line 207](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L207), MEASURED; [[concepts/psi-outer-correctors]] |
| The two-factor law | u_0 about h^3.5 (2D triple 3.62, 3.51, 3.69); G about h^-3.27 in 3D, h^-0.95 then h^-0.56 in 2D | [section 0d, line 266](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L266), MEASURED |
| The filter's damping flips sign | 1.61x worse at R/h=10.0, 5.86x better at 15.8 | [section 0d, line 326](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L326), MEASURED; [[decisions/psi-filter-none]] |
| BDF2 against Euler on a matched window | gain +11.1, +2.9, -3.0 % at N=64, 128, 256 (noise); the convective scheme inert to 0.06 % for N of at least 128 | [section 0d, line 351](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L351), MEASURED; [[decisions/momentum-schemes-bdf2-upwind]] |
| 2D volume error | orders 4.41, 4.01, 4.17 | [section 0d, line 385](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L385), MEASURED |
| The oscillation is the m=8 capillary wave | period 1.56 ms resolution-independent against 1.455 ms predicted; envelope -97 to +212 1/s against a physical -128 1/s; 83 % of the curvature error is non-absorbable variation | [section 0e, line 446](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L446), MEASURED |
| The force is at n+1 and the scheme is symplectic | det M = 1; centring at n+1/2 is spectrally identical | [section 0e, line 530](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530), DERIVED; [[concepts/force-time-centring]] |
| The operator swap does not bound | all eight `divGradAlphaSnGradAlpha` arms reach 0.14 to 0.84 m/s after one step and die within 148 to 493 steps | [section 0f, line 658](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L658), MEASURED |
| The clamp | growth persists with the curvature clamped to [-2000, 4000] 1/m, the band standard deviation held at 45 1/m against 2.1e4 unclamped | [section 0g, line 807](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L807), MEASURED |
| The freeze | no exponential departure; endpoint 2.20e-3 against 2.59e-2 (live) and 1.48e-2 (clamped) at 25 ms; shape 12x worse | [section 0g, line 830](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L830), MEASURED |
| The specification | g at most 4e-5 against the measured 4.4e-4 at N_L=120: a 10x reduction | [section 5, line 1010](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L1010), DERIVED |

## Decisions

- Score on g and on the non-absorbable content only; abandon t_blow and the end-state maximum as the metric ([section 1](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L891)) ([[retractions/t-blow-baseline]]).
- Keep BDF2 as the default on formal grounds ([section 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L351)) ([[decisions/momentum-schemes-bdf2-upwind]]).
- Compare arms on equal step counts; an unmatched comparison made BDF2 look 3.1x worse ([section 0d](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L351)) ([[concepts/error-vector-and-read-out-instants]]).
- The semi-implicit capillary force "buys nothing" after the gate of section 0c and is dropped from phase 3 unless revived; it was revived, deadlocked on the cluster on 2026-08-28 and is not needed with `projectedFlux` ([[models/semi-implicit-capillary-force]]).
- The pairing route (a force derived variationally from the advected state) is the only theorem-grade path ([section 0g](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L830)) ([[concepts/variational-capillary-force]]).
- Theta of the psi filter would have to vanish with the corrugation content; no filtering in production ([[decisions/psi-filter-none]]).

## Retracted or superseded inside it

1. "2D converges in every metric" is withdrawn: it was an artefact of the psi-filter seam bug ([section 0c, line 190](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L190)) ([[retractions/psi-filter-seam-bug]]).
2. "r dt constant across the ladder" was h and dt covarying; at fixed h the e-folds scale as dt^0.35 ([section 0c, line 145](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L145)).
3. The explicit-force reading and its numbers (+22 to 88 1/s, the n+1/2 remedy) are retracted ([section 0e, line 530](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L530)) ([[retractions/force-at-n-not-n-plus-1]]).
4. The tangential-structure diagnostic recommended earlier is deprioritised by section 3B ([line 945](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-shannon-parasitic-currents.md#L945)).
5. Later records that supersede parts of it: the translating droplet was a closed box at the time ([[retractions/closed-box-translating-droplet]]), so the SAAMPLE prediction of section 0g ("our translating droplet should grow slower than the stationary one") was never tested on a valid case; the parallel two-phase runs predate the seam fixes of 2026-09-27 ([[concepts/coupled-face-density-defect]]).

## What it does not cover

The curvature deliveries and the ellipse gate ([[studies/curvature-stabilization-campaign]]), the projectedFlux trace velocity that later became production ([[concepts/trace-velocity-projected-flux]]), the translating droplet after the closed-box repair, and any implementation of the variational force.

## Related

[[hubs/surface-tension]], [[hubs/verification]]. [[concepts/parasitic-current-mechanism]], [[concepts/force-time-centring]], [[concepts/capillary-time-step]], [[concepts/balanced-force-csf-flux]], [[concepts/variational-capillary-force]], [[concepts/psi-outer-correctors]], [[concepts/kang-gfm-and-sharp-heaviside]], [[concepts/cell-centre-inverse-curvature]], [[models/semi-implicit-capillary-force]], [[cases/stationary-droplet]]. Siblings: [[studies/curvature-stabilization-campaign]], [[studies/poly3d-roadmap]], [[studies/method-comparison]], [[studies/sl-quadratic-pre-print]].

## Log

### 2026-09-28
Created from the plan document at 8867581.
