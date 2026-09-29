---
title: "MASS_FLUX_BOUND_RHO: active with rhoLENT, evidence void"
description: "MASS_FLUX_BOUND_RHO true clips the rhoLENT density to [rho2, rho1]; its only measured basis (rhoDdtGate2D, rho to -72.28 without the bound) ran on the closed box and is void, so the decision is open"
aliases: [MASS_FLUX_BOUND_RHO true]
kind: decision
status: open
part: mass-flux
tags: [decision, part/mass-flux]
date: 2026-09-28
date_settled:
decided_by:
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, cases/default.parameter, config/gates/methodGate2D.yaml]
sources: ["METHOD 6 (L289-L299)", "METHOD 8.1 row MASS_FLUX_BOUND_RHO (L391)", "STATUS 11.13 (L3584)", "STATUS 11.14 (L3731-L3733)", "DP L768-L807", "G2 L68-L69"]
---
# MASS_FLUX_BOUND_RHO: active with rhoLENT, evidence void

> OPEN (2026-09-28). `MASS_FLUX_BOUND_RHO true` in the global default since 2026-09-02 ([DP L768-L807](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L768-L807)) and in the gate's coupling block ([methodGate2D.yaml L68-L69](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L68-L69)). The bound clips the rhoLENT auxiliary density to `[rho2, rho1]` after its solve. It is active in every study that does not override `MASS_FLUX`, because rhoLENT became the default the same day (82ca995) ([DP L800-L806](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L800-L806)). Its only measured basis, `config/rhoDdtGate2D` (density to -72.28 without the bound, clip L1 7.4e-04 with it), ran on the closed-box translating case and is void ([METHOD 6 L289-L299](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L289-L299), [METHOD 8.1 L391](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L391), [[retractions/closed-box-translating-droplet]]). No valid measurement of the bound exists after the fix 440107f; the record says it must be gated before it is called settled ([STATUS L3584](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3584)). The gate file cites the void study and carries the correction of 2026-09-28 ([methodGate2D.yaml L69](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L69)).

## The question

The rhoLENT density is the solution of an auxiliary transport equation, not the geometric density, so it is not bounded by construction. A density outside `[rho2, rho1]` is unphysical; a negative one flips the sign of `1/rho` in the pressure Laplacian and the run ends in a GAMG blow-up ([DP L772-L778](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L772-L778)). The clip injects a source `S_clip`, so `ddt(rho) + div(rhoPhi) = S_clip` and the momentum equation inherits `U S_clip`, the free-stream identity that rhoLENT exists to enforce; the solver reports `rhoClipL1` and `rhoClipFraction` so that the break is measured ([rhoLENTEqn.H L19-L61](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L19-L61)). The question is whether the bound is needed on the repaired case, and what it costs there.

## The measurement that decided it

None valid. The record:

| arm | metric | value | where |
|---|---|---|---|
| `rhoDdtGate2D`, N = 32, np 4, translating (closed box): `backward` without the bound | steps to FPE; min rho | 1074; -72.28 | VOID, [config L57-L60](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoDdtGate2D.yaml#L57-L60), [STATUS L55-L62](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L55-L62) |
| the same with `Euler` on the density equation | steps to FPE; min rho | 879; -23.28 | VOID, same |
| `matchedBDF2Translating2D`, `volumeCorrectionTranslating2D`: bound off | min rho | -27.7 and -1786 | VOID, [DP L943-L954](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L943-L954) |
| `rhoBoundGate2D` (4-rank collectives gate, closed box), 167 steps | clip activity; residual after the clip | 72 % of cells at 1.8e-05 against rho2 = 1.19; 2.3e-11 relative | VOID as physics; the collectives are safe (code), [config L32-L61](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoBoundGate2D.yaml#L32-L61) |
| `rhoLENTStationary2D`, N = 32 / 64 / 128 (valid, stationary) | clip L1 | 2.13e-12 / 2.56e-12 / 2.86e-12: at round-off | MEASURED, [config L44-L62](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L44-L62) |
| the 2026-09-01 GAMG blow-up to 6.6e+70 on the translating and oscillating droplets | — | the translating half is void; the oscillating half stands as an indicator only | [DP L800-L806](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L800-L806) |

Pre-registered read-out of the void gate ([rhoDdtGate2D.yaml L33-L45](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoDdtGate2D.yaml#L33-L45)): the convex-combination argument for Euler; both arms diverged, so the argument was falsified on that case, and the case was wrong.

What stands: the derivation that BDF2's homogeneous update `rho^{n+1} = (4/3) rho^n - (1/3) rho^{n-1} - ...` is an extrapolation and can undershoot ([rhoDdtGate2D.yaml L3-L18](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoDdtGate2D.yaml#L3-L18), DERIVED), and the matching argument that keeps `RHO_DDT_SCHEME backward` ([[decisions/momentum-schemes-bdf2-upwind]]).

## What it does not cover

1. The gate runs every droplet arm with the bound on ([methodGate2D.yaml L68](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L68)). A `boundRho false` control on the repaired translating and oscillating cases is the missing measurement; its read-out is the completion, `min rho`, `rhoClipL1` and the whole error vector.
2. `rhoClipFraction` counts round-off clips in pure cells, 37 to 48 % of the cells at a clip L1 of 1e-12, so as an error metric it measures round-off; it is reported, not scored ([STATUS L3731-L3733](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3731-L3733)).
3. The bound is inert for `geometricFaceDensity` and `interpolatedDensity`, which solve no density equation ([DP L768-L770](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L768-L770)).

## Related

[[hubs/mass-flux]] - [[models/mass-flux]] - [[concepts/bound-rho]] - [[concepts/rholent-mass-flux]] - [[concepts/ddt-scheme-pairing-bdf2]] - [[concepts/wrong-setup-voids]] - [[decisions/mass-flux-rholent]] - [[decisions/momentum-schemes-bdf2-upwind]] - [[retractions/closed-box-translating-droplet]] - [[cases/translating-droplet]] - [[decision-log]]

## Log

### 2026-09-28
OPEN. The basis was voided on 2026-09-27 (closed box); the bound stays active. Entered in [[decision-log#2026-09]] as REOPENED.
