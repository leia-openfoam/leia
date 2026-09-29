---
title: "boundRho: active, evidence void"
description: "The clip of the rhoLENT auxiliary density to [rho2, rho1] is on in every production run since 2026-09-02, its only measured basis ran on the closed box, and its clip fraction counts round-off; the decision is open"
aliases: []
kind: concept
status: open
part: mass-flux
tags: [concept, part/mass-flux]
date: 2026-09-28
date_settled:
decided_by:
code: [applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H, applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H, cases/default.parameter]
sources: [METHOD 6, METHOD 8.1 row MASS_FLUX_BOUND_RHO, STATUS 0, STATUS 11.13, STATUS 11.14, DP 765-803, G2 66-69]
---
# boundRho: active, evidence void

> Verdict (2026-09-28). `massFlux { boundRho true; }` clips the rhoLENT auxiliary density to `[rho2, rho1]` right after its solve ([`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L19-L61)). The token `MASS_FLUX_BOUND_RHO true` is the default since 28a1383 (2026-09-02), and because `MASS_FLUX` became `rhoLENT` the same day (82ca995) the clip is active in every study that does not override it ([STATUS 0](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L272-L275)). Its only measured basis, `rhoDdtGate2D` (rho to -72.28 without the bound, clipL1 7.4e-04 with it), ran on the closed-box translating case and is void ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L289-L299), [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L796-L802), [[retractions/closed-box-translating-droplet]]). No valid measurement of the bound exists after 440107f; METHOD 8.1 lists the row as "NONE valid" and says the bound must be gated before it is called settled ([METHOD 8.1](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L392), [[decisions/mass-flux-bound-rho]]).

## What it is

`rho1` and `rho2` are constants, so a density outside `[rho2, rho1]` is unphysical by construction. A negative density flips the sign of `1/rho` in the pressure Laplacian, the operator loses definiteness, and the run ends in a GAMG blow-up ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L769-L774)). The clip prevents that. It is not free: it injects a source `S_clip`, so the discrete identity becomes `ddt(rho) + div(rhoPhi) = S_clip` and the momentum equation inherits `U * S_clip`, the free-stream identity rhoLENT exists to enforce. The solver therefore writes `rhoClipL1` (volume mean of the clipped amount) and `rhoClipFraction` (fraction of cells clipped) to the metrics CSV, and the mass residual is computed after the clip ([`rhoLENTEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L28-L34)). The reductions of the clip are collective; every rank reaches them ([same](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/rhoLENTEqn.H#L54-L60)).

The code default is `false`; the token layer sets `true` ([`createMassFluxFields.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/createMassFluxFields.H#L117-L135), [token](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L803)). The Eulerian two-phase solver prints the two clip numbers to its log ([`alphaEqn.H`](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/alphaEqn.H#L120-L127)).

## Why it matters

The argument for the bound was written for matched BDF2: the homogeneous part of `backward` is the extrapolation `rho^{n+1} = (4/3) rho^n - (1/3) rho^{n-1}`, which undershoots by construction ([token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L791-L795), DERIVED). The falsified counter-hypothesis was that Euler on `ddt(rho)` restores boundedness by a convex-combination argument. The gate found that `rho_f` is a geometric face density, not an upwind value of `rho^n`, so the undershoot reaches about CFL x 997 in a cell that the interface cuts, and no time scheme rescues that ([`rhoDdtGate2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoDdtGate2D.yaml#L62-L74), DERIVED on a void run). The derivation stands as an argument; every number behind it is void.

## Evidence

| claim | number | where |
|---|---|---|
| without the bound, matched BDF2 drives rho negative | min rho -72.28 at step 1074 (`backward`); -23.28 at step 879 (`Euler`); 13.1 % and 6.4 % of the outer solves with rho below rho2 | [`rhoDdtGate2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoDdtGate2D.yaml#L50-L66), VOID (closed box) |
| with `boundRho false` over the full horizon | rho reached -27.7 and -1786 even with Euler on `ddt(rho)` | [token comment](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L938-L947), VOID |
| the clip on the translating droplet at 4 ranks | fires on all 167 steps, up to 72 % of the cells, magnitude 1.8e-05 against rho2 = 1.19 | [`rhoBoundGate2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoBoundGate2D.yaml#L32-L61), VOID |
| the clip on the stationary droplet | clipL1 2.13e-12 to 2.86e-12 at N = 32 to 128 | [`rhoLENTStationary2D`](https://github.com/leia-openfoam/leia/blob/8867581/config/rhoLENTStationary2D.yaml#L46-L52), MEASURED |
| `rhoClipFraction` counts round-off clips | 37 to 48 % of the cells at a clip L1 of 1e-12 in the pure phases | [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3712-L3715), [`make_gate_summary.py`](https://github.com/leia-openfoam/leia/blob/8867581/workflow/scripts/make_gate_summary.py#L57-L63), MEASURED |
| the gate config still cites the void gate | comment "rhoDdtGate2D: without it matched BDF2 drove rho to -72.28" | [`methodGate2D.yaml`](https://github.com/leia-openfoam/leia/blob/8867581/config/gates/methodGate2D.yaml#L69), stale at 8867581; corrected in the working tree on 2026-09-28 |

## Decisions

- `MASS_FLUX_BOUND_RHO true` (28a1383): [[decisions/mass-flux-bound-rho]], open. The method gate scores `rhoClipL1` and reports `rhoClipFraction` without a score (4bff922, [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3867-L3878)).

## Open questions

1. Re-establish the basis on the repaired case: run rhoLENT with and without the bound on `translatingDroplet2D` after 440107f, and read `rhoClipL1`, the mass residual and the error vector ([[hubs/mass-flux]], open item 1).
2. A clip fraction that counts round-off is not a metric. A threshold above round-off, or the clipped mass alone, is needed before the clip can be scored.
3. The two Eulerian and SL solvers share the clip; the Eulerian solver logs it but writes no droplet CSV yet ([[concepts/eulerian-solver-mass-flux-port]]).

## Related

- Hub: [[hubs/mass-flux]]. Model: [[models/mass-flux]].
- Siblings: [[concepts/rholent-mass-flux]], [[concepts/ddt-scheme-pairing-bdf2]], [[concepts/mass-flux-projection]], [[concepts/alphaf-source-donor-plane]].
- Decisions and retractions: [[decisions/mass-flux-bound-rho]], [[decisions/momentum-schemes-bdf2-upwind]], [[retractions/closed-box-translating-droplet]].
- Method: [[concepts/wrong-setup-voids]].

## Log

### 2026-09-28
Created.
