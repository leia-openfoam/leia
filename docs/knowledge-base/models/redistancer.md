---
title: "redistancer: noRedistancing, PDE, anchoredEikonal, planeFootWave"
description: "The redistancing line is closed: planeFootWave is second order on a static circle, but the frozen-band variant inflates the volume error from 0.017 to 1.77 at N = 256 and the PDE reinitialisation diverges; production runs noRedistancing."
aliases: [redistancer, REDISTANCER, reinitialisation]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-26
decided_by: [config/redistanceStatic2D.yaml, config/benchVortexGRLfrozen.yaml, "author decision 2026-08-07, PCS L323-L333"]
code: [src/leiaLevelSet/redistancer/redistancer.H, src/leiaLevelSet/redistancer/pdeRedistancer.H, src/leiaLevelSet/redistancer/anchoredEikonalRedistancer.H, src/leiaLevelSet/redistancer/planeFootWaveRedistancer.H]
sources: ["METHOD 9 item 4 (L777-L790)", "MC article sec:frozen (L304-L341)", "GRL article sec:static and sec:idempotency (L355-L411)", "GRL article sec:negative (L448-L468)", "PCT dead ends (L599-L604)", "PCS WP6 and item 7 (L302-L344, L378-L381)", "DP L42-L56 and L486-L498"]
---
# redistancer: noRedistancing, PDE, anchoredEikonal, planeFootWave

> Verdict (2026-09-28). The level set is initialised as a signed distance and is not maintained as one; no reinitialisation runs in the production method ([METHOD 1 L56-L57](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L56-L57), [DP L42-L45](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L42-L45)). The geometric redistancer `planeFootWave` is second order on the static gate (band error 4.48e-05 at h = 1/256, order 2.0) and exact on a plane to 1e-13 ([static_redistance_orders.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/static_redistance_orders.tex), [redistanceStatic2D_orders.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/redistanceStatic2D_orders.tex)). The line is closed anyway: a redistancer rebuilds distance to whatever zero set it is handed, including advection noise, so even the frozen-band variant with zero interface displacement inflates the volume error from 0.017 to 1.77 at N = 256 ([MC article L315-L341](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)). The PDE reinitialisation with the central Hamiltonian diverges (bulk gradient error 0.56 to 59.7) ([DP L495-L497](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L495-L497)). Redistancing of any kind, any cadence, is a listed dead end ([PCT L599-L604](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-combined-source-terms.md#L599-L604)).

## What it is

The key is `levelSet.redistancer.type`; the code default is `noRedistancing`, which is the base class ([redistancer.C L100](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.C#L100), [redistancer.H L170](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L170)). The base class owns the trigger: `interval` (every `redistanceInterval` steps), `gradPsiThreshold` (the band L1 mean of `abs(abs(grad psi) - 1)` over `{abs(psi) < 6h}` above `(h/L)^2`, after Abu-Al-Saud, Popinet and Tchelepi) and `signedDistanceBounds` (fire when the band `abs(grad psi)` leaves `[1/g, g]`) ([redistancer.H L33-L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L33-L49), [L88-L124](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L88-L124)). Tokens: `REDISTANCER`, `REDIST_TRIGGER`, `REDIST_THRESHOLD`, `REDIST_INTERVAL`, `REDIST_GRADBOUND`, `REDIST_REGION`, `REDIST_FREEZE`, `REDIST_ANCHOR_LAYERS`, `REDIST_HAMILTONIAN`, `REDIST_NITER` ([DP L42-L56](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L42-L56), [DP L486-L498](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L486-L498)). The solvers that build it are `leiaLevelSetFoam`, `leiaRedistancedLevelSetFoam` and `leiaLevelSetTwoPhaseFoam` ([createFields.H L64](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetFoam/createFields.H#L64), [L208](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L208)); the SL two-phase solver has none.

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| no redistancing | `noRedistancing` | settled, production | The base class; `correct()` is a no-op. | [redistancer.H L170](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L170), [DP L42-L45](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L42-L45) |
| PDE reinitialisation | `PDE` | retracted | Sussman-Smereka-Osher pseudo-time equation; a fixed pseudo-step is a mesh-dependent stability bomb; the central Hamiltonian injures saturated profiles; it increases the static band error at every resolution. | [pdeRedistancer.H L30-L38](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/pdeRedistancer.H#L30-L38), [GRL article L463-L467](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L448) |
| anchored Eikonal fill | `anchoredEikonal` | retracted (comparison model) | Plane anchors plus Tucker's advection-diffusion linearisation of the steady Eikonal equation; first order in the band on a plane (1.84e-03 at h = 1/256). | [anchoredEikonalRedistancer.H L30-L48](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/anchoredEikonalRedistancer.H#L30-L48), [redistanceStatic2D_orders.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/redistanceStatic2D_orders.tex) |
| plane-foot wave | `planeFootWave` | retracted (line closed) | Plane anchors plus a FaceCellWave nearest-donor fill that evaluates the continuous donor plane; second order static; one-signed O(h^2 kappa) displacement per event. | [planeFootWaveRedistancer.H L30-L44](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/planeFootWaveRedistancer.H#L30-L44), [GRL article L398-L411](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L385) |

The dictionary words are the `TypeName` strings ([PDE registration](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/pdeRedistancer.C#L55), [anchoredEikonal L94](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/anchoredEikonalRedistancer.H#L94), [planeFootWave L73](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/planeFootWaveRedistancer.H#L73)).

The frozen-band variant of `PDE` (`REDIST_FREEZE true`, `REDIST_ANCHOR_LAYERS 2`, `REDIST_HAMILTONIAN godunov`) keeps the transported sign-change band and two guard layers as Dirichlet anchors and refills the bulk only ([DP L489-L497](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L489-L497), [benchVortexGRLfrozen.yaml L1-L11](https://github.com/leia-openfoam/leia/blob/8867581/config/benchVortexGRLfrozen.yaml#L1-L11)).

## Why it matters

The zero level set and the band values are the same degrees of freedom, so any band rewrite displaces the interface by O(h^2 kappa) per event, and the displacement compounds ([MC article L238-L241](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L238-L241)). The quadratic value fit has no maximum principle, so the transported field carries advection noise that a redistancer entrenches ([METHOD 9 L786-L790](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L786-L790)). The signed-distance property is maintained instead by a source term that vanishes on the zero set ([[models/sdpls-source]], [[hubs/gradient-control]]).

## Where in the code

- Family: `src/leiaLevelSet/redistancer/`, library `libleiaRedistancer`.
- The plane anchors shared with the phase indicator and the geometric mass flux: `planeAnchors.{H,C}` ([anchoredEikonalRedistancer.H L31-L35](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/anchoredEikonalRedistancer.H#L31-L35)).
- The parallel-correct topological mask: [redistancer.H L133-L142](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L133-L142).
- Unit test: `applications/test/leiaTestRedistance`; cases `2DredistanceStatic`, `2DredistanceCircle`, `1DredistanceTest`.

## Evidence

| claim | number | where |
|---|---|---|
| Static gate, tanh profile, band L_inf after one event at h = 1/256 | planeFootWave 4.48e-05 (order 2.00), anchoredEikonal 3.53e-04 (1.98), PDE 4.07e-02 (the error grows from 2.03e-02) | MEASURED (L_inf, the record has no L2), [static_redistance_orders.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/static_redistance_orders.tex) |
| One event on an exact signed distance: plane | planeFootWave 1.61e-13 band error, volume change 2.36e-16; anchoredEikonal 1.84e-03, order 1.00 | MEASURED (L_inf), [redistanceStatic2D_orders.tex](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/data/tables/redistanceStatic2D_orders.tex) |
| One event on an exact circle: volume change per event | 1.11e-05 at h = 1/32 to 1.17e-07 at h = 1/256, ratio about 4 per level; `E_geom = E_vol` exactly (one-signed chord bias) | MEASURED, [GRL article L398-L407](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L385) |
| Frozen-band bulk-only PDE with Godunov fill on the Eulerian line, T = 8 | E_vol 0.017 to 1.77 at N = 256, 0.238 to 2.53 at N = 128; alpha integral inflated to 2.8x the droplet | MEASURED, [MC article L315-L334](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304) |
| Central Hamiltonian PDE fill | bulk gradient error 0.56 to 59.7, divergent; Godunov fill 0.56 to 0.009 per event | MEASURED, [DP L495-L497](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L495-L497), [MC article L247-L249](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L238) |
| Foot-point-cloud distances scallop the gradient | the point-cloud variant raised the band deviation an event must lower | MEASURED (qualitative), [GRL article L454-L457](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L448) |
| Plane-based band rewrite at any cadence | O(h^2 kappa) one-signed compounding; retired with measurements | MEASURED, [PCS L378-L381](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L378-L381) |

## Why it failed, or why we think so

Eulerian `linearUpwind` advection leaves bounded but nonzero over- and undershoots in the bulk: spurious sign changes far from the interface. The frozen anchors keep the true band, but the bulk fill rebuilds distance to the spurious zero crossings and the geometric indicator reads them as phase ([MC article L327-L341](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)). The injury does not come from moving the interface; it comes from rebuilding distance to advection artefacts. A smooth source that vanishes on the zero set cannot manufacture an interface, which is why the SDPLS source survives ([MC article L339-L341](https://github.com/leia-openfoam/leia/blob/8867581/docs/method-comparison/method-comparison-article/methodComparison.tex#L304)).

## Decisions

- `REDISTANCER noRedistancing` ([DP L42-L45](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L42-L45)); the gated frozen-band reset of PCS WP6 was promoted in v0.2 and never run in the coupled matrix ([PCS L302-L344](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L302-L344)).
- Binding rule (author, 2026-08-07): psi advection is never modified or suspended for a static case; any reset must pass the kinematic suite with the candidate active ([PCS L323-L333](https://github.com/leia-openfoam/leia/blob/8867581/docs/plan-curvature-stabilization.md#L323-L333)).

## Open questions

1. The GRL article's advected sections (reversed vortex, 3D shear and deformation, trigger ablation) are TODO ([GRL article L433-L445](https://github.com/leia-openfoam/leia/blob/8867581/docs/geometrically-redistanced-levelset/grl-level-set-article/geometricallyRedistancedLevelSet.tex#L433-L445)).
2. The static gate reports L_inf only; L2 and L1 are not recorded.

## Related

[[hubs/advection]] - [[hubs/gradient-control]] - [[concepts/redistancing-geometric-grl]] - [[models/phase-indicator]] - [[models/narrow-band]] - [[models/sdpls-source]] - [[models/level-set-advection]] - [[concepts/eulerian-fv-transport]] - [[studies/grl-pre-print]] - [[studies/method-comparison]] - [[decisions/psi-filter-none]]

## Log

### 2026-09-28
Created from the GRL article and tables, the method comparison article, the plan dead-end lists and the code.
