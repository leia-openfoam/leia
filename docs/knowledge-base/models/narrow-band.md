---
title: "narrowBand: none, empty, signChange, neighbours, distance, phaseIndicator"
description: "signChange is the default band everywhere; the neighbour dilation was decomposition-dependent until it exchanged once per layer; the band metrics use the unlimited gradient since 2026-09-26."
aliases: [narrowBand, NARROW_BAND, NARROWBAND]
kind: model
status: settled
part: advection
tags: [model, part/advection]
date: 2026-09-28
date_settled: 2026-09-26
decided_by: ["author decision, cases/default.parameter L562-L564", config/seamConsistency3Dserial.yaml, config/seamConsistency3Dpar4.yaml]
code: [src/leiaLevelSet/narrowBand/narrowBand.H, src/leiaLevelSet/narrowBand/signChangeNarrowBand.C, src/leiaLevelSet/narrowBand/neighbourNarrowBand.C, src/leiaLevelSet/narrowBand/distanceNarrowBand.C, src/leiaLevelSet/narrowBand/phaseIndicatorNarrowBand.H]
sources: ["CLAUDE 4-rank section (L260-L272)", "STATUS 4 psi-filter seam bug (L892-L917)", "STATUS 11.3 (L3261-L3269)", "DP L562-L564"]
---
# narrowBand: none, empty, signChange, neighbours, distance, phaseIndicator

> Verdict (2026-09-28). The narrow band marks the cells near the zero level set with a `{0, 1}` field named `NarrowBand` ([narrowBand.H L30-L31](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/narrowBand.H#L30-L31)). `signChange` is the historical default and the token default in every case ([DP L562-L564](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L562-L564)); the phase indicator, the redistancer anchors and the droplet band metrics read it. Two seam defects lived in this family or next to it: the neighbour dilation looped local topology only, so the band was torn at every processor boundary ([neighbourNarrowBand.C L54-L69](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/neighbourNarrowBand.C#L54-L69)), and the psi filter's own band dilation looped internal faces only, so the filtered cell set depended on the decomposition ([STATUS L903-L904](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L903-L904), [CLAUDE L265-L267](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L265-L267)). The band metrics of the droplet CSV use the unlimited `gradPsiMetric` since 2026-09-26 ([STATUS L3261-L3269](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3261-L3269)).

## What it is

The key is `levelSet.narrowBand.type`; the code default is `none` ([narrowBand.C L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/narrowBand.C#L49)). The kinematic solvers, the redistanced solver and the Eulerian two-phase solver build it with `narrowBand::New` ([leiaLevelSetFoam createFields.H L67](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetFoam/createFields.H#L67), [leiaSemiLagrangeLevelSetFoam L65](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangeLevelSetFoam/createFields.H#L65), [leiaLevelSetTwoPhaseFoam L217](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaLevelSetTwoPhaseFoam/createFields.H#L217)); the SL two-phase solver recomputes it on every outer corrector ([slAlphaEqn.H L172-L225](https://github.com/leia-openfoam/leia/blob/8867581/applications/solvers/leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H#L172)).

Two token spellings exist. `NARROW_BAND` renders `type` in `cases/1Dstretch` ([fvSolution.template L152](https://github.com/leia-openfoam/leia/blob/8867581/cases/1Dstretch/system/fvSolution.template#L152)); `NARROWBAND` with `NEIGHBOURS` renders `type` and `n` in `cases/3Ddeformation` and `cases/ellipsoidDroplet3D` ([3Ddeformation L75-L76](https://github.com/leia-openfoam/leia/blob/8867581/cases/3Ddeformation/system/fvSolution.template#L75-L76)); the droplet and vortex templates hardcode `signChange` ([stationaryDroplet2D L277-L280](https://github.com/leia-openfoam/leia/blob/8867581/cases/stationaryDroplet2D/system/fvSolution.template#L277-L280), [2Dvortex L126-L129](https://github.com/leia-openfoam/leia/blob/8867581/cases/2Dvortex/system/fvSolution.template#L126-L129)).

## Members

| member | dictionary word | status | verdict in one line | evidence |
|---|---|---|---|---|
| base, no band | `none` | settled | The base class; the code default of `narrowBand::New`. | [narrowBand.H L84](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/narrowBand.H#L84), [narrowBand.C L49](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/narrowBand.C#L49) |
| empty band | `empty` | settled | Owns the registered `NarrowBand` field, all zero; the base of the concrete bands. | [emptyNarrowBand.H L61](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/emptyNarrowBand.H#L61), [emptyNarrowBand.C L37-L50](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/emptyNarrowBand.C#L37-L50) |
| sign change | `signChange` | settled, production | Marks both cells of every face whose owner and neighbour psi values have opposite sign; reads the synced psi and the coupled-patch neighbour values. | [signChangeNarrowBand.C L42-L89](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/signChangeNarrowBand.C#L42-L89) |
| dilated band | `neighbours` (key `n`) | settled | The sign-change seed dilated `n` times through `cellCells`; one halo exchange per layer since the seam fix. | [neighbourNarrowBand.H L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/neighbourNarrowBand.H#L60), [neighbourNarrowBand.C L40-L90](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/neighbourNarrowBand.C#L40-L90) |
| distance band | `distance` (key `ncells`, default 5) | settled, assumes a signed distance | Marks `abs(psi)/deltaX_min < ncells`; a missing `mag()` once marked the whole interior phase. | [distanceNarrowBand.C L40-L73](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/distanceNarrowBand.C#L40-L73) |
| indicator-seeded band | `phaseIndicator` (keys `alpha`, `alphaTol`) | candidate | Seeded where alpha jumps across a face or the cell is partially filled, then dilated as `neighbours`; built for the extension band of the gradient-control line. | [phaseIndicatorNarrowBand.H L21-L84](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/phaseIndicatorNarrowBand.H#L21-L84), [phaseIndicatorNarrowBand.C L28-L29](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/phaseIndicatorNarrowBand.C#L28-L29) |

The dictionary words are the `TypeName` strings ([signChange L57](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/signChangeNarrowBand.H#L57), [distance L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/distanceNarrowBand.H#L60), [phaseIndicator L115](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/phaseIndicatorNarrowBand.H#L115)).

## Why it matters

The band decides where the phase indicator is evaluated ([geometricPhaseIndicator.C L123-L126](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/phaseIndicator/geometricPhaseIndicator.C#L123-L126)), where the GRL redistancer anchors its planes ([planeFootWaveRedistancer.H L40-L44](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/planeFootWaveRedistancer.H#L40-L44)) and where the droplet band metrics are read. A band that depends on the decomposition makes every band-restricted number decomposition-dependent, which is the defect class of [[concepts/seam-checks-and-decomposition-invariance]].

## Where in the code

- Family: `src/leiaLevelSet/narrowBand/`, part of `libleiaCore`.
- The seed exchange over coupled patches: [signChangeNarrowBand.C L67-L89](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/signChangeNarrowBand.C#L67-L89).
- The per-layer exchange of the dilation: [neighbourNarrowBand.C L54-L69](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/neighbourNarrowBand.C#L54-L69).
- The redistancer's own topological mask is parallel-correct, but its firing-criterion dilation is rank-local ([redistancer.H L133-L136](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.H#L133-L136), [redistancer.C L293-L294](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/redistancer/redistancer.C#L293-L294)).

## Evidence

| claim | number | where |
|---|---|---|
| The psi filter's band dilation and patch types made the filtered result decomposition-dependent | max abs U 2.43e-04 before the fix, 1.72e-06 after, 1.55e-06 with the filter off (3D, N_L = 60, 20 steps, serial against np 4) | MEASURED, [STATUS L906-L916](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L906-L916) |
| The neighbour dilation marked different cells on 4 and on 16 ranks before the per-layer exchange | stated in the fix comment; no table | MEASURED (qualitative), [neighbourNarrowBand.C L56-L62](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/neighbourNarrowBand.C#L56-L62) |
| A psi-based extension band is circular: the whole gradient error of 1Dstretch was the fade region | band mean abs(dpsi/dx) 0.942131, 0.999111, 1.000000, 1.000000 at nLayers 3, 8, 32, 200 | MEASURED, [phaseIndicatorNarrowBand.H L41-L51](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/phaseIndicatorNarrowBand.H#L41-L51) |
| Band metrics with the limited gradient against `gradPsiMetric` | all four CSVs identical at tolerance 0 over 30 steps on a smooth band | MEASURED, [STATUS L3265-L3269](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3265-L3269) |
| The `distance` band presupposes `abs(grad psi) = 1` | the marked width is not `ncells` cells where the gradient leaves 1 | DERIVED, [distanceNarrowBand.C L62-L65](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/distanceNarrowBand.C#L62-L65) |

## Decisions

- `NARROW_BAND signChange` stays the default so that every existing case is bit-unchanged ([DP L562-L564](https://github.com/leia-openfoam/leia/blob/8867581/cases/default.parameter#L562-L564)).
- Seam checks on every band change: the serial-against-np-4 pattern `config/seamConsistency3D{serial,par4}.yaml` ([CLAUDE L283-L286](https://github.com/leia-openfoam/leia/blob/d1e3414/CLAUDE.md#L283-L286)).

## Open questions

1. The sign-change band can miss a cell cut at a corner when no neighbouring centre changes sign ([signChangeNarrowBand.C L58-L60](https://github.com/leia-openfoam/leia/blob/8867581/src/leiaLevelSet/narrowBand/signChangeNarrowBand.C#L58-L60)); the `phaseIndicator` band catches it through the partial-fill test, and no ladder compares the two.
2. Two token spellings (`NARROW_BAND`, `NARROWBAND`) render the same dictionary entry in different templates.

## Related

[[hubs/advection]] - [[models/phase-indicator]] - [[models/redistancer]] - [[models/velocity-extension]] - [[concepts/seam-checks-and-decomposition-invariance]] - [[concepts/halo-limited-extension]] - [[retractions/psi-filter-seam-bug]] - [[decisions/psi-filter-none]]

## Log

### 2026-09-28
Created from the family sources, the seam-defect entries of STATUS and CLAUDE, and the templates.
