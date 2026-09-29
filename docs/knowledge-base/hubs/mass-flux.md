---
title: "Mass flux and density-ratio consistency"
description: "Hub: the rhoLENT mass-momentum-consistent flux, the closed-box void, the coupled-face density defect, the Eulerian port, and the late instability of the translating droplet"
kind: hub
status: settled
part: mass-flux
tags: [hub, part/mass-flux]
date: 2026-09-28
---
# Mass flux and density-ratio consistency

[[index]] <- back

## The question this part answers

How is the mass flux $\rho_f \phi$ of the momentum equation made consistent with the phase indicator at a density ratio of 840, so that the translating droplet keeps its velocity and the spurious current does not grow?

## Current verdict (2026-09-28)

The production flux is rhoLENT ([[concepts/rholent-mass-flux]], [[decisions/mass-flux-rholent]]): an auxiliary density transported with the volumetric flux, its face density $\rho_f$ built from the interface plane of the donor cell (`alphaFSource donorPlane`, [[concepts/alphaf-source-donor-plane]]), reset to the indicator density after the outer loop, with a relative mass residual of $10^{-13}$ ([METHOD 6](https://github.com/leia-openfoam/leia/blob/8867581/METHOD.md#L276-L320)). `boundRho` is active, but its evidence is VOID and the decision is open ([[concepts/bound-rho]]). The momentum time scheme is BDF2 everywhere and the pairing argument, not a table, decides it ([[concepts/ddt-scheme-pairing-bdf2]]). The density ratio is the amplifier of the parasitic current and the mass-momentum consistency is not the dominant term (dec002f, [[concepts/density-ratio-amplifier]]). Three defects shaped the record: the translating droplet case was a closed box until 2026-09-02, which voided every translating result before 440107f ([[retractions/closed-box-translating-droplet]], [[concepts/wrong-setup-voids]]); the face density on coupled faces was rank-local until 28d13f0 and the droplet metrics counted internal faces only until b1798c3, so every SL two-phase result on more than one rank before 2026-09-27 carries a decomposition error ([[concepts/coupled-face-density-defect]]); the Eulerian two-phase solver kept $\rho$, $\rho\phi$ and $\mu_f$ frozen until c094bd8 and moved the heavy droplet at 42 % of the stream ([[concepts/eulerian-solver-mass-flux-port]]). With the repaired case and binaries the translating droplet converges at order 1.7 to 3.1 in shape and travel over 0.1 s in a 20 mm box, and it shows a late instability: an outlet-triggered fast growth, and a slower interior growth that shrinks with $h$ between two rungs ([[cases/translating-droplet]], [STATUS 11.15](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3745-L4023)).

## Map

| kind | note | status | one line |
|---|---|---|---|
| model | [[models/mass-flux]] | settled | interpolatedDensity, geometricFaceDensity, rhoLENT and the sub-switches |
| concept | [[concepts/rholent-mass-flux]] | settled | the formulation and its measured basis |
| concept | [[concepts/bound-rho]] | open | active, evidence void |
| concept | [[concepts/alphaf-source-donor-plane]] | settled | donorPlane by author decision; the alternatives are void |
| concept | [[concepts/mass-flux-projection]] | voided | the curl-free correction, pre-fix |
| concept | [[concepts/ddt-scheme-pairing-bdf2]] | settled | BDF2 everywhere by the matching argument |
| concept | [[concepts/eulerian-solver-mass-flux-port]] | settled | frozen rho until c094bd8; the shared headers |
| concept | [[concepts/coupled-face-density-defect]] | settled | fixed 28d13f0 and b1798c3; the re-run question is open |
| concept | [[concepts/density-ratio-amplifier]] | settled | dec002f; ratio 1 completes |
| decision | [[decisions/mass-flux-rholent]], [[decisions/mass-flux-alphaf-donor-plane]], [[decisions/mass-flux-bound-rho]], [[decisions/momentum-schemes-bdf2-upwind]] | settled / open | the tokens |
| retraction | [[retractions/closed-box-translating-droplet]] | voided | every translating result before 440107f |
| retraction | [[retractions/mass-momentum-consistency-dominant-term]] | retracted | falsified by dec002f |
| retraction | [[retractions/late-translating-instability-is-the-outlet]] | retracted | the interior growth remains |
| case | [[cases/translating-droplet]], [[cases/popinet-translating-droplet]] | settled | the cases |

## Open, in order

1. `boundRho`: re-establish its basis on the repaired case ([[decisions/mass-flux-bound-rho]]).
2. Which earlier parallel SL two-phase studies to re-run or void after 28d13f0 and b1798c3 (author decision, [STATUS 11.14](https://github.com/leia-openfoam/leia/blob/8867581/STATUS.md#L3624-L3744)).
3. Void the five Eulerian two-phase studies made with the frozen density (author decision, plan item D-m).
4. The interior growth of the translating droplet: the third rung ($N = 200$, 40 mm box) on the cluster, then an order; a box-length token for the gates.
5. `rhoClipFraction` counts round-off clips in pure cells; it is reported, not scored.

## Log

### 2026-09-28
Created.
