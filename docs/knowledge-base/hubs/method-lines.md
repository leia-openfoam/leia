---
title: "The method lines and their status"
description: "Hub: the seven level-set method lines of leia (quadratic SL in production; linear SL, SDPLS, normal-projected SL, geometric redistancing, velocity extension, gradient control) with their status and documents"
kind: hub
status: settled
part: all
tags: [hub]
date: 2026-09-28
---
# The method lines and their status

[[index]] <- back

## The question this part answers

Which transport and gradient-control lines exist in leia, which one is production, which are closed, and where each one is documented?

## Current verdict (2026-09-28)

| line | status | one line | documents |
|---|---|---|---|
| Quadratic semi-Lagrangian level set | production | third-order transport on hexahedra, coupled two-phase solver; the reference of every gate | [[studies/sl-quadratic-pre-print]], [[hubs/advection]] |
| Linear semi-Lagrangian level set | measured | order about 1.1; the value fit must be quadratic | [[studies/sl-linear-pre-print]], [[concepts/linear-semi-lagrangian]] |
| Eulerian FV transport | reference | robust second at 20 times the error; the SDPLS sources live here | [[studies/method-comparison]], [[concepts/eulerian-fv-transport]] |
| SDPLS (source-term distance control, Eulerian) | kinematic positive, coupled negative | R converges kinematically; coupled it amplifies the curvature-seeded current 260 times | [[studies/sdpls-pre-print]], [[concepts/sdpls-source-eulerian]] |
| Normal-projected semi-Lagrangian | closed | the trace is clean; every write-back keyed on the fitted normals diverges | [[studies/npsl-design]], [[concepts/normal-projected-sl]] |
| Geometric redistancing (GRL) | closed | second order static; injurious under transport | [[studies/grl-pre-print]], [[concepts/redistancing-geometric-grl]] |
| Velocity extension (legacy models) | dominated | worst accuracy at 12 to 27 times the cost on pure advection | [[studies/velocity-extension-pre-print]], [[models/velocity-extension]] |
| Gradient control (SL source step, halo-limited extension) | first campaign failed | every candidate fails the 2D gate; the source discretisation is the first repair | [[studies/gcls-pre-print]], [[hubs/gradient-control]] |

## Map

| kind | note | status | one line |
|---|---|---|---|
| hub | [[hubs/advection]] | settled | transport and phase indicator |
| hub | [[hubs/viscosity]] | settled | face viscosity |
| hub | [[hubs/surface-tension]] | settled | capillary force and the parasitic current |
| hub | [[hubs/mass-flux]] | settled | rhoLENT and the density ratio |
| hub | [[hubs/gradient-control]] | open | sources and extensions |
| hub | [[hubs/verification]] | settled | the method |

## Open, in order

1. The gradient-control repairs ([[hubs/gradient-control]]).
2. The parasitic-current source (the curvature estimator) and its amplifiers ([[hubs/surface-tension]]).

## Log

### 2026-09-28
Created.
