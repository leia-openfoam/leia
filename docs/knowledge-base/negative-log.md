---
title: "Negative-findings log"
description: "One line per route that was tried and failed on a measurement, in time order: what failed, the number, the gate, and the condition under which it may be re-opened"
kind: index
status: settled
part: all
tags: [index]
date: 2026-10-13
---
# Negative-findings log

[[index]] <- back

One line per route that was tried and failed on a measurement:
`date · [[note]] · NEGATIVE route: number · gate or config · re-open only if: condition`.
A negative finding is a result. Before a route is proposed again, its line here is read; the route
is re-opened only when the stated condition holds, or with a new reason written next to it. Lines
are never edited; a correction is a new line (`CORRECTED`, `RE-OPENED`) that links the old one. A
note whose main finding is negative carries `outcome: negative` and must have a line here
(`check_kb.py` enforces it). "(recorded)" marks a line whose date is the date of the note, because
the record gives no date for the measurement. Query it with
`python3 docs/knowledge-base/.quartz/kb.py negative`.

## 2026-07

- 2026-07-10 · [[concepts/idec-defect-correction-failure]] · NEGATIVE the defect-corrected IDW reconstruction (iDEC) and its two iterative variants: the band gradient defect grows 39, 1.4e4, 4.7e9 under refinement (spectral radius above one) · `config/idwCompare2D.yaml` · re-open only if: an operator with spectral radius at most one
- 2026-07-30 · [[concepts/redistancing-geometric-grl]] · NEGATIVE geometric redistancing (GRL) under transport: the criterion fires every step and the volume error runs to about 6; the frozen-band variant raises the volume error from 0.017 to 1.77 at N = 256 · `config/bulkVortexGRL.yaml`, `config/benchVortexGRLfrozen.yaml` · re-open only if: an event that does not injure the field under transport

## 2026-08

- 2026-08-03 · [[concepts/normal-projected-sl]] · NEGATIVE the normal-projected SL write-back: every per-step write-back of a fit-derived offset diverges at x1.7 per 10 steps, the m > 4 corrugation is unchanged (0.223 h against 0.209 h) · `config/nslKinematicGate.yaml` · re-open only if: an offset that is not derived from the fitted normals
- 2026-08-25 · [[decisions/psi-filter-none]] · NEGATIVE the psi band filter as a stabiliser: 5.86x better at R/h = 15.8 and 1.61x worse at R/h = 10.0, a tuning knob and not a fix · `config/filterOffAmplifier3D.yaml` · re-open only if: never in production (CLAUDE.md, no filtering in the production method)
- 2026-08-27 · [[concepts/linear-semi-lagrangian]] · NEGATIVE the linear value fit as the production transport: shape order about 1.1, CFL 1 unstable at N >= 128, the linear Taylor extrapolation's gradient defect 1.4e9 at N = 128 · `config/linearConv2Dvortex.yaml` · re-open only if: no condition; the quadratic fit is load-bearing

## 2026-09

- 2026-09-05 · [[decisions/sl-fit-normal-equations]] · NEGATIVE householderQR as a cure of the polyhedral divergence: the same divergence step with the same phase volume to six figures; the amplifying cells are well conditioned (pivot 0.757) · `config/popinet3D_La12000_poly_dump4_qr.yaml`
- 2026-09-09 · [[decisions/mesh-family-hexahedral]] · NEGATIVE pMesh polyhedra for the quadratic SL fit: Lambda_max 1.2608 against 1.0527 (blockMesh), 1.0549 (cfMesh hexahedra) and 1.0566 (snappyHexMesh with layers); the growth exponent rho(B) - 1 is 2.3x larger · `workflow/scripts/fit_amplification_probe.sh` · re-open only if: a reconstruction whose amplification bound stays near the hexahedral value on polyhedra
- 2026-09-10 · [[concepts/value-bounds-and-clips]] · NEGATIVE bounds and clips of the reconstructed SL value as a fix: the clip costs +30.4 % volume error and, with the extremum exemption, fails gate G4 (step 506 against 527); the distance-cone bound is 14.1x worse on the one-way translation at N = 64; hard limiters collapse the order · `config/popinet3D_poly_sigma0_clipGate.yaml`, `config/advConv2Dtranslation.yaml`
- 2026-09-27 · [[concepts/halo-limited-extension]] · NEGATIVE the halo-limited extension at R = h as gradient control (HL0): the shear shape error 14x at order 1.39 against 3.07, the coupled translating and oscillating arms destroyed · `config/gates/methodGate2D.yaml` · re-open only if: E2.3 shows that the trace alone carried the order loss
- 2026-09-27 · [[concepts/closest-point-extension]] · NEGATIVE the closest-point extension (FP0) as production: the same order loss as HL0 (order 1.38), decomposition-dependent (kinematic seam 0.45, coupled seam 1.11) · `config/gates/methodGate2D.yaml`
- 2026-09-27 · [[concepts/why-the-candidates-failed]] · NEGATIVE the velocity-free laws with a centred q in the SL source step (S1, HL1q, HL1z, HL2): S1 shear shape 410x (order 0.06), HL1q 173x (order 0.02) and the stationary droplet destroyed, HL2 diverges in the translating arm at every rung · `config/gates/methodGate2D.yaml` · re-open only if: a one-sided q and a smooth band (E2.1, E2.2 of [[concepts/gradient-control-next-experiments]])
- 2026-09-28 (recorded) · [[concepts/sdpls-source-eulerian]] · NEGATIVE the SDPLS strain source R under coupling with the production curvature: the droplet drains at N = 32 and diverges at 64 and 128, the mode-4 current is amplified 260 times · the SDPLS article, two-way coupling · re-open only if: the curvature seed is removed (with exact curvature the current is 0.000 in 12 of 12 arms)
- 2026-09-28 (recorded) · [[concepts/material-form-transport]] · NEGATIVE the material form of the source (Rdiv): kinematic orders -3.8, -1.3, -1.5; withdrawn · the SDPLS article
- 2026-09-28 (recorded) · [[models/velocity-extension]] · NEGATIVE the legacy velocity extensions on pure advection: closestPoint 1.44e-03 at 1.61e4 s against the SL transport 1.49e-05 at 595 s (2D vortex, N = 512, T = 8) · `config/benchVortex*.yaml`
- 2026-09-28 (recorded) · [[concepts/face-curvature-deliveries]] · NEGATIVE the face-curvature delivery lever against the blow-up: every construction that lowered the gain lost an order, and no delivery changed the blow-up exponent · the curvature-stabilization plan
- 2026-09-29 · [[concepts/kang-gfm-and-sharp-heaviside]] · NEGATIVE the Kang face weights and the sharp Heaviside: the static balance improves 6x (58x with both) but both arms fail earlier than the default
- 2026-09-29 · [[concepts/integral-surface-tension-cst]] · NEGATIVE the 2D integral CST prototype: the only discrete equilibrium of the record at N = 64, but it diverges earlier than CSF at N = 128; the line is closed

## 2026-10

- 2026-10-02 · [[decisions/flux-projection-polyhedral]] · NEGATIVE the flux projection as a cure of `SL_CLIP false` on polyhedra: the volume error at T is 7.6 / 20.1 / 25.8 without it and +0.1 to +9.4 % higher with it · `config/gates/fluxProjectionGate3DshearPoly.yaml`
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE the centred q in the SL source step, at zero velocity (gate E1, OpenFOAM, N = 100, band 3): the interface cells flatten (median psi/d 0.18 at mu t = 80; the zero crossings move up to 0.74 h), a wider band is worse (band 12: psi/d 0.0000 by mu t = 40, numpy copy), and mu dt does not change it · `config/gates/sourceZeroFlow2D.yaml`, `workflow/scripts/source_step_replica.py` · re-open only if: never for q; a centred stencil carries the sawtooth envelope onto the zero set
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE a smooth band taper as the repair of the centred q (numpy copy, validated on E1): it delays the damage (sawtooth 4.1e-5 against 1.5e-4 at mu t = 40, N = 100), and the interface cells still flatten by mu t = 80 (median psi/d 0.95, zero crossings 0.68 h) · `workflow/scripts/source_step_replica.py table` · re-open only if: as a band-edge refinement on top of a one-sided q
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE the first-order one-sided q (Godunov; least squares over the upwind face or point neighbours; numpy copy): stable, but its fixed point is first order; with psi0 = d the zero crossings move 0.70 / 0.29 / 0.39 h at N = 50 / 100 / 200 by mu t = 80 · `workflow/scripts/source_step_replica.py table` · re-open only if: the interface cells are anchored (no lambda term there: zero shift 0 by construction, field error order 1.9 to 2.0), which leaves a wrong interface gradient uncorrected
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE the gradient of the cell's own centred quadratic fit as q (numpy copy): as unstable as the centred difference (median psi/d 0.31 / 0.52 / 0.92 and zero crossings 0.93 / 0.68 / 0.84 h at N = 50 / 100 / 200, mu t = 80) · `workflow/scripts/source_step_replica.py table` · re-open only if: never; a centred stencil reverses the sawtooth group velocity
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE the fit of the single most upwind face neighbour as q (numpy copy): the selection switches between neighbours, and the zero shift does not converge (0.22 / 0.058 / 0.20 h at N = 50 / 100 / 200) · `workflow/scripts/source_step_replica.py` · re-open only if: no condition; the upwind-weighted mean (`qScheme upwindFit`) replaces it
- 2026-10-13 · [[concepts/combined-source-note]] · NEGATIVE the combined source F = a + lambda (1 - q) and both single terms in the SL source step with the one-sided q, band 3, M_mu = 1, on the kinematic arms of methodGate2D: the gradient band error falls to 0.14 (CS1), 0.58 (a alone) and 0.30 (lambda (1 - q) alone) of the baseline, but the shape order at T falls from 3.25 to 0.88 / 1.24 / 0.26 and the volume error at T rises from 2.2e-4 to 1e-2 to 4e-2 (N = 136) · `config/candidates/CS1.yaml`, `CS0s.yaml`, `CS0r.yaml` · re-open only if: a form of the source whose shape order on the shear arm stays within the gate tolerance of the baseline
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE a band of 3 cells for the strain term: in the exact 1D stretch (a uniform) the band mean q at T is 0.534 at N = 256 (the closed form is 1); the source on the whole domain gives 1.000003 · `config/candidates/CS0s.yaml` with `SL_SOURCE_BAND_CELLS` 3 and 1000 · re-open only if: a band design that keeps q = 1 in the exact 1D arm
- 2026-10-13 · [[concepts/source-discretisation-defects]] · NEGATIVE anchoring the interface cells (F = 0 there) as the repair of the shape loss: with the anchor, the strain term on the whole domain still has 23x the baseline's shape error at T (1.29e-2 against 5.48e-4, N = 136); the direct zero-crossing shift is not the cause · `config/gates/sourceAnchorProbe2D.yaml` · re-open only if: no condition; it is a probe
- 2026-10-13 · [[concepts/combined-source-note]] · NEGATIVE four readings of the strain term's 2D shape loss (shear arm, whole-domain source, N = 136): the time step (the shape error 1.76e-2 / 1.80e-2 / 3.00e-2 at CFL 0.5 / 0.25 / 0.125), the deformation (ratio to the baseline 18.6 / 32.1 / 15.3 at T = 1 / 2 / 3), the direction n (centred: 34.3x the baseline), the unbounded quadratic fit (with the stencil clip: 20.2x) · `config/gates/sourceDtProbe2D.yaml`, `sourceDeformationProbe2D.yaml`, `sourceDirectionProbe2D.yaml`, `sourceFitProbe2D.yaml` · re-open only if: no condition; they are probes
