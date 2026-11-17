# One plane per face: the content-scaled flux (pre-registered 2026-11-18, before any run)

Tomislav chose option A (2026-11-18): ONE bounded experiment with a content-scaled one-plane flux; the track is closed
if it fails. leia branch `davof-ccs-oneplane` (head 260f25a8 at the start).

## Method change (one key, `sweptContentScale none|upwind`, default none; S1 `sweptZeroFluxFace` stays ON in every run of this experiment; S2 `sweptUpwindNoPlane none`)

For a face f whose ADJACENT upwind cell U (by the sign of phi_f) has NO own plane, and whose flux is computed by
clipping the swept volume with a donor plane D of another cell:

    V_D(U) = volume of cell U under the plane D (one cell clip of U by D, the same clip routine the positioning uses),
    L_U = alpha_U^n V_U (the liquid U holds), G_U = (1 - alpha_U^n) V_U.
    Liquid:  if V_D(U) > 0:        dVf_liquid <- dVf_clip * min(1, L_U / V_D(U))
    Gas (the mirrored statement): if V_U - V_D(U) > 0:  (phi_f dt - dVf) <- (phi_f dt - dVf_clip) * min(1, G_U / (V_U - V_D(U)))

Apply the liquid scaling when the plane claims more liquid than held, the gas scaling when it claims more gas than
held; at most one of the two acts (they cannot both be < 1). The face fraction alpha_f written from this face is scaled
by the same factor as its flux (as sweptBoundedFaceFraction scale does). No tolerance, no threshold, no other branch.
On a plane state V_D(U) = L_U to round-off, so the factor is 1 to round-off. Count per step: faces scaled, min factor,
sum of removed volume / V_k; print in the log.

## Pre-registered

- PA1 plane gate: every case of the 72-case gate that passes with the baseline one-plane setting (S1 only) passes with
  `upwind` too (L1(alpha_f) <= 1e-12, limiter silent). Falsifier: one baseline-passing case fails.
- PA2 curved, OF Co 0.2: shear2D N 64, shear3D N 32, deformation3D N 32 complete without abort or SIGFPE; the volume
  moved by the net-balance limiter per step, summed over the run, is at least 100 times below the baseline (S1 only)
  over the common steps. Falsifier: an abort, or less than 10 times.
- PA3 accuracy: E_shape within 1.10 of the pieces path (code 112c284a values in agents/diag2/stageC_tables.md) in the
  three cases. Falsifier: > 1.5.
- PA4 cost: wall clock at equal ranks below the pieces path in all three cases.

## Decision rule

PA1 and PA2 both hold -> run the 27-case matrix (three families x three N x OF Co 0.2/0.5/0.95, hex) once and report
the table with plicRDF and pieces columns (E_shape, orders, mass, min/max alpha before and after the limiter, limiter
volume, wall clock at equal ranks). PA1 or PA2 fails -> stop, no repair, no second variant; report the numbers and, if
cheap from the logs, where the first failing face is (cell states, factor, V_D(U), L_U).

Order: build -> unit check on one dumped plane step that the factor is 1 +- 1e-12 on every scaled face -> plane gate ->
three curved cases -> (matrix).

## Notes of the implementer, written before any run

- The baseline of PA1 (S1 only) is the Stage B onePlane row on the plane: S1 acted on 0 faces in the S1+S2 gate, so the
  S1-only gate equals the base gate (Stage B: 152 of 152 counted cells pass); that stored result is the baseline, no
  S1-only gate is re-run.
- V_D(U) is evaluated with the closed-form wet volume of the cell's fan-triangulated boundary against the plane (the
  routine of the prism clip: (1/3) sum (x_T - x_0) . a_T over the wet polygons, the cap contributes nothing).
- The unit check reports the raw ratio L_U / V_D(U) on the faces where both L_U and V_D(U) exceed 1e-12 V_U (a
  diagnostic filter only; the method has none).
- A negative alpha_U^n gives a negative factor as the statement is written; the number of negative factors is counted.
