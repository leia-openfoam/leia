# The swept-face flux: implementation and cluster verification (2026-10-28)

Agent report (Claude Fable 5.1) for Tomislav Maric. Fork branch `davof-swept`
of `tmaric/TwoPhaseFlow` (f7da47e, d33eafb, 4c7f3b5), leia branch
`davof-swept` (26a61898 and the tables' commit). Untracked file.

## What was built

- `isoAdvectionSwept.C`, `fluxIntegration sweptFace`: the vertices and the
  centre of a face traced backward along U by the midpoint rule (RK2, second
  order in dt; vertex velocities by volPointInterpolation, the centre's the
  cell-point interpolate in the donor, so for a uniform U the pre-image is
  the prism's f - d); the closed list of 4F triangles (face fan, pre-image
  fan reversed, two per lateral quad; outward for phi_f > 0, reversed for
  phi_f < 0; the signed volume as the check, flips counted); every triangle
  clipped by the donor's PLANE (davof::triNegativePolygon);
  V_f^- = (1/3) sum_T (p_T - x0).a_T^- and
  alpha_f^{n+1} = (a_{f'}^- . A_{f'}) / |A_{f'}|^2 from the same clips; dVf
  with the sign of phi_f. Scope in the header and etc/README-davof.md: planes
  only (a paraboloid above sweptPlaneTol 1e-6 in max|b| h, or a sphere, is a
  FatalError), one donor cell (out-of-donor faces counted), the Info line
  clause "swept faces, flips, out of donor, clipped triangles, closure max".
  prism stays the default and unchanged.
- `clippingSurface plane`: the paraboloid path with a zero tensor (the plane
  in the frame of the used normal positioned by the volume anchor; fitDof
  none or one). leia rows davof-planeSurf-exact, davof-planeSurf-ringArea.
- `davofOrbit.H`: the exact body turned about an axis (keys omega, axisPoint)
  in the exact normal model, the exactSphere surface and areaFractionError.
- leia: cases/davof/sphereRotation3D (the sphere turned once about the axis
  parallel to y through (0.75 . 0.5) by the fork's setFlow block, omega 2 pi,
  REFERENCE_SPEED 2.5), the arm rotation, the metrics parser's swept columns,
  make_davof_advection_table.py keys row+flux when a study sweeps the flux,
  three configs with the predictions P-S1..3 and the MEASURED paragraphs.

## Jobs (cluster /work/scratch/tm83tomy/davof/agents/swept, ledger .my_jobs)

55726034 build (f7da47e); 55726066/67/68 drivers cancelled by id before the
rebuild (the first cell lookup scanned every cell); 55726457 rebuild
(d33eafb); 55726458 plane, 55726460 sphere translation, 55726461 sphere
rotation: all COMPLETED (25, 9, 12 min).

## MEASURED

P-S1 (plane, 128 cases) HELD: every sweptFace case at round-off, run maxima
of L1(alphaf) 5.4e-15 (N 20, Co 0.2, normal arm) to 2.7e-13 (N 160, Co 0.5,
oblique), the prism rows of the same surface at the same magnitude (the end
values apart by 1e-15 to 2e-14 absolute); gate 128 PASS, 0 FAIL; closure
residual max 2.4e-17 to 7.8e-17; zero flips. The planeSurf rows pass too.

P-S2 (sphere translation, 8 cases): HELD for the exact normal (L1(alphaf)
identical to twelve digits in all 160 steps at N 32, end 1.152e-2; within
2.6e-12 at N 64, end 6.220e-3; L1(n), Eshape within 1e-11); FALSIFIED for
the ring normal (step 1 identical, step 2 apart by 2.4e-4 relative, the end
values within 2.0e-5 at N 32 and 1.4e-6 at N 64, the same interface-cell
counts). HYPOTHESIS: a per-face round-off difference that the L1 sum hides,
amplified by the ring normal (a cell at the wispTol threshold or a face that
the two clips put on different sides); open, needs a per-face comparison
(writeFaceError) after step 1. min(alpha) -6.5e-14 to -1.8e-12 both fluxes.

P-S3 (sphere rotation, 8 cases) FALSIFIED: at N 32 sweptFace ends ABOVE
prism, L1(alphaf) 2.384e-2 against 2.328e-2 (exact normal, +2.4 %) and
4.275e-2 against 4.233e-2 (ring, +1.0 %); at N 64 all four cases abort on
the interface-cell guard (steps 565/390 exact prism/swept, 268/270 ring),
min(alpha) -0.18 to -0.26: the plane-surface loop is not bounded under this
rotation at N 64 with either flux. Every swept face's pre-image leaves the
donor (3240 of 3240 faces every step at N 32): the one-donor clipping, shared
by both fluxes, is the binding assumption under rotation. Mass: sum(alpha V)
deviates by -2.2e-3 to -3.3e-3 % over the turn at N 32 for both fluxes
(1e-13 % under translation); cause not investigated. Zero flips; closure
below 4e-17.

Out-of-donor counts in the tables (fork d33eafb) include the non-sweeping
faces (2160 of 3240 on the translated sphere at N 32 are the faces parallel
to U); 4c7f3b5 counts the sweeping faces only (not re-run).

Cost per step, sweptFace / prism (Kromer-Bothe): 0.6 to 1.9 on the plane
(noisy, 15 ms to 7 s steps), 0.87 to 1.06 on the translated sphere, 1.12
to 1.24 on the rotation at N 32.

## Not done

No run on the laptop (build and template rendering only), so the four-rank
gate was not run; the swept routine's parallel path (a midpoint across a
processor boundary takes the Euler step) is untested. The multi-cell cut,
the curved cap, and the per-face comparison of P-S2's ring row are open.

## Rotation (2026-10-29): the plane-clipped loop under solid-body rotation

Case `cases/davof/sphereRotation3D`: the sphere turned once about the axis
parallel to y through (0.75 . 0.5), omega 2 pi, Co 0.2 at the sphere's far
point, clippingSurface plane (volume anchor), rows exact and ring normal.

### Fork additions (davof-swept 60c13fe .. 8e2f359)

- `sweptClassify`: per step the faces by the class of the cells their
  pre-image reaches; the first unbounded cell with its faces' fluxes; the
  mass through advection, bounding, clipping and the boundary patches.
- `fluxIntegration sweptFaceCells` (`davofSweptPieces.H`): the swept volume
  cut into the pieces inside the donor and its vertex-ring cells; full cells
  wet, empty dry, interface cells clipped by their own PLIC plane of t^n;
  piece volume from the clipped triangles and the k <= 3 cap areas solved
  from closure, no cut polygon. Unit check `apps/test/testSweptPieces`
  reproduces the note's `data/multicell.csv` to 1e-13 relative. Degenerate
  zero-thickness pieces (triangles lying in a cell plane) take no volume.
- `sweptScaleToPhi` (the kinematic scaling): the swept liquid volume times
  phi_f dt / V_swept. Under rotation the traced swept volume differs from
  phi_f dt by about (omega dt)^2/12 relative (4e-6 at N = 64), so a wholly
  wet swept volume into a full cell carries more than the upwind flux out of
  it, and the full cell overshoots (1 + 3.1e-7 at step 1). The prism's
  volume d.S_f is phi_f dt by construction.
- `sweptCellsFaceUpdate pieces | donor`: the face update from the pieces or
  from the donor's plane alone (the flux stays piecewise).

### Verdicts

| prediction | result |
|---|---|
| P-R0 comparators sound | HELD: plicRDF/isoAlpha bounded (min -8.9e-18 to -7.6e-13), mass 1.4e-8 % to 3e-8 % over the turn, Eshape order 1.9, L1(x) order 2.0 |
| step 2 classification | 67 % of swept faces reach an interface neighbour, 23 % empty, 8 % full, donor only 1.6 % (prism) / 0 % (traced) |
| P-R1 N = 64 abort removed | HELD: sweptFaceCells + scaling completes every turn at N = 32/64/128 |
| P-R2 bounding mass at round-off | HELD (1e-15 per step); total drift was the zeroGradient boundary inflow (P-R6) |
| P-R3 bounded / conservative | donor-face variant: excess <= 8.2e-11, min -5.2e-11 at N = 128, mass 4.8e-13 relative; pieces variant overshoots 2.9e-6 at N = 128 |
| P-R4 order | donor-face: ring 0.85 / 0.81 (held), exact 0.65 / 0.70; L1(x) 1.66 to 1.91, Eshape 1.74 to 1.88; pieces face update: 0.17 to 0.45 (falsified) |
| P-R5 at or below prism | donor-face at N = 32 below prism (2.288e-2 vs 2.328e-2); prism aborts at 64/128; pieces face update above prism (5.49e-2) |
| P-R6 boundary term | HELD: alpha fixedValue zero gives mass to 1e-12 % |
| H-F1 face update | HELD: the piecewise face update caused the face-fraction loss |

### Ladder (exact / ring normal; L1(alphaf) at the end of the turn)

| flux | N | L1(alphaf) | order | L1(x) | min alpha (after bound) | mass dev % | status | s/step |
|---|---|---|---|---|---|---|---|---|
| prism | 32 | 2.33e-2 / 4.23e-2 | | 1.53e-3 / 2.83e-3 | -8.5e-15 / -9.5e-13 | -3.8e-3 / -3.7e-3 (zG) | ok | 0.16 |
| prism | 64 | aborted | | | -0.28 / -0.18 | | step 565 / 268 | 1.04 |
| prism | 128 | aborted | | | -0.23 / -0.28 | | step 602 / 105 | 7.39 |
| sweptFaceCells, pieces face | 32 | 5.49e-2 / 6.70e-2 | | 1.50e-3 / 2.74e-3 | -7.2e-13 / -7.2e-13 | -1.1e-5 / -4.4e-5 (zG) | ok | 0.28 |
| sweptFaceCells, pieces face | 64 | 4.30e-2 / 4.90e-2 | 0.35 / 0.45 | 4.75e-4 / 6.74e-4 | -9.9e-14 / -3.9e-13 | -8.0e-6 / -3.2e-5 (zG) | ok | 1.44 |
| sweptFaceCells, pieces face | 128 | 3.83e-2 / 4.05e-2 | 0.17 / 0.28 | 1.39e-4 / 1.75e-4 | -7.8e-13 / -3.5e-12 (excess 2.9e-6 / 3.8e-7) | -7.5e-6 / -1.6e-5 (zG) | ok | 9.05 |
| sweptFaceCells, donor face | 32 | 2.29e-2 / 4.16e-2 | | 1.49e-3 / 2.81e-3 | -6.5e-13 / -3.0e-11 | 1.4e-12 / 1.5e-12 (fV) | ok | 0.28 |
| sweptFaceCells, donor face | 64 | 1.46e-2 / 2.31e-2 | 0.65 / 0.85 | 4.73e-4 / 7.59e-4 | -2.7e-13 / -1.2e-13 | 4.8e-12 / 5.8e-12 (fV) | ok | 1.45 |
| sweptFaceCells, donor face | 128 | 8.96e-3 / 1.32e-2 | 0.70 / 0.81 | 1.40e-4 / 2.02e-4 | -1.4e-11 / -5.2e-11 | 4.8e-11 / 4.8e-11 (fV) | ok | 9.08 |

zG: alpha zeroGradient on the patches; fV: alpha fixedValue zero (P-R6).

### Open

The donor-face variant is bounded and conservative but no longer
one-surface consistent (the flux uses the pieces' planes, the face update
the donor's); the first-step face error of the piecewise face update
(3e-2, independent of h) is not explained. inletOutlet for alpha needs phi
in initDavofFields. No four-rank run; the ring search is serial only.

## The first-step face error of the piecewise face update (2026-10-29/31)

Cause: under the rotation about y every y-normal face has phi_f = 0 and its
traced pre-image stays in the face's plane, which is a face plane of the
cells on both sides. The half-space test (x - x0).n <= 0 accepted that
coplanar pre-image with a zero signed distance in BOTH cells, so the face
fraction was the sum of two fractions (the donor's and the cell across).
Evidence (one-step dump, job 55759731): the 1080 y-normal faces carry 94 % of
the L1 sum, their pre-image coverage by the ring cells is 2.000 (1.000 on
sweeping faces); the worst 50 are all y-normal, phi = 0, 8 coplanar pieces,
pieces - donor = +0.50; face 87789: donor piece 0.439 + cell above 0.441 +
edge cells = 0.999 against exact 0.480. Python check with the note's
clip_polygon_wet: the ring's shares sum to 2.000 (s <= 0) and 1.000 (the
donor-side rule). Fix (fork c52a904, sweptCoplanar donorSide, default): a
coplanar part belongs to the cell on the donor's side of its plane.
One step at N = 32: 4.164e-2 before, 4.603e-3 after (donor plane 5.061e-3).
Ladder after the fix (job 55792290): first step below the donor update at
every N (P-F1 held), min(alpha) >= -3.8e-12 but the excess over 1 1e-11 to
8e-11 (P-F2 falsified on the excess), L1(alphaf) orders 0.67 / 0.70 exact
and 0.84 / 0.81 ring (P-F3 falsified); the end values equal the donor-plane
update's, so the per-turn order is set by the loop, not by the face update.

## Error accumulation under rotation (2026-11-01)

Short controls (N = 32 / 64, every step; config sphereRotation3D_hexHist*):
the base (anchored plane, sweptFaceCells) and the prism grow about linearly
over the turn (p 0.76 to 0.93, rate per unit time 1.6e-2 at N = 32, h-order
of the rate 0.45 -> per-turn order 0.67); the exact tangent plane and the
exact sphere (with the prism flux) do not accumulate; Co 0.1 gives the same
rate per unit time (1.18x); the pure spin and the translation saturate.
The accumulation needs the loop's own state-fed surface and an orbiting
interface; flux shape and step are not the cause. Hypothesis: the
orientation-dependent distortion the anchored plane imprints on alpha is
stationary under translation but has to be re-established while each
patch turns relative to the mesh under rotation. Next discriminator: a
circular translation without spin.

## Circular translation without spin (2026-11-01)

Setup (leia dd9f62fd, fork unchanged at 100eb6a, build 55838221; configs
`sphereRotation3D_hexHistOrbit` and `sphereRotation3D_hexHistOrbitExactPlane`,
predictions pre-registered in their headers before the runs). For a sphere
the rotation about the off-centre axis and the circular translation of the
sphere without spin have the same exact state at every time. Checked in the
fork: the exact normal model, the exactSphere surface and tangent plane
(isoAdvectionCurved.C) and areaFractionError take only the centre
davof::orbitOrigin(origin0, U0, omega, axisPoint, t) from davofOrbit.H, with
U0 from the dictionary (zero); nothing depends on an orientation. New token
TPF_SETFLOW_MODE (rotation | function): mode function of OpenFOAM's setFlow
(v2606 setFlow.C: U == velocity(t) on the cells and every patch, phi =
fvc::flux(U)) with velocity(t) = s omega_vec x (c(t) - axisPoint), a
Function1 sine whose level is a Function1 cosine, s = sin(theta/2)/(theta/2)
(theta = omega dt; 0.99998972 at N = 32), so that U(t_m) dt is the chord and
the step ends lie on the exact circle (the drift was removed, not
estimated). |U| = 1.571, interface Co 0.126 (the rotation: 0.05 to 0.2),
the same dt, 400 / 800 steps. The swept flux traces back with the
cell-point interpolate of U; for a uniform U the pre-image is x - U dt and
the kinematic scaling is 1 to round-off. The sphere stays 0.1 (3.2 h at
N = 32) from the walls.

Kinematic check (cluster, the written U and the alpha centroid at the four
writes): U uniform and equal to the formula at t - dt/2 to the 12 written
digits; the liquid centroid within 0.043 h (N = 32) and 0.023 h (N = 64) of
the exact centre over the turn (the rotation base: 0.084 h and 0.051 h).

Jobs 55841489 (anchor, N = 32 / 64, 17.7 min) and 55841490 (exact plane,
N = 32, 6.2 min), both COMPLETED; analysis davof_rotation_history.py on the
cluster, $A/rotation_history_orbit.csv.

L1(alphaf) at the first step / 25 / 50 / 75 / 100 % of the turn, the late
slope r over [0.25 T, T] per unit time, the growth exponents p over
[0.1 T, T] and [0.5 T, T], mass deviation at the end, min(alpha) and the
largest excess over 1 after the bounding:

| run | N | first | 25 % | 50 % | 75 % | 100 % | r | p_full / p_half | mass dev % | min alpha | excess |
|---|---|---|---|---|---|---|---|---|---|---|---|
| orbit, anchor | 32 | 5.13e-3 | 1.17e-2 | 1.34e-2 | 1.11e-2 | 9.52e-3 | -5.39e-3 | 0.02 / -1.44 | 3.7e-13 | -6.9e-13 | 3.3e-13 |
| orbit, anchor | 64 | 2.60e-3 | 6.48e-3 | 7.54e-3 | 6.12e-3 | 4.89e-3 | -3.64e-3 | -0.01 / -1.73 | 5.6e-12 | -9.9e-13 | 1.3e-12 |
| orbit, exact plane | 32 | 1.25e-2 | 9.06e-3 | 6.77e-3 | 5.28e-3 | 4.46e-3 | -5.58e-3 | decreasing | 9.2e-13 | -1.26e-2 | 1.6e-6 |
| rotation, anchor (base) | 32 | 4.60e-3 | 1.07e-2 | 1.58e-2 | 1.98e-2 | 2.31e-2 | 1.65e-2 | 0.82 / 0.80 | 1.3e-12 | -5.6e-13 | 6.3e-12 |
| rotation, anchor (base) | 64 | 2.35e-3 | 5.85e-3 | 9.13e-3 | 1.20e-2 | 1.45e-2 | 1.20e-2 | 0.93 / 0.91 | 5.6e-12 | -2.4e-13 | 1.1e-11 |
| rotation, exact plane | 32 | 1.20e-2 | 7.32e-3 | 5.13e-3 | 4.02e-3 | 3.45e-3 | -6.37e-3 | decreasing | 4.3e-13 | -1.9e-3 | 3.8e-6 |
| straight translation (x, own T = 1) | 32 | 5.79e-3 | 9.74e-3 | 1.07e-2 | 1.12e-2 | 1.15e-2 | 2.25e-3 | 0.30 / 0.21 | -6.7e-14 | -6.6e-14 | 2.6e-14 |
| straight translation (x, own T = 1) | 64 | 3.15e-3 | 5.52e-3 | 5.97e-3 | 6.15e-3 | 6.22e-3 | 8.0e-4 | 0.22 / 0.12 | -1.9e-12 | -1.1e-13 | 1.3e-14 |

The history at every eighth of the turn (anchor; orbit / rotation):
N = 32: T/8 7.71e-3 / 7.67e-3, T/4 1.17e-2 / 1.07e-2, 3T/8 1.12e-2 /
1.28e-2, T/2 1.34e-2 / 1.58e-2, 5T/8 1.05e-2 / 1.73e-2, 3T/4 1.11e-2 /
1.98e-2, 7T/8 7.95e-3 / 2.09e-2, T 9.52e-3 / 2.31e-2; the orbit's maximum
1.38e-2 at t = 0.52. N = 64: T/8 4.15e-3 / 4.08e-3, T/4 6.48e-3 / 5.85e-3,
T/2 7.54e-3 / 9.13e-3, 7T/8 4.45e-3 / 1.32e-2, T 4.89e-3 / 1.45e-2; the
orbit's maximum 7.77e-3 at t = 0.505. L1(x) at the end 7.51e-4 / 1.89e-4
(rotation 1.49e-3 / 4.72e-4).

Verdicts against the pre-registered thresholds:
- H-R HELD: r <= 0.2 r_rot at both N (-5.39e-3 against the threshold 3.3e-3
  at N = 32, -3.64e-3 against 2.4e-3 at N = 64) and p_half < 0.3 (-1.44 /
  -1.73). Secondary: L1(T)/L1(t1) = 1.86 / 1.88 (<= 2.5 predicted).
- H-D FALSIFIED: r is negative, below 0.5 r_rot at both N; p_full 0.02 /
  -0.01, outside [0.7, 1.3].
- H-O0 HELD: the exact-plane control falls over the turn (r = -5.58e-3), as
  under the rotation; its min(alpha) -1.26e-2 and L1(x) 0.15 at the end are
  the known mismatch of exact planes with the transported alpha (rotation:
  -1.9e-3, 0.28), and its interface cells grow from 416 to 1919.

What the numbers do not show as predicted: the orbit's history is not a
plateau like the straight translation's. Over the first quarter turn it
follows the rotation's to within 9 % (N = 32) and 11 % (N = 64), then
rises to a maximum near t = T/2 (2.7 / 3.0 times the first step) and falls
to 1.9 times the first step at the end; the eighth-turn samples are higher
at the quarter points (where the direction of travel is aligned with a mesh
axis) than at the neighbouring odd eighths (sampled at eighths only, so a
ripple is suggested, not resolved). Not explained.

Defect found in the error function objects (not fixed): volumeFractionError
(fork volumeFractionError.C, write()) takes U0 = U[0] whenever the velocity
field is uniform and adds U0 t to the orbit centre; under the circular
translation the field is uniform, so its exact sphere is displaced and the
orbit runs' Eshape (2.79e-2 at T/4, then 1.414e-2 = the sphere's whole
volume) is invalid. areaFractionError reads U0 from its dictionary (given,
zero) and is correct, so L1(alphaf), L1(x) and the verdicts are unaffected.
A fix (an optional U0 key as in areaFractionError) needs a fork commit and
a rebuild.

## Item 1 (2026-11-01): the U0 of volumeFractionError fixed

Agent: Claude Opus 5.5. Fork davof-swept 6315783, leia davof-swept 714f1c68.

Confirmed by reading the code: volumeFractionError::write() (fork 100eb6a)
set U0 = U[0] whenever max|U| - min|U| <= SMALL over the cells and patches
and passed it to davof::orbitOrigin, which adds U0 t to the orbit; under the
circular translation (uniform U(t)) the exact sphere was moved by
U(t) t on top of the orbit. Every other place was checked: the exact normal
model (exactNormal.C: U0 = the required key of the exactNormal dictionary),
the exactSphere surface and the planeSource exact tangent plane
(isoAdvectionCurved.C: U0 from the exactNormal dictionary, default zero),
the split-quadrature diagnostic of the debug path (the same dictionary),
areaFractionError (the key U0 when given, the field only without it) and
the apps (advectorVoF's U0 is the initial velocity field of its periodic
modes, not an exact body). volumeFractionError was the only consumer that
inferred U0 against a given orbit. The leia templates give U0 (zero on the
rotation arm) to areaFractionError and exactNormal; volumeFractionError had
none.

Fix: volumeFractionError reads the key U0 when given (precedence). Without
it, U0 is zero when omega is nonzero (the orbit already moves the body; a
uniform U under a rotation is the circular translation, whose exact body is
the orbit itself), and only without a rotation the former inference stays
as the fallback (valid for a constant uniform translation; the
sphereTranslation3D and planeTranslation3D templates give no U0 and rely on
it). The same omega guard in areaFractionError (no run changes: its
templates give U0). Not changed: the fallback's test on |U| only, and its
validity for a constant U only (a uniform U(t) without omega needs the key).
leia sphereRotation3D/system/controlDict.template gives
U0 (TPF_UX TPF_UY TPF_UZ) to volumeFractionError.

Verification (job 55844920 on the old binary 100eb6a, build 55844936 of
6315783, job 55845067 on the new one; 20-step copies of rendered cases in
$W/verify-u0, the cases' dictionaries unchanged, i.e. without the U0 key):
- straight translation (sphereTranslation3D_hexSwept case 0, no omega, the
  fallback path): volumeFractionError Time/Eshape/Emass/Ebound and
  areaFractionError.dat byte-identical before and after;
- rotation base (sphereRotation3D_hexHist case 0, omega, U not uniform):
  identical as well;
- circular translation (hexHistOrbit case 0, omega, U uniform): Eshape at
  step 20 1.068e-2 before, 3.700e-4 after (the rotation base's 3.698e-4 at
  the same step); Emass, Ebound, areaFractionError identical; the copy with
  U0 (0 0 0) added gives the same digits as the one without the key.

Reruns of the two circular-translation studies on 6315783 (jobs 55844938,
16.5 min, and 55844939, 6.2 min; the 100eb6a studies moved aside to
`*_preU0fix_100eb6a_2026-11-01`): L1(alphaf), L1(x), L1(n), devPercent,
min/max(alpha) and the interface-cell counts identical as strings to the
100eb6a runs in all 400 / 800 / 400 steps; only Eshape changed.
Eshape (now valid), at 25 / 50 / 75 / 100 % of the turn:

| run | N | Eshape 25 % | 50 % | 75 % | 100 % | L1(x) 100 % |
|---|---|---|---|---|---|---|
| orbit, anchor | 32 | 3.75e-4 | 3.81e-4 | 3.73e-4 | 3.72e-4 | 7.51e-4 |
| orbit, anchor | 64 | 9.67e-5 | 9.92e-5 | 9.56e-5 | 9.44e-5 | 1.89e-4 |
| orbit, exact plane | 32 | 6.90e-4 | 8.61e-4 | 9.91e-4 | 1.10e-3 | 1.50e-1 |
| rotation, anchor (base, 100eb6a) | 32 | 3.70e-4 | 3.88e-4 | 4.25e-4 | 4.68e-4 | 1.49e-3 |
| rotation, anchor (base, 100eb6a) | 64 | 9.50e-5 | 1.02e-4 | 1.20e-4 | 1.40e-4 | 4.72e-4 |

The orbit's Eshape stays within 1.5 % (N = 32) and 3.7 % (N = 64) of its
T/4 value at the eighth-turn samples from T/4 to T; its order between
N = 32 and 64 at the end is 1.98; the rotation's Eshape grows by 26 % /
47 % from T/4 to T.

## Item 2 (2026-11-02): davof-swept-par merged into davof-swept

Fork: merge commit 0f7bf4d (davof-swept-par fc08178 onto 6315783), no
conflicts (the branches touch disjoint hunks of isoAdvection.C/H and
isoAdvectionCurved.C; checked by reading the parallel branch's changes to
buildCellSurfaces and curvedFlux against 100eb6a's planeSource exact path:
the exact tangent planes are built before the exchanges, which carry them
like any other surface). leia: merge commit 3a09bbff (davof-swept-par
95b8fcc8), two conflicts, both resolved by keeping both sides
(sphereRotation3D.parameter: TPF_PLANE_SOURCE next to TPF_SWEPT_FACE_LOG and
TPF_WRITE_PRECISION; fvSolution.template: planeSource next to
sweptFaceLog). The davof-swept-par branches are kept.

Verification on the merged binary (build 55845299):
- serial, unchanged (job 55845300, the 20-step copies of item 1): the
  translation, the rotation base and the orbit give volumeFractionError
  (Time, Eshape, Emass, Ebound) and areaFractionError.dat byte-identical to
  6315783, and the translation and the rotation base byte-identical to
  100eb6a as well; the swept counters of the log identical.
- decomposition (studies sphereRotation3D_hexSweptParN32np{1,2,4,8,16},
  jobs 55845301-05, 4.0 to 5.5 min; davof_decomposition_compare.py against
  np 1, $W/cmp-merge): every face of the serial run has a value at every
  step in every decomposition (nMissing 0, nDup 0). Rows exact / ring:

| np | step 1 max dalpha | step 1 max dalphaf | steps 1-19 max dalpha | steps 1-19 max dalphaf | step 20 max dalpha | step 20 max dalphaf |
|---|---|---|---|---|---|---|
| 2 | 2.9e-14 / 6.8e-14 | 4.1e-14 / 5.7e-14 | 4.2e-11 / 2.5e-9 | 1.1e-10 / 5.2e-9 | 6.9e-12 / 4.3e-10 | 1.8e-11 / 1.2e-9 |
| 4 | 3.5e-14 / 6.8e-14 | 8.5e-14 / 5.7e-14 | 6.0e-11 / 2.1e-9 | 8.0e-10 / 4.3e-9 | 2.1e-11 / 3.5e-10 | 4.2e-11 / 9.7e-10 |
| 8 | 8.8e-14 / 6.7e-14 | 8.5e-14 / 5.7e-14 | 6.0e-11 / 2.1e-9 | 9.9e-10 / 4.4e-9 | 2.1e-11 / 4.3e-10 | 4.3e-11 / 1.2e-9 |
| 16 | 9.8e-14 / 6.8e-14 | 5.7e-14 / 5.7e-14 | 5.5e-11 / 1.5e-9 | 2.4e-10 / 4.4e-9 | 1.1e-6 / 6.4e-7 | 3.6e-3 / 2.7e-3 |

  The parallel report's numbers are reproduced, including the sixteen-rank
  face at step 20 (the cap-rank fallback, item 3). The merged serial np 1
  runs equal the parallel branch's np 1 runs (fork fc08178) to the last
  written digit at all 20 steps in both rows (alpha, alphaf, dVf: 0).

## Item 3 (2026-11-02): the cap set of a piece, independent of the face order

Fork davof-swept de7aea9 (build 55845642); leia 4291efb9 (token
TPF_SWEPT_CAP_SET, config sphereRotation3D_hexSweptParN32np1Crossing with
P-C1..P-C4 pre-registered before the runs).

The defect, read from davofSweptPieces.H and the parallel run's debug face
(job 55812211): sweptPieceVolume flagged a cell plane as a cap when it
crossed the polygon at its turn in the clipping order (the cell's face
list), so (i) a plane whose new edge a later plane removed again still
counted, and (ii) the clipped remnant of zero area of a triangle touching a
cell plane in one edge set a flag as well. The face at step 20 is an
x-normal face whose lateral quad of the edge y = 0.375 lies IN that plane
(the rotation about y keeps y); the edge neighbour below it
(alpha 0.063) received those triangles with a zero signed distance, and its
cap set was {x, y, z+, PLIC} (four: rank-deficient, fallback) in the serial
face order and three planes in the decomposed one (volume 8.8e-24, no
fallback).

Decision (commit message and file header of davofSweptPieces.H): make the
cap set a property of the piece, i.e. of the swept volume intersected with
the half-spaces, rather than only muting zero-volume pieces (which leaves
the order dependence for pieces of positive volume):
(a) a boundary triangle lying in a cell plane (every vertex within 1e-10 of
    its longest edge) is a contact, dropped from that cell's piece, when the
    donor lies on the other side of the plane; with the donor on the cell's
    side it bounds the piece and is not clipped by that plane (the side
    rule of the face update's coplanar rule donorSide; the triangle's own
    orientation is not used because a lateral quad of a twisted swept
    volume is a bow tie with triangles pointing to both sides);
(b) a plane carries a cap when it produced an edge of the FINAL wet polygon
    (edge labels carried through the Sutherland-Hodgman clip), longer than
    1e-12 of the triangle's longest edge, on a polygon of area above 1e-14
    of that edge squared.
Key sweptCapSet final (default) | crossing (the earlier rule, kept for the
comparison); sweptFaceLog gains the per-face fallback flag; the debug face
prints the caps, the rank status and the contacts of every piece.

Verification (jobs 55845643-48, 55845649, 55846097; $W/cmp-caps):

- P-C1 HELD. N = 32, np 2/4/8/16 against np 1 (final), 20 steps, rows
  exact / ring: max over all 20 steps of max|dalpha| 1.2e-13 to 2.5e-12 and
  of max|dalphaf| 1.4e-12 to 6.1e-12 (threshold 1e-8); step 20 at np 16:
  2.2e-12 / 2.2e-12 exact, 1.9e-13 / 1.0e-12 ring (merge binary 1.1e-6 /
  3.6e-3 and 6.4e-7 / 2.7e-3). The fallback flags are 0 on every face at
  every step in every decomposition; nMissing 0, nDup 0. The steps 1-19
  floor of the merge binary (up to 1e-9 exact, 5.2e-9 ring in alphaf) is
  gone as well: 1.4e-12 to 6.1e-12.
- P-C2 HELD. crossing on de7aea9 equals the merge binary's np 1 to the last
  written digit (alpha, alphaf, dVf; 20 steps, both rows).
- P-C3 HELD as stated, with a finding. final against crossing on the same
  binary: at step 1 no face differs above 1e-12 (max dVf 4.5e-14 h^3); the
  crossing run falls back from step 12 on (2 to 19 / 4 to 22 faces per
  step), final never (0 at every step). Before the first fallback, from
  step 5 on, faces without a fallback differ as well: 22 / 15 faces at step
  5, dVf up to 4.8e-11 / 1.8e-9 h^3, alphaf at round-off; 88 to 330 faces at
  steps 6-11, alphaf up to 1.0e-10 / 4.3e-9 (the loop's amplification). The
  debug face of the largest one (ring row, step 5, z-normal face at
  (0.390625 0.578125 0.46875), job 55846097): under crossing the interface
  cell 31948 (alpha 0.279) held a piece of contact triangles only (coverage
  1.75e-4 m^2, three caps) with the volume -5.57e-14 m^3, exactly the
  difference of the face's liquid volume (5.2045871742e-6 final,
  5.2045871185e-6 crossing); under final that cell has 2 contacts, no
  coverage, no volume. So the earlier rule also gave contact pieces a
  spurious nonzero volume where their cap system was solvable. From step 12
  the faces with a changed fallback flag carry the large differences
  (alphaf up to 2.3e-2), the others up to 3.8e-4 (propagated).
- P-C4 HELD. Full turn sphereRotation3D_hexHist (N = 32 / 64; the 100eb6a
  study moved to *_100eb6a_2026-11-02):

| N | rule | rank fallbacks (face-steps / steps with one) | L1(alphaf) end | L1(x) end | Eshape end | min alpha | excess | mass dev % |
|---|---|---|---|---|---|---|---|---|
| 32 | crossing (100eb6a) | 7552 / 308 of 400 | 2.30545e-2 | 1.49130e-3 | 4.6807e-4 | -5.6e-13 | 6.3e-12 | 1.3e-12 |
| 32 | final (de7aea9) | 0 / 0 | 2.30541e-2 | 1.49131e-3 | 4.6805e-4 | -5.8e-13 | 6.3e-12 | 4.8e-13 |
| 64 | crossing (100eb6a) | 32717 / 655 of 800 | 1.452023e-2 | 4.72148e-4 | 1.39877e-4 | -2.4e-13 | 1.1e-11 | 5.6e-12 |
| 64 | final (de7aea9) | 0 / 0 | 1.452021e-2 | 4.72147e-4 | 1.39876e-4 | -2.7e-13 | 1.1e-11 | 5.1e-12 |

  The end values change by -1.7e-5 (N = 32) and -2.0e-6 (N = 64) relative.
  Degenerate pieces (closure residual above 1e-6, volume zero) 62 / 1270
  face-steps (crossing 60 / 1082).

Not changed and open: the closure residual of the ACCEPTED pieces is not at
round-off under either rule: per step its maximum lies between 1e-11 and
1e-6 (crossing, steps without a fallback: up to 6.4e-7 / 9.9e-7; final: up to
9.99e-7 / 1.0e-6, i.e. pieces just below the 1e-6 degenerate threshold are
taken); the 1.0 of the crossing runs' column were the rank-deficient pieces.
Which pieces close only to 1e-7..1e-6 is not investigated.

## Item 4 (2026-11-02): the finer circular-translation test

Fork de7aea9 (build 55845642), leia 2caa4293 (configs, analysis script and
the predictions F-A..F-D, committed before the runs) and 1927347f; jobs
55846276 (sphereRotation3D_hexFineOrbit, 20.6 min), 55846277
(_hexFineOrbit2, 20.6 min), 55846278 (_hexFineLineDiag, 11.1 min), 55846279
(_hexFineLineX, 10.2 min); the rotation base is sphereRotation3D_hexHist on
de7aea9 (job 55845649, item 3). All serial, every step sampled (400 / 800
samples per turn of the orbit and the rotation, 200 / 400 per turn of
orbit2, 180 / 360 on the lines), the same dt everywhere (0.0025 / 0.00125),
|U| = pi/2. Analysis: davof_rotation_angle.py on the cluster
($W/angle_steps_fine.csv: per step t, omega t, the direction of travel phi,
its angle psi to the nearest mesh axis, L1(alphaf), L1(x), Eshape;
$W/angle_summary_fine.csv). N = 128 was not run.

The setups as rendered (checked in the case dictionaries): orbit2 omega
4 pi about (0.625 . 0.5); the lines with both orbital scales zero and the
level (1.1107 0 1.1107) resp. (1.5708 0 0), the exact body's U0 the same
and axisPoint the sphere's centre; the written U of lineX uniform
(1.57079632679 0 0). OpenFOAM's Courant number (sum of the components) is
0.178 on the diagonal and 0.126 on lineX and the orbit (the displacement
per step 0.089 h along each of x and z on the diagonal, 0.126 h along x).
All eight runs bounded (min(alpha) >= -9.9e-13, excess <= 4.5e-12),
mass deviation <= 6.3e-12 %, no rank fallback. The orbit on de7aea9 against
its 6315783 run: L1(alphaf) apart by at most 2.6e-7 (N = 32) and 1.7e-8
(N = 64) absolute (the contact rule of item 3), the quarter values equal to
the printed digits.

L1(alphaf) (orbit, orbit2 and rotation at t = 1/8 .. 1 of the end time;
orbit2 makes two turns in that time), the largest value and its time, the
time of the maximum of the quarter-turn mean S, the ripple fitted at
4 omega and 8 omega (r = amplitude / mean L1, delta the phase: 0 = maxima at
psi = 0, travel along a mesh axis), L1(x) and Eshape at the end:

| run | N | first | 1/8 | 1/4 | 3/8 | 1/2 | 5/8 | 3/4 | 7/8 | 1 | max (t) | S max at | r4 | delta4 | r8 | delta8 | L1(x) | Eshape |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| orbit | 32 | 5.13e-3 | 7.71e-3 | 1.17e-2 | 1.12e-2 | 1.34e-2 | 1.05e-2 | 1.11e-2 | 7.95e-3 | 9.52e-3 | 1.38e-2 (0.518) | 0.480 | 0.070 | 21.9 | 0.014 | 76.2 | 7.51e-4 | 3.72e-4 |
| orbit | 64 | 2.60e-3 | 4.15e-3 | 6.48e-3 | 6.51e-3 | 7.54e-3 | 6.24e-3 | 6.12e-3 | 4.45e-3 | 4.89e-3 | 7.77e-3 (0.505) | 0.485 | 0.063 | 16.9 | 0.015 | 50.7 | 1.89e-4 | 9.44e-5 |
| orbit2 (2 omega) | 32 | 5.11e-3 | 8.56e-3 | 9.65e-3 | 8.28e-3 | 7.52e-3 | 8.19e-3 | 8.79e-3 | 8.14e-3 | 7.71e-3 | 1.01e-2 (0.255) | 0.243, 0.755 | 0.065 | 32.4 | 0.016 | 51.2 | 6.32e-4 | 3.61e-4 |
| orbit2 (2 omega) | 64 | 2.60e-3 | 4.79e-3 | 5.25e-3 | 4.61e-3 | 4.18e-3 | 4.62e-3 | 4.95e-3 | 4.54e-3 | 4.24e-3 | 5.40e-3 (0.256) | 0.248, 0.751 | 0.073 | 22.2 | 0.024 | 36.7 | 1.70e-4 | 9.23e-5 |
| rotation (base) | 32 | 4.60e-3 | 7.67e-3 | 1.07e-2 | 1.28e-2 | 1.58e-2 | 1.73e-2 | 1.98e-2 | 2.09e-2 | 2.31e-2 | 2.33e-2 (0.995) | 0.875 (window end) | 0.044 | -64.5 | 0.010 | -40.2 | 1.49e-3 | 4.68e-4 |
| rotation (base) | 64 | 2.35e-3 | 4.08e-3 | 5.85e-3 | 7.23e-3 | 9.13e-3 | 1.04e-2 | 1.20e-2 | 1.32e-2 | 1.45e-2 | 1.48e-2 (0.995) | 0.875 (window end) | 0.030 | -71.9 | 0.004 | -53.2 | 4.72e-4 | 1.40e-4 |

The straight translations (t = 0.1125, 0.225, 0.3375, 0.45), the mean over
[0.1, 0.45] and L1(end) / max L1:

| run | N | first | 0.1125 | 0.225 | 0.3375 | 0.45 | mean [0.1, 0.45] | end / max | L1(x) end | Eshape end |
|---|---|---|---|---|---|---|---|---|---|---|
| lineDiag (psi 45) | 32 | 4.31e-3 | 7.98e-3 | 1.15e-2 | 1.47e-2 | 1.75e-2 | 1.263e-2 | 1.00 | 1.20e-3 | 4.13e-4 |
| lineDiag (psi 45) | 64 | 2.19e-3 | 4.35e-3 | 6.35e-3 | 8.39e-3 | 1.03e-2 | 7.17e-3 | 1.00 | 3.41e-4 | 1.14e-4 |
| lineX (psi 0) | 32 | 5.14e-3 | 8.64e-3 | 1.03e-2 | 1.13e-2 | 1.14e-2 | 1.054e-2 | 0.964 | 8.66e-4 | 3.83e-4 |
| lineX (psi 0) | 64 | 2.61e-3 | 5.12e-3 | 5.86e-3 | 6.19e-3 | 6.28e-3 | 5.88e-3 | 0.989 | 2.29e-4 | 9.99e-5 |

Mean L1(alphaf) for t >= T/8 binned by psi ([0, 15), [15, 30), [30, 45]
degrees): orbit 1.140e-2 / 1.065e-2 / 1.002e-2 (N = 32), 6.27e-3 / 5.85e-3 /
5.58e-3 (N = 64); orbit2 8.27e-3 / 7.85e-3 / 7.58e-3, 4.56e-3 / 4.22e-3 /
4.09e-3; rotation 1.635e-2 / 1.590e-2 / 1.564e-2, 9.74e-3 / 9.61e-3 / 9.49e-3.

Verdicts against the pre-registered thresholds:
- F-A HELD: the orbit's ripple at period T/4 has r4 = 0.070 / 0.063 (>= 0.02)
  and delta4 = 21.9 / 16.9 degrees (<= 45): maxima 0.015 T / 0.012 T after
  the quarter points, where the travel is along a mesh axis. Secondary:
  orbit2 r4 0.065 / 0.073, delta4 32.4 / 22.2 degrees (HELD); the rotation
  base's ripple does NOT share the orbit's phase (delta4 -64.5 / -71.9,
  apart by 86 / 89 degrees > 45: FALSIFIED; its maxima lie 0.045 T / 0.050 T
  before the quarter points; r4 0.044 / 0.030).
- F-B HELD (the turned angle): orbit2's S peaks at t = 0.243 / 0.248 and
  0.755 / 0.751, half a turn into each of its two turns, exceeding the
  minimum of S on [0.4, 0.6] by 0.23 / 0.25 and 0.13 / 0.18 of mean L1
  (>= 0.10); the elapsed-time reading is not held (the global maximum of S
  is not in [0.4, 0.6]). Check held: the orbit's S peaks at 0.480 / 0.485.
  The second hump of orbit2 is lower than the first (S excess 0.13 against
  0.23 at N = 32).
- F-C FALSIFIED: the mean of lineX over that of lineDiag is 0.834 / 0.821
  (<= 1.00): the axis-aligned straight translation has the LOWER face error.
  The saturation clause holds for lineX (0.964 / 0.989) and fails for
  lineDiag: its L1(alphaf) grows over the whole run, about linearly (the
  increments per quarter of the run 3.5e-3, 3.2e-3, 2.8e-3 at N = 32),
  slope over the last three quarters 2.77e-2 / 1.80e-2 per unit time
  (rotation base: 1.65e-2 / 1.20e-2 over its turn), h-order of the slope
  0.62.
- F-D HELD: the orbit's Eshape within 1.4 % / 3.7 % of its T/4 value at the
  eighth-turn samples after T/4.

What the numbers show without further interpretation: (i) the rise and fall
of the circular translation's face error recurs once per turn of the
direction of travel (the 2 omega run has two humps at half of each turn,
the omega run one at half its turn), not once per unit of elapsed time;
(ii) on top of it a ripple of period T/4 with 6 to 7 % of the mean, largest
shortly after the travel aligns with a mesh axis, in both orbits; (iii) a
straight translation along the x-z diagonal accumulates face error without
any turning, at a slope above the rotation's, while the one along x
saturates. Surprising against the earlier picture (the straight
translation along x saturates, so translation was taken as
non-accumulating): the diagonal direction accumulates. Not examined: the
diagonal line at another Co (its per-axis displacement per step is 0.089 h
against 0.126 h on lineX), other diagonal directions, N = 128.

## The diagonal line's confound (2026-11-02)

Tomislav's approval via the coordinating session. At the same speed and dt
the x-z diagonal moved 0.0889 h per step along each axis and the x-line
0.1257 h. The fixed step alone crosses the two variables (dt = Co h / 2.5;
speed pi/2, window END_TIME 0.45, fork de7aea9 / build 55845642, serial):
LineDiagM = the x-z diagonal at Co 0.2 sqrt(2) (0.1257 h per axis);
LineXs = the x-line at Co 0.2/sqrt(2) (0.0889 h per axis); Line3D = the
space diagonal (pi/2)(1 1 1)/sqrt(3) from (0.25 0.25 0.25) at Co 0.2
(0.0726 h per axis) and Co 0.2 sqrt(3) (0.1257 h), N = 32 only. leia
fc97d8bd: configs, davof_line_confound.py and the predictions P-L1..P-L4
(in sphereRotation3D_hexFineLineDiagM.yaml), committed before the runs.
Jobs 55854681 (LineDiagM, 8.3 min), 55854682 (LineXs, 12.1 min), 55854683
(Line3D, 4.6 min). The rendered dictionaries were checked (deltaT
0.003535533906 / 0.001767766953 / 0.0025 / 0.004330127019 at N = 32, the
arm velocity as setFlow's level and as U0, axisPoint the sphere's centre).

Measure: the late slope s, the least-squares slope of L1(alphaf) against t
over [0.25 T, T] (as reported for lineX and lineDiag in item 4), and its
log-position x = ln(s/s_X)/ln(s_D/s_X) between the existing x-line and
diagonal of the same N: D-like x >= 0.75, X-like x <= 0.25. "per-axis" is
the displacement per step along each moving mesh axis in cells, "|U|dt"
along the travel; mean L1 over t in [0.1, 0.45]; bounds and mass over all
steps, the mass deviation at the end; no rank fallback in any run:

| run | N | per-axis (h) | U dt (h) | steps | late slope | x | class | mean L1 | L1 first | L1 max | L1 end | L1(x) end | Eshape end | min alpha | excess | mass dev % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineX | 32 | 0.1257 | 0.1257 | 180 | 8.059e-3 | 0 | ref | 1.054e-2 | 5.14e-3 | 1.177e-2 | 1.135e-2 | 8.66e-4 | 3.83e-4 | -3.4e-14 | 9.9e-13 | -2.0e-13 |
| lineX | 64 | 0.1257 | 0.1257 | 360 | 2.988e-3 | 0 | ref | 5.884e-3 | 2.61e-3 | 6.346e-3 | 6.278e-3 | 2.29e-4 | 9.99e-5 | -5.4e-14 | 2.5e-12 | 3.9e-13 |
| lineDiag | 32 | 0.0889 | 0.1257 | 180 | 2.773e-2 | 1 | ref | 1.263e-2 | 4.31e-3 | 1.752e-2 | 1.752e-2 | 1.20e-3 | 4.13e-4 | -2.4e-14 | 1.2e-12 | 6.7e-13 |
| lineDiag | 64 | 0.0889 | 0.1257 | 360 | 1.804e-2 | 1 | ref | 7.170e-3 | 2.19e-3 | 1.028e-2 | 1.028e-2 | 3.41e-4 | 1.14e-4 | -6.9e-14 | 4.5e-12 | 3.2e-12 |
| LineDiagM | 32 | 0.1257 | 0.1777 | 127 | 2.405e-2 | 0.885 | D-like | 1.123e-2 | 3.64e-3 | 1.542e-2 | 1.542e-2 | 1.07e-3 | 3.94e-4 | -2.3e-13 | 7.2e-13 | 2.4e-13 |
| LineDiagM | 64 | 0.1257 | 0.1777 | 255 | 1.550e-2 | 0.916 | D-like | 6.325e-3 | 1.95e-3 | 9.049e-3 | 9.005e-3 | 3.02e-4 | 1.07e-4 | -4.2e-13 | 4.4e-12 | 2.5e-12 |
| LineXs | 32 | 0.0889 | 0.0889 | 255 | 8.914e-3 | 0.082 | X-like | 1.138e-2 | 5.43e-3 | 1.281e-2 | 1.237e-2 | 9.10e-4 | 3.91e-4 | -1.7e-14 | 1.2e-12 | 6.7e-14 |
| LineXs | 64 | 0.0889 | 0.0889 | 509 | 3.192e-3 | 0.037 | X-like | 6.346e-3 | 2.72e-3 | 6.825e-3 | 6.776e-3 | 2.45e-4 | 1.02e-4 | -4.8e-14 | 3.6e-12 | 5.0e-13 |
| Line3D Co 0.2 | 32 | 0.0726 | 0.1257 | 180 | 2.356e-2 | 0.868 | D-like | 1.496e-2 | 3.60e-3 | 1.981e-2 (t 0.348) | 1.763e-2 | 1.85e-3 | 4.31e-4 | -2.5e-14 | 1.3e-12 | 1.6e-12 |
| Line3D Co 0.346 | 32 | 0.1257 | 0.2177 | 104 | 1.067e-2 | 0.227 | X-like | 1.924e-2 | 2.21e-3 | 2.577e-2 (t 0.368) | 1.854e-2 | 2.22e-3 | 4.56e-4 | -3.2e-14 | 1.0e-12 | 1.0e-12 |

L1(alphaf) at t = 0.05, 0.10, ..., 0.45: LineDiagM N = 32 5.72e-3 6.86e-3
7.97e-3 9.54e-3 1.09e-2 1.19e-2 1.30e-2 1.42e-2 1.54e-2 (monotone, max at
the end, as lineDiag); LineXs N = 32 8.01e-3 9.16e-3 1.02e-2 1.09e-2 1.15e-2
1.17e-2 1.22e-2 1.24e-2 1.24e-2; Line3D Co 0.2 8.48e-3 1.04e-2 1.02e-2
1.33e-2 1.70e-2 1.50e-2 1.91e-2 1.76e-2 1.76e-2; Line3D Co 0.346 1.42e-2
1.27e-2 1.96e-2 2.08e-2 2.20e-2 2.22e-2 2.29e-2 2.05e-2 1.85e-2.

Verdicts:
- P-L1 HELD (the oblique travel): LineDiagM is D-like at both N (x 0.885 /
  0.916; slope 87 % / 86 % of lineDiag's) and LineXs X-like at both N (x
  0.082 / 0.037; slope 111 % / 107 % of lineX's). Matching the per-axis
  displacement to the x-line's leaves the diagonal accumulating, and
  reducing the x-line's to the diagonal's leaves it levelling off.
- P-L2 NOT HELD.
- P-L3 (secondary, N = 32): by the threshold the space diagonal is D-like
  at Co 0.2 (x 0.868) and X-like at Co 0.2 sqrt(3) (x 0.227), the pattern
  the per-axis reading predicts, not the oblique one. Both 3D histories
  are non-monotone (maxima 1.98e-2 at t 0.348 and 2.58e-2 at t 0.368,
  then falling), unlike every planar line; the Co 0.346 run's slope is
  small because it falls after t 0.37, while its mean L1 (1.92e-2) is the
  largest of all runs. A late slope does not separate accumulation from
  levelling off on these two histories, so the 3D pair does not decide
  between the readings.
- P-L4 HELD: min(alpha) >= -4.2e-13, excess <= 4.4e-12, |mass deviation|
  <= 2.5e-12 % in all six new runs.

Converter note: with an irrational fixed dt the time of a step can print
differently in two postProcessing files, and tpf_advection_metrics.py then
writes that step as two rows, each with the other's columns empty (8 / 10 /
10 / 24 / 0 / 6 such pairs in the six runs). davof_line_confound.py reads
every column from its non-empty entries; the converter is not changed.
(Superseded: the converter fixed in leia 24ebf8df, next section.)

## The converter's join on the time (2026-11-02)

Cause: the sources print a step's time with different precisions
(areaFractionError.dat 0.03181980515, the log's `t = ` 0.031819805154); the
join key, the time rounded to 9 significant digits, made 0.0318198051 of the
first and 0.0318198052 of the second, and the outer join wrote two rows.
Fix (leia 24ebf8df, tpf_advection_metrics.py): the times are matched
numerically, one step when within tol = 1e-3 dt_min (dt_min the smallest gap
between consecutive times within one source). The printed times of one step
differ by about 5e-10 |t| at most (10 or more significant digits), so tol
lies three orders of magnitude above that and below dt for |t| < 1e6 dt. Not
the step index: the sources do not share a row set (reconstructionError at
the end only unless every step, a t = 0 row in some .dat files, restarts).
TIME is the 9-digit key of the log's time where the log has the step, else
of the first source with it; unchanged for every unsplit step.

Verification (job 55855815, 2.1 min): all 272 cases under leia/studies with
advectorVoF.csv and log.advectorVoF re-converted by the previous (2aa541ee)
and the new converter on the same inputs: 267 byte-identical; the 5 others
are exactly the split confound cases (LineDiagM 00000 / 00001, LineXs 00000
/ 00001, Line3D 00001: 135 / 265 / 265 / 533 / 110 rows -> 127 / 255 / 255 /
509 / 104, one per step), in each the sequence of non-empty values of every
column identical, the constant columns identical, and no row left with the
area-fraction columns filled and the log's empty or the reverse. The six
confound cases were re-converted in place (the old files kept as
advectorVoF.preJoinFix.csv; tpfCommit and tpfVersion carried over from them,
de7aea9). The workaround of davof_line_confound.py was then removed (leia
0940279b): the plain row-wise reading on the re-converted files prints the
confound table identical, byte for byte, to the workaround on the split
files. The studies of the next two sections were converted by the new
converter as they ran.

## The exact plane on the straight lines (2026-11-02)

Does the oblique build-up need the volume-anchored plane? The x-z diagonal
and the x-line at Co 0.2 with planeSource exact (the tangent plane of the
exact sphere per interface cell). Fork f974845 (de7aea9 plus the opt-in
faceErrorCsv, off here), build 55855814; leia 4a4311b6 (configs and P-E1..P-E3
with the exact plane's expected effects, committed before the runs); jobs
55855816 (diagonal, 12.9 min), 55855817 (x-line, 10.1 min). Late slope s and
x as in the confound section; "error sum" = L1 nFaces (the sum of
|alphaf - alphaf_exact| over the counted faces), its late slope over the
same window; interface cells = the davof count of the transported state at
the first step and the end; no rank fallback in any run:

| run | N | late slope s | mean L1 | L1 first | L1 end | L1(x) end | Eshape end | nFaces first -> end | interface cells | error-sum late slope | min alpha | excess | mass dev % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineX, anchor | 32 | 8.059e-3 | 1.054e-2 | 5.14e-3 | 1.135e-2 | 8.66e-4 | 3.83e-4 | 1668 -> 1765 | 416 -> 444 | 14.13 | -3.4e-14 | 9.9e-13 | -2.0e-13 |
| lineX, anchor | 64 | 2.988e-3 | 5.884e-3 | 2.61e-3 | 6.278e-3 | 2.29e-4 | 9.99e-5 | 6990 -> 6954 | 1760 -> 1756 | 21.02 | -5.4e-14 | 2.5e-12 | 3.9e-13 |
| lineDiag, anchor | 32 | 2.773e-2 | 1.263e-2 | 4.31e-3 | 1.752e-2 | 1.19e-3 | 4.13e-4 | 1656 -> 1668 | 416 -> 414 | 48.71 | -2.4e-14 | 1.2e-12 | 6.7e-13 |
| lineDiag, anchor | 64 | 1.804e-2 | 7.170e-3 | 2.19e-3 | 1.028e-2 | 3.41e-4 | 1.14e-4 | 6984 -> 7056 | 1760 -> 1732 | 126.3 | -6.9e-14 | 4.5e-12 | 3.2e-12 |
| lineX, exact plane | 32 | 4.723e-3 | 1.218e-2 | 1.25e-2 | 1.270e-2 | 9.15e-4 | 4.90e-4 | 1668 -> 1756 | 416 -> 432 | -1.62 | -2.3e-17 | 3.3e-13 | 1.9e-13 |
| lineX, exact plane | 64 | 3.749e-3 | 6.138e-3 | 6.33e-3 | 6.440e-3 | 2.96e-4 | 1.46e-4 | 7046 -> 6942 | 1760 -> 1748 | -5.22 | -4.0e-17 | 2.4e-12 | 1.1e-13 |
| lineDiag, exact plane | 32 | -1.163e-2 | 7.755e-3 | 1.18e-2 | 6.295e-3 (max 1.28e-2) | 1.60e-1 | 9.52e-4 | 1665 -> 4028 | 416 -> 1491 | 19.74 | -4.1e-13 | 6.0e-8 | 5.4e-13 |
| lineDiag, exact plane | 64 | -2.743e-3 | 4.340e-3 | 6.06e-3 | 3.482e-3 (max 6.42e-3) | 1.21e-1 | 3.67e-4 | 7008 -> 29256 | 1760 -> 3310 | 89.03 | -2.6e-4 | 2.5e-7 | 2.8e-12 |

L1(alphaf) of the exact-plane diagonal at t = 0.1125 / 0.225 / 0.3375:
9.35e-3 / 8.88e-3 / 6.18e-3 (N = 32), 4.42e-3 / 4.13e-3 / 3.78e-3 (N = 64);
of the exact-plane x-line 1.20e-2 / 1.13e-2 / 1.29e-2, 5.15e-3 / 5.78e-3 /
6.43e-3.

Verdicts against the pre-registered thresholds:
- P-E1 HELD by its measure: D_E = s(diagonal) - s(x-line) = -1.635e-2 /
  -6.49e-3 (<= 4.92e-3 / 3.76e-3). P-E2 NOT HELD.
- The bias named in the header decides this verdict. The exact-plane
  diagonal smears: the transported state's interface cells grow 416 ->
  1491 / 1760 -> 3310, the counted faces 2.4 / 4.2 times, L1(x) reaches
  0.16 / 0.12 (anchored 1.2e-3 / 3.4e-4), min(alpha) -2.6e-4 at N = 64; the
  mean L1 falls because the added faces enter the denominator, while the
  error SUM rises with a late slope of 19.7 / 89.0 against -1.6 / -5.2 of
  the exact-plane x-line (anchored: 48.7 / 126.3 against 14.1 / 21.0). The
  exact-plane x-line does not smear (432 / 1748 interface cells at the end,
  L1(x) 9.2e-4 / 3.0e-4, bounded to 1e-16, level L1(alphaf)). So with the
  exact plane the diagonal still differs from the x-line, now as a growing
  error sum on a smearing field, and the registered measure (the mean L1)
  does not show it; the run does not show that the build-up needs the
  anchor.
- P-E3: mass deviation <= 2.8e-12 % (HELD). Bounds: the x-line stays within
  [-4e-17, 1 + 2.4e-12]; the diagonal reaches -2.6e-4 and 1 + 2.5e-7 at
  N = 64 (the header expected both lines outside; only the diagonal is).

Surprising: an exact tangent plane per cell, independent of alpha, keeps the
transported alpha on the sphere along x (L1(x) at the anchor's level) and
lets it drift and smear along the x-z diagonal at the same speed and step.

## Where the anchored diagonal's error builds up (2026-11-02)

The anchored lines rerun with the per-face output (fork f974845:
areaFractionError faceErrorCsv, the counted faces of every step with centre,
S_f, alphaf, alphaf_exact; leia 4a4311b6: configs, davof_face_bins.py and
B-0..B-3 committed before the runs); jobs 55855818 (diagonal, 10.8 min) and
55855819 (x-line, 12.3 min). Per step e = |alphaf - alphaf_exact| per
counted face; the contribution of a bin C_b = sum_bin e / nFaces (the C_b of
a binning add up to L1), its late slope over [0.25 T, T] and its share of
the late slope of L1 (the shares add up to 1). Bins: the face direction;
the side, theta = angle of the exact outward normal to the travel: front
[0, 60), tangent [60, 120), rear [120, 180] degrees; theta in 30-degree
bins; the angle of the normal to the nearest mesh axis. "faces" = the bin's
share of the counted faces at the end.

Checks (B-0 HELD): areaFractionError.dat of all four runs byte-identical to
the original lineDiag / lineX runs (fork de7aea9); the row sums of the CSVs
reproduce the L1 column to 4.8e-13 relative at every step (180 / 360 steps).
By symmetry the x- and z-normal faces of the diagonal and the y- and
z-normal faces of the x-line carry identical numbers.

Late slope per bin (share of the late slope of L1 in brackets) and the
bin's share of the faces, N = 32 / N = 64:

| bin | lineDiag N 32 | lineDiag N 64 | lineX N 32 | lineX N 64 | faces (diag / x) |
|---|---|---|---|---|---|
| total L1 | 2.773e-2 | 1.804e-2 | 8.059e-3 | 2.988e-3 | |
| x-normal | 9.92e-3 (0.358) | 6.32e-3 (0.350) | 3.25e-3 (0.403) | 1.17e-3 (0.391) | 0.33 / 0.33 |
| y-normal | 7.90e-3 (0.285) | 5.39e-3 (0.299) | 2.41e-3 (0.299) | 9.10e-4 (0.304) | 0.33 / 0.33 |
| z-normal | 9.92e-3 (0.358) | 6.32e-3 (0.350) | 2.41e-3 (0.299) | 9.10e-4 (0.304) | 0.33 / 0.33 |
| front (0-60) | 1.152e-2 (0.416) | 6.77e-3 (0.375) | 1.14e-3 (0.141) | -7.3e-6 (-0.002) | 0.25 / 0.25-0.27 |
| tangent (60-120) | 6.37e-3 (0.230) | 5.37e-3 (0.297) | 5.79e-3 (0.718) | 3.01e-3 (1.006) | 0.50-0.51 / 0.46-0.50 |
| rear (120-180) | 9.84e-3 (0.355) | 5.90e-3 (0.327) | 1.13e-3 (0.140) | -1.0e-5 (-0.003) | 0.24-0.25 / 0.25-0.27 |
| theta 0-30 | 3.44e-3 (0.124) | 2.18e-3 (0.121) | -1.9e-5 (-0.002) | 8.5e-7 (0.000) | 0.07 / 0.06-0.07 |
| theta 30-60 | 8.09e-3 (0.292) | 4.59e-3 (0.254) | 1.16e-3 (0.143) | -8.1e-6 (-0.003) | 0.17-0.18 / 0.19-0.20 |
| theta 60-90 | 3.24e-3 (0.117) | 2.82e-3 (0.156) | 2.95e-3 (0.366) | 1.45e-3 (0.485) | 0.26-0.27 / 0.21-0.26 |
| theta 90-120 | 3.13e-3 (0.113) | 2.54e-3 (0.141) | 2.84e-3 (0.352) | 1.56e-3 (0.521) | 0.24-0.25 / 0.24-0.25 |
| theta 120-150 | 6.80e-3 (0.245) | 3.89e-3 (0.216) | 1.08e-3 (0.134) | -1.1e-5 (-0.004) | 0.17-0.18 / 0.19-0.20 |
| theta 150-180 | 3.05e-3 (0.110) | 2.01e-3 (0.111) | 4.8e-5 (0.006) | 9.4e-7 (0.000) | 0.07 / 0.06 |
| axis 0-15 | 1.49e-3 (0.054) | 7.85e-4 (0.044) | 2.99e-4 (0.037) | 2.61e-4 (0.087) | 0.07 / 0.07 |
| axis 15-30 | 7.05e-3 (0.254) | 5.30e-3 (0.294) | 2.52e-3 (0.313) | 1.07e-3 (0.359) | 0.26-0.29 / 0.27-0.29 |
| axis 30-45 | 1.596e-2 (0.575) | 9.60e-3 (0.532) | 3.53e-3 (0.438) | 1.27e-3 (0.426) | 0.51-0.54 / 0.50-0.51 |
| axis 45-55 | 3.23e-3 (0.117) | 2.35e-3 (0.130) | 1.71e-3 (0.212) | 3.83e-4 (0.128) | 0.13 / 0.13-0.16 |

Contribution C_b at the first step -> the end, N = 32 (N = 64): diagonal
front 1.15e-3 -> 6.45e-3 (5.84e-4 -> 3.57e-3), tangent 2.11e-3 -> 6.14e-3
(1.05e-3 -> 3.90e-3), rear 1.05e-3 -> 4.94e-3 (5.55e-4 -> 2.82e-3); x-line
front 1.44e-3 -> 2.59e-3 (7.01e-4 -> 9.50e-4), tangent 2.27e-3 -> 6.71e-3
(1.24e-3 -> 4.43e-3), rear 1.43e-3 -> 2.05e-3 (6.64e-4 -> 9.00e-4). Mean
error per face at the end: diagonal front 2.63e-2 / 1.42e-2, rear 2.03e-2 /
1.14e-2, tangent 1.20e-2 / 7.77e-3; x-line front 9.48e-3 / 3.83e-3, rear
7.66e-3 / 3.60e-3, tangent 1.46e-2 / 8.82e-3. Within the diagonal's front
and rear, the x-, y- and z-normal faces carry 0.13 to 0.14 of the slope each
at N = 32 (y/front 0.142, y/rear 0.091) and 0.09 to 0.13 at N = 64.

Verdicts:
- B-1 NOT HELD: the y-normal faces carry 0.285 / 0.299 of the diagonal's
  late slope (< 0.50), their share of the faces 0.33.
- B-2 NOT HELD: the x- and z-normal faces carry 0.716 / 0.701 (< 0.75).
  Hence MIXED: the build-up is not tied to a face direction; each direction
  carries about its share of the faces in both lines.
- B-3 FALSIFIED: the tangent band carries 0.230 / 0.297 of the diagonal's
  late slope (< 0.50) with 0.51 / 0.50 of the faces.

Where the build-up concentrates: the diagonal's excess late slope over the
x-line (1.967e-2 / 1.505e-2) lies on the front and rear sides of the
sphere: there the diagonal's slope is 2.136e-2 / 1.267e-2 against 2.27e-3 /
-1.7e-5 on the x-line, i.e. 1.909e-2 / 1.269e-2 = 97 % / 84 % of the excess;
the tangent band's slope is close in both lines (6.37e-3 / 5.37e-3 against
5.79e-3 / 3.01e-3; 3 % / 16 % of the excess). Within front and rear the
30-60 and 120-150 degree bins carry the most (0.292 + 0.245 / 0.254 + 0.216
of the diagonal's slope with 0.35 of the faces). On the x-line the front
and rear sides stop growing at N = 64 (slopes -7e-6 and -1e-5) and its whole
late slope sits in the tangent band. Against the mesh axes the diagonal's
slope follows the face counts of the angle bins (0.575 / 0.532 in 30-45
degrees with 0.54 / 0.51 of the faces), the x-line's likewise.

## Setup check of the straight lines (2026-11-02, no new runs)

Asked because the exact-plane diagonal smears where the exact-plane x-line
does not. Read from the fork's code (f974845) and the rendered cases on the
cluster; no defect found, no run repeated.

Code: every object that defines the exact body reads origin, U0, omega and
axisPoint from its own dictionary and places it with davof::orbitOrigin,
which turns x0 about axisPoint and then adds U0 t (davofOrbit.H): the
exact normal model (exactNormal.C, U0 required), the exactSphere surface and
the planeSource exact plane (isoAdvectionCurved.C, the exactNormal subdict
of the advection's model dictionary, at tState = t - dt, i.e. t^n),
areaFractionError and volumeFractionError. With axisPoint = x0 (the lines'
setting) the turn is the identity and the centre is x0 + U0 t in every
direction. Inference of U0 from the field is left only without a U0 key and
with omega zero (areaFractionError, volumeFractionError); every rendered
dictionary of these cases gives U0, so nothing is inferred.

Rendered dictionaries (LineDiagExactPlane, LineXExactPlane, LineDiag, LineX,
N = 32 and 64; Line3D both Co): in fvSolution's exactNormal, in
areaFractionError and in volumeFractionError, origin, radius 0.15, U0,
omega (0 -2 pi 0) and axisPoint are identical, U0 equal to the arm's
velocity in all three components ((1.1107207345395915 0 1.1107207345395915),
(1.5707963267948966 0 0), (0.9068996821171089 x3)) and axisPoint equal to
origin; setFlow's two orbital scales are zero and its level is the same
vector; the written U at the first write is uniform and equal to it to the
12 written digits; the log's exact-normal line repeats the values. The
areaFractionError centre written by the Faces runs at t = 0.45 is (0.749824
0.5 0.749824) and (0.956858 0.5 0.5), x0 + U0 t.

The alpha centroid of the written fields (uniform hex mesh, 2N x N x N cells
on [0,2] x [0,1]^2) against x0 + U0 t, in cells h, at t = 0.15 / 0.30 /
0.45: LineDiagExactPlane 0.054 / 0.041 / 0.0015 (N = 32), 0.024 / 0.008 /
0.036 (N = 64); LineXExactPlane 0.041 / 0.046 / 0.045, 0.018 / 0.020 /
0.021; LineDiag 0.027 / 0.048 / 0.064, 0.014 / 0.026 / 0.036; LineX 0.013 /
0.017 / 0.019, 0.004 / 0.006 / 0.006; Line3D Co 0.2 0.027 / 0.051 / 0.070,
Co 0.346 0.006 / 0.016 / 0.026 (the y component moving with x and z). The
volume is conserved to the printed digits in every run. The exact body and
the transported body move together in all of them.

What the smearing is (cells by alpha at the end, N = 32 unless stated):

| run | [1e-12, 1e-6) | [1e-6, 1e-3) | [1e-3, 0.01) | [0.01, 0.99) | [0.99, 1 - 1e-6) | volume below 1e-3 (cells) |
|---|---|---|---|---|---|---|
| lineX, anchor | 260 | 57 | 40 | 336 | 70 | 7.5e-4 |
| lineDiag, anchor | 5076 | 33 | 30 | 356 | 55 | 3.4e-4 |
| lineX, exact plane | 0 | 8 | 16 | 320 | 180 | 3.9e-3 |
| lineDiag, exact plane | 5538 | 3480 | 1040 | 464 | 188 | 0.46 |
| lineDiag, exact plane, N 64 | 30662 | 20237 | 6426 | 1666 | 2436 | 3.29 |
| Line3D Co 0.2 | 20800 | 2281 | 45 | 361 | 32 | 3.4e-2 |
| Line3D Co 0.346 | 19210 | 2029 | 18 | 361 | 41 | 4.7e-2 |

(of 463.25 cell volumes at N = 32, 3706 at N = 64). The exact-plane
diagonal's growth of the interface cells is a film of cells with alpha below
1e-2 around the body, holding 0.1 % of the volume at N = 32; the band
[0.01, 0.99) grows only 398 -> 464. The oblique anchored runs leave a film
as well, at smaller values: 5076 cells below 1e-6 on the x-z diagonal
against 260 on the x-line, and 23 000 cells below 1e-3 on the space
diagonal. Not a setup defect; its cause is not examined here.

## The consistency of the swept flux under oblique travel (2026-11-02)

Tomislav's approval via the coordinating session. The problem: under travel
oblique to the mesh a region of a cell lies in the swept volumes of two
faces (a corner of C in those of C|N and N|NE under north-east travel), and
the cell that receives it balances only if both faces give that region the
same liquid. Two rules of the swept flux can differ there:
- (1) mixed faces: only the faces of candidate cells (the interface cells
  and their face neighbours) are computed; every other face keeps the
  algebraic upwind flux alpha_donor phi_f dt, also where its swept volume
  reaches a cell in another state;
- (2) wisps: a cell with 1e-8 < alpha < 1 - 1e-8 and no surface is wet or
  dry in the pieces by alpha >= 1/2, while the algebraic faces release
  alpha phi_f dt.

The hypothesis under test: these leave residues in the cells beside the
body (the low-alpha film) and drive the diagonal's front/rear build-up.

Fork davof-swept c635611 (both keys opt-in, default false):
- sweptConsistencyLog (serial) probes every algebraic face whose vertices
  touch a cell with a surface or another alpha. It cuts the face's swept
  volume into the ring pieces and compares the consistent rule (a surface
  clips, any other cell gives alpha_c V_piece) with the algebraic flux. It
  also sums the wisp pieces' 0/1-against-alpha mismatch of the computed
  faces, and writes the per-cell residue R (the used change of a cell's
  liquid minus the consistent rule's) per step.
- sweptConsistent extends the candidates to the vertex ring of every active
  cell: a surface, an interface cell, or a wisp. The vertex marks are synced
  over the processor patches. Such faces are computed with or without a
  surface, and a piece without a surface carries alpha_c V_piece (alpha_c
  times its pre-image area in the face update).
- The stencil: a face's swept volume reaches only cells sharing a vertex of
  the face while the displacement per step stays below h per axis (0.089 h
  and 0.126 h here). Every such cell shares that vertex with the donor.

leia:
- 126761fc: tokens, configs, davof_consistency.py and the predictions,
  committed before the runs.
- 2f240340: the converter fix below.

Build 55858735. Serial line studies: jobs 55859067 (DiagConsLog) and
55859068 (XConsLog), the A runs; 55859069 (DiagCons) and 55859070 (XCons),
the C runs. Parallel: 55859071-75.

Verification (jobs 55858715 before / 55858737 after, 20-step copies):
- With the keys off, the x-z diagonal, the x-line and the rotation base
  (N = 32) reproduce f974845 byte for byte: fields at steps 10 and 20, every
  .dat, and the log apart from timings.
- With the log on, the flux, the fields and the .dat files are unchanged.
- The full ConsLog runs give areaFractionError.dat and volumeFractionError.dat
  byte-identical to the original lineDiag / lineX runs in all eight cases.

Disclosed in the headers: the 20-step numbers were seen before the
predictions were written. The agent's own expectation (C-1 unlikely to hold)
is stated there separately.

Parallel (sphereRotation3D_hexSweptConsParN32np{1,2,4,8,16}: rotation, rows
exact and ring, 20 steps; davof_decomposition_compare.py, $W/cmp-cons):
- Step 1: max|dalpha| 1.3e-14 to 8.8e-14 and max|dalphaf| 4.1e-14 to
  8.5e-14.
- Steps 1 to 20: max|dalpha| at most 2.6e-12 (exact) and 4.8e-13 (ring),
  below the default family's 6.0e-11 / 2.5e-9.
- nMissing 0 and nDup 0 at every step.
- At step 1, 368 candidates added and 1092 faces without a surface for
  every np.
- No rank fallback.
- P-C1, P-C2 and P-C3 HELD.

A: the diagnostic on the default flux (sums over the run in h^3; per-step
means; residue cells: |R| > 1e-15 V; film: alpha in [1e-12, 1e-3) at the
end):

| | lineDiag N 32 | lineDiag N 64 | lineX N 32 | lineX N 64 |
|---|---|---|---|---|
| algebraic faces probed per step | 28010 | 218100 | 1672 | 11650 |
| reaching another state per step (surface / wisp / pure) | 1144 (0 / 981 / 163) | 10020 (0 / 9103 / 918) | 0 | 0 |
| algebraic mismatch: sum of abs / largest (h^3) | 6.2e-5 / 3.8e-8 | 1.2e-3 / 4.3e-8 | 1.1e-10 / 4.3e-14 | 2.4e-8 / 1.9e-13 |
| wisp pieces per step | 500 | 2872 | 61 | 487 |
| wisp mismatch: sum / sum of abs / largest (h^3) | 3.9e-4 / 2.1e-3 / 1.3e-6 | 1.1e-2 / 2.9e-2 / 1.5e-6 | 1.0e-3 / 1.9e-3 / 2.3e-6 | 2.7e-2 / 3.7e-2 / 5.7e-6 |
| residue: sum of abs (h^3) / largest abs(R)/V | 3.8e-3 / 2.6e-6 | 5.2e-2 / 2.9e-6 | 2.9e-3 / 2.3e-6 | 5.0e-2 / 5.7e-6 |
| film cells at the end that got a residue | 99.8 % | 99.5 % | 95.9 % | 99.3 % |
| residue cells that end in the film | 66 % | 56 % | 15 % | 23 % |
| late abs(R): front / tangent / rear | 0.255 / 0.460 / 0.285 | 0.286 / 0.529 / 0.185 | 0.229 / 0.541 / 0.230 | 0.238 / 0.639 / 0.122 |
| abs(R) by distance: inside / 0-1 h / 1-2 h / >2 h | 0.58 / 0.14 / 0.23 / 0.05 | 0.68 / 0.09 / 0.18 / 0.05 | 0.72 / 0.10 / 0.15 / 0.03 | 0.81 / 0.05 / 0.12 / 0.02 |

Verdicts A:
- A-0 HELD.
- A-1 HELD: on the x-line no algebraic face reaches another state, and its
  summed mismatch is 2e-6 / 2e-5 of the diagonal's.
- A-2 HELD by its measure, but the measure does not separate. Almost every
  film cell receives a residue above 1e-15 V on the x-line too, and in the
  consistent runs, whose residue is round-off (below).
- A-3 NOT HELD at N = 32 (front + rear 0.540) and FALSIFIED at N = 64
  (0.471). The residue sits mostly in the tangent band and inside the
  sphere.

Of the two mechanisms, the wisp rule moves 34 / 25 times more than the mixed
faces, and it moves as much on the x-line as on the diagonal. No algebraic
face reaches an interface cell: the mixed faces reach wisps and pure cells
only.

C: the consistent variant (sweptConsistent, with the log and faceErrorCsv).
Late slope over [0.25 T, T] and x as in the confound section; front + rear
from davof_face_bins.py; bands at the end; the default values in brackets:

| run | N | late slope | x | L1 end | max abs(dL1) over the steps | front + rear slope | film [1e-12, 1e-3) | [1e-12, 1e-6) / [1e-6, 1e-3) / [0.99, 1 - 1e-6) | mass dev % end | min alpha / excess |
|---|---|---|---|---|---|---|---|---|---|---|
| lineDiag cons | 32 | 2.7731e-2 (2.7733e-2) | 1.000 | 1.7517e-2 (1.7518e-2) | 7.0e-7 | 2.136e-2 (2.136e-2) | 7239 (5109) | 7232 / 7 / 24 (5076 / 33 / 55) | 5.6e-7 (6.7e-13) | -9.7e-15 / 1.0e-13 |
| lineDiag cons | 64 | 1.8032e-2 (1.8037e-2) | 1.000 | 1.0279e-2 (1.0280e-2) | 1.7e-6 | 1.267e-2 (1.267e-2) | 39223 (30711) | 39187 / 36 / 167 (30404 / 307 / 638) | 8.4e-8 (3.2e-12) | -1.8e-13 / 6.9e-13 |
| lineX cons | 32 | 8.0591e-3 (8.0589e-3) | 0.000 | 1.1354e-2 (1.1354e-2) | 1.0e-6 | 2.27e-3 (2.27e-3) | 1925 (317) | 1902 / 23 / 45 (260 / 57 / 70) | 3.0e-13 (-2.0e-13) | -1.5e-14 / 1.4e-13 |
| lineX cons | 64 | 2.9861e-3 (2.9884e-3) | -0.000 | 6.2768e-3 (6.2776e-3) | 4.7e-6 | -2.0e-5 (-1.7e-5) | 13257 (4182) | 13147 / 110 / 233 (3817 / 365 / 753) | 2.1e-12 (3.9e-13) | -6.6e-14 / 8.6e-13 |

Per step in the consistent runs:
- Candidates added: 1883 / 9874 (diagonal), 1279 / 7039 (x-line).
- Faces without a surface: 5970 / 30220 and 4023 / 21470.
- Fallbacks of those faces: 0. Algebraic faces reaching another state: 0.
  Wisp mismatch: 0.
- The residue left per step is at most 5e-11 V: algebraic faces between
  cells below 1e-8 of different alpha.
- Interface cells at the end: 416 / 1762 (diagonal) and 433 / 1724 (x-line),
  as in the default runs.
- Eshape and L1(x) at the end are within 0.03 % of the default.

Verdicts C:
- C-0: the consistency clauses HELD. The mass clause NOT HELD on the
  diagonal: the whole deviation is liquid leaving through the walls (summed
  boundary outflow 5.626e-9 / 8.42e-10 of the mass, from step 146 / 332).
  - The film reaches the outflow walls: 370 / 1071 wall cells with alpha
    above 1e-12, against 173 / 585 in the default.
  - The swept flux of a computed wall face carries alpha_c V_piece, where
    the algebraic face takes the wall's fixedValue 0. In the default runs no
    computed face reaches a wall.
  - By the registered rule the verdict on the hypothesis is therefore
    INCONCLUSIVE.
- C-1 FALSIFIED by its threshold: the consistent diagonal stays D-like
  (x = 1.000 at both N). Removing both inconsistencies changes its
  L1(alphaf) history by at most 1.7e-6 absolute and its late slope by
  0.006 % / 0.02 %.
- C-2 FALSIFIED: the film grows, by 1.42 / 1.28 times on the diagonal and
  6.1 / 3.2 times on the x-line, while the near-full band shrinks.
- C-3 FALSIFIED: the front + rear slope is unchanged to 4 digits.
- C-4 HELD.
- R-1: the rotation base was not run (conditional on C-1).

What the runs show: the mixed-face and wisp mismatches are real but small
(1e-8 to 1e-6 of a cell per face and step). Removing them leaves the face
error history of both lines unchanged to 1e-6, so they do not drive the
oblique build-up. The film is not their residue either. The 0/1 wisp rule
holds the liquid of the sub-1/2 wisps back, and the alpha_c V_piece rule
transports it upwind like the algebraic faces do, so the consistent variant
spreads a larger film. On the diagonal that film reaches the outflow walls.

Surprising:
- The 0/1 wisp rule, the inconsistency, keeps the film smaller.
- The swept flux does not apply the alpha boundary condition on the wall
  faces it computes.

Not examined: the cause of the film itself (a candidate: the bounding of
limitFluxes); the consistent variant without the log (cost). With the log
on, N = 64 took 2048 s (default flux) and 2322 s (consistent) on the
diagonal, and 355 s and 773 s on the x-line.

## The converter's time directory (2026-11-02)

tpf_advection_metrics.py chose the newest subdirectory of a function
object's output by float(name). With faceErrorCsv true, areaFractionError
writes a directory faceError next to 0. float('faceError') failed, and the
fallback sorted()[-1] took faceError, which holds no .dat, so every
area-fraction column (L1_alphaf ...) was empty. Affected: the faceErrorCsv
runs, sphereRotation3D_hexFineLine{Diag,X}Faces and _Cons. The earlier
per-face analysis read the .dat itself and is unaffected.

Fix (leia 2f240340): only directories whose name is a number are candidates.
Job 55860362:
- All 298 study cases converted by the previous and the fixed converter:
  290 byte-identical, and the 8 others exactly the faceErrorCsv cases.
- Those 8 re-converted in place (old files kept as
  advectorVoF.preDirFix.csv, tpfCommit carried over): 6 columns filled, no
  non-empty value changed.

Open defect of the opt-in variant (recorded, not fixed; 2026-11-02):
sweptConsistent computes the walls' faces as well. On such a face a piece of
a cell without a surface carries alpha_c V_piece, so liquid leaves through a
wall where the algebraic flux takes the boundary value (fixedValue 0). This
is the whole mass deviation of the consistent diagonal above: 5.6e-9 / 8.4e-10
of the mass. The default path (sweptConsistent false) never computes a wall
face in these runs.

## Is the oblique build-up in alpha or in the face update? (2026-11-02)

Tomislav's approval via the coordinating session. The question: does the
oblique build-up of L1(alphaf) belong to the transported volume fractions,
which any plane-based geometric transport would show, or to the DAVOF face
update? The test reconstructs face fractions from alpha alone and scores
them as the loop's are.

Fork davof-swept, pushed to tmaric only:
- f568faa: the post-processing app staticFaceFractions
  (apps/benchmark/staticFaceFractions; serial; no library source changes).
- be237e4: Allwmake's executable bit, lost in f568faa.
- 394afb6: a diagnostic line of the app.

leia 401f0bc6: configs and davof_static_faces.py, with the predictions in the
header of sphereRotation3D_hexFineLineDiagW.yaml. They were committed before
any run. The numbers seen before are listed there: the loop's slopes,
Eshape, the film counts, and plicRDF on the rotation. My own expectation
from the Eshape ratios was H-face.

Jobs:
- 55863223, 55863347: builds. 55863196 failed on the mode bit and its four
  dependent studies were cancelled.
- 55863224-27: the four studies.
- 55863587: the app over 8 cases x 37 written times. 55863228 failed on an
  unknown option.
- 55863306, 55863348: a smoke test at t = 0.

The static rule:
- Interface cells: 1e-8 < alpha < 1 - 1e-8.
- The plane: the exact sphere normal, positioned by alpha to 1e-11
  (vofCutCell), then repositioned with the normal at the first plane's
  centroid.
- A face takes its donor cell's clip, the owner for phi_f >= 0
  (cutFacePLIC). This is the dt -> 0 limit of the sweptFaceCells pieces with
  donorSide.
- A cell without a plane gives 0 or 1.
- The mean of the two cells' clips is a sensitivity check.

Scoring: areaFractionError's L1(alphaf), its counted faces and its
exactSphere state. Ec = mean |alpha - alpha_exact| over the exact state's
interface cells (416 / 1762 cells). Late slope: over the 28 written times in
[0.25 T, T] (writes every 0.0125).

Note on Eshape: volumeFractionError builds its alpha_exact with cutCellIso
(the cut by the isosurface of the vertex values; corrected 2026-11-02, this
note first said cutCellImpFunc, which is only a commented-out line there), not
the exactSphere state. The app's sum |alpha - alpha_exact| V is 0.09 to 0.88
of Eshape over the run (0 at t = 0), so Eshape carries a construction floor
that dilutes ratios between the lines. Fixed on 2026-11-02 (next section). The
Eshape column of the table below is the former one.

| run | N | loop L1 end | loop late slope (per step / written) | static L1 end | static late slope | floor (alpha_exact) L1 end / slope | Ec end | Ec late slope | Eshape end | L1(x) end |
|---|---|---|---|---|---|---|---|---|---|---|
| ours, x-z diagonal | 32 | 1.752e-2 | 2.773e-2 / 2.783e-2 | 1.626e-2 | 2.524e-2 | 6.67e-3 / -2.1e-4 | 2.373e-2 | 4.224e-2 | 4.126e-4 | 1.195e-3 |
| ours, x-z diagonal | 64 | 1.028e-2 | 1.804e-2 / 1.790e-2 | 1.005e-2 | 1.613e-2 | 3.29e-3 / -4.9e-5 | 1.330e-2 | 2.630e-2 | 1.142e-4 | 3.412e-4 |
| ours, x-line | 32 | 1.135e-2 | 8.059e-3 / 8.082e-3 | 1.096e-2 | 4.918e-3 | 6.71e-3 / 7.7e-5 | 1.227e-2 | 1.281e-2 | 3.829e-4 | 8.661e-4 |
| ours, x-line | 64 | 6.278e-3 | 2.988e-3 / 3.045e-3 | 6.239e-3 | 2.837e-3 | 3.40e-3 / -4.8e-5 | 6.374e-3 | 4.605e-3 | 9.986e-5 | 2.290e-4 |
| plicRDF, x-z diagonal | 32 | - | - | 2.986e-2 | 5.493e-2 | (as ours) | 4.156e-2 | 8.494e-2 | 6.130e-4 | 1.982e-3 |
| plicRDF, x-z diagonal | 64 | - | - | 1.674e-2 | 3.266e-2 | | 2.221e-2 | 4.804e-2 | 1.691e-4 | 5.308e-4 |
| plicRDF, x-line | 32 | - | - | 2.022e-2 | 3.005e-2 | | 2.527e-2 | 4.604e-2 | 5.048e-4 | 1.410e-3 |
| plicRDF, x-line | 64 | - | - | 1.108e-2 | 1.487e-2 | | 1.294e-2 | 2.215e-2 | 1.408e-4 | 3.764e-4 |

All eight runs completed, bounded and conservative:
- min alpha at or above -3.1e-13; excess at most 4.5e-12.
- Mass deviation at the end at most 5.3e-12 %.
- plicRDF has no alphaf, so there is no loop column.
- The N = 64 rows' tpfCommit reads 394afb6: the converter ran after that
  build. The solver and libraries were unchanged, and K-0 holds.

Verdicts:
- K-0 HELD:
  - The W re-runs' areaFractionError.dat and volumeFractionError.dat are
    byte-identical to the originals in all four cases.
  - The app's L1 of the loop's alphaf matches the .dat at all 36 written
    times to at most 5.4e-12 relative.
  - 0 failed plane positionings.
- H-alpha HELD: phi = (s_S,D - s_S,X)/(s_L,D - s_L,X) = 1.029 (N = 32) and
  0.895 (N = 64). The static face fractions of the transported alpha build up
  on the diagonal like the loop's (static s_D/s_X 5.13 / 5.69 against the
  loop's 3.44 / 5.88). The mean rule gives phi 0.998 / 0.963.
- H-face FALSIFIED (it needs phi <= 0.25).
- K-1 HELD: the floor's slopes are at most 2.1e-4 in magnitude, against the
  thresholds 0.1 s_L,X = 8.1e-4 / 3.0e-4. Floor-corrected phi: 1.043 / 0.895.
- E-1 HELD: s_Ec,D/s_Ec,X = 3.30 / 5.71. The excess sits in the exact
  interface cells: our Eout at the end is at most 2.2e-6 per exact interface
  cell, so the low-alpha film carries none of it.
- R-1: phi_R = 1.260 / 1.198, so plicRDF's alpha carries the oblique excess
  too.
  - Its own ratio is smaller: static s_D/s_X 1.83 / 2.20, Ec 1.85 / 2.17.
  - But its absolute excess is at least ours: static s_D - s_X 2.49e-2 /
    1.78e-2 against ours 2.03e-2 / 1.33e-2.

What the runs show:
- The oblique build-up is in the transported volume fractions. The face
  update reproduces what its alpha implies: at the end the loop's L1(alphaf)
  is within 0.6 % to 8 % of the static L1 of the same alpha. It starts from
  0 against the reconstruction floor's 6.7e-3 / 3.3e-3, and the floor stays
  flat.
- plicRDF + isoAdvector, from the same initial alpha with the same dt, has
  the oblique excess as well, with larger errors on both lines.
- Our Ec late slope is 3.6 / 4.8 times smaller than plicRDF's on the x-line,
  but only 2.0 / 1.8 times smaller on the diagonal.

Surprising:
- The registered expectation (H-face) was wrong. Eshape's diagonal/x-line
  ratio at the end (1.08 / 1.14) hides a factor of about 2 in Ec (1.94 /
  2.09), because of its construction floor.
- From N = 32 to 64, Ec's late slope falls by 1.61 on our diagonal against
  2.78 on our x-line, and by 1.77 / 2.08 for plicRDF.

Not examined: the mechanism in alpha. One candidate is the volume flux of the
oblique swept volumes.

## The shape error fixed, the preprint re-scored, the fair rows (2026-11-02)

Tomislav's decision via the coordinating session: fix the shape-error
measure, re-score the preprint's shape errors, and test the oblique excess on
fair rows (a computed normal for our loop too). The preprint itself is edited
by a separate writing agent.

Commits:
- Fork davof-swept 1f3f238, pushed to tmaric only.
- leia 16f59466: the fair-row configs with the predictions, committed before
  the runs.
- leia 2a865522: davof_rescore_eshape.py.
- This section and the MEASURED paragraph.

Jobs:
- 55866604: the build.
- 55866605: V-2.
- 55866606 and 55866618: the re-score (N <= 128, and N = 256).
- 55866607-10: the four fair-row studies.
- 55866611: the app over them.

### The fix

volumeFractionError took alpha_exact from cutCellIso: each cell cut by the
isosurface 0 of the implicit function's vertex values, linear along the
edges. That is not the exact sphere.

The new key exactState:
- exactSphere: the cell fractions of davofInit::exactSphere at origin(t).
  This is the state of initDavofFields and the reference of
  areaFractionError and staticFaceFractions.
- cutCellIso: the former construction, kept for comparison.
- auto, the default: exactSphere when functionType is sphere and every cell
  is an axis-aligned box (the condition of davofInit::exactSphere, which is
  fatal elsewhere), cutCellIso otherwise.

Why this default: it is the only construction here that is the exact sphere,
and it is the one every other measure of these studies uses. Bodies and
meshes where it does not apply keep their output byte for byte.

What changes, exactly:
- On spheres on axis-aligned boxes, the Eshape column.
- The first .dat header line, which now reads "volumeFractionError
  (exactState exactSphere)".
- Eshape keeps its normalisation: none. It is sum_c |alpha_c -
  alpha_exact,c| V_c in m^3, against the sphere volume 1.414e-2.
- Emass, Ebound, the solver and every other output are unchanged.

V-2 HELD (zero-velocity copy of the x-line N = 32 case, 3 steps, so the state
stays the initial exact one):
- exactState auto: Eshape 1.567e-19 m^3 at every step.
- exactState cutCellIso: 3.580e-4 m^3, the former floor (2.5 % of the
  sphere's volume).
- The app at t = 0: 0 exactly.

V-1 HELD: in all eight fair-row cases (fields written at 17 digits), the
fixed object's Eshape equals the app's sum |alpha - alpha_exact| V at every
one of the 36 written times to at most 6.2e-12 relative.

### The re-score of the preprint (arXiv v1, tag davof-v1.0)

The written end fields of the 114 completed cases were re-scored without
re-solving:
- Studies: sphereTranslation3D_hex, _hexExactSurface, _hexParaboloid,
  _hexParaboloidRingFit, _hexParaboloidN128par, _hexParaboloidN256.
- Method: staticFaceFractions -entry volumeFractionError -davofMethod
  exactSphere -evolOnly, at t = 0 and t = 1.
- The aborted runs have no end field and are tabulated as aborted.

Checks:
- The v1 tables' old Eshape_end equals the runs' .dat last rows in all 114
  cases (difference 0).
- The new E_shape at t = 0 is 0, or 3e-16 to 7e-16 m^3 for the parallel runs,
  in every study but sphereTranslation3D_hex.
- sphereTranslation3D_hex starts from its own quadraticFaces state, 1.064e-4
  / 2.848e-5 / 7.277e-6 m^3 from the exact sphere at N = 32 / 64 / 128, and
  its new E_shape includes that.

Full mapping: $W/preprint_eshape_map.csv (every cited number) and
$W/rescore_eshape.csv (every case), $W = /work/scratch/tm83tomy/davof/agents/swept.

N = 128, old -> new (order over N = 32, 64, 128: old -> new).

table_advection_sphere (the plane-clipped loop):

| row | Co 0.1 | Co 0.2 | Co 0.5 |
|---|---|---|---|
| DAVOF, ring normal | 3.902e-5 -> 3.338e-5 (1.82 -> 1.82) | 3.047e-5 -> 2.441e-5 (1.83 -> 1.84) | 1.576e-5 -> 7.394e-6 (2.02 -> 2.02) |
| DAVOF, exact normal | 2.190e-5 -> 1.585e-5 (1.91 -> 1.90) | 1.935e-5 -> 1.288e-5 (1.92 -> 1.90) | 1.575e-5 -> 7.293e-6 (2.00 -> 1.94) |
| DAVOF, plicRDF normal | 6.750e-5 -> 6.346e-5 (1.51 -> 1.48) | 6.264e-5 -> 5.830e-5 (1.39 -> 1.32) | 3.403e-5 -> 2.844e-5 (1.44 -> 1.02) |
| DAVOF, ring quadric | 3.980e-5 -> 3.447e-5 (1.86 -> 1.85) | 3.096e-5 -> 2.500e-5 (2.61 -> 2.74) | 1.575e-5 -> 7.316e-6 (2.00 -> 1.93) |
| plicRDF | 3.928e-5 -> 3.432e-5 (1.90 -> 1.89) | 3.078e-5 -> 2.478e-5 (1.91 -> 1.93) | 1.614e-5 -> 7.960e-6 (1.98 -> 1.94) |
| isoAlpha | 4.582e-5 -> 4.136e-5 (1.86 -> 1.85) | 3.526e-5 -> 2.999e-5 (1.86 -> 1.86) | 1.710e-5 -> 9.549e-6 (1.95 -> 1.76) |

table_advection_curved (the dual advection; old 2.302e-5 to 2.304e-5 in every
row, order 1.98):

| row | Co 0.1 | Co 0.2 | Co 0.5 | new order (0.1 / 0.2 / 0.5) |
|---|---|---|---|---|
| exact sphere, exact / own / ring (identical) | 2.127e-12 | 1.202e-12 | 1.792e-13 | none (quadrature level) |
| paraboloid, exact | 3.223e-8 | 2.397e-8 | 2.978e-9 | 4.04 / 4.05 / 3.64 |
| paraboloid, ring | 1.868e-6 | 1.751e-6 | 1.655e-6 | 2.63 / 2.67 / 2.74 |
| paraboloid, ring quadric | 1.631e-7 | 1.473e-7 | 1.293e-7 | 3.42 / 3.25 / 2.30 |

The old values were the construction floor itself: 3.580e-4 / 9.078e-5 /
2.302e-5 at N = 32 / 64 / 128, identical in every row whose transport is
exact or nearly so.

tab:tpf-n256 (Co 0.2; old "all rows 2.30e-05 5.74e-06 2.00"):

| row | N = 128 | N = 256 | pairwise order |
|---|---|---|---|
| exact (control) | 2.397e-8 (parallel run 2.436e-8) | 2.294e-9 | 3.39 |
| ring quadric | 1.473e-7 (parallel 1.476e-7) | 3.206e-8 | 2.20 |
| ring mean | 1.751e-6 | 3.298e-7 | 2.41 |

Preprint statements that change in substance (the edits are the writing
agent's):
- "The centroid and the shape error are second order (1.98) in every row"
  (the exact-sphere surface rows): the shape error is 1.8e-13 to 2.1e-12 m^3.
  The 1.98 was the floor's order.
- Paraboloid, exact normal, "the centroid (1.99) and the shape error (1.98)
  as in the control": 2.40e-8 at N = 128, order 4.05, not the control's.
- Paraboloid, ring normal, "the shape error (1.98) is 2.30e-5 against
  3.08e-5": 1.75e-6 (order 2.67) against plicRDF's 2.48e-5, 14 times below.
- Ring quadric, "the shape error is the control's": 1.47e-7 (order 3.25),
  6.1 times the exact-normal control.
- tab:tpf-n256's single E_shape row: three rows, orders 3.39 / 2.20 / 2.41.
- Plane-clipped loop against plicRDF ("the shape error is plicRDF's", 3.05e-5
  against 3.08e-5): 2.441e-5 against 2.478e-5. It holds, within 1.5 %. Both
  include the initial state's 7.28e-6; at Co 0.5 both are within 2 % and 9 %
  of it.
- Not re-scorable here: the perturbed-mesh statements ("converging shape
  errors", "the shape error below it"). davofInit::exactSphere needs
  axis-aligned boxes, so these stay on cutCellIso. Their comparison with a
  hexahedral value must not mix constructions. An exact sphere-polyhedron
  cell volume would be needed; it is not implemented.

### The fair rows

Row (b), RingSwept: the swept flux with the ring normal. Row (c), RingPlane:
the preprint's plane-clipped loop with the ring normal. Both use LineDiagW /
LineXW settings.

All eight runs completed, bounded (min alpha >= -6.4e-13) and conservative
(mass deviation <= 3.8e-12 %).

E_shape is the fixed one at the end. For (a) and plicRDF it is the app's sum,
equal to the fixed object by V-1.

| row | N | line | loop L1(alphaf) end | loop late slope | static L1 end | static late slope | Ec end | Ec late slope | E_shape end | L1(x) end |
|---|---|---|---|---|---|---|---|---|---|---|
| (a) ours, exact normal, swept | 32 | diag | 1.752e-2 | 2.783e-2 | 1.626e-2 | 2.524e-2 | 2.373e-2 | 4.224e-2 | 3.013e-4 | 1.195e-3 |
| | 32 | x | 1.135e-2 | 8.082e-3 | 1.096e-2 | 4.918e-3 | 1.227e-2 | 1.281e-2 | 1.632e-4 | 8.661e-4 |
| | 64 | diag | 1.028e-2 | 1.790e-2 | 1.005e-2 | 1.613e-2 | 1.330e-2 | 2.630e-2 | 8.928e-5 | 3.412e-4 |
| | 64 | x | 6.278e-3 | 3.045e-3 | 6.239e-3 | 2.837e-3 | 6.374e-3 | 4.605e-3 | 4.242e-5 | 2.290e-4 |
| plicRDF + isoAdvector | 32 | diag | - | - | 2.986e-2 | 5.493e-2 | 4.156e-2 | 8.494e-2 | 5.280e-4 | 1.982e-3 |
| | 32 | x | - | - | 2.022e-2 | 3.005e-2 | 2.527e-2 | 4.604e-2 | 3.364e-4 | 1.410e-3 |
| | 64 | diag | - | - | 1.674e-2 | 3.266e-2 | 2.221e-2 | 4.804e-2 | 1.492e-4 | 5.308e-4 |
| | 64 | x | - | - | 1.108e-2 | 1.487e-2 | 1.294e-2 | 2.215e-2 | 8.615e-5 | 3.764e-4 |
| (b) ring normal, swept | 32 | diag | 3.279e-2 | 5.816e-2 | 2.845e-2 | 5.264e-2 | 4.651e-2 | 8.610e-2 | 5.906e-4 | 2.071e-3 |
| | 32 | x | 1.839e-2 | 2.314e-2 | 1.846e-2 | 2.299e-2 | 2.240e-2 | 3.624e-2 | 2.982e-4 | 1.286e-3 |
| | 64 | diag | 1.831e-2 | 3.659e-2 | 1.770e-2 | 3.420e-2 | 2.415e-2 | 5.178e-2 | 1.621e-4 | 5.696e-4 |
| | 64 | x | 1.033e-2 | 1.276e-2 | 1.039e-2 | 1.243e-2 | 1.195e-2 | 1.922e-2 | 7.953e-5 | 3.495e-4 |
| (c) ring normal, plane-clipped (preprint) | 32 | diag | 3.363e-2 | 5.909e-2 | 2.317e-2 | 3.588e-2 | 4.821e-2 | 8.965e-2 | 6.123e-4 | 2.141e-3 |
| | 32 | x | 1.941e-2 | 2.313e-2 | 1.854e-2 | 2.184e-2 | 2.235e-2 | 3.620e-2 | 2.976e-4 | 1.280e-3 |
| | 64 | diag | 1.890e-2 | 3.745e-2 | 1.362e-2 | 2.403e-2 | 2.491e-2 | 5.354e-2 | 1.674e-4 | 5.896e-4 |
| | 64 | x | 1.069e-2 | 1.298e-2 | 1.047e-2 | 1.264e-2 | 1.192e-2 | 1.922e-2 | 7.936e-5 | 3.486e-4 |

The oblique excess. Delta = s_Ec,D - s_Ec,X; "share" is Delta / Delta(a);
phi = (s_S,D - s_S,X)/(s_L,D - s_L,X), each row against its own loop.

| row | Delta N 32 / 64 | share | s_Ec,D/s_Ec,X | phi |
|---|---|---|---|---|
| (a) | 2.943e-2 / 2.169e-2 | 1 / 1 | 3.30 / 5.71 | 1.03 / 0.89 |
| plicRDF | 3.891e-2 / 2.589e-2 | 1.32 / 1.19 | 1.85 / 2.17 | - |
| (b) | 4.985e-2 / 3.256e-2 | 1.69 / 1.50 | 2.38 / 2.69 | 0.85 / 0.91 |
| (c) | 5.345e-2 / 3.431e-2 | 1.82 / 1.58 | 2.48 / 2.79 | 0.39 / 0.47 |

Verdicts:
- O-b PRESENT and O-c PRESENT: shares >= 0.5 and ratios >= 1.5 at both N.
- A-b HELD (phi 0.85 / 0.91): with the swept face update, the face fractions
  follow the alpha.
- A-c PARTIAL (phi 0.39 / 0.47): the plane-clipped loop's face fractions
  build up on the diagonal beyond what its alpha implies. The loop's late
  slope is 5.9e-2 / 3.7e-2 against its static 3.6e-2 / 2.4e-2.

What the runs show:
- With a computed normal, our transport's alpha is plicRDF's on both lines.
  - Ec late slope, (b) against plicRDF, x-line: 3.62e-2 vs 4.60e-2 (N = 32),
    1.92e-2 vs 2.22e-2 (N = 64).
  - Diagonal: 8.61e-2 vs 8.49e-2, 5.18e-2 vs 4.80e-2.
  - The exact-normal row's advantage over plicRDF (1.8 to 4.8 times in Ec's
    slope) is the normal's.
- The swept flux and the plane flux give nearly the same alpha with the ring
  normal (Ec at the end within 4 %).
- They differ in the face fractions: the plane-clipped loop's face update
  adds an oblique excess of its own, and the swept face update does not.

Surprising:
- The former E_shape of every exact or nearly exact transport in the
  preprint was the construction floor of the reference.
- With the fixed measure, the paraboloid rows' shape error at N = 128 is 14
  (ring) to 1030 (exact normal) times below plicRDF's, at orders 2.3 to 4.05.
  That comparison uses plicRDF from another initial state; the same-state
  ratios are in the next section.

## The shape error on the same state; the paraboloid loop on oblique translation (2026-11-02)

Tomislav's decision via the coordinating session. The dates of the previous
section's headers were corrected from 2026-11-03 to 2026-11-02; the commit
messages of fork 1f3f238 and leia 16f59466 / 2a865522 / dd575392 keep the
wrong date.

Commits:
- leia 1a3335d7: configs and predictions, committed before the runs.
- This section, the MEASURED paragraph and the date corrections.

Jobs:
- 55868371: K-0c.
- 55868372: the re-score of the timing study.
- 55868373-74: the two paraboloid line studies (6 cases each).
- 55868375: the app over them.

CSVs: $W/same_state_eshape.csv and $W/oblique_rows.csv, with $W =
/work/scratch/tm83tomy/davof/agents/swept.

### The shape error on the same state

The preprint's cost paragraph points to plicRDF on the exact initial state of
the dual advection: sphereTranslation3D_hexParaboloidTiming, rows plicRDF and
ring mean.
- What it holds: N = 32 / 64 / 128 at Co 0.2 only, three repetitions each.
  The repetitions agree to every digit.
- Its end fields were re-scored with the fixed measure; nothing was re-run.
- E(t = 0) is 0.
- Its ring-mean row reproduces the table study's E_shape and L1(x) to every
  digit.
- Co 0.1 and 0.5 do not exist for plicRDF on this state.

Same state, Co 0.2, end of the run. Ratios are row / plicRDF.

| row | N | E_shape | ratio | L1(n) | ratio | L1(x) | ratio |
|---|---|---|---|---|---|---|---|
| plicRDF + isoAdvector | 32 | 3.231e-4 | 1 | 3.852e-2 | 1 | 1.385e-3 | 1 |
| paraboloid, exact normal | 32 | 6.540e-6 | 0.0202 | 1.926e-4 | 0.0050 | 4.081e-4 | 0.295 |
| paraboloid, ring quadric | 32 | 1.324e-5 | 0.0410 | 1.338e-3 | 0.0347 | 4.100e-4 | 0.296 |
| paraboloid, ring mean | 32 | 7.053e-5 | 0.218 | 3.273e-2 | 0.850 | 5.178e-4 | 0.374 |
| plicRDF + isoAdvector | 64 | 8.109e-5 | 1 | 1.938e-2 | 1 | 3.411e-4 | 1 |
| paraboloid, exact normal | 64 | 3.649e-7 | 0.0045 | 2.840e-5 | 0.0015 | 9.586e-5 | 0.281 |
| paraboloid, ring quadric | 64 | 9.896e-7 | 0.0122 | 3.613e-4 | 0.0186 | 9.606e-5 | 0.282 |
| paraboloid, ring mean | 64 | 9.638e-6 | 0.119 | 1.760e-2 | 0.908 | 1.169e-4 | 0.343 |
| plicRDF + isoAdvector | 128 | 2.157e-5 | 1 | 9.125e-3 | 1 | 9.069e-5 | 1 |
| paraboloid, exact normal | 128 | 2.397e-8 | 0.00111 | 1.120e-5 | 0.00123 | 2.574e-5 | 0.284 |
| paraboloid, ring quadric | 128 | 1.473e-7 | 0.00683 | 1.632e-4 | 0.0179 | 2.572e-5 | 0.284 |
| paraboloid, ring mean | 128 | 1.751e-6 | 0.0812 | 8.191e-3 | 0.898 | 2.851e-5 | 0.314 |

Orders over N = 32, 64, 128 (E_shape / L1(n) / L1(x)):
- plicRDF: 1.95 / 1.04 / 1.97.
- Paraboloid, exact normal: 4.05 / 2.05 / 1.99.
- Paraboloid, ring quadric: 3.25 / 1.52 / 2.00.
- Paraboloid, ring mean: 2.67 / 1.00 / 2.09.

plicRDF on the exact state against plicRDF on the plane-loop study's
quadraticFaces state, N = 128, Co 0.2:
- E_shape 2.157e-5 against 2.478e-5.
- L1(n) 9.125e-3 against 9.20e-3 (the preprint's number).
- L1(x) 9.069e-5 against 1.04e-4.

### The paraboloid loop on oblique translation

The preprint's dual advection (rows exact normal, ring quadric, ring mean;
Kromer-Bothe, prism flux) at the LineDiagW / LineXW settings.

K-0c HELD. Copies of the preprint's N = 32, Co 0.2 cases, re-run with fork
1f3f238, end at the original values:
- Ring mean and ring quadric: bit for bit in areaFractionError and
  reconstructionError.
- Exact normal: within 4.7e-7 relative.

B-0, all 12 runs:
- Completed, min alpha >= -6.1e-13, excess <= 2.7e-13.
- Mass deviation <= 4.8e-12 %.
- No fit fallback, no anchor failure.
- CPU about 45 s (N = 32) and 465 to 485 s (N = 64) per run, serial.

N = 128 was not run. Serial it would take about 2 h per case. With ranks,
the 36 written 4.2-million-cell states would have to be reconstructed for
the serial app.

| row | N | line | loop L1 end | loop late slope | Ec end | Ec late slope | E_shape end | L1(x) end |
|---|---|---|---|---|---|---|---|---|
| paraboloid, exact normal | 32 | diag | 8.362e-4 | 1.280e-3 | 1.115e-3 | 1.997e-3 | 1.428e-5 | 4.087e-4 |
| | 32 | x | 5.768e-4 | 3.824e-4 | 5.464e-4 | 5.835e-4 | 7.275e-6 | 4.025e-4 |
| | 64 | diag | 1.218e-4 | 1.912e-4 | 1.487e-4 | 2.940e-4 | 1.009e-6 | 9.592e-5 |
| | 64 | x | 7.138e-5 | 2.691e-5 | 6.640e-5 | 4.717e-5 | 4.427e-7 | 1.006e-4 |
| paraboloid, ring quadric | 32 | diag | 1.264e-3 | 2.397e-3 | 1.734e-3 | 3.491e-3 | 2.214e-5 | 4.151e-4 |
| | 32 | x | 9.128e-4 | 1.247e-3 | 1.058e-3 | 1.774e-3 | 1.408e-5 | 4.069e-4 |
| | 64 | diag | 2.316e-4 | 3.659e-4 | 2.883e-4 | 5.268e-4 | 1.951e-6 | 9.619e-5 |
| | 64 | x | 1.415e-4 | 1.195e-4 | 1.552e-4 | 1.991e-4 | 1.033e-6 | 1.008e-4 |
| paraboloid, ring mean | 32 | diag | 1.594e-2 | 2.259e-2 | 2.265e-2 | 3.168e-2 | 2.878e-4 | 1.058e-3 |
| | 32 | x | 5.838e-3 | 4.210e-4 | 5.001e-3 | 5.423e-4 | 6.658e-5 | 5.389e-4 |
| | 64 | diag | 8.528e-3 | 1.562e-2 | 1.103e-2 | 2.176e-2 | 7.454e-5 | 2.867e-4 |
| | 64 | x | 1.929e-3 | -3.714e-5 | 1.457e-3 | 7.813e-5 | 9.699e-6 | 1.200e-4 |

The excess, next to the plane rows of the previous section. Delta = s_Ec,D -
s_Ec,X, the same N; "share" is Delta / Delta(a); p_D / p_X is the order
log2(Ec(32)/Ec(64)) of Ec at the end on each line.

| row | Delta N 32 / 64 | share | s_Ec,D/s_Ec,X | p_D / p_X |
|---|---|---|---|---|
| (a) plane, exact normal, swept | 2.943e-2 / 2.169e-2 | 1 / 1 | 3.30 / 5.71 | 0.84 / 0.94 |
| plicRDF | 3.891e-2 / 2.589e-2 | 1.32 / 1.19 | 1.85 / 2.17 | 0.90 / 0.97 |
| (b) plane, ring, swept | 4.985e-2 / 3.256e-2 | 1.69 / 1.50 | 2.38 / 2.69 | 0.95 / 0.91 |
| (c) plane, ring, plane-clipped | 5.345e-2 / 3.431e-2 | 1.82 / 1.58 | 2.48 / 2.79 | 0.95 / 0.91 |
| paraboloid, exact normal | 1.414e-3 / 2.468e-4 | 0.048 / 0.011 | 3.42 / 6.23 | 2.91 / 3.04 |
| paraboloid, ring quadric | 1.717e-3 / 3.277e-4 | 0.058 / 0.015 | 1.97 / 2.65 | 2.59 / 2.77 |
| paraboloid, ring mean | 3.114e-2 / 2.168e-2 | 1.06 / 1.00 | 58 / 278 | 1.04 / 1.78 |

Verdicts:
- H-curved-free NOT HELD for any row: it needs both Delta < 0.2 Delta(a)
  and ratio < 1.5.
- H-curved-excess HELD for all three rows by the registered measure (ratio
  >= 1.5 at both N). The ratio measure and the absolute excess disagree for
  the exact normal and the ring quadric:
  - Their Delta is 0.011 to 0.058 of the plane row's.
  - Their Delta falls 5.7 / 5.2 times from N = 32 to 64.
  - Their Ec converges at third order on both lines (2.6 to 3.0).
  - So the diagonal keeps a relative excess at a much smaller error level.
- The ring mean keeps the full absolute excess (Delta 1.06 / 1.00 of the
  plane row's). Its x-line levels off (Ec late slope 5.4e-4 / 7.8e-5) while
  its diagonal grows at 3.2e-2 / 2.2e-2. Its Ec order is 1.04 on the
  diagonal against 1.78 on the x-line.
- O-N: not applicable, since H-curved-free did not hold. The orders are
  reported above.
- phi is not meaningful for these rows: the static construction is plane
  based, the paraboloid loop is not.

What the runs show:
- With the exact normal or the ring quadric, the paraboloid loop reduces the
  oblique excess of Ec 17 to 88 times in absolute terms against the plane
  rows. Ec, E_shape and L1(alphaf) on the diagonal stay 2 to 3 times the
  x-line's, at third order in N on both lines.
- With the ring mean, the excess is as large as the plane rows'.
- L1(x) shows no oblique excess for the exact normal and the ring quadric
  (diagonal 4.09e-4 / 9.59e-5 against x-line 4.03e-4 / 1.01e-4).

## The swept-face advection in the preprint's framing: A to D (2026-11-02)

Via the coordinating session. The preprint presents sweptFaceCells as its
advection scheme.

Commits (leia):
- f7ea4b7f: davof_rescore_cases.py.
- bc2a21e8: the TPF_SWEPT_SCALE token in the translation templates; the B, C
  and D configs; the predictions, committed before the runs.
- cececa3f: davof_swept_ladders.py.
- abccf02d: G-2, pre-registered after G-1 failed.
- This section and the MEASURED paragraphs.

The binary is the fork 1f3f238 build (55866604); a7e6853 changes comments only.

Jobs:
- A: 55868997 and 55868998.
- B: 55869109; G-2: 55870149.
- C: 55869110-13.
- D: 55869114-15.

CSVs ($W = /work/scratch/tm83tomy/davof/agents/swept):

| item | file |
|---|---|
| A | $W/rescoreA_eshape_end.csv (end values, orders) and $W/rescoreA_eshape_times.csv (every written time) |
| B | $W/itemB_plane_gate.csv and $W/itemB_plane_gate_L1.csv |
| C | $W/itemC_rotation_ladder.csv |
| D | $W/itemD_sphere_ladder.csv (all cases) and $W/itemD_comparison.csv (N = 128 with orders, next to the preprint's plane-clipped rows and plicRDF) |

The token TPF_SWEPT_SCALE (sweptScaleToPhi) was added to the plane and
sphere translation templates. Its default is false, the fork's default, so
every earlier study renders the same dictionary values. The other swept keys
of the final configuration are the fork's defaults: pieces, donorSide, cap set
final, planeSource anchor.

### A: the fixed E_shape of the rotation ladder and the line runs (no re-solve)

110 current cases of sphereRotation3D_* were re-scored at every written time;
superseded copies (VOID, pre-fix) were excluded. E(t = 0) is 0 in every case.
For an aborted case the last written time is the abort time, not the end.

Rotation ladder, sphereRotation3D_hexSweptCellsFixed (pieces, donorSide), end
of the turn:

| row | N 32 | N 64 | N 128 | order |
|---|---|---|---|---|
| exact normal | 3.997e-4 | 1.285e-4 | 3.869e-5 | 1.68 (old 1.76) |
| ring normal | 7.602e-4 | 2.114e-4 | 5.738e-5 | 1.86 (old 1.87) |

Other studies (old -> new at the end):
- The donor-face variant (hexSweptCellsDonorFace) agrees to 3 digits.
- plicRDF on the rotation (hexComparators, zeroGradient patches): 9.081e-4 /
  2.437e-4 -> 8.999e-4 / 2.416e-4.
- The lines (N 32 / 64), diagonal / x-line:
  - Default swept flux: 3.013e-4 / 1.632e-4 and 8.928e-5 / 4.242e-5 (old
    4.126e-4 / 3.829e-4 and 1.142e-4 / 9.986e-5).
  - plicRDF: 5.280e-4 / 3.364e-4 and 1.492e-4 / 8.615e-5.

The full list is in the CSVs.

The fix moves the rotation values by at most 15 %; the floor (3.58e-4 at
N = 32) is small against these errors. It moves the lines' x-line values by
a factor of 2.3: there the floor was most of the old number.

### B: the translated-plane gate with the swept-face advection

The preprint's planeTranslation3D matrix: Co 0.2, N = 20 / 40 / 80 / 160, both
arms, rows exact and ring normal (planeTranslation3D_hexSwept).

G-1 FALSIFIED: no case passes the 1e-12 gate. Run maxima:

| row, arm | L1(alphaf) N 20 / 40 / 80 / 160 | L1(n) |
|---|---|---|
| exact, normal | 1.52e-7 / 3.04e-8 / 1.60e-8 / 2.06e-7 | <= 1.7e-14 |
| exact, oblique | 5.75e-6 / 2.46e-6 / 3.20e-6 / 2.84e-6 | <= 1.5e-14 |
| ring, normal | 1.62e-7 / 3.37e-8 / 1.61e-8 / 2.05e-7 | 1.8e-7 to 7.3e-9 |
| ring, oblique | 6.04e-6 / 2.65e-6 / 3.85e-6 / 3.37e-6 | 9.4e-6 to 4.9e-6 |

For comparison, the preprint's plane-clipped loop on the same cases has run
maxima 2.9e-15 to 4.2e-14.

The rest of the measurement:
- Bounded: min alpha >= -9.7e-13, excess <= 5.2e-11.
- The mass deviation at the end, 97.52 % (normal) and 51.16 % (oblique), is
  the case's own: the plane leaves through the boundary. The preprint's
  plane-clipped runs give the same values (97.516 / 51.157), so G-1's mass clause
  (1e-10 %) was mis-set for this open case.
- The error appears at the second step. At N = 20 normal, step 1 is 6.1e-15
  and step 2 is 1.0e-8.
- Per step, 8 % (oblique) to 44 % (normal) of the swept faces take the
  one-donor fallback (nSweptRankFallback, in serial).
- The pieces' closure residual reaches 0.62 to 0.72.
- The partial cells without a surface grow from 147 to 507 over the run
  (N = 20, normal).

G-2 (the case scaled to a unit box, N = 20 / 40) NOT HELD: it fails as at L =
0.004 (maxima within a factor 1.15, the same fallback shares 0.10 to 0.44).
The failure is not a scale effect; its cause is open.

### C: the rotation ladder with the final code, and plicRDF on the same rotation

The sphere is turned once, Co 0.2, from the exactSphere state, with alpha
fixedValue. N = 128 runs on 64 ranks. Values are at the end of the turn.

| row | N | L1(alphaf) | L1(n) | L1(x) | E_shape | min alpha / excess | mass dev % |
|---|---|---|---|---|---|---|---|
| swept, exact normal | 32 | 2.305e-2 | 2.469e-3 | 1.491e-3 | 3.997e-4 | -5.6e-13 / 6.3e-12 | 4.2e-13 |
| | 64 | 1.452e-2 | 8.221e-4 | 4.722e-4 | 1.285e-4 | -2.7e-13 / 1.1e-11 | 4.9e-12 |
| | 128 | 8.943e-3 | 2.419e-4 | 1.392e-4 | 3.869e-5 | -3.3e-13 / 8.0e-11 | 1.1e-12 |
| swept, ring normal | 32 | 4.158e-2 | 4.938e-2 | 2.809e-3 | 7.603e-4 | -8.0e-13 / 4.0e-12 | 1.0e-12 |
| | 64 | 2.323e-2 | 2.277e-2 | 7.587e-4 | 2.114e-4 | -1.0e-14 / 1.6e-11 | 6.3e-12 |
| | 128 | 1.322e-2 | 1.022e-2 | 2.024e-4 | 5.738e-5 | -1.1e-13 / 8.2e-11 | 9.3e-13 |
| plicRDF + isoAdvector | 32 | - | 5.350e-2 | 3.483e-3 | 8.999e-4 | -8.9e-18 / 2.2e-16 | 1.3e-12 |
| | 64 | - | 2.132e-2 | 8.642e-4 | 2.416e-4 | -7.6e-13 / 4.4e-16 | 1.0e-11 |
| | 128 | - | 9.302e-3 | 2.310e-4 | 6.435e-5 | -3.1e-5 / 7.1e-7 | 2.0e-12 |

Orders over N = 32, 64, 128:

| row | L1(alphaf) | L1(n) | L1(x) | E_shape |
|---|---|---|---|---|
| exact | 0.68 | 1.68 | 1.71 | 1.68 |
| ring | 0.83 | 1.14 | 1.90 | 1.86 |
| plicRDF | - | 1.26 | 1.96 | 1.90 |

Verdicts:
- R-0 HELD (no rank fallback).
- R-1 HELD (0.83 / 0.68).
- R-2 HELD.
- R-3 HELD.
- R-4 HELD. At N = 128:
  - The exact row's L1(x) is 1.392e-4 and its E_shape 3.869e-5, against
    plicRDF's 2.310e-4 and 6.435e-5.
  - The ring row's are 2.024e-4 and 5.738e-5, below plicRDF's.

The final code reproduces the earlier ladder (sphereRotation3D_hexSweptCells
Fixed): E_shape to 4 digits, L1(alphaf) 8.94e-3 against 8.96e-3 (exact,
N = 128).

Cost at N = 128 on 64 ranks: rank 0 CPU 1143 s (exact) and 2105 s (ring).

plicRDF on 64 ranks at N = 128 is not bounded: min alpha -3.1e-5, excess
7.1e-7. Serial at N = 32 / 64 it stays within -7.6e-13.

### D: the translated-sphere ladder with the swept-face advection

The preprint's tab:tpf-sphere: the table's quadraticFaces state, Co 0.1 / 0.2 /
0.5; N = 128 on 64 ranks. The plane-clipped rows and plicRDF are the
preprint's own runs on the same state, with E_shape re-scored.

| Co | row | L1(alphaf) N128 (p) | L1(n) N128 (p) | L1(x) N128 (p) | E_shape N128 (p) |
|---|---|---|---|---|---|
| 0.1 | swept, exact | 3.723e-3 (0.87) | 5.365e-4 (1.45) | 7.928e-5 (1.92) | 1.587e-5 (1.90) |
| 0.1 | swept, ring | 7.666e-3 (0.79) | 1.261e-2 (0.94) | 1.341e-4 (1.86) | 3.340e-5 (1.82) |
| 0.1 | plane-clipped, exact | 3.842e-3 (0.85) | 5.553e-4 (1.43) | 7.926e-5 (1.92) | 1.585e-5 (1.90) |
| 0.1 | plane-clipped, ring | 7.723e-3 (0.80) | 1.261e-2 (0.94) | 1.335e-4 (1.86) | 3.338e-5 (1.82) |
| 0.1 | plicRDF | - | 1.142e-2 (1.03) | 1.397e-4 (1.90) | 3.432e-5 (1.89) |
| 0.2 | swept, exact | 2.963e-3 (0.86) | 4.140e-4 (1.54) | 6.753e-5 (1.93) | 1.289e-5 (1.90) |
| 0.2 | swept, ring | 5.556e-3 (0.83) | 1.091e-2 (0.95) | 1.024e-4 (1.86) | 2.443e-5 (1.84) |
| 0.2 | plane-clipped, exact | 3.322e-3 (0.82) | 4.704e-4 (1.48) | 6.758e-5 (1.93) | 1.288e-5 (1.90) |
| 0.2 | plane-clipped, ring | 5.763e-3 (0.81) | 1.092e-2 (0.95) | 1.022e-4 (1.86) | 2.441e-5 (1.84) |
| 0.2 | plicRDF | - | 9.203e-3 (1.04) | 1.038e-4 (1.94) | 2.478e-5 (1.93) |
| 0.5 | swept, exact | 1.284e-3 (0.98) | 2.137e-4 (1.81) | 4.748e-5 (1.98) | 7.301e-6 (1.93) |
| 0.5 | swept, ring | 1.622e-3 (1.12) | 8.170e-3 (1.03) | 4.876e-5 (2.00) | 7.342e-6 (2.04) |
| 0.5 | plane-clipped, exact | 3.462e-3 (0.99) | 5.543e-4 (1.90) | 4.759e-5 (1.98) | 7.293e-6 (1.94) |
| 0.5 | plane-clipped, ring | 2.958e-3 (1.00) | 8.120e-3 (1.02) | 4.880e-5 (1.99) | 7.394e-6 (2.02) |
| 0.5 | plicRDF | - | 6.030e-3 (0.97) | 4.871e-5 (1.99) | 7.960e-6 (1.94) |

All 18 swept cases are bounded (min alpha >= -1.0e-12, excess <= 5.0e-11)
and conservative (mass dev <= 2.9e-12 %), with no rank fallback.

Verdicts:
- T-0 HELD.
- T-1 HELD (Co 0.2: L1(x) 1.93 / 1.86, E_shape 1.90 / 1.84).
- T-2 HELD (0.86 / 0.83).
- T-3 HELD: the ring row's E_shape is 0.986 and its L1(x) 0.987 of
  plicRDF's; the exact row is below plicRDF in both.
- T-4 HELD at every Co.

What the runs show:
- On the axis-aligned translation the swept-face advection gives the plane
  loop's L1(x) and E_shape to 3 digits.
- Its L1(alphaf) is at or below the plane loop's at every Co. At Co 0.5 it
  is 2.7 times lower with the exact normal and 1.8 times lower with the ring
  normal.

Surprising:
- B: the swept-face advection does not keep a translated plane exact (errors
  up to 6e-6), while it keeps a translated sphere within plicRDF's errors.
- C: plicRDF on 64 ranks is not bounded at N = 128.

## The plane gate of the swept-face advection: the wisp pieces (2026-11-02/03)

The problem: B above. The piecewise swept flux (fluxIntegration
sweptFaceCells) does not keep a translated plane exact, although every piece
of a plane case is cut by the same plane. Fork davof-swept 0ea9caa .. fab068c,
leia davof-swept 0b4fdc20 .. edfddf2b. Predictions F-0..F-2, A-1..A-3, M-1,
M-2 and M9-1, M9-2 are in the config headers (planeTranslation3D_hexSweptWisp,
..WispA, ..WispM, ..WispM9), each written before its runs.

Jobs (ledger .my_jobs):
- diagnosis: 55870504, 55870604, 55870680, 55872695;
- F-0 before / after: 55871135, 55871147;
- neighbour studies: 55871149-51;
- variants: 55880564 (pg-var), 55880908 (pg-var2), 55883022 (pg-var3),
  55883256 (pg-var4);
- anchored studies: 55880912, 55880915, 55880916;
- matched studies: 55881463-65 (M), 55883569-71 (M9);
- builds: 55871146, 55880563, 55880907, 55881462, 55883255.

### The cause

Diagnosis on N = 20, normal arm, exact normal (jobs 55870504, 55870604,
55870680):
- After step 1, alpha is off by up to 1.1e-6 in interface cells.
- Of the 1550 faces computed in step 1, the 872 one-donor fallbacks and every
  face fraction are exact. Five faces computed by pieces are short by up to
  8.0e-7 V.
- On the worst face, one ring cell has alpha 1.64e-6: the plane cuts a small
  corner of it, so |m_c| is below the interface threshold wispTol V^(2/3)
  and the cell has no plane. The piece in it (6.2e-14 m^3) was taken dry by
  the wisp's 0/1 rule (wet for alpha >= 1/2).
- sweptConsistent (alpha_c V_piece) errs the same way: it spreads the
  corner's liquid over the whole cell.

The mechanism: a cell holding liquid but no plane gives its pieces a liquid
volume that is not the plane's.

### The fix and its three rounds

Keys (all opt-in; the defaults reproduce the build before, F-0 HELD bit for
bit on 20-step copies of the plane, x-line and rotation cases, job 55871147):
- sweptWispSurface: none (default) | neighbour | anchored | matched.
- sweptWispTol: the band of the wisp class (default 1e-8).
- sweptWispMatchTol: the tolerance of matched (default 1e-11).

The rounds:
- neighbour: the plane of the closest face neighbour with a surface clips
  the wisp's pieces and its part of the face update.
  - Exact on the plane.
  - On the rotation, accuracy within 1e-3 but NOT bounded: min alpha -1.7e-5
    (N = 64, exact normal), -9.8e-6 (ring). The neighbour's plane extended
    into a wisp of a curved interface holds more liquid there than the wisp
    has. F-2 FALSIFIED on boundedness.
- anchored: the neighbour's normal, the plane positioned in the wisp by
  alpha^n V.
  - Bounded on the rotation (A-3 HELD: min alpha >= -9.5e-13).
  - NOT exact on the plane with sweptWispTol 0 (A-1 FALSIFIED: oblique arm
    3.5e-11 to 3.6e-9). The cause is conditioning: a wisp of round-off alpha
    gets a positioned corner, and the face fraction of a corner goes as
    alpha^(2/3), so (1e-15)^(2/3) gives about 1e-10.
- matched: the neighbour's plane where the volume it puts in the wisp is
  alpha^n V within sweptWispMatchTol V, else anchored.
  - With the tolerance 1e-11 (M), the oblique arm at N = 160 fails (see
    below).
  - With 1e-9 (M9, the final rule), it is restored.

Two further causes surfaced in the plane matrix (job 55880564; N = 20 unless
stated):

| case | study keys | + sweptWispTol 0 | + sweptConsistent | both |
|---|---|---|---|---|
| oblique Co 0.2, max L1(alphaf) | 1.2e-8 | 8.6e-14 | 1.2e-8 | 6.0e-14 |
| oblique Co 0.95 | 1.3e-8 | 3.9e-14 | 1.3e-8 | 3.9e-14 |
| normal Co 0.95 | 0.115 | 0.116 | 7.8e-12 | 1.6e-13 |

- The oblique arm: one face fraction of step 1 was short by 1.05e-5. A
  corner below 1e-8 V of a ring cell was counted empty. Cured by
  sweptWispTol 0.
- The normal arm at Co 0.95 (OpenFOAM Co 1.53): the candidate band limits
  Co. The faces outside the face-neighbour band of the interface cells keep
  the algebraic flux and a stale face fraction, and at a displacement of
  0.77 h per axis their swept volumes reach liquid. In step 1, E_shape is
  2.5e-11 while L1(alphaf) is 0.115. Cured by the vertex-ring candidates
  (sweptConsistent true).

### The plane matrix with the final rule

The final rule: sweptWispSurface matched, sweptWispMatchTol 1e-9,
sweptWispTol 0, sweptConsistent true.

planeTranslation3D_hexSweptWispM9, job 55883569 (M with 1e-11: job
55881463). The gate: run maxima of
L1(alphaf) and L1(n) below 1e-12, min alpha >= -1e-12, excess <= 1e-11, and
Eshape / sum(alpha V) below 1e-12.

The mass clause is corrected for this open case. Eshape = int |alpha -
alpha_exact| dV bounds |M - M_exact| at whatever time the run ends. At
Co 0.95, END_TIME 8e-4 is not a multiple of dt and the run stops at 7.6e-4,
so the earlier comparison with the plane-clipped loop's mass at 8e-4 was not
valid there.

The Courant number: TARGET_ADVECTIVE_CO is of |U| (dt = Co h). OpenFOAM's
0.5 sum|phi_f| dt / V equals sum_k |U_k| dt / h, which is 1.612 Co on the
normal arm and Co on the oblique arm. The largest per-axis displacement is
0.8124 Co h and Co h respectively: below h at every Co, and the "beyond
layer" count is 0 from N = 40 on.

| class | passed | max L1(alphaf) / L1(n) | min alpha | max Eshape/M |
|---|---|---|---|---|
| OpenFOAM Co <= 1, N 20 / 40 | 23 / 24 | 4.9e-13 / 1.1e-12 | -1.5e-13 | 3.4e-13 |
| OpenFOAM Co <= 1, N 80 / 160 | 6 / 24 | 9.5e-12 / 8.3e-12 | -8.2e-13 | 7.1e-12 |
| OpenFOAM Co > 1 (normal arm, Co 0.8 / 0.95) | 2 / 16 | 3.1e-7 / 5.4e-7 | -4.3e-4 | 1.1e-4 |

31 / 64 pass, against 15 / 64 with neighbour (sweptWispTol 1e-8, the
candidate band), 9 / 64 with anchored, and 31 / 64 with M (1e-11).

The failures, by class:
- With OpenFOAM Co <= 1, every case is below 1e-11 in both norms, bounded
  to -8.2e-13, and conserves mass to 7.1e-12 (Eshape/M).
- At N <= 40 the only failure is the ring row, normal arm, Co 0.5, N = 40:
  L1(n) 1.14e-12 (the ring normal's own error; its L1(alphaf) is 4.9e-13).
- At N = 80 every oblique case at Co >= 0.5 passes. The normal arm and the
  oblique Co 0.2 sit at 1.3e-12 to 1.9e-12; at N = 160 everything sits at
  3.5e-12 to 9.5e-12. This is the floor (iii)-b below.
- With M (1e-11), the oblique arm at N = 160 had reached 7.6e-11 to 1.4e-9
  (M-1 FALSIFIED there). The error was there at step 2, not accumulated:
  - Three steps at N = 160 (job 55883022): 1.1e-11 with matched against
    1.4e-12 with neighbour. Excluding round-off cells (sweptWispTol 1e-14)
    changes nothing.
  - matched sent wisps that the plane does cut to the anchored branch, since
    their mismatch exceeded 1e-11 V at this resolution.
  - With sweptWispMatchTol 1e-9 or 1e-7 (job 55883256), step 3 is 1.378e-12,
    neighbour's value. The N = 32 rotation with 1e-7 stays bounded (min alpha
    -7.9e-13) and within 4e-7 of 1e-11.
- Verdicts: M9-1 HELD on the N = 160 clause (oblique 3.5e-12 to 9.5e-12,
  below 1e-11). M9-1 FALSIFIED on the count (31, predicted >= 33) and on N = 80
  (oblique Co 0.2 at 1.5e-12). "Within a factor 2 of M" is FALSIFIED in the
  OpenFOAM Co > 1 class only, where min alpha changes by up to 6.5x (the growth
  (iii)-a is sensitive to any change).

Open, with their discriminators:
- (iii)-a, OpenFOAM Co > 1 on the normal arm: min alpha grows by about x1.4
  per step. With the final rule it reaches -1.9e-10 at N = 80 and -1.9e-6 at
  N = 160 for Co 0.8, and -4.3e-4 at N = 160 for Co 0.95. It is
  unchanged by every wisp rule and by sweptConsistent. Thousands of one-donor
  rank fallbacks occur per step, and the bounding is extended to a second
  layer in hundreds of cells.
  - The hypothesis that the fallback, a donor-cell flux, is unstable above
    OpenFOAM Co 1 is FALSIFIED at 1.05: the normal arm at Co 0.65 (OpenFOAM
    1.05, N = 80, 24 steps) stays at -3e-14 (job 55880908).
  - The threshold lies between OpenFOAM Co 1.05 and 1.29. Cause open.
- (iii)-b, a floor that grows with N on the normal arm even at Co 0.2:
  L1(alphaf) 1.2e-13 / 1.3e-12 / 8.7e-12 at N = 40 / 80 / 160, Eshape/M up
  to 7e-12.
  - It is not the anchor's tolerance: anchorTolRel and curvedTolRel 1e-14
    change nothing, bit for bit (job 55880908).
  - The prism and the one-donor sweptFace keep the same cases at round-off
    through N = 160 (L1(alphaf) <= 2.4e-13, Eshape/M <= 1e-12, study
    planeTranslation3D_hexParaboloidSwept), so the floor is specific to the
    pieces: about 4e-14 per step against 5e-16. Cause open.

### The curved cases with the final rule

sphereRotation3D_hexRotFinalWispM9 and
sphereTranslation3D_hexSweptLadderWispM9 (jobs 55883570, 55883571; Co 0.2),
against the default:

| case | L1(alphaf) | L1(x) | E_shape | min alpha / excess | mass dev % |
|---|---|---|---|---|---|
| rotation exact N 32 | -6.8e-5 | -8.8e-5 | -1.3e-4 | -6.0e-13 / 3.9e-12 | 2.0e-12 |
| rotation exact N 64 | -2.9e-4 | -3.1e-4 | -5.4e-4 | -9.5e-13 / 1.1e-11 | 1.4e-11 |
| rotation ring N 32 | -6.4e-5 | -7.8e-4 | -1.1e-4 | -7.9e-13 / 5.6e-12 | 1.7e-12 |
| rotation ring N 64 | -3.2e-4 | -7.2e-4 | -4.8e-4 | -8.1e-13 / 1.6e-11 | 1.4e-11 |
| sphere x-line exact N 32 | -7.8e-8 | -4.5e-7 | -2.8e-6 | -2.3e-14 / 4.1e-13 | 5.0e-13 |
| sphere x-line exact N 64 | -2.5e-4 | +1.6e-4 | -2.0e-4 | -4.0e-14 / 9.1e-12 | 3.0e-12 |
| sphere x-line ring N 32 | -1.0e-5 | -2.3e-5 | -5.7e-5 | -2.2e-14 / 3.7e-13 | 5.1e-13 |
| sphere x-line ring N 64 | -3.8e-4 | -2.4e-4 | -5.6e-4 | -3.3e-14 / 9.1e-12 | 2.7e-12 |

The error columns are relative changes against the default; every change is
below 1e-3 and slightly favourable (M-2 and M9-2 HELD; the M and A rounds
agree to two digits). Bounded and conservative as the default (the
default's min alpha is -5.6e-13 to -1.0e-14; with neighbour it was -1.7e-5).

### Verdicts and the recommendation

The verdicts:
- F-0 HELD.
- F-1 FALSIFIED (the first rule, neighbour with the band 1e-8: 15 / 64).
- F-2 FALSIFIED on boundedness.
- A-1 FALSIFIED, A-2 FALSIFIED on both hypotheses, A-3 HELD.
- M-1 FALSIFIED (oblique N = 160), M-2 HELD.
- M9-1 HELD on N = 160 and FALSIFIED on the count and on N = 80, M9-2 HELD.

What the runs demonstrate:
- The plane gate's failure at Co 0.2 had two causes, both in the classes of
  ring cells without a plane:
  - the 0/1 rule of a wisp (the normal arm);
  - the empty band 1e-8 (the oblique arm).
- A third cause bounds Co: the face-neighbour candidate band.
- With the final rule, every plane case with OpenFOAM Co <= 1 (both arms up
  to the target Co 0.95 on the oblique one) is exact to 1e-11 in L1(alphaf)
  and L1(n), bounded to -8.2e-13, and conserves mass to 7.1e-12 relative.
  The 1e-12 gate is met at N <= 40 (23 / 24).
- The curved cases do not change beyond 1e-3 and stay bounded.

Limitations and open questions:
- (iii)-b: the pieces' round-off floor grows with N, to 9.5e-12 at N = 160.
  The prism stays at 2.4e-13 there.
- (iii)-a: min alpha grows on the normal arm above an OpenFOAM Co between
  1.05 and 1.29. Both causes are open.
- The cost of sweptConsistent (more computed faces) was not measured.

The recommendation: the final rule as the sweptFaceCells setting of the
studies (sweptWispSurface matched, sweptWispMatchTol 1e-9, sweptWispTol 0,
sweptConsistent true). The fork's defaults stay at the earlier behaviour
until Tomislav decides.

## The wall faces of the consistent variant (2026-11-03)

The problem: sweptConsistent (part of the final wisp rule) computed every
face of a candidate, also a wall face without a surface on either side or in
the second ring. A piece of a cell without a surface then carried
alpha_c V_piece through a wall whose alpha is fixedValue 0. This was the whole
mass loss of the consistent diagonal: 5.6e-7 % / 8.4e-8 % (N = 32 / 64), the
summed boundary outflow 5.626e-9 / 8.42e-10 of the mass.

The fix (fork eb17c0a, under sweptConsistent only): such a face keeps the
boundary condition's flux, the value of the default path, which skips every
face without a surface. This applies only where the patch's alpha condition
fixes the value. On a zeroGradient outflow patch (the plane gate's rim) the
algebraic alpha_owner phi_f dt is the approximation the swept flux replaces,
not a condition, and the swept flux stays. Counted per step ("consistent
boundary faces kept at the boundary condition").

Jobs: 55884523 / 55884531 (W-3, before / after), 55884530 (the build),
55884532-34 (W-1, W-2), 55885589 (the E_shape re-score).

| run | N | mass dev % | summed boundary outflow | late slope L1(alphaf) | L1(x) end | E_shape end (re-scored) | min alpha |
|---|---|---|---|---|---|---|---|
| default diagonal | 32 | 6.7e-13 | 0 | 2.77329e-2 | 1.19481e-3 | 3.01304e-4 | -1.4e-14 |
| consistent, before | 32 | 5.6e-7 | 5.63e-9 | 2.77312e-2 | 1.19474e-3 | 3.01286e-4 | -1.3e-18 |
| consistent, fixed | 32 | 5.7e-13 | 0 | 2.77312e-2 | 1.19474e-3 | 3.01286e-4 | -1.2e-18 |
| default diagonal | 64 | 3.2e-12 | 0 | 1.80365e-2 | 3.41190e-4 | 8.92760e-5 | -1.0e-14 |
| consistent, before | 64 | 8.4e-8 | 8.42e-10 | 1.80325e-2 | 3.41134e-4 | 8.92518e-5 | -5.0e-22 |
| consistent, fixed | 64 | 3.3e-12 | 0 | 1.80325e-2 | 3.41134e-4 | 8.92518e-5 | -5.1e-22 |

How the E_shape column was obtained: the consistent diagonal of 2026-11-02 ran
before the E_shape fix (its own Eshape column used the cutCellIso reference:
4.126e-4 / 1.141e-4). All three are re-scored here from the end fields with
staticFaceFractions -entry volumeFractionError -davofMethod exactSphere
-evolOnly, as in the preprint's re-score. The fixed run's own Eshape
(3.0129e-4 / 8.9252e-5) equals its re-score.

Verdicts:
- W-1 HELD: mass deviation at round-off; the counter positive from step 146 /
  332 (4764 / 7437 faces in total); slope, L1(x) and the re-scored E_shape
  within 3e-7 of the consistent run before the fix.
- W-2 HELD: the curved ladders of the final rule (rotation and sphere
  translation, N = 32 / 64, both normals) are identical to the last digit of
  the CSV; the counter is 0 in every step.
- W-3 HELD: two 20-step plane copies are bit for bit identical; the counter is
  0.

No number of the plane matrix (zeroGradient patch, the branch cannot trigger)
or of the curved ladders changes.

## The UFVFC pyramid volume correction (2026-11-03)

The problem: under a non-uniform velocity the traced swept volume V_g is not
the flux volume V_f = |phi_f| dt.
- V_g = (1/3) sum_T x_T . a_T over the closed list of 4F triangles.
- V_f uses phi of the half step (setFlow at t + dt/2), the temporal pair of
  the explicit midpoint trace.
- The scaling (sweptScaleToPhi) multiplies the liquid by V_f / V_g and
  divides by V_g, which vanishes for a face twisted near a rotation axis.

The mechanism of the pyramid (Tomislav's design; fork 720cfc6, key
sweptVolumeCorrection auto | none | scale | pyramid, default auto = the mode
of sweptScaleToPhi):
- Only the face points are traced.
- The apex of the pre-image fan starts at the mean of the traced vertices.
- The apex moves along the pre-image's outward area vector A_pre by
  delta = 3 (V_f - V_g) / |A_pre|. This is exact because the volume is linear
  in the apex and the fan's area vector does not depend on it.
- The corrected closed list is clipped and capped as before, and the face
  update uses the corrected fan.
- No tetrahedra or prisms are formed.

The two references:
- The correction is that of the unsplit face-vertex flux calculation:
  T. Marić, H. Marschall, D. Bothe, An enhanced un-split face-vertex
  flux-based VoF method, J. Comput. Phys. 371 (2018) 967-993,
  doi:10.1016/j.jcp.2018.03.048. No arXiv version was found.
- The volume sum is that of T. Tolle, D. Gründing, D. Bothe, T. Marić,
  triSurfaceImmersion: Computing volume fractions and signed distances from
  triangulated surfaces immersed in unstructured meshes, Comput. Phys.
  Commun. 273 (2022) 108249, doi:10.1016/j.cpc.2021.108249,
  arXiv:2101.08511.

The review against the design: the saved first version started the apex at
the traced face centre. It now starts at the vertex mean, so nothing but the
face points is traced.

Logged per step: max |V_f - V_g|/V_f, max |delta|/h, the twisted faces, the
pyramid skips and the largest per-axis displacement over h.

The default path is unchanged bit for bit (jobs 55884662 / 55884669: four
20-step copies). A sanity run (55884670) gates the studies.

Studies (jobs 55885449-57, the final wisp rule in every case). Tomislav's
change before the runs: Co 0.8 / 0.95 dropped, so Co 0.2 / 0.5 only; the
first submission was cancelled before it started:
- planeTranslation3D_hexSweptPyr: 32 cases against WispM9.
- sphereRotation3D_hexLineXCo{Scale,Pyr}: N 32 / 64.
- sphereRotation3D_hexRotCo{Scale,Pyr}: N 32 / 64, both normals.
- sphereRotation3D_hexSpinCo{Scale,Pyr}: N 32 / 64.
- sphereRotation3D_hexRotN128Co{Scale,Pyr}: 64 ranks.
- CSV: csv/pyramid_item1.csv on the cluster.

Measured:
- Uniform translation.
  - Plane: every run maximum within 2.9e-13 of the scaling's, with the same
    17 / 17 gate passes.
  - x-line: L1(alphaf), L1(x) and E_shape within 1.3e-10 relative.
  - max |V_f - V_g|/V_f <= 5.4e-13 and max |delta|/h <= 1.2e-13 (round-off).
- The rotation at target Co 0.2 (OpenFOAM Co 0.86 / 0.87 / 0.88 at
  N = 32 / 64 / 128):
  - max |V_f - V_g|/V_f is 1.4e-5 / 3.3e-6 / 8.0e-7, and max |delta|/h is
    3.2e-7 / 3.9e-8 / 4.7e-9.
  - The pyramid's L1(alphaf), L1(x) and E_shape lie 2e-8 to 8e-7 below the
    scaling's, both normals.
- The pure spin at target Co 0.2 (OpenFOAM Co 0.99): defect 5.8e-6 / 1.3e-6,
  E_shape 6e-6 / 2e-6 lower. min alpha 0 for both modes.
- No face is twisted (V_g below 0.1 of the unsigned tetrahedra) in any run,
  also not at the poles of the spin.
- Target Co 0.5 is OpenFOAM Co 2.16 to 2.19 on the rotation and 2.47 on the
  spin:
  - On the rotation the per-axis displacement reaches 1.57 h. Parts of the
    swept volumes leave the one-layer ring (largest uncovered share 0.82 to
    1.0), there are up to 1.1e8 rank fallbacks per run, and min alpha is
    -7e-8 (N = 32), -8e-7 (N = 64) and -1.7e-4 (N = 128) for both modes.
  - The pyramid lowers the excess there (to 4e-13 to 2.7e-12, from up to
    1e-8 with the scaling) but not the undershoot.
  - The spin at Co 0.5 (displacement 0.24 h at the interface) stays bounded.

Verdicts:
- Y-0 HELD.
- Y-1 FALSIFIED, narrowly and for both modes, by the rotation at N = 128 and
  OpenFOAM Co 0.88 (64 ranks): excess 2.8e-10 (scaling) and 3.3e-10
  (pyramid) against 1e-10. The mass stays within 1e-11 %. The ladder before
  the wisp rule had 8e-11 there. Every other case with OpenFOAM Co <= 1 holds.
- Y-2 HELD on accuracy (pyramid within 7e-5 of the scaling everywhere).
  FALSIFIED on the twisted faces: there are none.
- Y-3 HELD in direction: the E_shape gain grows from 8e-8 to 8e-7 at Co 0.2
  to 1e-5 to 5e-5 at Co 0.5, both N and normals. It is negligible in size.

By the rule registered for the comparison (pyramid only if Y-0 and Y-1
held), the comparison with plicRDF runs the scaling. The two differ by less
than 1e-6 relative at OpenFOAM Co <= 1.

What this demonstrates: the pyramid is exact on uniform translation, removes
the division by V_g, and is at least as accurate and bounded as the scaling
wherever the swept volumes stay in the one-layer ring. On these cases the two
are indistinguishable at OpenFOAM Co <= 1.

The limitation: beyond OpenFOAM Co 1 on the rotation the one-layer ring, not
the volume correction, sets the error.

## The comparison with plicRDF + isoAdvector at OpenFOAM Co 1, timed (2026-11-03)

The question (Tomislav, via the coordinating session): at OpenFOAM Co 1 (and
Co 0.2 for reference), from the same exact initial sphere, on the same mesh,
with the same dt and the same rank count, how do the swept-face advection and
plicRDF + isoAdvector compare in error norms and in total wall-clock time?
The standing rule applies: error norms are compared with the total wall-clock
time at the same number of MPI ranks.

Predictions C-1 .. C-4: config/davof/sphereRotation3D_hexCmpLineXSweptN32.yaml,
written before the runs.

Commits:
- leia 736b93d6: the tokens TPF_N_ALPHA_BOUNDS, TPF_SNAP_TOL, TPF_CLIP
  (defaults = the template before) and 18 preparation configs.
- Fork 720cfc6, built, unchanged during the timings.

Jobs:
- 55885736-53: the preparations.
- 55887088-106 and 55889018-22, 55889179-89, 55889532-40, 55889642-45,
  55889869-71: the timings.
- Some were cancelled before their start and resubmitted; the ledger
  .my_jobs says which.

CSVs for the writing agent, on the cluster:
- /work/scratch/tm83tomy/davof/agents/swept/csv/cmp_runs.csv: every run.
- csv/cmp_table.csv: per case and method.
- Raw data: cmp-time/<flow>_N<N>_co<0|1>/.

### The setup

Cases and ranks:
- The translated sphere along x and along the x-z diagonal, END_TIME 0.45.
- The rotation, one turn, END_TIME 1.
- alpha fixedValue walls; N = 32 / 64 / 128 on 1 / 8 / 64 ranks, for every
  method.

The Courant number is OpenFOAM's: the maximum over the cells of
0.5 sum |phi_f| dt / V, logged.
- Measured ratio to the target Co of the template: pi/5 (x),
  sqrt(2) pi/5 (diagonal), 4.320 / 4.359 / 4.379 (rotation at N = 32 / 64 /
  128, at the domain's corners).
- The runs are set to OpenFOAM Co 1.0000016 and 0.2000003 (dt rounded to six
  digits).
- They end at n dt with n = round(END_TIME / dt): for example 23 steps and
  t = 0.4576 on the x-line at N = 32 and Co 1.
- Note: the earlier rotation ladders "at Co 0.2" ran at OpenFOAM Co 0.86 to
  0.88.

The rows:
- swept, ring: the swept-face advection with the final wisp rule (matched,
  match tolerance 1e-9, band 0, sweptConsistent, the wall faces kept at the
  boundary condition), the ring normal, and the volume correction scale. By
  the rule registered before, the pyramid only if Y-0 and Y-1 held; Y-1 was
  falsified, and the two differ by less than 1e-6 here.
- swept, exact: the same with the exact normal, as a control.
- plicRDF: isoAdvection with reconstructionScheme plicRDF, set as the fork's
  advection test suite sets it:
  - snapTol 0, clip false, nAlphaCorr 1, nAlphaSubCycles 1;
  - nAlphaBounds 3 on the translations (test-discUniFlow);
  - nAlphaBounds 5 on the rotation (test-deformationSphere, test-Leveque,
    test-vortexShearedDisc, the fork's rotational tests).
  - The swept rows keep nAlphaBounds 3.

Diagnostics at the end only for every method: reconstructionError,
areaFractionError and volumeFractionError at the end, no profiling, the
fields written at the end.

The timing:
- One exclusive i01 node per job: Intel Xeon Platinum 9242, 2 x 48 cores,
  384 GB.
- Three repetitions of every method, the order rotated per repetition.
- Every solve is one srun step (--ntasks np --cpu-bind cores) on the same
  launcher and node.
- Reported: wall = the step's wall time (date around srun; sacct's Elapsed
  agrees to the second where recorded), median [min, max]; clock = the
  solver's ClockTime at the end, median; step = the median of timeStepTime.
- Where 30-minute jobs were needed, repetitions ran on different i01 nodes.
  The deflt partition stopped starting jobs at about 14:00.
- The three repetitions give identical error norms in every case.

Scored at the end: L1(x), E_shape (the exact sphere, m^3), L1(n) against the
exact sphere, min alpha / excess, mass deviation, and L1(alphaf) for ours.

### The translated sphere along x

| N, ranks | OF Co | row | L1(x) | E_shape | L1(n) | L1(alphaf) | min / excess | mass % | wall s | clock s | step s | wall / plicRDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32, 1 | 0.2 | swept ring | 1.025e-3 | 2.329e-4 | 4.204e-2 | 1.494e-2 | -2.0e-14 / 2.1e-13 | 7.8e-14 | 22.5 [22.5, 22.5] | 23 | 0.17 | 4.5 |
| | | swept exact | 7.351e-4 | 1.283e-4 | 5.539e-3 | 9.595e-3 | -1.8e-14 / 2.1e-13 | 3.3e-14 | 23.4 [23.4, 23.4] | 24 | 0.18 | 4.7 |
| | | plicRDF | 1.134e-3 | 2.562e-4 | 3.758e-2 | - | -1.9e-24 / 0 | 1.1e-14 | 5.0 [4.9, 5.0] | 5 | 0.02 | 1 |
| 32, 1 | 1 | swept ring | 3.998e-4 | 2.364e-7 | 3.152e-2 | 6.990e-3 | -1.1e-13 / 2.9e-13 | 1.4e-13 | 7.5 [7.5, 8.9] | 8 | 0.20 | 2.6 |
| | | swept exact | 4.024e-4 | 2.364e-7 | 1.595e-3 | 6.647e-3 | -6.4e-14 / 5.3e-13 | 2.0e-13 | 7.7 [7.7, 7.8] | 8 | 0.22 | 2.7 |
| | | plicRDF | 4.012e-4 | 2.363e-7 | 2.072e-2 | - | -1.0e-14 / 1.0e-14 | 6.7e-14 | 2.9 [2.9, 2.9] | 3 | 0.02 | 1 |
| 64, 8 | 0.2 | swept ring | 2.792e-4 | 6.185e-5 | 2.066e-2 | 8.394e-3 | -3.8e-14 / 5.0e-12 | 8.2e-13 | 84.2 [84.2, 85.9] | 84 | 0.37 | 5.8 |
| | | swept exact | 1.932e-4 | 3.371e-5 | 1.286e-3 | 5.334e-3 | -4.7e-14 / 4.9e-12 | 9.7e-13 | 87.7 [87.5, 88.1] | 88 | 0.38 | 6.0 |
| | | plicRDF | 2.919e-4 | 6.578e-5 | 1.759e-2 | - | -9.9e-15 / 0 | 2.8e-13 | 14.7 [14.6, 14.7] | 15 | 0.05 | 1 |
| 64, 8 | 1 | swept ring | 9.531e-5 | 2.325e-7 | 1.774e-2 | 3.485e-3 | -1.2e-13 / 2.4e-12 | 2.2e-14 | 20.3 [20.2, 20.8] | 20 | 0.38 | 3.5 |
| | | swept exact | 9.563e-5 | 2.325e-7 | 5.732e-4 | 3.271e-3 | -1.7e-13 / 2.4e-12 | 7.8e-14 | 21.0 [21.0, 21.1] | 21 | 0.39 | 3.6 |
| | | plicRDF | 9.549e-5 | 2.326e-7 | 1.232e-2 | - | -1.0e-14 / 1.0e-14 | -2.2e-14 | 5.9 [5.8, 5.9] | 6 | 0.05 | 1 |
| 128, 64 | 0.2 | swept ring | 7.666e-5 | 1.729e-5 | 1.008e-2 | 4.781e-3 | -9.7e-14 / 9.7e-12 | 2.4e-12 | 434.3 [430.0, 443.9] | 434 | 0.94 | 6.7 |
| | | swept exact | 5.272e-5 | 9.251e-6 | 3.595e-4 | 3.004e-3 | -8.8e-14 / 9.7e-12 | 2.5e-12 | 448.7 [447.9, 448.8] | 448 | 0.98 | 6.9 |
| | | plicRDF | 7.725e-5 | 1.721e-5 | 8.262e-3 | - | -1.0e-14 / 1.0e-14 | 1.6e-13 | 64.7 [64.6, 65.2] | 64 | 0.12 | 1 |
| 128, 64 | 1 | swept ring | 2.577e-5 | 2.325e-7 | 8.031e-3 | 1.882e-3 | -8.1e-14 / 2.9e-11 | -2.0e-13 | 94.9 [94.6, 96.4] | 95 | 0.98 | 4.8 |
| | | swept exact | 2.576e-5 | 2.325e-7 | 1.912e-4 | 1.791e-3 | -3.6e-13 / 2.9e-11 | -2.9e-13 | 98.7 [98.4, 99.2] | 99 | 1.02 | 5.0 |
| | | plicRDF | 2.500e-5 | 2.325e-7 | 5.835e-3 | - | -1.0e-14 / 1.0e-14 | 7.8e-14 | 19.9 [19.5, 20.0] | 20 | 0.13 | 1 |

At OpenFOAM Co 1 on the x-line the displacement is exactly h per step. Every
step shifts the volume by one cell, which any upwind scheme does exactly: all
methods reach E_shape 2.3e-7, the discretisation of the initial sphere. These
rows do not discriminate.

### The translated sphere along the x-z diagonal

| N, ranks | OF Co | row | L1(x) | E_shape | L1(n) | L1(alphaf) | min / excess | mass % | wall s | clock s | step s | wall / plicRDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32, 1 | 0.2 | swept ring | 2.009e-3 | 5.729e-4 | 5.356e-2 | 3.175e-2 | -1.7e-14 / 2.5e-13 | 1.2e-13 | 42.4 [42.3, 43.9] | 42 | 0.25 | 7.0 |
| | | swept exact | 1.155e-3 | 2.894e-4 | 2.368e-3 | 1.687e-2 | -3.6e-13 / 2.0e-13 | 1.2e-13 | 43.3 [43.3, 43.4] | 43 | 0.25 | 7.1 |
| | | plicRDF | 1.930e-3 | 5.119e-4 | 4.230e-2 | - | -2.5e-13 / 2.2e-16 | 8.5e-13 | 6.1 [6.0, 6.1] | 6 | 0.02 | 1 |
| 32, 1 | 1 | swept ring | 7.864e-4 | 1.809e-4 | 3.870e-2 | 1.118e-2 | -3.9e-11 / 1.4e-13 | 1.4e-13 | 10.5 [10.4, 11.0] | 10 | 0.24 | 3.4 |
| | | swept exact | 4.187e-4 | 1.182e-5 | 2.392e-3 | 4.846e-3 | -3.0e-14 / 2.0e-13 | -4.4e-14 | 10.7 [10.7, 10.7] | 11 | 0.25 | 3.4 |
| | | plicRDF | 6.469e-3 | 1.295e-3 | 1.246e-1 | - | -1.0e-14 / 9.3e-15 | 1.3e-13 | 3.1 [3.1, 3.2] | 3 | 0.03 | 1 |
| 64, 8 | 0.2 | swept ring | 5.517e-4 | 1.568e-4 | 2.408e-2 | 1.770e-2 | -8.9e-13 / 5.2e-12 | 4.4e-13 | 164.5 [164.1, 166.8] | 164 | 0.54 | 8.3 |
| | | swept exact | 3.285e-4 | 8.579e-5 | 8.397e-4 | 9.891e-3 | -5.9e-13 / 5.3e-12 | 4.3e-13 | 169.1 [169.1, 169.2] | 169 | 0.56 | 8.5 |
| | | plicRDF | 5.127e-4 | 1.436e-4 | 1.922e-2 | - | -4.5e-10 / 4.7e-13 | 5.3e-13 | 19.9 [19.8, 20.2] | 20 | 0.05 | 1 |
| 64, 8 | 1 | swept ring | 1.764e-4 | 3.887e-5 | 1.907e-2 | 5.045e-3 | -1.5e-13 / 2.0e-12 | 3.1e-13 | 36.4 [36.4, 36.7] | 37 | 0.54 | 4.0 |
| | | swept exact | 1.007e-4 | 3.784e-6 | 7.224e-4 | 2.357e-3 | -3.2e-11 / 1.5e-12 | 3.8e-13 | 37.3 [37.2, 37.3] | 37 | 0.56 | 4.1 |
| | | plicRDF | 5.077e-3 | 6.306e-4 | 1.392e-1 | - | **-1.5e-3 / 1.2e-3** | -4.4e-14 | 9.1 [9.0, 9.3] | 9 | 0.08 | 1 |
| 128, 64 | 0.2 | swept ring | 1.283e-4 | 3.689e-5 | 9.731e-3 | 8.820e-3 | -2.8e-11 / 1.2e-11 | 1.2e-12 | 982.5 [979.5, 1010.8] | 982 | 1.51 | 9.9 |
| | | swept exact | 9.017e-5 | 2.443e-5 | 2.419e-4 | 5.713e-3 | -3.7e-11 / 1.2e-11 | 1.5e-12 | 1016.6 [1009.4, 1022.9] | 1011 | 1.57 | 10.2 |
| | | plicRDF | 1.231e-4 | 3.437e-5 | 7.793e-3 | - | -1.2e-7 / 1.7e-13 | 8.1e-12 | 99.5 [98.8, 100.0] | 94 | 0.14 | 1 |
| 128, 64 | 1 | swept ring | 3.549e-5 | 7.184e-6 | 8.292e-3 | 2.178e-3 | -3.6e-11 / 1.3e-11 | 1.5e-12 | 204.3 [204.2, 204.8] | 204 | 1.53 | 5.0 |
| | | swept exact | 2.640e-5 | 8.771e-7 | 1.866e-4 | 1.242e-3 | -1.3e-13 / 1.3e-11 | 1.6e-12 | 211.2 [210.6, 211.2] | 211 | 1.59 | 5.2 |
| | | plicRDF | 2.999e-3 | 3.221e-4 | 1.484e-1 | - | **-6.6e-3 / 4.1e-3** | 4.3e-13 | 40.5 [40.5, 40.6] | 41 | 0.25 | 1 |

OpenFOAM Co 1 breaks plicRDF + isoAdvector on the diagonal in parallel:
- At N = 64 / 128 on 8 / 64 ranks, alpha leaves [0, 1] by up to 1.5e-3 /
  6.6e-3 (min alpha and excess).
- Its L1(n) is 0.14 to 0.15, and its L1(x) is 29 / 85 times the swept ring
  row's.
- Serially (N = 32) it stays bounded, at an L1(x) 8 times ours.
- At OpenFOAM Co 0.2 on the diagonal it undershoots slightly in parallel:
  -4.5e-10 (N = 64) and -1.2e-7 (N = 128).
- Not tuned (the settings are the fork's tests').

### The rotation, one turn

| N, ranks | OF Co | row | L1(x) | E_shape | L1(n) | L1(alphaf) | min / excess | mass % | wall s | clock s | step s | wall / plicRDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32, 1 | 0.2 | swept ring | 3.958e-3 | 1.059e-3 | 5.728e-2 | 5.704e-2 | -9.8e-13 / 1.1e-11 | 2.1e-12 | 438.1 [437.0, 439.1] | 438 | 0.25 | 9.7 |
| | | swept exact | 1.959e-3 | 5.261e-4 | 3.181e-3 | 2.981e-2 | -9.8e-13 / 2.5e-12 | 2.5e-12 | 443.5 [442.1, 445.0] | 444 | 0.25 | 9.8 |
| | | plicRDF | 4.389e-3 | 1.155e-3 | 5.752e-2 | - | -5.1e-19 / 2.2e-16 | 1.2e-12 | 45.2 [45.0, 45.4] | 45 | 0.02 | 1 |
| 32, 1 | 1 | swept ring | 2.629e-3 | 7.071e-4 | 4.784e-2 | 3.823e-2 | -8.7e-13 / 9.5e-13 | 2.0e-12 | 90.5 [90.5, 91.6] | 90 | 0.25 | 8.0 |
| | | swept exact | 1.399e-3 | 3.758e-4 | 2.511e-3 | 2.159e-2 | -4.4e-13 / 3.2e-12 | 1.4e-12 | 91.4 [91.3, 91.7] | 92 | 0.25 | 8.1 |
| | | plicRDF | 3.328e-3 | 8.477e-4 | 5.443e-2 | - | -5.1e-18 / 2.2e-16 | 1.2e-12 | 11.3 [11.3, 11.3] | 11 | 0.03 | 1 |
| 64, 8 | 0.2 | swept ring | 1.057e-3 | 2.956e-4 | 2.514e-2 | 3.220e-2 | -4.4e-14 / 6.0e-10 | 4.7e-12 | 1875.0 [1871.7, 1875.4] | 1875 | 0.52 | 9.8 |
| | | swept exact | 6.274e-4 | 1.718e-4 | 8.168e-4 | 1.929e-2 | -5.1e-14 / 2.3e-10 | 5.0e-12 | 1895.0 [1894.8, 1896.7] | 1895 | 0.52 | 9.9 |
| | | plicRDF | 1.139e-3 | 3.127e-4 | 2.407e-2 | - | -2.6e-19 / 2.2e-16 | 1.1e-12 | 192.2 [189.9, 192.3] | 192 | 0.05 | 1 |
| 64, 8 | 1 | swept ring | 7.093e-4 | 1.974e-4 | 2.250e-2 | 2.191e-2 | -7.8e-13 / 3.5e-11 | 3.3e-12 | 380.0 [379.9, 380.6] | 380 | 0.52 | 9.0 |
| | | swept exact | 4.449e-4 | 1.208e-4 | 7.539e-4 | 1.379e-2 | -5.6e-13 / 1.6e-11 | 3.2e-12 | 383.7 [383.3, 384.4] | 384 | 0.52 | 9.1 |
| | | plicRDF | 8.467e-4 | 2.297e-4 | 2.233e-2 | - | -2.0e-17 / 4.4e-16 | 9.4e-13 | 42.2 [41.3, 44.0] | 42 | 0.05 | 1 |
| 128, 64 | 0.2 | swept ring | 2.829e-4 | 8.042e-5 | 1.127e-2 | 1.830e-2 | -1.4e-13 / 7.5e-10 | 7.5e-12 | 9370.9 [9330.0, 9409.6] | 9370 | 1.21 | 9.8 |
| | | swept exact | 1.874e-4 | 5.219e-5 | 2.555e-4 | 1.192e-2 | -1.2e-12 / 5.1e-10 | 8.5e-12 | 9473.1 [9456.8, 9473.1] | 9473 | 1.22 | 9.9 |
| | | plicRDF | 2.949e-4 | 8.374e-5 | 9.560e-3 | - | -5.0e-13 / 4.4e-16 | 1.8e-12 | 959.4 [958.6, 970.1] | 959 | 0.13 | 1 |
| 128, 64 | 1 | swept ring | 1.891e-4 | 5.342e-5 | 9.880e-3 | 1.227e-2 | -4.9e-11 / 6.5e-11 | 5.9e-12 | 1882.7 [1881.8, 1882.7] | 1882 | 1.21 | 9.2 |
| | | swept exact | 1.303e-4 | 3.620e-5 | 2.392e-4 | 8.362e-3 | -4.9e-11 / 1.3e-10 | 6.4e-12 | 1905.8 [1903.4, 1909.7] | 1906 | 1.22 | 9.3 |
| | | plicRDF | 2.234e-4 | 6.096e-5 | 9.432e-3 | - | -9.8e-13 / 2.9e-13 | 1.9e-12 | 204.3 [201.6, 207.6] | 204 | 0.13 | 1 |

The deflt partition started no job from about 14:00 to 16:15. Every run of
the rotation cells exceeds the 30-minute partition, so they waited. They ran
one job per repetition (each with the three methods interleaved).

### Verdicts

- C-1 FALSIFIED, narrowly and only by the exact-normal control: the rotation
  at N = 128, OpenFOAM Co 1 has excess 1.3e-10 against 1e-10.
  - The ring row there: min alpha -4.9e-11, excess 6.5e-11.
  - Every other case, N and normal at OpenFOAM Co 1: min alpha >= -4.9e-11,
    excess <= 6.5e-11, mass <= 6.4e-12 %.
  - At OpenFOAM Co 0.2 the rotation has excess 6.0e-10 / 7.5e-10 (ring, N =
    64 / 128) and 2.3e-10 / 5.1e-10 (exact).
  - The mass is conserved everywhere.
- C-2 FALSIFIED.
  - Holds on the rotation: the ring row's L1(x) is 0.79 to 0.96 and its
    E_shape 0.83 to 0.96 of plicRDF's.
  - Holds on the diagonal at OpenFOAM Co 1 (0.012 to 0.12; plicRDF breaks).
  - Fails on the diagonal at OpenFOAM Co 0.2: L1(x) 1.04 to 1.08 and E_shape
    1.07 to 1.12 of plicRDF's, all three N.
  - Fails at the edge on the x-line at N = 128: E_shape 1.005 (Co 0.2),
    L1(x) 1.031 (Co 1, the degenerate shift).
  - The exact normal (the control) is below plicRDF everywhere (0.009 to
    0.73) except on the degenerate x-line shift at Co 1 (1.00 to 1.03).
- C-3 FALSIFIED (the expectation was too low):
  - The swept-face advection's total wall clock at the same ranks is 2.6 to
    6.7 times plicRDF's on the x-line, 3.4 to 9.9 on the diagonal, and 8.0 to
    9.8 on the rotation, with the ring normal (expected 1.5 to 4 and 2 to 6).
  - The exact normal is not cheaper: 2.7 to 10.2.
  - Per step the median is 0.17 to 1.53 s against 0.02 to 0.25 s. The ratio
    grows with N on the translations (Co 0.2: 4.5 / 5.8 / 6.7 on x, 7.0 / 8.3
    / 9.9 on the diagonal).
- C-4:
  - HELD serially: plicRDF bounded at both Co on all three cases.
  - In parallel the diagonal is not bounded at either Co (above).
  - The x-line and the rotation (every N and Co, nAlphaBounds 5) are bounded.
    The registered expectation of an unbounded N = 128 rotation on 64 ranks
    (seen with nAlphaBounds 3 at OpenFOAM Co 0.88) is FALSIFIED: with
    nAlphaBounds 5 its min alpha is -5.0e-13 (Co 0.2) and -9.8e-13 (Co 1).

What this demonstrates:
- At OpenFOAM Co 1 the swept-face advection stays bounded and conservative
  where plicRDF + isoAdvector, on the same mesh, ranks and dt, loses
  boundedness on oblique translation.
- On the rotation it is 4 to 21 % more accurate in L1(x) with the ring
  normal; the gain shrinks with N at Co 0.2 (0.90 / 0.93 / 0.96 of plicRDF's
  at N = 32 / 64 / 128).
- It costs 3 to 10 times the wall clock at the same ranks.
- On the diagonal at Co 0.2, with the ring normal, it is 4 to 12 % less
  accurate than plicRDF; there the normal, not the flux, sets the error (the
  exact normal is 0.6 to 0.7 of plicRDF's).

## The cost of the swept-face advection: HPCToolkit profiles (2026-11-04)

Tomislav's item: profile the swept-face advection with HPCToolkit on
Lichtenberg and isolate the serial and the parallel bottlenecks. Nothing was
changed in the fork; the recommendations at the end are not implemented.

### The setup

- HPCToolkit 2025.0.1 (spack, ~/spack/opt/spack/linux-x86_64/hpctoolkit-
  2025.0.1-ky3i4lsw; the 2017.06 copy in ~/Software/hpctoolkit-install lacks
  its libraries). hpcrun samples CPUTIME every 1000 us and the LLC misses
  (libpfm ix86arch::LLC_MISSES, PAPI-less) every 100000 misses, on every
  rank; hpcstruct (cached) and hpcprof in the job. The installation has no
  hpcviewer and no Python bindings, so the databases (format v4: meta.db,
  profile.db) are read by a reader written from the installed FORMATS.md
  (workflow/scripts/davof_prof/hpdb.py) and analysed by hpan.py: exclusive
  time = HPCToolkit's "function" scope (the lexical subtree of a function,
  not crossing a call or an inlined call: an inlined helper is its own
  entry), inclusive = the "execution" scope. The samples are attributed to
  categories by the nearest ancestor context matching a rule (a function
  name, a source-line range of fork 720cfc6, or a load module;
  rules_swept.py); libc's malloc/free count for the caller and are listed
  apart. Scripts: workflow/scripts/davof_prof/ (README there).
- The build: a detached git worktree of the fork at davof-swept HEAD 720cfc6
  (prof/TwoPhaseFlow on the cluster), WM_COMPILER=Gcc (11.5.0),
  WM_COMPILE_OPTION=Opt (c++OPT = -O3, the default), WM_OPTIONS
  linux64GccDPInt32Opt, plus FOAM_EXTRA_CXXFLAGS="-g -fno-omit-frame-pointer"
  (appended to c++FLAGS by wmake/rules/linux64Gcc/c++), its own platforms/.
  One compile line: g++ -std=c++17 -m64 -pthread -DOPENFOAM=2606 -DWM_DP
  -DWM_LABEL_SIZE=32 -Wall ... -O3 -DNoRepository -ftemplate-depth=1000 -g
  -fno-omit-frame-pointer. The default build is untouched; OpenFOAM's own
  libraries stay as installed (-O3, no -g: attributed by function only).
- The cases: the comparison's rotation, OpenFOAM Co 1, prepared from the
  comparison's study cases exactly as davof-cmp-time.sbatch does (ring
  normal, sweptVolumeCorrection scale, the final wisp rule; plicRDF with
  nAlphaBounds 5), on exclusive i01 nodes (Xeon Platinum 9242, 2 x 48
  cores), --cpu-bind=cores. Serial: N = 32 and 64 full turns (346 and 697
  steps), N = 128 20 steps (2M cells). Parallel: N = 64 on 8 ranks and N =
  128 on 64 ranks, full turns (697 and 1401 steps), the comparison's scotch
  decompositions. plicRDF + isoAdvector profiled the same way on every case.
- The window: the full turn rather than 20 steps, because the comparison's
  timeStepTime is NOT stationary in parallel: at 8 ranks the mean over 50
  steps goes 0.34 / 0.49 / 0.73 / 0.82 / 0.82 / 0.70 / 0.45 / 0.35 ... s
  over the turn, at 64 ranks 1.03 to 2.11 s, while serially it is 1.05 to
  1.18 s (N = 64). The variation is the imbalance moving through the
  partitions as the sphere turns (below), not the algorithm; the first 20
  steps (0.31 s at 8 ranks) would under-represent the turn (0.52 s).
- hpcrun's overhead: 0.4 to 1.3 % of the wall clock (serial N = 32: 91.1 vs
  91.7 s; serial N = 64, 20 steps: 37.8 vs 38.2 s; 8 ranks: 0.538 vs 0.544
  s/step; 64 ranks: 1883 vs 1902 s per turn).
- Jobs: 55914082 (build, deflt_short), 55924284 (serial), 55924287 (8
  ranks), 55924288 (64 ranks); 55912929-32 were the first submission
  (cancelled / a spec-parsing failure). CSVs on the cluster:
  csv/hp_categories.csv (run x category), hp_ranks.csv (per rank),
  hp_hotspots.csv (top 25 exclusive per run), hp_timing.csv; the full
  analyses in prof/an/ (flat, lines, tree, callers, categories per rank).

### (a) The serial breakdown per step, N = 64, one turn (697 steps)

Seconds per step and the share of the step; exclusive by category (the
categories partition the step; the inclusive time of a routine is the sum
of its categories). The swept-face advection (ring normal) against plicRDF +
isoAdvector on the same case.

| category (fork routine) | swept s/step | % | plicRDF s/step | % |
|---|---|---|---|---|
| pieces: clip by the cell planes and the PLIC plane (sweptPieceVolume, davofSweptPieces.H 296-394, with the inlined clipByHalfSpace) | 0.4685 | 42.2 | - | - |
| pieces: cap closure and volume (395-496) | 0.0025 | 0.2 | - | - |
| face update: the pre-image fan clipped (sweptPieceArea) | 0.1819 | 16.4 | - | - |
| swept cells: the ring search (vertex ring, bounding boxes; isoAdvectionSwept.C 694-731) | 0.0412 | 3.7 | - | - |
| swept cells: the cell planes, state and wisp rule (774-846) | 0.0302 | 2.7 | - | - |
| swept cells: the pre-image fan, the rest of sweptFaceCellsFlux | 0.0079 | 0.7 | - | - |
| trace: the velocity at the traced points (interpolationCellPoint, cellPointWeight::findTetrahedron) | 0.0742 | 6.7 | 0.0019 | 0.9 |
| trace: the cell holding a point (ring, octree) | 0.0174 | 1.6 | - | - |
| trace: the midpoint trace, the triangle list | 0.0104 | 0.9 | - | - |
| the point velocity (volPointInterpolation of U, once per step) | 0.0222 | 2.0 | 0.0113 | 5.4 |
| curvedFlux: the candidate band | 0.0056 | 0.5 | - | - |
| curvedFlux: the wisp surfaces (neighbour plane, anchoring) | 0.0036 | 0.3 | - | - |
| curvedFlux: the second ring / surfaceOf | 0.0014 | 0.1 | - | - |
| curvedFlux: the face loop and the rest | 0.0106 | 1.0 | - | - |
| the remote exchange (serial: the zoneDistribute stencil setup) | 0.0030 | 0.3 | 0.0113 | 5.4 |
| the surface build (buildCellSurfaces: fit, anchor, kromerBothe) | 0.0451 | 4.1 | - | - |
| the reconstruction (davof::reconstruct + ring normal; plicRDF::reconstruct) | 0.0365 | 3.3 | 0.0684 | 32.7 |
| the bounding (limitFluxes, boundCellContent, second layer) | 0.0102 | 0.9 | 0.0077 | 3.7 |
| isoAdvector's timeIntegratedFlux (the faces no plane set) | 0.0147 | 1.3 | 0.0051 | 2.4 |
| isoAdvector's advect, the rest (alpha update, counters) | 0.0293 | 2.6 | 0.0248 | 11.8 |
| setFlow (the function object) | 0.0309 | 2.8 | 0.0306 | 14.6 |
| the solver loop's field algebra (main) | 0.0355 | 3.2 | 0.0352 | 16.8 |
| I/O and start-up | 0.0243 | 2.2 | 0.0105 | 5.0 |
| the Courant number | 0.0015 | 0.1 | 0.0015 | 0.7 |
| total | 1.1105 | 100 | 0.2090 | 100 |

- The swept advection's own parts (pieces, face update, swept cells, trace,
  point velocity, curvedFlux, surface build) are 0.926 s/step = 83 % of the
  step; the two clips alone (pieces + face update) 0.65 s = 59 %.
- The ratio to plicRDF is 5.3 serially (the comparison measured 9.0 at 8
  ranks: the difference is the parallel imbalance, (c)). Function by
  function the swept row adds: pieces +0.472, face update +0.182, trace
  +0.100, swept cells +0.079, surface build +0.045, curvedFlux +0.021, point
  velocity +0.011, timeIntegratedFlux +0.010 s/step; its reconstruction is
  cheaper than plicRDF's (-0.032; the ring normal against plicRDF's
  iterations); setFlow, the solver loop and the alpha update are the same.
- Allocation is 2.6 % of the step (20.3 s of 774: the trace's tetrahedron
  lists 13.8 s, the solver loop 2.2, the DynamicLists of sweptFaceCellsFlux
  1.2, the surface build 0.9). The LLC misses (12.2e9 over the turn) sit in
  the mesh-wide field operations (setFlow 2.5e9, the point velocity 2.5e9,
  the alpha update 2.3e9, the solver loop 2.0e9, the reconstruction 0.8e9);
  the two clips have 0.8 % of them: the clipping is compute-bound.
- The in-code timers of the comparison logs agree: advect is 97 to 99 % of
  the step, curvedFlux 95 %, the exchange 1.7 to 1.8 % (its wait 0.1 to 0.3
  %).

### (b) The scaling with N, serial

The full turns at N = 32 and 64 (346 and 697 steps); the growth factor
against the interface-cell count (431 -> 1723, x 4.00; swept faces 4500 ->
17574, x 3.91; pieces 6261 -> 23939, x 3.82; cells x 8). The 20-step windows
at N = 64 and 128 (the serial N = 128 turn is not affordable); there the
start-up and the field write of step 20 are excluded.

| category | N = 32 s/step | N = 64 s/step | x (32 -> 64) | N = 64, 20 steps | N = 128, 20 steps | x (64 -> 128) |
|---|---|---|---|---|---|---|
| pieces: clip | 0.1205 | 0.4685 | 3.89 | 0.464 | 1.830 | 3.94 |
| face update | 0.0472 | 0.1819 | 3.85 | 0.182 | 0.691 | 3.80 |
| ring search | 0.0106 | 0.0412 | 3.89 | 0.043 | 0.167 | 3.9 |
| cell planes, wisp rule | 0.0075 | 0.0302 | 4.03 | 0.027 | 0.114 | 4.2 |
| trace: velocity at the points | 0.0191 | 0.0742 | 3.88 | 0.085 | 0.403 | 4.7 |
| trace: cell search | 0.0045 | 0.0174 | 3.87 | 0.022 | 0.097 | 4.4 |
| surface build | 0.0130 | 0.0451 | 3.47 | 0.043 | 0.179 | 4.2 |
| reconstruction | 0.0064 | 0.0365 | 5.70 | 0.028 | 0.148 | 5.3 |
| point velocity (mesh-wide) | 0.0025 | 0.0222 | 8.9 | 0.029 | 0.235 | 8.1 |
| setFlow (mesh-wide) | 0.0034 | 0.0309 | 9.1 | 0.033 | 0.366 | 11 |
| alpha update (mesh-wide) | 0.0028 | 0.0293 | 10.5 | 0.030 | 0.334 | 11 |
| bounding (mesh-wide) | 0.0009 | 0.0102 | 11 | 0.010 | 0.110 | 11 |
| total | 0.2606 | 1.1105 | 4.26 | | ~5.5 (steady) | |

- Every part of the swept advection grows with the interface cells (x 3.5 to
  4.0 from N = 32 to 64; the trace's cell search and the velocity lookup
  somewhat faster, x 4.4 to 4.7, from 64 to 128: the octree and the
  tetrahedron search are logarithmic in the mesh), the reconstruction x 5.3
  to 5.7 (the ring stencil), and the mesh-wide operations shared with
  plicRDF (setFlow, the alpha update, the point velocity, the bounding, the
  Courant number) x 8 to 11. No part of the swept flux grows faster than the
  interface-cell count.
- The ratio to plicRDF therefore FALLS with N serially (9.6 at N = 32, 5.3 at
  N = 64): the mesh-wide part that both share grows faster than the swept
  part. The growth of the ratio with N in the comparison (4.5 / 5.8 / 6.7 on
  the x-line) comes from the ranks (1 / 8 / 64) changing with N, not from
  the algorithm: (c).

### (c) Parallel: the time in MPI against compute, the imbalance

CPUTIME summed over the turn per rank (its helper threads included; a rank
busy-waits in MPI, so CPUTIME is the wall clock). The MPI time is also split
by the call site (the nearest fork routine above the MPI call).

| | N = 64, 8 ranks (697 steps) | N = 128, 64 ranks (1401 steps) |
|---|---|---|
| step (CPUTIME per rank / steps; the wall clock) | 0.545 s (0.551) | 1.341 s (1.358) |
| MPI library, share of all CPUTIME | 72.1 % | 91.2 % |
| MPI per rank: mean / max / min | 274 / 359 / 72 s | 1715 / 1841 / 929 s |
| Pstream wrappers | 0.1 % | 0.1 % |
| compute per rank: max / mean (max/mean) | 303 / 104 s (2.92) | 931 / 156 s (5.98) |
| ranks with compute < 10 % of the max | 2 of 8 | 42 of 64 |
| compute per rank, s | 235 59 303 58 43 27 85 20 | 75 ... 931 (ranks 27, 7, 15, 28, 29, 26, 14, 6: 931 758 692 691 483 419 385 359; 42 ranks at 37-80) |
| partial cells per rank at t = 0 (1e-8 < alpha < 1 - 1e-8) | 532 438 394 396 0 0 0 0; max/mean 2.42 | 55 ranks with 0; max 1746 of 6920 (25 %), mean 108, max/mean 16.2 |
| MPI by call site, share of the MPI time | syncProcPatches 90.8 %, sweptExchange (mapDistribute) 6.4 %, the reduce in davof::reconstruct 1.5 %, the alpha update 0.4 %, the second exchange 0.2 %, the surface exchanges 0.2 % | syncProcPatches 90.1 %, sweptExchange 6.2 %, davof::reconstruct 1.6 %, the alpha update 0.7 %, the second exchange 0.6 %, setFlow 0.2 %, the surface exchanges 0.2 %, the bounding 0.1 % |
| the exchange proper (pack, map, unpack: compute) | 0.0049 s/step/rank | 0.0058 s/step/rank |
| the busiest rank: compute, of which | rank 2: 0.435 s/step: pieces 0.204, face update 0.083, trace velocity 0.031, surface build 0.019, ring 0.018; MPI 0.104 | rank 27: 0.665 s/step: pieces 0.299, face update 0.129, trace velocity 0.046, surface build 0.027, ring 0.027; MPI 0.663 |
| an idle rank | rank 7: compute 0.029 s/step (the solver loop 0.0046, setFlow 0.0039, the alpha update 0.0037, the point velocity 0.0033); MPI 0.514 | rank 61: compute 0.026 s/step (0.0044, 0.0043, 0.0044, 0.0037); MPI 1.314 |
| the balanced step (the sum of compute / ranks / steps) | 0.151 s: the step is 3.6 x | 0.111 s: the step is 12 x |
| cap-closure fallback / candidate band, on the busiest rank | 0.0003 / 0.0008 s/step | 0.0016 / 0.0011 s/step |
| plicRDF on the same case: MPI share; compute max/mean; its interface work of its compute | 47 %; 2.4; - | 65.6 %; 1.94; 27 % (reconstruction 908 s + bounding 196 + timeIntegratedFlux 70 of 4373 s) |

- The parallel bottleneck is the load imbalance of the interface band under
  the comparison's scotch decomposition of the whole mesh, not the
  communication: the exchange proper is 0.5 % of the step; 90 % of the MPI
  time is the all-to-all of syncProcPatches inside curvedFlux, where the
  ranks without interface cells wait for the rank with the most. The swept
  advection's compute is ~ 95 % interface work (plicRDF's 27 %), so the
  interface band's distribution sets its parallel efficiency directly: at 64
  ranks 55 of 64 ranks hold no partial cell and the busiest holds 25 %.
- The busiest rank's compute is only half the step at 64 ranks (0.665 of
  1.34 s) because the band moves through the partitions as the sphere turns:
  eight ranks take turns as the busiest (931 to 359 s over the turn) and the
  step is each step's maximum; hence the 2 x variation of the per-step time
  over the turn in parallel (1.03 to 2.11 s per 50 steps at 64 ranks) against
  a constant serial step.
- The speed-up against the serial step: N = 64, 1.11 -> 0.545 s at 8 ranks =
  2.0 (efficiency 25 %); N = 128, ~5.5 s (the serial 20-step profile without
  the start-up and the write) -> 1.017 s (the unprofiled 20 steps at 64
  ranks) = 5.4 (8.5 %). Amdahl's law would read a serial fraction f = (64/5.4
  - 1)/63 = 0.17 into the 64-rank figure (0.42 into the 8-rank one), but
  there is no serial fraction: the balanced step would be 0.111 s, and the
  18-19 x of the earlier note is not reached by the measured 5.4. (The
  origin of the 18-19 x is not in this report; the measured values stand.)
- Neither the cap-closure fallback nor the candidate band dominates on any
  rank (0.1 to 0.3 % of the busiest rank's compute); on every busy rank the
  piece clip is first and the face-update clip second.
- The remote exchange's wait (6 % of the MPI time, 0.076 s/step/rank at 64
  ranks) is also an imbalance wait: the idle ranks reach the mapDistribute of
  sweptExchange before the busy ones leave the surface build.

### (d) The top ten exclusive hotspots (serial, N = 64, one turn; 774 s)

| # | % | s | function (file:line; fork 720cfc6) | kind |
|---|---|---|---|---|
| 1 | 15.9 | 123.1 | Foam::minusOp (ops.H:209): the vector subtraction inlined into clipByHalfSpace (davofSweptPieces.H:184, 205) and the contact test (341) | algorithmic (vertices x planes) |
| 2 | 10.6 | 82.1 | davof::clipByHalfSpace (davofSweptPieces.H:181-212; own lines 195, 185, 196, 192, 181, 205); 64.0 s under sweptPieceVolume, 18.1 s under sweptPieceArea | algorithmic |
| 3 | 7.2 | 56.0 | Foam::Vector::inner (VectorI.H:157): the signed distances of the clip (184) and of the contact / coplanar tests (341, 545) | algorithmic; half redundant: the contact test (336-352) computes (t_j - x0).n for every vertex and plane and clipByHalfSpace recomputes it at 184 |
| 4 | 5.8 | 45.2 | davof::sweptPieceVolume, own lines (362, 343, 336, 346, 371, 360, 296): the contact rule's loop over the planes and the per-triangle bookkeeping | algorithmic |
| 5 | 4.3 | 33.6 | std::__copy_move (stl_algobase.h:430): the FixedList<point, 16> copy poly.p = out, poly.e = lab at the end of every clip (davofSweptPieces.H:210-211) | redundant (a 384-byte copy per clip; ping-pong buffers or an in-place clip remove it) |
| 6 | 2.9 | 22.6 | Foam::Vector::cross (VectorI.H:169): areaVector (154) and the triangle normals | algorithmic |
| 7 | 2.6 (5.7 with its callees) | 20.3 (44) | cellPointWeight::findTetrahedron (libfiniteVolume) + polyMeshTetDecomposition::cellTetIndices 4.8 s + faceTetIndices 5.4 s + the allocation of their lists 13.8 s: the velocity at a traced point (sweptTriangles' lambda, isoAdvectionSwept.C 329-344) | redundant + allocation: the cell's tetrahedra are rebuilt and the list allocated per traced point, and a mesh point is traced by every face that uses it |
| 8 | 2.4 (5.6 with magSqr, mag) | 18.8 (43.6) | plusEqMagSqrOp2 (ops.H:79) + magSqr (scalarImpl.H:339) 14.8 s + mag (doubleScalar.H:78) 10.0 s: the mag() of the triangle edges (300) and of sweptPieceArea's coplanar test (541-546: per plane and vertex) | redundant in part: the polygon size L is recomputed for every plane (541-546), the signed distances again at 557 |
| 9 | 2.3 | 17.5 | Foam::boundBox::boundBox(points, cellPoints) (libOpenFOAM; the touches lambda, isoAdvectionSwept.C 705-709): a cell bounding box per ring candidate per face, every step | redundant (a static mesh: once) |
| 10 | 1.8 | 13.7 | volPointInterpolation::interpolateInternalField (libfiniteVolume): the point velocity of the whole mesh, once per step | algorithmic, mesh-wide (x 8.1 from N = 64 to 128); could be restricted to the band's points |

Then: sweptPieceArea's own lines 12.5 s (1.6 %), sweptFaceCellsFlux's own
lines 10.5 s (1.4 %), operator*(scalar, vector) 10.1 s, FixedList::operator[]
8.4 s, kbIntegrateFace 7.4 s (the surface build), fvc::surfaceIntegrate 6.4
s, malloc + free 11.6 s (1.5 %). At 8 and 64 ranks the same functions lead
the compute (the clip's lines, the tetrahedron search), behind the MPI
progress functions (mca_btl_vader_component_progress 51 % / 74 % of all
samples, opal_progress 10 % / 8 %: the busy-wait).

### The changes with the largest expected gain (not implemented)

1. An exact pre-test of the ring cells before the two clips (serial; the
   largest). Evidence: the ring of a face holds 18 to 27 cells (the donor's
   vertex ring whose bounding boxes overlap the swept volume's: the box of a
   face's swept volume touches 3 x 3 cell columns and 2 to 3 layers, and
   boundBox::overlaps is inclusive), while the pieces per face are 1.36
   (23939 pieces / 17574 faces per step at N = 64; the same at N = 128). So
   ~ 93 % of the sweptPieceVolume and sweptPieceArea calls clip 16 (and 4)
   triangles by 6 to 7 planes and find nothing; the two clips are 0.65 of
   the 1.11 s serial step (59 %) and 0.43 of the busiest rank's 0.665 s at 64
   ranks. A separating-plane test of the swept volume's vertices (the face's
   points, their pre-images and the two centres) against the cell's planes
   (all outside one plane -> no piece) costs ~ 60 dot products per candidate
   instead of up to 16 x 7 clips. Expected: the clips drop to ~ 0.15 s/step
   (the true pieces and the near misses) -> the serial step from 1.11 to ~
   0.6 s (factor ~ 1.8), and the busiest rank's compute in parallel alike
   (the 64-rank step from 1.34 to ~ 0.9 s).
2. Balance the interface band across the ranks (parallel; the dominant
   loss). Evidence: (c): 91 % of all CPU time at 64 ranks is the wait at
   syncProcPatches, the busiest rank holds 25 % of the partial cells and 55
   ranks none, the balanced step would be 0.111 s against 1.34 (12 x), at 8
   ranks 0.151 against 0.545 (3.6 x); the exchange proper is 0.5 % of the
   step. A decomposition weighted by the interface band (decomposePar with
   cell weights from the alpha band) addresses the start; since the band
   moves through the partitions over the turn (the per-step time varies 2 x),
   a periodic re-decomposition (redistributePar) or a dynamic weighting is
   needed to hold the balance. Expected factor 3 to 6 at 64 ranks and 2 to 3
   at 8 ranks (the bound 12 and 3.6; the band's width against the rank
   count and its motion limit it). With 1 and 2 together the 64-rank step
   would be ~ 0.2 to 0.3 s against plicRDF's 0.14: a ratio of ~ 2 instead of
   9.
3. Trace each mesh point once per step (serial; redundant + allocation).
   Evidence: hotspot 7: the velocity at a traced point rebuilds the cell's
   tetrahedra (an allocated list) per lookup, and a mesh point is traced by
   every face that uses it; the trace's velocity lookup and cell search are
   0.092 s/step at N = 64 (8.3 %) and grow x 4.5 to N = 128. A per-step cache
   of the traced points keyed by the mesh point (the midpoint trace depends
   on the point only, not on the face) and of the face centres removes the
   repeats, and the tet decomposition can be reused per cell. Expected ~ 6 %
   of the serial step (factor 1.06), ~ 7 % of the busiest rank's compute.
4. Within the clip (micro, after 1 it matters little): the polygon copy at
   the end of clipByHalfSpace (4.3 %), the per-plane recomputation of L and
   of the signed distances in the contact and coplanar tests (~ 4 % of the
   step: hotspots 3 and 8), the cell bounding boxes once (hotspot 9, 2.3 %).
   Together ~ 10 % of the serial step today; ~ 1 % after 1.

What this demonstrates: serially the swept-face advection's cost is the
clipping of ~ 20 ring cells per face of which ~ 1.4 hold a piece (59 % of
the step in the two clips), every part of it scales with the interface-cell
count, allocation is 2.6 % and the work is compute-bound; in parallel its
cost is set by the rank holding the most interface at each step, which under
the mesh-wide decomposition leaves 42 of 64 ranks idle and makes the
comparison's 9 x against plicRDF (5.3 x serially) a decomposition effect.

## The volume consistency of the PLIC plane positioning (2026-11-04)

The question (Tomislav, 2026-11-04): is surfaceIteratorPLIC::vofCutCell, as
davof::positionPlane uses it, volume-consistent, i.e. after the positioning,
does clipping every interface cell by its plane give back alpha_c V_c, and
does the sum over the interface cells close to machine epsilon times their
number? Measured on the sphere.

Pre-registered (the header of the app, fork fdd649d, before the runs):
- E1: per cell |V_cut - alpha_c V_c| < 1e-14 V_c.
- E2: |sum over the interface cells of (V_cut - alpha_c V_c)| < N_interface
  x 1e-15 x V_c.
- E3: every status 0 (converged) except the cells with alpha outside
  [tol, 1 - tol].
- The agent's reading of the code before the runs, also in the header: the
  iteration stops at |alpha(p) - alpha| <= isoFaceTol = 1e-8, so residuals
  between 1e-14 and 1e-8 are admitted by construction; E1 and E2 hold only
  if the first Newton exit sits at round-off in every cell.
- Before the tolerance sweep (commit 9394337): at 1e-14 E1 HELD everywhere
  with no secant-bad cell; at 1e-15 some cells reach the secant stage and a
  few exhaust it.

Commits (fork TwoPhaseFlow, branch davof-swept, remote tmaric):
- fdd649d: apps/benchmark/planeVolumeConsistency (the app, Allwmake entry).
- 9394337: the cap-area check with the sign of faceArea_ (into the liquid);
  the worst-cell trace with the normal actually used.

Jobs (ledger .my_jobs): 55926460 build; 55926498 / 55926499 / 55926500 the
first run (N = 32 serial, 64 on 8 ranks, 128 on 64 ranks, positionTol 1e-8);
55926744 build; 55926769 / 55926770 / 55926771 the sweep (N = 32 and 64 at
1e-8, 1e-10, 1e-12, 1e-14, 1e-15; N = 128 at 1e-8 and 1e-14). Each
invocation took 5 to 12 s wall at N = 128 on 64 ranks.

Raw data: /work/scratch/tm83tomy/davof/agents/swept/vol-consistency/ (the
first run: <case>_<t0|end>.dat and .log) and vol-consistency/run2/ (the
sweep: <case>_<t0|end>_tol<tol>.dat and .log; the .log holds the ten worst
cells per rank with the branch and the iterates, every non-converged cell
with its full iterate history, and every status-0 cell whose residual is
above the tolerance).

### The code (read 2026-11-04, fork 720cfc6)

surfaceIteratorPLIC::vofCutCell(celli, alpha, tol, maxIter = 100, normal):
- The normal is normalised in place; a zero normal returns sign(alpha -
  0.5) without a cut (a warning).
- The vertex values f_i = (p_i - x_c).n are sorted; a1 = 1 is assumed at the
  lowest level, a2 = 0 at the highest (alpha(f) decreases in f).
- alpha < tol returns -1, 1 - alpha < tol returns +1, before any cut: the
  cut object keeps the previous cell's cut, which davof::positionPlane
  handles by reading the RETURNED status (the phantom-plane fix of
  2026-10-13).
- Bisection over the sorted vertex levels (a cut per step) until two
  adjacent levels bracket alpha; an exact hit returns that cut.
- Adjacent levels closer than 1e-14: the cut at the lower level is
  returned. Bracket flatter than tol in alpha: the midpoint cut is returned.
- Otherwise two more cuts at the thirds, the cubic through the four
  (Vandermonde, LU), Newton on the cubic from f = 1/3 until the model
  residual is <= tol (at most 1000 steps), one cut at the root, and if its
  residual is <= tol the cut is returned. Between adjacent vertex levels
  alpha(f) is exactly a cubic (its derivative, the section area, is
  quadratic), so the model is exact up to round-off.
- Else the secant from that point and a second point x1 = clamp(max(1e-3
  (f2 - f1), 1e-13), f1, f2): an absolute offset near the cell centre, not
  a perturbation of the Newton point. The loop runs while the residual is
  above tol, fewer than maxIter steps, and g1 != g2. Whatever ends it, the
  status of the LAST cut is returned, so a loop that stalls (g1 == g2) or
  exhausts maxIter returns a plane with a residual above tol as status 0.
- The returned plane is in every branch the last cut; cutValue(),
  VolumeOfFluid(), surfaceArea() and surfaceCentre() belong to it.
- The liquid side: cutFacePLIC counts a vertex as liquid for (p - base).n >
  0 after snapping |.| < 1e-15 to 0; status -1 is the fully liquid cell,
  +1 the fully gas cell. cutCell::calcGeomDataCutFace turns faceArea_ into
  the sub-cell (into the liquid). cutCell::calcCellData sums pyramids from
  the mean of the sub-cell's face centres with |.| per pyramid, exact for
  the convex sub-cell of a convex cell.
- Production: isoFaceTol 1e-8 (the dictionary key, the fvSolution of every
  comparison run), surfCellTol 1e-8, maxIter 100.

The swept pieces: a SECOND positioning. isoAdvection::buildCellSurfaces
(isoAdvectionCurved.C) takes vofCutCell's plane only as the start, S.o =
the cut polygon's centroid and the used normal, and with fitAnchor volume
(every comparison run: clippingSurface plane, planeSource anchor, fitDof
none, anchorMethod ccs, anchorTolRel 1e-11) re-solves the offset c so that
the Kromer-Bothe volume under the plane equals alpha_c V_c, to max(1e-11
V, 1e-14 V) in the volume and 1e-11 h in the offset (the bracket closed, the
exit codes 2 to 9), the full and empty cells placed at the edge of the flat
interval. planeOf in isoAdvectionSwept.C then takes (S.o + S.c e3, e3) of
that anchored surface as the plane the pieces are clipped by. So the
1e-8 tolerance of vofCutCell enters the swept flux only through the anchor's
start; the clipped plane is the anchored one.

### The app

apps/benchmark/planeVolumeConsistency (fork 9394337): for a time directory's
alpha and alphaf, the interface cells as davof::reconstruct selects them
(m_c = sum_f alpha_f S_f^out, the all-equal face rule, |m_c| > wispTol
V^(2/3)); per cell the production call vofCutCell(celli, alpha,
positionTol, 100, n) with n one of areaNormal (m_c), exact (the sphere's
normal at the cell centre, re-positioned with the normal at the first
plane's centroid, as the exact model) or field (the written
interfaceNormal, the normal the run used); the returned plane clipped by a
fresh cutCellPLIC (V_plic) and by the divergence-theorem volume of the cut
polyhedron, V_ind = (1/3)[sum_F x_F . S_F + base . S_cap] over the wet
polygons of the cell faces and the closing cap (independent of cutFacePLIC
and cutCellPLIC); d_c = (V_ind - alpha_c V_c)/V_c, its histogram, the sums;
a verbatim copy of vofCutCell with the branch and the iterates recorded,
checked bitwise against the production call; the orientation (a flip is
V_ind closer to (1 - alpha) V than to alpha V; on the sphere the liquid
sub-cell centroid must be nearer the centre than the gas one). Serial or
parallel; -out writes the table outside the case.

### Results at the production tolerance 1e-8

30 files: the exact initial sphere (t = 0, initDavofFields exactSphere) at
N = 32 / 64 / 128 on the rotation case (centre (0.5 0.5 0.5)) and the
diagonal case (centre (0.25 0.5 0.25)); the transported spheres at the end of
the comparison runs of 2026-11-03 (rotation one turn and the diagonal line,
OpenFOAM Co 0.2 and 1, the ring row and the exact row, repetition 1); the
normals areaNormal and exact at t = 0, areaNormal, exact and field at the
end. V_c = 3.05e-5 / 3.81e-6 / 4.77e-7, h = 1/32 / 1/64 / 1/128; the liquid
volume 0.0141372 everywhere (mass conserved).

Checks that held in every file:
- The traced copy returned the production cutValue bitwise (0 mismatches of
  117 000 positionings).
- Every inside cell converged in the Newton branch: 0 secant, 0 stalled, 0
  flat brackets, 0 coincident levels, 0 bisection hits, 0 planes outside the
  vertex range; at most 8 to 11 Newton steps on the cubic.
- The returned plane is the last cut: max |V_plic/V - VolumeOfFluid()| = 0.
- The two clips agree: max |V_plic - V_ind|/V_c 3.6e-15 (N = 32) to 4.0e-14
  (N = 128); max |S_cap + A_plic|/h^2 1.2e-15 to 6.6e-15.
- The orientation: 0 flips and 0 sphere-side violations (5232 cells tested
  at N = 128).
- The outside cells (alpha < 1e-8 or > 1 - 1e-8 while |m_c| passes the wisp
  threshold; 0 to 13 per file) all returned the designed -1 / +1; their
  alpha residual is at most 1.2e-9.
- E3 HELD in all 30 files.

The per-cell mismatch, areaNormal (the exact and field normals give the same
numbers to within 10 %):

| case, time | N_int (inside) | max \|d\| | mean \|d\| | cells >= 1e-14 | sum signed | sum abs | sum abs / (N eps V) |
|---|---|---|---|---|---|---|---|
| sphere t = 0, N 32 (both centres) | 416 | 2.0e-10 | 1.4e-11 | 344 | +1.3e-13 | 1.8e-13 | 6.2e4 |
| sphere t = 0, N 64 | 1760 | 8.4e-9 | 7.0e-10 | 1328 | -2.2e-13 | 4.7e-12 | 3.2e6 |
| sphere t = 0, N 128 | 6680 | 9.7e-9 | 6.4e-10 | 4915 | -2.2e-13 | 2.0e-12 | 2.9e6 |
| rotation end, Co 0.2, ring, N 32 | 421 (415) | 9.7e-9 | 6.7e-10 | 291 | +1.2e-12 | 8.5e-12 | 3.0e6 |
| rotation end, Co 1, ring, N 32 | 421 (416) | 9.8e-9 | 8.6e-10 | 289 | +3.6e-12 | 1.1e-11 | 3.8e6 |
| diagonal end, Co 0.2, ring, N 32 | 414 (412) | 9.7e-9 | 1.1e-9 | 290 | +5.3e-12 | 1.3e-11 | 4.7e6 |
| diagonal end, Co 1, ring, N 32 | 416 | 9.8e-9 | 5.9e-10 | 264 | -1.0e-12 | 7.5e-12 | 2.7e6 |
| rotation end, Co 0.2, ring, N 64 | 1745 (1737) | 1.0e-8 | 7.6e-10 | 1260 | +7.8e-13 | 5.0e-12 | 3.4e6 |
| rotation end, Co 1, ring, N 64 | 1758 (1753) | 9.8e-9 | 7.1e-10 | 1276 | +9.8e-13 | 4.7e-12 | 3.2e6 |
| diagonal end, Co 0.2, ring, N 64 | 1752 (1742) | 1.0e-8 | 8.8e-10 | 1356 | +1.5e-13 | 5.8e-12 | 3.9e6 |
| diagonal end, Co 1, ring, N 64 | 1760 | 8.9e-9 | 7.6e-10 | 1291 | -9.8e-13 | 5.1e-12 | 3.4e6 |
| rotation end, Co 0.2, ring, N 128 | 6921 (6908) | 1.0e-8 | 7.2e-10 | 5333 | +3.3e-13 | 2.4e-12 | 3.2e6 |
| rotation end, Co 1, ring, N 128 | 6773 (6766) | 1.0e-8 | 6.6e-10 | 5229 | +3.9e-14 | 2.1e-12 | 3.0e6 |
| diagonal end, Co 0.2, ring, N 128 | 6800 | 1.0e-8 | 7.3e-10 | 5407 | -8.1e-14 | 2.4e-12 | 3.3e6 |
| diagonal end, Co 1, ring, N 128 | 6736 | 1.0e-8 | 8.1e-10 | 5473 | -4.4e-13 | 2.6e-12 | 3.7e6 |

(d = (V_ind - alpha V)/V; the sums in volume units; eps = 2.22e-16; V the
mean interface-cell volume. The exact rows at the end: the same ranges.)

The histogram of |d| at the end of the rotation, N = 64, Co 1, ring,
areaNormal (1753 cells; the bins < 1e-16, 1e-16..1e-15, ..., 1e-9..1e-8, >=
1e-8): 6, 91, 380, 197, 177, 169, 207, 243, 283, 0. The residual is spread
over every decade between 1e-16 and 1e-8: Newton on the cubic stops at the
first iterate below 1e-8, wherever its quadratic sequence lands. No cell
exceeds the tolerance.

Verdicts at 1e-8:
- E1 FALSIFIED in all 30 files: max |d| 9.7e-9 to 1.0e-8 (2.0e-10 on the
  exact N = 32 sphere, whose cuts repeat by symmetry); 60 to 82 % of the
  cells are at or above 1e-14.
- E2 FALSIFIED in all 30 files: |sum signed| 3.9e-14 to 5.3e-12 against the
  bound 3.2e-18 to 1.3e-17 (10^4 to 10^6 times); relative to the liquid
  volume 3e-12 to 4e-10. The absolute sum is 2e-12 to 1.3e-11, 6e4 to 5e6
  times N eps V.
- E3 HELD.
- The orientation convention HELD: alpha is the volume on the side (x -
  base).n > 0 of the plane whose normal is m_c, the exact inward normal or
  the written normal; no cell shape flipped it.

What this demonstrates: at its production tolerance the positioning is
exactly as consistent as that tolerance, |V_cut - alpha_c V_c| <= 1e-8 V_c
in every cell with the whole range used, and the per-cell errors cancel to
a signed sum of order 1e-11 of the liquid volume. There is no wrong plane,
no silently non-converged plane, no bracket failure and no orientation
defect at 1e-8. Machine-epsilon consistency is not what the loop asks for.

### The opt-in: the tolerance sweep (isoFaceTol is the key; the default stays 1e-8)

The sweep's pre-registration: 1e-14 HELD everywhere without secant-bad
cells; 1e-15 with some secant cells, a few exhausted.

| tol | N | max \|d\| over the files | cells in the secant | secant-bad (stalled / exhausted) | \|sum signed\| | E1 | E2 | E3 |
|---|---|---|---|---|---|---|---|---|
| 1e-10 | 32, 64 | 3.3e-11 to 1.0e-10 | 0 | 0 | 4e-17 to 6e-14 | F | F | H |
| 1e-12 | 32, 64 | 1.6e-13 to 1.0e-12 | 0 to 1 | 0 | 1e-18 to 7e-16 | F | F (H in 1) | H |
| 1e-14 | 32 | 3.4e-15 to 1.25e-14 | 0 to 6 | 0 | 2e-20 to 9e-18 | H in 12 of 20, F by 1 cell (1.0 to 1.25e-14) in 8 | H | H |
| 1e-14 | 64 | 1.0e-14 to 2.0e-14 | 3 to 65 | 0 to 2 | 5e-20 to 2e-18 | F (2 to 42 cells, by at most 2x) | H | H in 10 of 20, F in 10 (1 to 2 cells) |
| 1e-14 | 128 | 2.2e-14 to 3.1e-14, with outliers 1.3e-13, 2.5e-13, 6.2e-13, 1.0e-12, 2.7e-12 | 707 to 1218 (11 to 18 %) | 70 to 209 | 3e-19 to 2e-18 | F (526 to 1039 cells) | H | F in all 15 |
| 1e-15 | 32 | 4.1e-15 to 1.9e-14 | 75 to 119 (20 to 28 %) | 20 to 123 | 7e-19 to 2e-18 | H in 18 of 20 | H | F in all 20 |
| 1e-15 | 64 | 9.8e-15 to 7.2e-14 | 248 to 471 (20 to 27 %) | 400 to 818 (23 to 46 %) | 2e-19 to 1.5e-18 | F (1 to 18 cells) | H | F in all 20 |

The mean |d| at 1e-14 is 1.2e-15 (N = 32), 2.0e-15 to 3.0e-15 (N = 64) and
4.1e-15 to 4.4e-15 (N = 128): the round-off floor of the cut volume relative
to V_c grows with N (the cells are smaller, the cut's intermediate
quantities are not), and at N = 128 it sits above 1e-14 in 8 to 15 % of the
cells. The absolute sum at 1e-14 is 3 to 7 times N eps V at N = 32, 7 to 13
at N = 64 and 19 to 20 at N = 128; the signed sum 0.3 to 2 times N eps V.

Effect of the opt-in on the sums: tightening isoFaceTol from 1e-8 to 1e-14
takes the signed sum from 1e-13..5e-12 to 1e-19..2e-18 (E2 HELD in all 55
files at 1e-14) and the absolute sum from 2e-12..1.3e-11 to 1e-17..2e-17,
at the price of 0 to 18 % of the cells passing through the secant stage.

### The defect found (latent at 1e-8, exposed below the round-off floor)

The secant stage returns its LAST iterate with status 0 when it stalls (g1
== g2, the two latest cuts equal to the bit) or exhausts maxIter. Below the
floor this happens in 1.2 % of the cells at N = 128 (tol 1e-14) and in 23 to
46 % at N = 64 (tol 1e-15); across the sweep 22 986 stalled, 16 exhausted, 1
coincident-level cut above tol. Most such planes sit at the floor (the
residual 1e-15 to 3e-14), but the last iterate is not the best one: the
worst returned planes are off by 2.67e-12 V (LineDiag N = 128, Co 0.2, exact
row, exact normal), 1.01e-12 V (rotation N = 128, Co 0.2, ring row, m_c; the
same file has 3 inside cells returned as full/empty: alpha 1.3e-14 to 5e-15,
m_c of 4e-7 to 4e-5 h^2, the secant jumped to the vertex level f = +h/2
twice and stalled there) and 6.1e-13 V, i.e. 60 to 270 times the requested
tolerance, reported as converged. At the production tolerance 1e-8 this
branch is never entered (0 of 117 000 positionings).

Contributing cause: the secant's second point x1 = clamp(1e-3 (f2 - f1)) is
an absolute cut value near the cell centre, not a perturbation of the Newton
root, so the first secant step leaves the neighbourhood of the root; the
stalls come after 2 to 7 steps.

Proposed fix, NOT applied to production: (i) keep the best iterate (the cut
with the smallest |residual| among those evaluated) and return it when the
secant loop ends without reaching tol; (ii) take the second secant point as
the Newton root shifted by a small fraction of the bracket (or use the
bracket ends f1, f2 with a regula-falsi step, which always exists); (iii)
return a distinct code, or at least count, when the loop ends above tol, so
that davof::positionPlane can see it. Each is opt-in by a new key; the
effect to show is the max |d| at tol 1e-14 on N = 128 falling from 2.7e-12
to the floor (~ 3e-14) with the same cell count in the secant. Not done in
this task: the production default is unchanged, and isoFaceTol 1e-8 never
reaches the branch.

What this demonstrates: the positioning is volume-consistent to its
tolerance, which is 1e-8 and not machine epsilon; the dictionary key
isoFaceTol takes it to the round-off floor (1e-14 at N <= 64) with the
Newton branch alone; below the floor the secant's exit returns unconverged
planes as converged, up to 270 times the tolerance off. The swept flux clips
the anchored plane (anchorTolRel 1e-11 in volume, 1e-11 h in offset), not
vofCutCell's.

## The Consecutive Cubic Spline positioning in the reconstruction (2026-11-04)

The question (Tomislav, via the coordinating session, after the section
above): position the plane to machine epsilon with Consecutive Cubic Spline
positioning from the closed-form explicit position, which should converge
within one iteration, so that the flux integration clips properly
positioned planes; configurable and replaceable with the previous algorithm.

Pre-registered before the runs:
- C1..C5 in the header of apps/benchmark/planeVolumeConsistency (fork
  de416fa): C1 per cell |V_cut - alpha V| < 1e-14 V; C2 |sum| of the order
  N eps V; C3 one to two cuts per cell [the agent: two to three, one on a
  planar state]; C4 the end-of-interval cells get a plane; C5 the cost per
  positioning not above vofCutCell's.
- P-1, P-2 in config/davof/planeTranslation3D_hexSweptWispM9Ccs.yaml:
  identical or better errors; the wall clock not worse than 5 %.
- The anchor's start with the plane at round-off: 0 to 1 iterations (the
  task).

Commits:
- Fork davof-swept (remote tmaric): de416fa (the key positioning
  vofCutCell | ccs with positioningTol and positioningMaxIter; the shared
  kernel davofCcsPositioning.H; the anchor's start statistics and the key
  anchorSkipConverged; the app's -positioning ccs), 264778b (the CCS step
  proper, the target clamped to [0, V]), b6a7dc9 (the floor exit only within
  a hundred tolerances; positioningMaxIter 30). The isoAdvection edits of
  the first commit travelled in 3aaf49d of the other session, which commits
  in the same worktree (see "Surprises").
- leia davof-swept (origin): f12875ee (the tokens TPF_POSITIONING,
  TPF_POSITIONING_TOL, TPF_ANCHOR_SKIP_CONVERGED, defaults = before; the
  gate config), 7519fae0 and a605efec (the gate config of the later
  builds), d809acc0 (the tokens TPF_SWEPT_PRE_TEST none and
  TPF_SWEPT_STATE_SHORTCUT false, see "Surprises").

Jobs (ledger .my_jobs): builds 55928069 (de416fa), 55929209 (rc 126: the
worktree's Allwmake not executable) and 55929401 (264778b), 55929937
(b6a7dc9, in the worktrees TwoPhaseFlow-v2 and -v3 so that the running jobs
kept their library); the verification 55928070-72, 55929402-04,
55929938-40; the anchor pairs 55929615, 55929941; the plane gate 55928074
(de416fa) and 55929942 (b6a7dc9); the timed comparison 55928075-105
(de416fa, the fork defaults of 3aaf49d active) and 55929943-67 (b6a7dc9,
pinned); the pre-test check 55929616. Raw data: vol-consistency/run3,
run4, run5 (the verification of the three builds), vol-consistency/
anchor-skip22 and anchor-skip3, vol-consistency/pretest-check, cmp-pos and
cmp-pos3, leia/studies/davof/planeTranslation3D_hexSweptWispM9Ccs and
..M9Ccs2; the tables csv/pos_table.csv, csv/pos_table_pos3.csv,
csv/gate_pos_Ccs.csv, csv/gate_pos_Ccs2.csv (vol-consistency/pos_an.py,
gate_an.py).

### What was built

The reconstruction key positioning (default vofCutCell: the path of
720cfc6 bit for bit; the face moment is now formed before the first-pass
positioning, which does not read it). With ccs, davof::positionPlane
positions the plane along the selected normal by davof::ccsPositionPlane
(davofCcsPositioning.H): the start p0 = pPlane_c, the Gauss-theorem
position; the volume and the slope dV/dp = -A (A the cut polygon's area)
from cutCellPLIC::calcSubCell, the same cut as vofCutCell and the
diagnostics; the bracket from the vertex extrema (V = V_c and 0, A = 0,
no cut); the first step Newton's from the single evaluated point; every
later step the root of the cubic Hermite interpolant through the two
consecutive evaluated points with their slopes (T. Maric, J. Comput.
Phys.: X 11 (2021) 100093, doi:10.1016/j.jcpx.2021.100093,
arXiv:2012.02227; exact once both lie between the same vertex levels as
the root, since V(p) is a cubic there), through the bracket ends once both
are evaluated, else through the previous and the latest point,
extrapolated, taken where the root lies inside the bracket; bisection
where an iterate leaves the bracket or the bracket was not halved over two
steps; the stop at |V - alpha V| <= positioningTol V (default 1e-14, the
floor measured above), at a non-decreasing residual at a cubic step within
a hundred tolerances (the floor; the best cut is returned), at a closed
bracket, or at positioningMaxIter. alpha outside [0, 1] is clamped; alpha
<= 0 or >= 1 takes the vertex extreme (no plane, as vofCutCell); the cells
with alpha within isoFaceTol of 0 or 1 are positioned. normal_ and
centre_ come from the final cut as before; the kernel hermiteRoot01 is the
anchor's, moved unchanged and shared. The Info line reports the cuts per
cell (histogram), the residuals, the floor and maxIter exits and the
end-of-interval cells; the per-cell signed residual is kept for the
anchor.

The anchor: the start residual's maximum and mean, the cells within the
tolerance, and |V_KromerBothe - V_cut|/V for the same plane (positioning
ccs) in the curved report; the key anchorSkipConverged (default false)
keeps a plane ccs left within the anchor's tolerance without the start
evaluation.

The three builds: de416fa took the cubic through the bracket ends only,
which the one-sided Newton iterates rarely completed (3.8 cuts per cell,
13 000 Newton against 3 300 cubic steps at N = 128) and ran a round-off
wisp of alpha -3e-15 to maxIter (1.1e-12 V off); 264778b the CCS step
proper and the clamp; b6a7dc9 the floor exit armed only near the
tolerance, after the sphere's exact-normal row had returned planes up to
3.7e-3 V off from a non-decrease far from the root (the anchor's start
residual mean 1e-4 against 1e-14 on the ring row). The verification
statistics of 264778b and b6a7dc9 are identical on the 30 files.

### The verification (planeVolumeConsistency -positioning ccs, b6a7dc9, jobs 55929938-40)

The same 30 files as the section above (the exact initial sphere and the
transported ends, N = 32 / 64 / 128, the normals areaNormal, exact,
field), each positioned by vofCutCell (isoFaceTol 1e-8) and by ccs in the
same job on the same node.

| N | max \|d\| (ccs) | mean \|d\| | cuts per cell, mean (histogram 2 / 3 / 4+) | floor exits | \|sum signed\| / (N eps V) | ns per positioning ccs / vofCutCell |
|---|---|---|---|---|---|---|
| 32 | 4.4e-15 to 1.2e-14 | 1.0e-15 to 1.6e-15 | 2.81 to 2.98 (t = 0: 72 / 320 / 24 of 416) | 0 | 0.1 to 1.8 | 1659 to 1959 / 3221 to 4061 |
| 64 | 8.9e-15 to 1.8e-14 | 2.1e-15 to 2.9e-15 | 2.84 to 2.97 | 0 to 7 of 1760 | 0.1 to 1.3 | 1669 to 2037 / 3566 to 4300 |
| 128 | 2.5e-14 to 4.0e-14 | 4.1e-15 to 5.6e-15 | 2.91 to 3.22 (777 / 4567 / 1445 of 6789) | 180 to 420 of 6800 (3 to 6 %) | 0.3 to 2.9 | 2039 to 2578 / 3677 to 4248 |

(d = (V_ind - alpha V)/V by the independent clip; eps = 2.22e-16; the
means over the 84 file-normal pairs: ccs 1991 ns, vofCutCell 3779 ns.)
The histogram's mode is 3 cuts: the start, Newton's step, the consecutive
cubic, which lands at the floor; 2 cuts where V(p) is linear (an
axis-aligned normal); 4 or more where a vertex level lies between the
iterates. The closed-form start is 1 to 8 % of V off on the sphere (the
realizability defect of a curved interface in a cell), exact on a plane:
on the plane case the first step gives 416 cells with 1 cut and the
residual max 5.7e-15.

Verdicts:
- C1 HELD at N = 32 in 25 of 28 file-normal pairs (the others by one cell
  at 1.0 to 1.2e-14); FALSIFIED at N = 64 (max 1.8e-14) and N = 128 (max
  4.0e-14) by the round-off floor of the cut volume, as pre-registered for
  N = 128.
- C2 HELD in all 84: the signed sum is 0.1 to 2.9 N eps V.
- C3 FALSIFIED (1 to 2 cuts): 2.8 to 3.2; the agent's 2 to 3 holds at N
  <= 64 on the mean and is exceeded by 0.2 at N = 128.
- C4 HELD except the cells with alpha <= 0 or below the floor (1 to 7 per
  file), which take the vertex extreme.
- C5 HELD: 47 % less per positioning.
- Unchanged from the section above: the orientation (0 flips, 0
  sphere-side violations), the two clips within 4e-14 V.

### The anchor (job 55929941, rotation N = 32, Co 0.2, the ring row, pre-test keys pinned)

| run | anchor evaluations per cell (step 1 / last) | start \|r\|/V max, mean (step 1) | within tol (step 1 / last) | \|V_KB - V_cut\|/V max over the turn | wall s | L1(alphaf), E_shape at the end |
|---|---|---|---|---|---|---|
| vofCutCell | 1.37 / 1.55 | 2.1e-10, 1.6e-11 | 336 of 416 / 230 of 423 | - | 158.8 | 5.703502513999e-02, 1.058669062821e-03 |
| ccs | 1.00 / 1.00 | 1.1e-14, 1.7e-15 | 416 of 416 / 420 of 423 | 1.6e-14 | 148.5 | 5.703502513976e-02, 1.058669062804e-03 |
| ccs + anchorSkipConverged | 0 / 0.014 (skipped 416 / 417) | - | - | - | 144.3 | 5.703502512518e-02, 1.058669062810e-03 |

The three cells not within tolerance at the end are the ones without a
plane (alpha at 0 or 1, the start residual 0.5). The Kromer-Bothe volume
and the cutCellPLIC cut of the same plane agree to 1.6e-14 V on this case,
3.7e-14 at N = 64 and 8.4e-14 at N = 128 (cmp-pos). The end alpha fields:
max |alpha_ccs - alpha_vofCutCell| = 4.0e-11, max |alpha_skip - alpha_ccs|
= 5.4e-10, the sums equal to 2e-14 relative. The skip changes the pieces'
fluxes within the anchor's own offset tolerance (1e-11 h), not within
machine epsilon: with the anchor the plane is still shifted by up to
1e-11 h from the reconstruction's; min alpha -8.8e-13 to -9.8e-13 and the
excess 1.1e-11 in all three.

### The timed comparison (b6a7dc9, pre-test keys pinned, jobs 55929943-67; cmp-pos3)

ring and exact normal x vofCutCell and ccs, one exclusive i01 node per
(flow, N, Co, repetition), the four solves in rotated order, everything
else as the Co 1 comparison of 2026-11-03. The rotation at N = 64, Co 0.2
exceeds the 30-minute jobs (four solves of 16 min) and has one repetition
or none; N = 128 ran in the first batch only (fork defaults of 3aaf49d
active: wall ratio 0.956 to 0.986, the norms as below).

| case | row | E_shape, L1(alphaf) ratio ccs / vofCutCell | L1(x) ratio | L1(n) ratio | wall ratio (median of 3) | cuts per cell | anchor evaluations per cell ccs / vofCutCell |
|---|---|---|---|---|---|---|---|
| diagonal N 32, Co 0.2 | ring / exact | 1.000000 / 1.000000 | 1.0017 / 1.0000 | 1.0000 / 1.0000 | 0.978 / 0.977 | 2.98 / 2.93 | 1.00 / 1.52 |
| diagonal N 32, Co 1 | ring / exact | 1.000000 / 1.000000 | 1.0000 / 1.0000 | 1.0000 / 1.0000 | 0.983 / 0.980 | 3.08 / 2.91 | 1.00 / 1.49 |
| diagonal N 64, Co 0.2 | ring / exact | 1.000000 / 1.000000 | 1.0046 / 1.0021 | 1.0088 / 1.0031 | 0.982 / 0.972 | 2.94 / 2.91 | 1.00 / 1.53 |
| diagonal N 64, Co 1 | ring / exact | 1.000000 / 1.000000 | 1.0000 / 1.0000 | 1.0000 / 1.0000 | 0.982 / 0.988 | 2.93 / 2.87 | 1.00 / 1.50 |
| rotation N 32, Co 0.2 | ring / exact | 1.000000 / 1.000000 | 1.0064 / 1.0034 | 1.0072 / 1.0821 | 0.977 / 0.973 | 2.96 / 2.93 | 1.00 / 1.55 |
| rotation N 32, Co 1 | ring / exact | 1.000000 / 1.000000 | 1.0060 / 1.0027 | 1.0048 / 1.1286 | 0.978 / 0.973 | 2.98 / 2.94 | 1.00 / 1.57 |
| rotation N 64, Co 1 | ring / exact | 1.000000 / 1.000000 | 1.0022 / 1.0000 | 1.0023 / 1.0000 | 0.980 / 0.972 | 2.97 / 2.93 | 1.00 / 1.55 |

Bounded and conservative alike: min alpha -1.7e-14 to -8.3e-11, the
excess 1.3e-13 to 3.5e-11, the mass within 3e-12 % in both. The ccs
residual's maximum over all steps and cells is 1.0e-14 to 6.6e-14 V.
- P-2 HELD: the wall clock is 1.2 to 2.8 % lower with ccs in every case.
- P-1 HELD for the transport: E_shape and L1(alphaf) are identical to
  seven digits in every case (the anchor re-solves the offset either way,
  so the clipped plane is the same within 1e-11 h). FALSIFIED for the
  reconstruction norms of the end state, L1(x) up to +0.6 % and L1(n) up
  to +13 % (the exact row at N = 32): with ccs the cells with alpha within
  1e-8 of 0 or 1 carry a plane, a corner sliver whose centroid and normal
  enter the norm; vofCutCell leaves them without.

### The plane gate (b6a7dc9, config planeTranslation3D_hexSweptWispM9Ccs2, job 55929942; 57 of 64 cases at the time of writing)

Against planeTranslation3D_hexSweptWispM9, case by case: the gate verdicts
are identical in all 57 (52 pass both, the 5 failing M9 cases fail with
ccs too); the run maxima of L1(alphaf) are round-off numbers (4e-14 to
1e-11) and differ by -18 % to +24 % in both directions (30 identical or
better, 27 worse); min alpha and the mass clause alike. The ClockTimes are
integer seconds on shared nodes (median ratio 1.2 between 0.6 and 1.6):
not a timing measurement; the timing is the comparison above. At N = 20
the ccs run maxima are 3.0e-14 to 1.7e-13 against M9's 2.8e-14 to 2.0e-13.

### Surprises

- The fork defaults of 3aaf49d (the other session's sweptPreTest
  touching and sweptStateShortcut true, in every build from de416fa on)
  are not inert: the N = 20 plane case loses the plane at step 1
  (L1(alphaf) 4.8e-3 against 6e-15, the step-1 swept line "pre-test: ring
  cells skipped 15", the beyond-layer pieces 0 against 16). The 2 x 2 check
  (job 55929616, the same case): touching / true gives L1(alphaf) max
  9.0e-2 (vofCutCell) and 9.1e-2 (ccs), min alpha -2.2e-2; none / false
  gives 4.23e-14 (vofCutCell, the M9 reference exactly) and 4.08e-14
  (ccs). On the sphere cases the norms agree with the 720cfc6 reference to
  four digits with either setting, and the step is 37 % faster with the
  pre-test (rotation N = 32, Co 0.2, ring: 217 against 344 s). The first
  gate (job 55928074) and the first comparison ran with it; the tokens
  TPF_SWEPT_PRE_TEST (none) and TPF_SWEPT_STATE_SHORTCUT (false) now pin
  the 720cfc6 path in the templates; the fork default is theirs to decide.
- Two sessions commit in the same worktree: 3aaf49d carried this task's
  isoAdvection edits two minutes before de416fa.
- "Within one iteration" holds on a planar state only (1 cut); on the
  sphere the closed-form start is 1 to 8 % of V off and the third cut is
  the answer.

What this demonstrates: positioning ccs is a drop-in replacement, opt-in
by the key: the planes sit at the round-off floor (1e-15 V mean, 4e-14
max at N = 128) with three cuts per cell, 47 % less per positioning and
1 to 3 % less per step; the anchor then starts within its tolerance in
every cell with a plane and reduces to its start evaluation, or none with
anchorSkipConverged, the fluxes changed by at most 5e-10 in alpha; the
transported fields are identical in E_shape and L1(alphaf). The default
stays vofCutCell.

## The efficiency of the piecewise swept-face advection (2026-11-04)

Tomislav's decision after the profile: improve the efficiency of the
original piecewise algorithm, guided by the HPCToolkit breakdown, every
change giving the same numbers bit for bit (or round-off, stated and
explained per change) on 20-step copies, keyed when it alters any number,
the serial step and the breakdown measured after each; then the
decomposition. His addition: the state shortcut (a face whose touched
cells are all full or all empty takes the upwind flux without a clip). The
plane gate entered the verification set after the first shortcut proved
wrong there (below). The fork's defaults are unchanged (sweptPreTest none,
sweptStateShortcut false, sweptTraceCache false, since 158f486): every
build reproduces 720cfc6 bit for bit unless a key is set; leia's templates
carry TPF_SWEPT_PRE_TEST and TPF_SWEPT_STATE_SHORTCUT (d809acc0, the CCS
session), not yet a token for sweptTraceCache.

### The verification and the measure

- The identity set, 20 steps each, serial, writePrecision 17, every written
  field compared bit for bit (cmp) and alpha and alphaf by max |difference|,
  the log's counters diffed: the rotation N 32 and 64 (ring row, OF Co 1),
  the diagonal line N 32 (ring row, OF Co 1), the plane N 40 normal (exact
  row, hexSweptWispM9, Co 0.95) and the plane gate (hexSweptPyr: N 20 and
  40, both arms, ring row, Co 0.5, the final wisp rule matched + consistent;
  sweptVolumeCorrection scale everywhere). The reference is the code with
  the keys off (ref3, fork 158f486), itself bit for bit 720cfc6 / 9394337
  (ref2 at de416fa against "before") on every case: the CCS positioning
  commits in between preserve the defaults.
- The serial timing: the rotation N 64, OF Co 1, final wisp rule, 200
  steps, three repetitions on an exclusive i01 node, the mean of the
  per-step time over steps 2 to 199, the median and spread of the three
  (the serial step is stationary over the turn: 1.05 to 1.18 s per 50-step
  bin). The breakdown: 100 steps under hpcrun (CPUTIME 1 kHz) as in the
  profile section, the categories by routine or source-line range of the
  current sources (rules_swept2.py; rules_swept.py for the original); a
  100-step profile carries the start-up and the final write (0.18 s/step
  of "I/O" and "main" that the 200-step timing excludes).
- The cycle: workflow/scripts/davof_prof/opt-cycle.sbatch (the build of the
  profiling worktree at a commit, the identity set against a reference tag,
  the timing, the profile; the keys per run inserted into fvSolution); the
  profiling build throughout (-O3 -g -fno-omit-frame-pointer).
- The hazard met: a second session (CCS positioning) committed in the same
  worktree; 3aaf49d carried its uncommitted include of a file it added two
  minutes later (de416fa), so 3aaf49d alone does not build, and the first
  cycles on it ran a binary with one stale object. From 158f486 on, my
  edits live in a worktree of their own (TwoPhaseFlow-opt) and are
  committed through the index (git apply --cached) so neither session's
  pending edits enter the other's commits; the shared clone's working tree
  then lags HEAD for the edited files until reset.

### The changes (fork davof-swept; every commit names the keys it adds)

| commit | change | key (default) | identity |
|---|---|---|---|
| 3aaf49d, 5e4d0c0, 85cc209 | the pre-test of a ring cell before the two clips: the cell is skipped when one of its planes has every vertex of the swept volume (x_f, the face points, the apex after the pyramid, the traced points) on its outer side. strict: beyond 1e-10 L (L the longest triangle edge, the bound of every polygon size the clips meet and the tolerance of the contact and coplanar tests): the clips of such a cell remove every triangle. touching: signed distance >= 0 (the cells sharing only a face, an edge or a point with the swept volume) unless a pre-image triangle lies in the plane (the coplanar rule of the face update keeps it). Two passes over the ring: the planes, state and counters of every cell as before, then the clips of the touched cells in the same order with the same sums. | sweptPreTest none / strict / touching (none) | strict and touching both bit for bit on all eight cases |
| 1e279d4 | the clip internals: the signed distances of the vertices computed once (planeDistances) and shared between the contact / coplanar tests and the clip; the clip into a second polygon, the two buffers alternated (no copy of 16 points and labels back after every clip, hotspot 5); the triangle's three distances per plane kept from the contact test and reused while the polygon is still the triangle (hotspot 3), a plane the triangle lies beyond by 1e-12 L drops it before any clip, a plane it lies within by 1e-12 L is not clipped after a crossing (no crossing point's round-off reaches the margin); the polygon size L measured for the coplanar test only when the largest distance is within 1e-10 of the triangle's longest edge (hotspot 8); the cell bounding boxes of the ring search cached per step (hotspot 9) | none (unconditional) | bit for bit (c2 against ref2) |
| 4ec4bb7, 557fab7 | the trace cache: per mesh point the cell that held its midpoint and the velocity there; another face holding the point repeats only the cell search from its own donor (the ring order of that donor, as the per-face trace) and takes the kept velocity when the same cell results, else interpolates anew (hotspot 7: the tetrahedra of the cell rebuilt per lookup). The first version (4ec4bb7) reused the traced point itself: a midpoint within the tolerance of a cell face is claimed by the first cell of the donor's search order, the velocity there differs by the extrapolation over the tolerance, and the wisp rules amplified it (rotation N 64, 20 steps: alpha 1.1e-9, alphaf 1.2e-8) | sweptTraceCache (false) | 557fab7 bit for bit on all eight cases (c4); 4ec4bb7 round-off, amplified |
| 3aaf49d, 5e4d0c0, 85cc209 | the state shortcut: a face whose touched ring cells are all full (alpha^n exactly 1, no surface) takes V_liq = \|phi_f\| dt and alpha_f = 1, all empty 0 and 0, without a clip; exact 0 and 1 rather than the wisp tolerance since sweptConsistent gives a piece of a cell at alpha = 1 - 1e-10 alpha V_piece (the faces the tolerance would take in addition are counted). The condition that the pieces cover the swept volume: 3aaf49d none (wrong, below), 5e4d0c0 every vertex inside the donor (never true under a velocity with a component along the face: 0 of 17528 faces on the rotation), 85cc209 every vertex inside a touched ring cell | sweptStateShortcut (false) | round-off, amplified: rotation N 64 alpha 1.2e-10, alphaf 1.0e-6 after 20 steps; plane N 40 oblique 5.7e-13 / 1.7e-12; three of the eight cases bit for bit |
| 5e4d0c0, 0703ee6, 85cc209, 09a6c07 | the debug key: every skipped ring cell clipped and counted when it holds a piece by the pieces' own gate (coverage above 1e-12 of the boundary) and not a degenerate slab, or a nonzero wetted area; a shortcut face run through the pieces and counted when a sum differs by more than 1e-10; the first forty printed with the face, the cells, the coverages, the separating plane, the piece's volume, caps and contacts | sweptPreTestCheck (false) | - |
| 158f486 | the defaults back to none / false / false (the code of 720cfc6 by default); dc75317 before it changed only the executable bits of two Allwmake scripts under the same message | - | - |

The first shortcut (3aaf49d) and why it was wrong: it assumed the touched
cells cover the swept volume. On the plane gate the inflow boundary faces
have their pre-image outside the mesh; the pieces give V_liq = 0 and
alpha_f = 0 there (the uncovered pre-image, the "beyond layer" faces), the
shortcut gave \|phi_f\| dt and 1: L1(alphaf) 4.8e-3 at step 1, min alpha
-2.2e-2 over the run (job 55929616, the CCS session's 2 x 2 check; the
sphere norms were inert to four digits, so the sphere copies had not
caught it: no face of theirs has a swept volume leaving its ring). The
defaults were reverted (158f486) before anything else.

What the check found (sweptPreTestCheck, cycles chkT, chkS, chkS2 and the
debug runs 55933853 / 55934488):

- The pre-test: on the sphere cases, the plane N 40 and the plane N 20
  normal no skipped cell holds a piece. On the plane N 20 the four inflow
  boundary faces per step whose donor the touching rule skips: the swept
  volume lies in the donor's boundary plane (the face points on it, the
  pre-image outside the mesh), the donor's piece is the face triangles kept
  as lying in the plane (the contact rule, the donor inside): a
  zero-thickness slab, degenerate, volume 0, wetted area 0, 25 to 29 % of
  the boundary as coverage, nothing taken from it. On the rotation N 32 a
  few skipped cells per step return a wetted area of 1.4e-14 of A_pre from
  zero-area polygons (coverage 0, two contacts): the round-off the clips
  of the unskipped code add to the wetted area vector; the fields stay bit
  for bit over the 20 steps on every case, so these round-off areas did not
  reach a written digit. The strict rule has no such case by construction;
  the touching rule's identity is empirical on the eight cases.
- The shortcut with the coverage condition of 85cc209: the faces the check
  flags are faces whose donor piece has an unsolvable cap system (ok 0,
  three or four caps), i.e. faces the pieces hand to the rank fallback
  (the one-donor clip), which for a full donor without a surface returns
  exactly the shortcut's values; the coverage of those faces is 1.000.
  The shortcut is therefore the pieces' result up to round-off on every
  face, but that round-off is amplified by the discrete wisp rules
  (matching at 1e-11) to 1e-6 in alphaf within 20 steps on the rotation N
  64, as any round-off change is on that case (the first trace cache
  showed the same). And it gains nothing: the candidate band is one to two
  layers around the interface, so nearly every candidate face has an
  interface cell in its ring; on the rotation N 64 the shortcut takes 540
  of 17571 faces (3 full, 536 empty at the last step) and the coverage
  test costs what it saves: 0.542 s/step against 0.541 without it (the
  unsafe first version, no coverage test: 0.478). Not achieved as asked;
  the key stays opt-in, default false.

### The serial step after each change (rotation N 64, OF Co 1, final wisp rule)

| configuration | commit (keys) | s/step, median of 3 (spread) | per step: ring cells skipped / clipped; trace points reused / traced | against 1.110 |
|---|---|---|---|---|
| the original | 9394337 (none); de416fa keys off | 1.110 (0.1 %); 1.122 | 0 / 316278 | 1 |
| clip internals | 1e279d4, keys off | 0.937 (0.2 %) | 0 / 316278 | 0.84 |
| + pre-test strict | b6a7dc9, strict | 0.978 (0.3 %) | 0 / 316278 | 0.88 |
| + pre-test touching | b6a7dc9, touching | 0.541 (0.2 %) | 244156 / 72122 | 0.49 |
| touching + the first trace cache (round-off) | b6a7dc9 (4ec4bb7), touching + cache | 0.488 (0.6 %) | 244156 / 72122; 63282 / 7070 | 0.44 |
| touching + the exact trace cache (THE RESULT) | 09a6c07 (557fab7), touching + cache | 0.519 (0.5 %) | 244156 / 72122; 45273 / 24799 | 0.47 |
| touching + the unsafe shortcut (3aaf49d, wrong on the plane) | b6a7dc9 | 0.478 (0.5 %) | faces full / empty / pieces 4 / 536 / 17031 | 0.43 |
| touching + the covered shortcut (round-off, amplified) | 85cc209 | 0.542 (0.2 %) | 3 / 536 / 17032 | 0.49 |

- The result: 1.110 -> 0.519 s/step, 2.14 times, bit for bit; against
  plicRDF + isoAdvector serial on the same case (0.209 s/step) the ratio
  is 2.5 instead of 5.3.
- The strict pre-test skips nothing on the rotation: no ring cell lies
  beyond 1e-10 L of the swept volume, the touching cells ARE the ring; it
  costs its own test (4 %). The touching rule skips 77 % of the ring cells
  (4.1 clipped of 18 per face against 1.36 pieces per face) and halves the
  step.
- The exact trace cache reuses 65 % of the vertex traces (45273 of 70072
  per step; the first version 90 %, the rest are midpoints whose cell
  search from another donor finds another cell within the tolerance).

### The breakdown before and after (100-step profiles, s/step, the same case)

| category | original (ref2) | after (c4) |
|---|---|---|
| pieces: the clip by the cell planes and the PLIC plane | 0.465 | 0.061 |
| face update: the pre-image fan clipped | 0.189 | 0.034 |
| swept cells: the pre-test | - | 0.039 |
| swept cells: the ring search | 0.041 | 0.018 |
| swept cells: the cell planes, state, wisp rule | 0.029 | 0.027 |
| swept cells: the rest of sweptFaceCellsFlux | 0.017 | 0.011 |
| trace: the velocity at the traced points (findTetrahedron) | 0.077 | 0.054 |
| trace: the cell search (ring, octree) | 0.025 | 0.023 |
| trace: the triangle list, the midpoint trace | 0.011 | 0.010 |
| the point velocity (volPointInterpolation of U) | 0.023 | 0.024 |
| the surface build | 0.043 | 0.043 |
| the reconstruction | 0.028 | 0.027 |
| curvedFlux (candidate band, wisp surfaces, face loop) | 0.018 | 0.023 |
| the bounding, timeIntegratedFlux, the alpha update | 0.049 | 0.049 |
| setFlow, the Courant number, the solver loop | 0.116 | 0.117 |
| I/O and start-up (the 100-step window's write and reads) | 0.110 | 0.109 |
| total of the window | 1.264 | 0.677 |

- The two clips fell from 0.654 to 0.095 s/step (the pre-test's 0.039
  included: 0.134), the swept advection's own parts from 0.93 to 0.37.
  What is left above 5 % of the step: the trace's velocity lookup (8 %,
  OpenFOAM's cellPointWeight rebuilding the cell's tetrahedra per lookup:
  a per-cell cache of the decomposition would take it further), the pre-test
  itself (6 %: 60 dot products per ring cell, 18 cells per face; a cheaper
  first test by the cell's bounding box against each vertex would halve it),
  the surface build and the reconstruction (6 and 4 %), the point velocity
  of the whole mesh (3.5 %, mesh-wide), and the parts shared with plicRDF
  (setFlow, the alpha update, the Courant number, 17 %).
- The top exclusive functions after: the vector subtraction and dot product
  of the clips (3.5 and 3.8 %), sweptFaceCellsFlux's own lines (the
  pre-test), findTetrahedron, volPointInterpolation, the per-step rebuild
  of the cell bounding boxes (1.5 %: once per mesh would do), malloc / free
  (3.4 %).

### The decomposition (item 5): the optimised code in parallel

The runs: the comparison's rotation (ring row, OF Co 1, final wisp rule),
sweptPreTest touching and sweptTraceCache on, full turns on exclusive i01
nodes, the undecomposed study case decomposed in the job
(workflow/scripts/davof_prof/opt-par.sbatch): scotch over the whole mesh
(the comparison's decomposition) against scotch with the cell weights of
cellweight.py (W for the cells within two face layers of a partial cell
at t = 0, 1 elsewhere; decomposeParDict weightField). N = 64 on 8 ranks
three turns, N = 128 on 64 ranks two turns; the per-step mean over the
turn, the 50-step-bin history, the partial cells per rank at t = 0.

| N, ranks | decomposition | cells per rank | partial cells per rank at t = 0 | s/step (turn means) | 50-step bins | before the optimisation (scotch, profile section) |
|---|---|---|---|---|---|---|
| 64, 8 | scotch | 64885 to 66139 | 532 438 394 396 0 0 0 0 (max/mean 2.4) | 0.2355 / 0.2366 / 0.2360 | 0.16 0.21 0.30 0.32 0.32 0.28 0.19 0.16 0.17 0.20 0.26 0.28 0.25 0.20 | 0.545 |
| 64, 8 | weighted W 100 | 1422 to 137284 | 0 0 0 0 399 466 459 436 (2.1) | 0.2933 / 0.2935 / 0.2928 | 0.17 0.23 0.36 0.41 0.39 0.29 0.22 0.20 0.23 0.33 0.42 0.38 0.27 0.18 | - |
| 64, 8 | weighted W 1000 | 805 to 517098 | 0 220 203 239 278 310 258 252 (1.4) | 0.4991 / 0.5000 | 0.35 0.43 0.51 0.55 0.54 0.53 0.53 0.54 0.54 0.55 0.56 0.52 0.45 0.39 | - |
| 64, 8 | W 100, redistributePar every 100 steps (7 cycles) | 2k to 134k after each renewal | at the renewals: max/mean 1.76 to 2.03, 0 to 2 ranks without | 0.2305 solve + 6.9 s per redistribution (41 s of 235 s per turn) | 0.17 0.26 0.21 0.33 0.19 0.26 0.21 0.29 0.18 0.24 0.18 0.23 0.23 0.24 | - |
| 128, 64 | scotch | 64881 to 66191 | 55 ranks with 0, max 1746 (16.2) | 0.5834 / 0.5824 | 0.47 0.50 0.48 0.50 0.51 0.52 0.50 0.49 0.58 0.69 0.79 0.85 0.85 0.69 0.53 0.51 0.58 0.63 0.66 0.62 0.59 0.52 0.49 0.53 0.60 0.61 0.54 0.50 0.67 | 1.341 |
| 128, 64 | weighted W 100 | 992 to 101162 | 35 ranks with 0, max 373 (3.5) | 0.7102 / 0.7097 | 0.35 0.42 0.44 0.54 0.66 0.84 0.95 0.88 0.82 0.83 0.75 0.66 0.77 0.85 0.83 0.77 0.69 0.77 0.80 0.71 0.68 0.86 1.03 0.96 0.75 0.55 0.41 0.30 0.94 | - |
| 128, 64 | W 100, redistributePar every 100 steps (7 cycles = half a turn, then redistributePar killed, rc 137, at the seventh renewal) | 1k to 100k after each renewal | at the renewals: 36 to 39 ranks without, max/mean 3.6 to 4.0 | 0.4195 solve over the half turn + 13.9 s per redistribution (84 s of 439 s) | 0.35 0.43 0.42 0.53 0.34 0.45 0.39 0.50 0.31 0.43 0.39 0.45 0.36 0.53 | - |

- The optimised code in parallel with the comparison's decomposition: 0.545
  -> 0.236 s/step at 8 ranks (2.3 times), 1.341 -> 0.583 at 64 ranks (2.3
  times); the parallel step follows the serial one (the busiest rank's
  compute shrank with it), the imbalance is unchanged: the history over the
  turn still varies 2 x (0.16 to 0.32, 0.47 to 0.85 s per 50 steps).
- The speed-up against the serial step: 0.519 -> 0.236 at 8 ranks = 2.2
  (before 2.0); against the serial N = 128 steady step the 64-rank step is
  now 0.58 s (the serial N = 128 step was not re-measured after the
  optimisation; with the serial parts scaling as at N = 64, ~ 2.6 s, the
  speed-up stays ~ 4.5, the balanced step would be ~ 0.05 s).
- The static weighting does not help; it hurts. With W 100 scotch gives the
  band at t = 0 half of the weight and spreads it over four ranks of 1.4k
  to 10k cells, the other four ranks take 100k to 137k cells each: the
  band is as concentrated as before (four ranks hold it, max/mean 2.1
  against 2.4) and as the sphere turns the band leaves the small ranks and
  crosses the big ones (bins 0.17 to 0.42 against 0.16 to 0.32), while the
  mesh-wide work (setFlow, the alpha update, the point velocity: 0.14
  s/step serial) is now 2 x unbalanced. With W 1000 one rank holds 517k of
  524k cells and the step doubles. At N = 128 W 100 halves the band's
  concentration at t = 0 (35 ranks without a partial cell against 55,
  max/mean 3.5 against 16) with ranks of 1k to 101k cells, and the turn is
  slower (0.710 against 0.583 s/step): the first 100 steps are faster (0.35
  against 0.47, the band on its weighted ranks), then the band leaves them
  and the step rises to 0.8 to 1.0 s on the big ranks.
- The renewed weighted decomposition (5 b, N 64, every 100 steps) does not
  help either: the solve step is scotch's (0.2305 against 0.236 s) and the
  six redistributions add 41 s to a 193 s turn. At every renewal scotch with
  the weights keeps the band compact (max/mean 1.8 to 2.0 over the band
  ranks, two ranks without a partial cell at four of the six renewals) and
  the rest of the mesh on ranks of 2k to 134k cells. The reason is scotch
  itself: it balances the total weight per rank while minimising the edge
  cut, so with a heavy compact band it gives each rank either a piece of
  the band or a block of the rest, never a slice of both; the band's share
  of the ranks does not grow with the weight. Balancing the swept work
  needs a decomposition that slices the band across every rank (the band's
  cells partitioned on their own, each rank one slice, the rest filled in
  around it: a custom decomposition, or scotch's multi-constraint mode
  which OpenFOAM's interface does not expose), and then a renewal as the
  band moves. Not done here.

### The restart and the redistribution (item 5 b)

The fork restarts bit for bit: the rotation N 32, 40 steps in one go
against 20 steps, a stop and 20 more from the written state (writePrecision
17), every field at step 40 identical (job 55931978; the last time value
differs in its 17th digit from the accumulated dt, which changes nothing
written). OpenFOAM v2606's decomposePar and redistributePar read a
weightField from decomposeParDict (a volScalarField of cell weights,
scotch weighted; applications/utilities/parallelProcessing), so a
decomposed run can be renewed between steps without any change to the fork:
opt-redist.sbatch runs K steps, writes, computes the weights per processor
from the written alpha (cellweight_par.py, the layers grown inside each
processor), redistributePar -parallel -overwrite -latestTime, and continues
from the written state; the solve and the redistribution wall clocks are
kept apart. The weight field needs the processor patches typed processor
(OpenFOAM names them procBoundary<i>to<j> and tries the last listed
pattern first). Results: in the table above. N 64: the solve step is
scotch's, the renewals cost 20 % of the turn. N 128: over the half turn it
completed the solve step is 0.42 s against scotch's 0.47 to 0.52 over the
same window (the band on the weighted ranks after each renewal, 36 to 39
ranks still without a partial cell), the renewals cost 13.9 s each (0.12
s per step, which cancels the gain), and the seventh renewal was killed
(rc 137, job 55936389; the cause is not in the log, the memory of 64
ranks at 3000 MB per rank redistributing 2M cells is the suspect). The
same conclusion as at N 64: with scotch's partition of a weighted graph
the band's share of the ranks does not grow, so a renewal moves the
compact band along with the sphere and cannot balance it.

### The 18-19 x of the parallel branch, reconciled

The parallel branch's report (config sphereRotation3D_hexSweptParN128np*,
the MEASURED block, 2026-10-31 / 11-01, fork fc08178): N = 128, 10 steps,
target Co 0.2 (OF Co 0.88), the median of steps 2 to 9, rows exact / ring:
np 1 8.90 / 8.86 s, np 64 0.490 / 0.470 s, 18.2 / 18.9 times; the advect
call 3.33 / 3.30 s serial and 0.34 to 0.36 s at np 64. The same two
configs rerun today on the current code with the keys off (jobs 55934645
/ 55934646, the default build at 09a6c07, the study's own scotch
decomposeParDict, the study's tokens: no wisp rule, no consistent band,
anchorMethod ccs): np 1 11.59 / 11.55 s per step with the advect call
3.11 / 3.05 s; np 64 0.595 / 0.530 s with the advect call 0.356 / 0.316 s;
19.5 / 21.8 times. So the branch's numbers reproduce in kind, and the
reconciliation with the 5.4 of the profile section is:

- What the two serial figures measure. The branch's "serial 8.9 s" is the
  whole step of a 10-step run whose per-step diagnostics (areaFractionError
  and volumeFractionError every step on 2M cells, the Courant number, the
  field algebra) cost 5.6 s of the 8.9 (8.5 of 11.6 today); the swept
  advection (the advect call) is 3.3 s. That part parallelises perfectly,
  the advect call does not: 3.3 -> 0.34 s at 64 ranks is 9.6 times
  (today 3.05 -> 0.32, 9.6 times), the diagnostics 5.6 -> 0.1. The 18-19
  times is the ratio of two totals dominated by the diagnostics. The
  profile section's 5.4 divides the steady serial step (5.5 s at OF Co 1,
  the comparison's settings, no per-step diagnostics) by the 64-rank
  20-step mean (1.02 s) of the comparison's runs.
- The code: fc08178 has no wisp rule and no consistent band; the
  comparison's final rule computes every face whose donor shares a vertex
  with an active cell (4752 faces without a surface per step at N 64) and
  the wisp surfaces, which doubles the swept work per step (the comparison
  64-rank step 1.02 s against 0.53 here, both scotch).
- Co: OF 0.88 against 1 changes the pieces little (the clipped triangles
  per face are the same); the step count: 10 steps from the start against
  a turn, over which the busiest rank changes and the per-step time varies
  2 x (the 10-step window sits at 0.47 to 0.50 of the N 128 history above
  while the turn's mean is 0.58).
- The decomposition: both scotch over the whole mesh; the branch's run
  shows the same concentration ("the interface band sits in a few ranks").

The two figures are consistent: the advection itself scales 9 to 10 times
at 64 ranks in both, and neither better, because the band sits on 9 of 64
ranks.

What this demonstrates: the piecewise swept-face advection runs 2.14 times
faster serially (1.110 -> 0.519 s/step on the rotation N 64) with the same
numbers bit for bit on eight 20-step cases, by not clipping the ring cells
the swept volume only touches (the touching pre-test, 77 % of the ring),
by the clip internals, and by tracing each mesh point once per step with
an exact reuse; in parallel the step follows (0.545 -> 0.236 s at 8 ranks,
1.341 -> 0.583 s at 64 ranks, 2.3 times) and the imbalance of the
interface band is untouched, because scotch with cell weights, static or
renewed, keeps the band compact; the state shortcut is equivalent to the
pieces only up to a round-off that the wisp rules amplify and gains
nothing on the rotation, so it stays opt-in false; the 18-19 x of the
parallel branch is the scaling of a step dominated by per-step diagnostics,
its advection scaled 9.6 times then and now.



## The reference setting of the swept-face advection (2026-11-04)

Tomislav's instruction (via the coordinating session): make the measured
best configuration the fork's default, so that a dictionary without the keys
runs it, and verify the defaults by a rerun with no keys in the dictionaries.

### The defaults (fork davof-swept 407f969)

| key | old default | new default | provenance |
|---|---|---|---|
| sweptWispSurface | none | matched | the plane gate (2026-11-02/03) |
| sweptWispMatchTol | 1e-11 | 1e-9 | the plane gate |
| sweptWispTol | 1e-8 | 0 | the plane gate |
| sweptConsistent | false | true | the consistency under oblique travel (2026-11-02); the wall faces kept at the boundary condition (eb17c0a) |
| sweptVolumeCorrection | auto (= none) | pyramid | the UFVFC pyramid (2026-11-03): exact on uniform translation, within 1e-6 of the scaling at OpenFOAM Co <= 1 |
| positioning | vofCutCell | ccs (positioningTol 1e-14, positioningMaxIter 30 unchanged) | the Consecutive Cubic Spline positioning (2026-11-04) |
| anchorSkipConverged | false | true | the same: the anchor skipped where ccs is within its tolerance |
| sweptPreTest | none | touching | the efficiency work (2026-11-04), bit for bit |
| sweptTraceCache | false | true | the same, bit for bit |

Every old value stays selectable by writing the key. The fork's
etc/README-davof.md carries the same table. leia f74c3e29: the token defaults
of the three davof cases are these values, the token TPF_SWEPT_TRACE_CACHE
added; the configs that pin an old value are named in its message (the Cmp*
and *CoScale configs pin sweptVolumeCorrection scale, the Wisp / WispA
configs the surfaces neighbour / anchored, the *ConsLog configs sweptConsistent
false, planeTranslation3D_hexSweptWispM9Ccs2 sweptPreTest none); the
sweptFaceCells configs before the wisp rule (Fine*, Hist*, RotFinal*,
SweptPar*, SweptCons*, SweptLadder*, the first plane configs) pin none of
the nine tokens and would run the reference setting if re-rendered.

### The rerun with no keys (PRE-REGISTERED 2026-11-04, 23:40, before the submission)

The runner: workflow/scripts/davof_prof/ref-run.sbatch copies the rendered
study cases and deletes the thirteen lines of the reference-setting keys
(the nine above, positioningTol, positioningMaxIter, sweptStateShortcut,
sweptScaleToPhi) from alpha.water, so that the fork's defaults act; the
converter scores every case; the main clone's default build (407f969) runs.
The cells:
- the plane gate matrix: planeTranslation3D_hexSweptWispM9Ccs2's cases at
  target Co 0.2 and 0.95, N 20 / 40 / 80 / 160, both arms, exact and ring (32
  cases, serial);
- the rotation (one turn), the sphere along x and along the x-z diagonal
  (END 0.45): N 32 / 64 / 128 on 1 / 8 / 64 ranks at OpenFOAM Co 1 and 0.2,
  the comparison's cases; the ring row three repetitions, the exact row one
  (the control's norms; its wall clock was within 4 % of the ring row's in
  the comparison); one exclusive i01 node per (flow, N, Co), wall = the
  srun step's wall time, as the comparison;
- the pure spin: sphereRotation3D_hexSpinCoPyr's cases (N 32 / 64, target
  Co 0.2 / 0.5, exact row), serial;
- the parallel identity: sphereRotation3D_hexSweptParN32np{1,2,4,8,16} (N
  32, 20 steps, both rows, every step written at precision 17, the face
  log), scored by davof_decomposition_compare.py against np 1.
Predictions:
- D-1 (the errors equal the opt-in runs'): the plane gate's verdicts are
  those of M9Ccs2 case by case (the cases that fail M9Ccs2's gate at target
  Co 0.2 / 0.95 fail here too, none else); every run
  maximum of L1(alphaf) and L1(n) of a passing case is a round-off number
  (<= 1e-11); the rotation, the lines and the spin have E_shape and
  L1(alphaf) within 1e-6 relative of the comparison's (ring and exact; the
  pyramid against the scaling: 1e-6 at OpenFOAM Co <= 1; the positioning:
  identical to seven digits) and within 1e-6 of hexSpinCoPyr's; L1(x) within
  1 %; L1(n) at the end within 1 % on the ring row and up to +13 % on the
  exact row at N 32 (the corner slivers that ccs gives a plane). FALSIFIED
  by any case outside these bands.
- D-2 (bounded, conservative): every case min alpha >= -1e-10, excess <=
  1.5e-10 (the comparison's largest, 1.3e-10, was the exact row at N 128,
  OpenFOAM Co 1), |mass deviation| <= 1e-9 %.
- D-3 (the wall clock): at equal ranks the ring row's median wall clock is
  2.0 to 2.5 times lower than the comparison's in every one of the 18
  timed cells (the serial step 2.14 x, the scotch cells 2.3 x in the
  efficiency work; ccs and the skipped anchor 1 to 3 %; the pyramid at the
  cost of the scaling). FALSIFIED by a cell below 2.0 or above 2.5.
- D-4 (the parallel identity): at step 1 max |alpha - alpha_np1| <= 1e-12
  and max |alphaf - alphaf_np1| <= 1e-12 for np 2 / 4 / 8 / 16, both rows;
  nMissing 0, nDup 0 at every step; over the 20 steps max |dalpha| <= 6e-10
  (exact) / 2.5e-8 (ring), the default family's bands; a step above it is
  reported with its cause (the rank-fallback discontinuity of P-P2, if it
  returns).
- D-5 (the profile): the serial N = 64 step of the profiling build falls
  from 1.11 to 0.50 to 0.55 s (the 100-step cycle measured 0.519 s at
  touching + trace cache; ccs and the skipped anchor take 1 to 3 % more
  off, the pyramid nothing), the two clips to 0.09 to 0.10 s/step; at 8
  ranks the step from 0.545 to 0.23 to 0.26 s, at 64 ranks from 1.34 to
  0.55 to 0.62 s (the comparison's scotch decompositions as before); the
  MPI share stays above 70 / 90 %.


### MEASURED (2026-11-05, 00:37 to 06:30; every cell finished)

Jobs: the main clone's build 55939665, the profiling build 55939666; the
plane gate 55944952 (N 20 / 40 / 80) and 55940563 (N 160); the spin
55944116; the parallel identity 55944115; the timed cells 55944119-30,
55945423-24, and, one 30-minute job per repetition on separate i01 nodes,
55946942-48 (Rot N 64 Co 0.2), 55946950-56 (Rot N 128 Co 1), 55946958 /
55946960 (LineDiag N 128 Co 0.2); Rot N 128 Co 0.2 55940570 (four runs of
an hour on 64 ranks, in deflt from 02:19); the profiles
55944117 (serial), 55944118 (8 ranks), 55946962 / 55946963 (64 ranks: the
measurement, then hpcstruct + hpcprof). The cluster's deflt partition
started no job of this set between 23:45 and 02:18; the cells that fit 30
minutes were moved to deflt_short.
Every run's startup line reads "sweptPreTest touching, sweptStateShortcut
0, sweptTraceCache 1, sweptWispSurface matched, sweptWispTol 0,
sweptWispMatchTol 1e-09, sweptConsistent true" with no key in the
dictionary (the strip left 0 key lines).

CSVs on the cluster (/work/scratch/tm83tomy/davof/agents/swept):
csv/ref_setting.csv (every case and method of the rerun: the norms, their
relative difference to the opt-in runs, the wall clocks and the ratios);
csv/hp2_categories.csv, hp2_ranks.csv, hp2_hotspots.csv, hp2_timing.csv
(the two profiles side by side); the databases prof/runs2/serial/ring32_full,
ring64_20, ring64_full, prof/runs2/par8/ring64_full, prof/runs2/par64/
ring128_full (each <run>/hpctoolkit-database, the first profile's under
prof/runs/); the analyses prof/an2/<mode>_<run>.{summary.txt,flat.csv,
lines.csv,tree.txt,callers.txt,cat.csv,catranks.csv,ranks.csv} (the first
profile's under prof/an/). Scripts: workflow/scripts/davof_prof/ref-run.sbatch,
ref_an.py, hp-run2.sbatch, rules_swept3.py, hp2_an.sh, hp_side.py.

#### The plane gate (32 cases)

The rule, applied alike to both sides: the run maxima of L1(alphaf) and
L1(n) below 1e-12, min alpha >= -1e-12, the excess <= 1e-11, E_shape / V
< 1e-12.

| N | Co | arm | row | gate M9Ccs2 / no keys | L1(alphaf) max M9Ccs2 / no keys | L1(n) max M9Ccs2 / no keys | min alpha | excess | steps |
|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.2 | normal | exact | pass / pass | 4.08e-14 / 4.23e-14 | 3.4e-15 / 2.8e-15 | -2.7e-14 | 1.5e-13 | 20 |
| 20 | 0.2 | oblique | exact | pass / pass | 6.49e-14 / 6.34e-14 | 1.9e-15 / 1.8e-15 | -1.1e-14 | 3.6e-13 | 20 |
| 20 | 0.2 | normal | ring | pass / pass | 4.31e-14 / 3.85e-14 | 5.2e-14 / 5.9e-14 | -2.0e-14 | 1.5e-13 | 20 |
| 20 | 0.2 | oblique | ring | pass / pass | 6.50e-14 / 6.93e-14 | 8.2e-14 / 7.5e-14 | -1.3e-14 | 3.6e-13 | 20 |
| 20 | 0.95 | normal | exact | pass / pass | 1.43e-13 / 1.49e-13 | 2.3e-15 / 2.3e-15 | -1.3e-13 | 5.9e-14 | 4 |
| 20 | 0.95 | oblique | exact | pass / pass | 3.94e-14 / 4.21e-14 | 1.7e-15 / 1.8e-15 | -1.3e-14 | 4.7e-14 | 4 |
| 20 | 0.95 | normal | ring | fail / fail | 1.69e-13 / 1.57e-13 | 2.6e-13 / 2.5e-13 | -1.6e-13 | 5.9e-14 | 4 |
| 20 | 0.95 | oblique | ring | pass / pass | 4.37e-14 / 4.07e-14 | 4.7e-14 / 4.6e-14 | -1.2e-14 | 4.7e-14 | 4 |
| 40 | 0.2 | normal | exact | pass / pass | 1.16e-13 / 1.18e-13 | 4.5e-15 / 4.5e-15 | -1.1e-13 | 8.6e-13 | 40 |
| 40 | 0.2 | oblique | exact | pass / pass | 4.40e-13 / 4.39e-13 | 3.8e-15 / 4.4e-15 | -2.3e-14 | 1.9e-12 | 40 |
| 40 | 0.2 | normal | ring | pass / pass | 1.24e-13 / 1.24e-13 | 1.4e-13 / 1.4e-13 | -1.1e-13 | 8.6e-13 | 40 |
| 40 | 0.2 | oblique | ring | pass / pass | 4.72e-13 / 4.66e-13 | 6.5e-13 / 6.4e-13 | -1.8e-14 | 1.9e-12 | 40 |
| 40 | 0.95 | normal | exact | fail / fail | 3.04e-12 / 2.41e-12 | 4.2e-15 / 3.7e-15 | -1.9e-11 | 2.4e-13 | 8 |
| 40 | 0.95 | oblique | exact | pass / pass | 1.62e-13 / 1.60e-13 | 3.9e-15 / 4.0e-15 | -1.1e-13 | 4.9e-13 | 8 |
| 40 | 0.95 | normal | ring | fail / fail | 3.18e-12 / 3.82e-12 | 5.5e-12 / 6.4e-12 | -1.9e-11 | 2.4e-13 | 8 |
| 40 | 0.95 | oblique | ring | pass / pass | 1.66e-13 / 1.61e-13 | 1.7e-13 / 1.7e-13 | -8.7e-14 | 4.9e-13 | 8 |
| 80 | 0.2 | normal | exact | fail / fail | 1.28e-12 / 1.26e-12 | 8.7e-15 / 9.7e-15 | -1.1e-14 | 4.5e-12 | 80 |
| 80 | 0.2 | oblique | exact | fail / fail | 1.49e-12 / 1.51e-12 | 7.5e-15 / 7.5e-15 | -4.9e-14 | 1.1e-11 | 80 |
| 80 | 0.2 | normal | ring | fail / fail | 1.27e-12 / 1.27e-12 | 5.3e-13 / 5.2e-13 | -1.2e-14 | 4.5e-12 | 80 |
| 80 | 0.2 | oblique | ring | fail / fail | 1.53e-12 / 1.54e-12 | 1.6e-12 / 1.6e-12 | -4.9e-14 | 1.1e-11 | 80 |
| 80 | 0.95 | normal | exact | fail / fail | 1.38e-11 / 1.34e-11 | 1.1e-14 / 9.5e-15 | -1.2e-08 | 5.5e-12 | 17 |
| 80 | 0.95 | oblique | exact | pass / pass | 8.89e-13 / 8.87e-13 | 7.4e-15 / 7.3e-15 | -2.6e-14 | 6.1e-12 | 17 |
| 80 | 0.95 | normal | ring | fail / fail | 1.34e-11 / 1.39e-11 | 1.8e-11 / 1.5e-11 | -1.2e-08 | 5.5e-12 | 17 |
| 80 | 0.95 | oblique | ring | pass / pass | 8.92e-13 / 9.05e-13 | 7.6e-13 / 7.7e-13 | -1.8e-14 | 6.1e-12 | 17 |
| 160 | 0.2 | normal | exact | fail / fail | 8.71e-12 / 8.69e-12 | 1.8e-14 / 1.8e-14 | -5.2e-14 | 1.9e-11 | 160 |
| 160 | 0.2 | oblique | exact | fail / fail | 9.18e-12 / 9.22e-12 | 1.5e-14 / 1.5e-14 | -1.3e-14 | 6.3e-11 | 160 |
| 160 | 0.2 | normal | ring | fail / fail | 8.66e-12 / 8.68e-12 | 1.9e-12 / 1.9e-12 | -5.1e-14 | 1.9e-11 | 160 |
| 160 | 0.2 | oblique | ring | fail / fail | 9.51e-12 / 9.26e-12 | 8.3e-12 / 7.9e-12 | -1.6e-14 | 6.3e-11 | 160 |
| 160 | 0.95 | normal | exact | fail / fail | 2.87e-07 / 2.87e-07 | 1.8e-14 / 1.8e-14 | -6.1e-04 | 3.7e-07 | 34 |
| 160 | 0.95 | oblique | exact | fail / fail | 5.13e-12 / 4.68e-12 | 1.5e-14 / 1.5e-14 | -2.3e-14 | 3.9e-11 | 34 |
| 160 | 0.95 | normal | ring | fail / fail | 3.12e-07 / 3.12e-07 | 5.4e-07 / 5.4e-07 | -5.7e-04 | 3.7e-07 | 34 |
| 160 | 0.95 | oblique | ring | fail / fail | 5.16e-12 / 5.04e-12 | 4.7e-12 / 4.8e-12 | -1.8e-14 | 3.9e-11 | 34 |

32 of 32 verdicts equal to M9Ccs2's (15 pass, 17 fail on both sides: the
normal arm at target Co 0.95, OpenFOAM Co 1.53, the open growth (iii)-a,
which at N = 160 reaches L1(alphaf) 3e-7 and min alpha -6e-4 on both sides
alike; and the N = 80 / 160 cases at Co 0.2 and the oblique arm at N = 160,
Co 0.95, whose run maxima sit at 1.3e-12 to 9.5e-12, above the 1e-12 of
this rule on both sides). Every run maximum is a number of the same size as
M9Ccs2's (within -21 % to +21 %; 3e-15 to 1.4e-11 outside the (iii)-a
cases); the mass is conserved to 1e-12 %.

#### The timed cells (ring three repetitions, exact one; 14 of 18)

| cell | OF Co | row | E_shape (rel. to the comparison) | L1(alphaf) (rel.) | L1(x) (rel.) | L1(n) (rel.) | min / excess | wall no keys, median [min, max] s | wall comparison | comparison / no keys | step s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| x-line N 32, 1 rank | 0.2 | ring | 2.3293e-4 (+1.4e-10) | 1.4939e-2 (-1.5e-10) | 1.0264e-3 (+1.7e-3) | 4.2658e-2 (+1.5e-2) | -1.3e-14 / 2.1e-13 | 10.2 [10.2, 10.5] | 22.5 | 2.19 | 0.07 |
| | | exact | 1.2829e-4 (-2.9e-10) | 9.5951e-3 (+5.0e-10) | 7.3486e-4 (-3.8e-4) | 6.3177e-3 (+1.4e-1) | -1.7e-14 / 2.1e-13 | 10.8 | 23.4 | 2.17 | 0.07 |
| x-line N 32, 1 rank | 1 | ring | 2.3635e-7 (-1.2e-10) | 6.9899e-3 (-5.7e-13) | 3.9975e-4 (+8.6e-8) | 3.1516e-2 (+6.9e-9) | -7.3e-14 / 3.0e-13 | 4.3 [4.3, 4.5] | 7.5 | 1.73 | 0.07 |
| | | exact | 2.3641e-7 (+1.6e-10) | 6.6466e-3 (+1.4e-9) | 4.0244e-4 (+5.5e-10) | 1.5948e-3 (-8.9e-8) | -8.9e-14 / 3.5e-13 | 4.4 | 7.7 | 1.73 | 0.07 |
| x-line N 64, 8 ranks | 0.2 | ring | 6.1850e-5 (+2.8e-10) | 8.3941e-3 (+1.8e-10) | 2.7945e-4 (+8.5e-4) | 2.0956e-2 (+1.4e-2) | -4.2e-14 / 5.0e-12 | 36.7 [36.2, 37.2] | 84.2 | 2.29 | 0.15 |
| | | exact | 3.3711e-5 (-2.3e-10) | 5.3345e-3 (-7.0e-10) | 1.9291e-4 (-1.7e-3) | 1.7480e-3 (+3.6e-1) | -4.6e-14 / 5.0e-12 | 38.4 | 87.7 | 2.29 | 0.16 |
| x-line N 64, 8 ranks | 1 | ring | 2.3252e-7 (-4.2e-10) | 3.4847e-3 (-3.4e-12) | 9.5311e-5 (-1.1e-7) | 1.7740e-2 (+1.5e-8) | -1.5e-13 / 2.4e-12 | 10.6 [10.5, 11.2] | 20.3 | 1.91 | 0.15 |
| | | exact | 2.3252e-7 (-1.6e-10) | 3.2711e-3 (+3.7e-9) | 9.5627e-5 (-2.2e-7) | 5.7317e-4 (-1.9e-7) | -1.4e-13 / 2.4e-12 | 10.9 | 21.0 | 1.92 | 0.16 |
| x-line N 128, 64 ranks | 0.2 | ring | 1.7291e-5 (+2.9e-10) | 4.7810e-3 (+1.5e-9) | 7.6630e-5 (-3.4e-4) | 1.0102e-2 (+1.9e-3) | -8.8e-14 / 9.6e-12 | 179.4 [179.1, 179.6] | 434.3 | 2.42 | 0.38 |
| | | exact | 9.2505e-6 (-8.5e-10) | 3.0040e-3 (+1.1e-8) | 5.2724e-5 (+1.6e-4) | 3.8696e-4 (+7.6e-2) | -1.1e-13 / 9.6e-12 | 191.2 | 448.7 | 2.35 | 0.40 |
| x-line N 128, 64 ranks | 1 | ring | 2.3246e-7 (-2.8e-10) | 1.8815e-3 (+2.0e-9) | 2.5787e-5 (+5.5e-4) | 8.0368e-3 (+6.8e-4) | -8.6e-14 / 2.9e-11 | 44.5 [44.2, 44.9] | 94.9 | 2.13 | 0.41 |
| | | exact | 2.3246e-7 (-3.0e-10) | 1.7913e-3 (-6.7e-10) | 2.5761e-5 (+4.1e-7) | 1.9116e-4 (+4.8e-6) | -2.5e-13 / 2.9e-11 | 47.4 | 98.7 | 2.08 | 0.43 |
| diagonal N 32, 1 rank | 0.2 | ring | 5.7286e-4 (-2.5e-11) | 3.1749e-2 (-2.0e-9) | 2.0126e-3 (+1.7e-3) | 5.3562e-2 (+2.6e-6) | -2.0e-14 / 2.8e-13 | 16.1 [16.1, 16.6] | 42.4 | 2.64 | 0.08 |
| | | exact | 2.8940e-4 (-8.0e-10) | 1.6873e-2 (-1.1e-9) | 1.1550e-3 (+4.2e-7) | 2.3675e-3 (+6.8e-8) | -3.2e-13 / 2.0e-13 | 18.1 | 43.3 | 2.40 | 0.09 |
| diagonal N 32, 1 rank | 1 | ring | 1.8090e-4 (-5.5e-13) | 1.1178e-2 (-6.7e-11) | 7.8639e-4 (-8.6e-9) | 3.8703e-2 (+2.1e-8) | -6.0e-13 / 2.6e-13 | 5.7 [5.3, 6.3] | 10.5 | 1.84 | 0.08 |
| | | exact | 1.1818e-5 (-1.8e-9) | 4.8457e-3 (+4.0e-10) | 4.1866e-4 (+5.9e-10) | 2.3924e-3 (-3.0e-8) | -2.3e-14 / 1.4e-13 | 5.4 | 10.7 | 1.97 | 0.08 |
| diagonal N 64, 8 ranks | 0.2 | ring | 1.5684e-4 (+3.5e-11) | 1.7696e-2 (+4.0e-11) | 5.5372e-4 (+3.7e-3) | 2.4281e-2 (+8.5e-3) | -2.0e-11 / 5.3e-12 | 60.4 [60.4, 60.5] | 164.5 | 2.72 | 0.18 |
| | | exact | 8.5795e-5 (-1.0e-9) | 9.8911e-3 (-1.3e-9) | 3.2914e-4 (+2.1e-3) | 8.4231e-4 (+3.1e-3) | -1.2e-11 / 5.3e-12 | 63.4 | 169.1 | 2.67 | 0.19 |
| diagonal N 64, 8 ranks | 1 | ring | 3.8868e-5 (-1.8e-12) | 5.0451e-3 (-7.2e-11) | 1.7639e-4 (+9.3e-8) | 1.9075e-2 (+2.0e-8) | -1.2e-13 / 2.0e-12 | 15.9 [15.4, 17.1] | 36.4 | 2.29 | 0.18 |
| | | exact | 3.7843e-6 (-3.1e-9) | 2.3569e-3 (+6.4e-9) | 1.0070e-4 (+2.0e-7) | 7.2236e-4 (-2.2e-7) | -3.3e-11 / 1.5e-12 | 16.1 | 37.3 | 2.32 | 0.19 |
| diagonal N 128, 64 ranks | 0.2 | ring | 3.6888e-5 (-3.7e-11) | 8.8200e-3 (+4.6e-11) | 1.2827e-4 (-4.7e-7) | 9.7307e-3 (+2.9e-7) | -1.8e-11 / 1.2e-11 | 348.9 [347.6, 350.8] | 982.5 | 2.82 | 0.53 |
| | | exact | 2.4429e-5 (-8.7e-10) | 5.7134e-3 (-5.9e-10) | 9.0167e-5 (-6.8e-7) | 2.4190e-4 (+4.0e-7) | -2.4e-11 / 1.2e-11 | 369.6 | 1016.6 | 2.75 | 0.56 |
| diagonal N 128, 64 ranks | 1 | ring | 7.1839e-6 (+3.3e-9) | 2.1775e-3 (-9.4e-9) | 3.5490e-5 (-1.6e-7) | 8.2923e-3 (+8.7e-8) | -3.6e-11 / 1.3e-11 | 77.5 [77.3, 77.7] | 204.3 | 2.64 | 0.54 |
| | | exact | 8.7709e-7 (+4.9e-9) | 1.2420e-3 (+2.5e-9) | 2.6397e-5 (-2.5e-6) | 1.8660e-4 (+1.2e-5) | -1.4e-13 / 1.3e-11 | 82.2 | 211.2 | 2.57 | 0.58 |
| rotation N 32, 1 rank | 0.2 | ring | 1.0587e-3 (-1.6e-9) | 5.7035e-2 (-1.6e-9) | 3.9776e-3 (+5.1e-3) | 5.7649e-2 (+6.5e-3) | -6.9e-13 / 1.9e-11 | 164.4 [164.2, 164.4] | 438.1 | 2.67 | 0.09 |
| | | exact | 5.2607e-4 (-4.1e-9) | 2.9807e-2 (-3.0e-9) | 1.9624e-3 (+1.7e-3) | 3.1766e-3 (-1.5e-3) | -2.1e-13 / 2.1e-12 | 168.3 | 443.5 | 2.64 | 0.10 |
| rotation N 32, 1 rank | 1 | ring | 7.0707e-4 (-6.0e-7) | 3.8232e-2 (-5.2e-7) | 2.6449e-3 (+6.0e-3) | 4.8100e-2 (+5.5e-3) | -9.2e-13 / 6.0e-13 | 34.8 [34.8, 35.1] | 90.5 | 2.60 | 0.09 |
| | | exact | 3.7576e-4 (-1.4e-6) | 2.1590e-2 (-1.3e-6) | 1.4072e-3 (+6.2e-3) | 3.5835e-3 (+4.3e-1) | -7.1e-13 / 6.1e-13 | 35.6 | 91.4 | 2.57 | 0.09 |
| rotation N 64, 8 ranks | 0.2 | ring | 2.9561e-4 (-3.7e-10) | 3.2196e-2 (-2.0e-10) | 1.0577e-3 (+3.5e-4) | 2.5174e-2 (+1.4e-3) | -6.2e-14 / 6.2e-10 | 705.1 [704.2, 706.1] | 1875.0 | 2.66 | 0.20 |
| | | exact | 1.7175e-4 (-2.0e-9) | 1.9290e-2 (-1.9e-9) | 6.2757e-4 (+3.1e-4) | 8.4615e-4 (+3.6e-2) | -7.7e-13 / 1.0e-10 | 720.4 | 1895.0 | 2.63 | 0.20 |
| rotation N 64, 8 ranks | 1 | ring | 1.9738e-4 (-1.3e-7) | 2.1912e-2 (-1.4e-7) | 7.1023e-4 (+1.4e-3) | 2.2567e-2 (+2.8e-3) | -5.3e-13 / 6.1e-12 | 145.3 [144.7, 145.5] | 380.0 | 2.62 | 0.20 |
| | | exact | 1.2085e-4 (-3.6e-7) | 1.3792e-2 (-3.5e-7) | 4.4486e-4 (-5.1e-7) | 7.5390e-4 (-1.6e-6) | -5.6e-12 / 5.3e-12 | 147.9 | 383.7 | 2.60 | 0.20 |
| rotation N 128, 64 ranks | 0.2 | ring | 8.0418e-5 (+6.1e-10) | 1.8297e-2 (+5.6e-10) | 2.8343e-4 (+1.8e-3) | 1.1299e-2 (+2.7e-3) | -1.8e-13 / 7.4e-10 | 3541.6 [3538.6, 3551.4] | 9370.9 | 2.65 | 0.47 |
| | | exact | 5.2188e-5 (-3.8e-10) | 1.1923e-2 (+1.3e-10) | 1.8746e-4 (+1.2e-4) | 2.6495e-4 (+3.7e-2) | -7.9e-14 / 5.1e-10 | 3622.8 | 9473.1 | 2.61 | 0.47 |
| rotation N 128, 64 ranks | 1 | ring | 5.3417e-5 (-5.2e-8) | 1.2274e-2 (-3.2e-8) | 1.8942e-4 (+1.7e-3) | 9.9009e-3 (+2.1e-3) | -5.4e-11 / 7.2e-11 | 727.8 [717.6, 851.4] | 1882.7 | 2.59 | 0.47 |
| | | exact | 3.6197e-5 (-1.1e-7) | 8.3619e-3 (-6.9e-8) | 1.3028e-4 (-3.0e-7) | 2.3922e-4 (-5.1e-6) | -6.6e-11 / 1.4e-10 | 739.1 | 1905.8 | 2.58 | 0.48 |

The mass deviation is within 9.6e-12 % in every cell. The three
repetitions give the same norms to every printed digit. The split cells'
repetitions ran on separate i01 nodes (the rotation N 128 at Co 1: 718,
728 and 851 s, the last on a node of its own at the same time as the other
three; the comparison's spread at that cell was 1882 to 1883 s on one node).

#### The pure spin (sphereRotation3D_hexSpinCoPyr, exact row) and the parallel identity

| N | target Co | OF Co | E_shape (rel. to hexSpinCoPyr) | L1(alphaf) (rel.) | L1(x) (rel.) | L1(n) (rel.) | min / excess | mass % | step s |
|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.2 | 0.990 | 2.5486e-5 (-2.9e-9) | 6.1976e-3 (+2.1e-10) | 4.7471e-4 (+1.4e-7) | 7.4966e-4 (+3.3e-7) | 0 / 9.3e-13 | 1.0e-12 | 0.26 |
| 32 | 0.5 | 2.474 | 2.1231e-5 (+1.1e-9) | 4.5104e-3 (+8.2e-10) | 4.5696e-4 (-2.0e-7) | 1.1874e-3 (-2.0e-7) | 0 / 5.5e-13 | 3.7e-13 | 0.26 |
| 64 | 0.2 | 0.997 | 7.0677e-6 (+2.6e-9) | 3.2924e-3 (+1.8e-9) | 1.2361e-4 (+3.0e-6) | 2.8416e-4 (-1.4e-5) | 0 / 1.8e-11 | 2.1e-13 | 1.52 |
| 64 | 0.5 | 2.494 | 6.3340e-6 (+2.1e-9) | 2.4365e-3 (+2.9e-9) | 1.1542e-4 (+4.6e-7) | 4.3037e-4 (+1.9e-6) | 0 / 4.4e-12 | 1.1e-12 | 1.52 |

| row | np | step 1 max dalpha / dalphaf | max over 20 steps dalpha / dalphaf | nMissing / nDup | rank fallbacks |
|---|---|---|---|---|---|
| exact | 2 | 0 / 0 | 0 / 0 | 0 / 0 | 0 |
| exact | 4 | 2.1e-14 / 5.7e-14 | 1.3e-12 / 4.4e-9 | 0 / 0 | 0 |
| exact | 8 | 2.0e-14 / 9.9e-14 | 6.6e-13 / 4.2e-10 | 0 / 0 | 0 |
| exact | 16 | 3.4e-14 / 1.1e-13 | 2.3e-12 / 1.6e-9 | 0 / 0 | 0 |
| ring | 2 | 2.9e-14 / 5.7e-14 | 6.2e-12 / 8.3e-10 | 0 / 0 | 0 |
| ring | 4 | 3.2e-14 / 5.7e-14 | 9.4e-12 / 1.4e-9 | 0 / 0 | 0 |
| ring | 8 | 3.0e-14 / 1.1e-13 | 4.7e-12 / 7.4e-10 | 0 / 0 | 0 |
| ring | 16 | 2.9e-14 / 5.8e-14 | 9.5e-12 / 1.4e-9 | 0 / 0 | 0 |

The two-rank exact row is bit for bit the serial run at every step (the
exact normal needs no stencil across the processor boundary and no face of
this decomposition takes a rank fallback); every other decomposition stays
within 1e-11 in alpha over the 20 steps, 50 to 2500 times below the bands
of the default family (6e-10 / 2.5e-8), and the np 16 step-20 discontinuity
of P-P2 does not return (no rank fallback on any rank).

#### Verdicts

- D-1 HELD for the plane gate (32 of 32 verdicts equal, the maxima
  numbers of the same size), for E_shape and L1(alphaf) within 1e-8
  everywhere except the rotation at OpenFOAM Co 1 (ring -6.0e-7 / -5.2e-7
  at N 32, -1.3e-7 at N 64, -5.2e-8 at N 128; exact -1.4e-6 / -1.3e-6 at N
  32, -3.6e-7 at N 64, -1.1e-7 at N 128): the pyramid against the scaling
  at OpenFOAM Co 1, the band of 1e-6 exceeded by the exact row at N 32
  (1.4e-6). FALSIFIED there, narrowly. L1(x) within 0.62 % everywhere
  (HELD). L1(n) at the end: FALSIFIED on both rows: the ring row +1.4 to
  +1.5 % on the x-line at OpenFOAM Co 0.2 (N 32 / 64), within 0.9 %
  elsewhere; the exact row +43 % on the rotation at N 32, Co 1, +36 % on
  the x-line at N 64, Co 0.2, +14 % at N 32, +7.6 % at N 128, +3.6 % on the
  rotation at N 64, Co 0.2, within 0.3 % at Co 1 on the lines. The
  positioning section measured +8 to +13 % at N 32 with the anchor active;
  with the anchor skipped (anchorSkipConverged) the corner slivers that ccs
  gives a plane carry their own plane into the end-state norm. The
  transported fields do not see it (E_shape, L1(alphaf) within 1.4e-6).
  The spin HELD (E_shape and L1(alphaf) within 3e-9, L1(x) within 3e-6,
  L1(n) within 1.4e-5).
- D-2: min alpha >= -6.6e-11 and the mass within 9.6e-12 % in every
  finished cell, the spin's min alpha 0; the excess <= 2.9e-11 at OpenFOAM
  Co 1 except the rotation N 128 exact (1.4e-10, the comparison's 1.3e-10),
  and 6.2e-10 (ring) / 1.0e-10 (exact) on the rotation N 64 at Co 0.2
  (the comparison's 6.0e-10 / 2.3e-10), 7.4e-10 / 5.1e-10 at N 128 (the
  comparison's 7.5e-10 / 5.1e-10). FALSIFIED as written (the band 1.5e-10
  was taken from the Co 1 cells); the values equal the opt-in runs' within
  round-off.
- D-3 FALSIFIED on both sides. The ratio of the comparison's wall clock to
  the rerun's (the ring row, medians) is 1.73 to 2.72: below 2.0 in the
  cells whose runs last 4 to 11 s (the x-line at Co 1 at N 32 / 64, the
  diagonal at Co 1 at N 32: 1.73 / 1.91 / 1.84, where the start-up, the
  mesh read and the end write are the same seconds on both sides), 2.1 to
  2.4 on the x-line otherwise, and 2.6 to 2.7 on the rotation and the
  diagonal (the rotation N 32: 2.67 / 2.60; N 64 on 8 ranks: 2.66 / 2.62;
  N 128 on 64 ranks: 2.65 / 2.59; the diagonal N 64 at Co 0.2: 2.72, N
  128: 2.82 / 2.64). Per step: the rotation N 64 on 8 ranks 0.52 -> 0.20 s
  (2.6 x; the efficiency work's scotch cell gave 0.545 -> 0.236 with
  touching + trace cache alone), N 128 on 64 ranks 1.21 -> 0.47 s (one
  turn 9371 -> 3542 s), N 32 serial 0.25 -> 0.09 s. The exact row's wall clock is within 7 % of the
  ring row's, as in the comparison.
- D-4 HELD (above).
- D-5: the serial N = 64 step 0.491 s (predicted 0.50 to 0.55: below the
  band by 2 %), the two clips 0.101 s/step (HELD); 8 ranks 0.207 s/step
  (predicted 0.23 to 0.26) and 64 ranks 0.508 s/step (predicted 0.55 to
  0.62): both below the band; the MPI share 64.5 % at 8 ranks and 84.6 % at
  64 (predicted above 70 / 90 %): FALSIFIED, the imbalance eased (below).

#### The HPCToolkit profile of the reference setting, side by side with the first

The same build flags (Opt, -O3, plus -g -fno-omit-frame-pointer in the
profiling worktree prof/TwoPhaseFlow at 407f969), the same events (CPUTIME
every 1000 us, LLC misses every 100000), the same cases (the comparison's
rotation at OpenFOAM Co 1, the ring row, the comparison's scotch
decompositions) with the thirteen key lines deleted, the same readers
(hpdb.py, hpan.py; the categories rules_swept3.py, the line ranges of
407f969). hpcrun's overhead: 1.1 % (serial N 32: 35.2 against 35.5 s; N 64,
20 steps: 25.6 against 25.8 s). plicRDF was not profiled again: none of the
nine keys is read by plicRDF + isoAdvector (its numbers stand).

(a) The serial breakdown per step, N = 64, one turn (697 steps):

| category | before (720cfc6) s/step | % | after (407f969) s/step | % | after / before |
|---|---|---|---|---|---|
| pieces: the clip by the cell planes and the PLIC plane | 0.4685 | 42.2 | 0.0660 | 13.4 | 0.14 |
| pieces: the cap closure and the volume | 0.0032 | 0.3 | 0.0012 | 0.2 | 0.38 |
| face update: the pre-image fan clipped (sweptPieceArea) | 0.1819 | 16.4 | 0.0353 | 7.2 | 0.19 |
| swept cells: the ring search | 0.0412 | 3.7 | 0.0127 | 2.6 | 0.31 |
| swept cells: the cell planes, state, wisp rule | 0.0302 | 2.7 | 0.0291 | 5.9 | 0.96 |
| swept cells: the pre-test | - | - | 0.0372 | 7.6 | - |
| swept cells: the shortcut and the check (both off) | - | - | 0.0013 | 0.3 | - |
| swept cells: the pre-image fan and the rest | 0.0079 | 0.7 | 0.0122 | 2.5 | 1.54 |
| trace: the velocity at the traced points (findTetrahedron) | 0.0743 | 6.7 | 0.0302 | 6.2 | 0.41 |
| trace: the cell search (ring, octree) | 0.0174 | 1.6 | 0.0153 | 3.1 | 0.88 |
| trace: the midpoint trace, the triangle list | 0.0104 | 0.9 | 0.0107 | 2.2 | 1.03 |
| the point velocity (volPointInterpolation of U, the whole mesh) | 0.0222 | 2.0 | 0.0223 | 4.5 | 1.00 |
| curvedFlux: the band, the wisp surfaces, the exchanges, the face loop | 0.0212 | 1.9 | 0.0226 | 4.6 | 1.07 |
| the remote exchange (serial: the stencil setup) | 0.0030 | 0.3 | 0.0031 | 0.6 | 1.03 |
| the surface build (fit, anchor, volume rule) | 0.0451 | 4.1 | 0.0155 | 3.2 | 0.34 |
| the reconstruction (davof, ring normal, positioning) | 0.0365 | 3.3 | 0.0299 | 6.1 | 0.82 |
| the bounding | 0.0102 | 0.9 | 0.0113 | 2.3 | 1.11 |
| isoAdvector: timeIntegratedFlux | 0.0147 | 1.3 | 0.0150 | 3.0 | 1.02 |
| isoAdvector: advect, the rest (alpha update, counters) | 0.0293 | 2.6 | 0.0291 | 5.9 | 0.99 |
| setFlow | 0.0309 | 2.8 | 0.0309 | 6.3 | 1.00 |
| the Courant number, other function objects | 0.0015 | 0.1 | 0.0016 | 0.3 | 1.07 |
| I/O and start-up | 0.0243 | 2.2 | 0.0218 | 4.4 | 0.90 |
| the solver loop (main) | 0.0355 | 3.2 | 0.0355 | 7.2 | 1.00 |
| total | 1.1105 | 100 | 0.4907 | 100 | 0.44 |

- The swept advection's own parts fell from 0.926 to 0.315 s/step (64 % of
  the step, 83 % before); the two clips from 0.650 to 0.101 (0.139 with the
  pre-test that replaces them); the trace's velocity lookup from 0.074 to
  0.030 (the cache: each mesh point once); the surface build from 0.045 to
  0.016 and the reconstruction from 0.037 to 0.030 (ccs with the anchor
  skipped: 0.036 s/step, 3 % of the old step, 7 % of the new). The parts
  shared with plicRDF (the reconstruction, the bounding, timeIntegratedFlux,
  the alpha update, setFlow, the Courant number, I/O, the solver loop) are
  unchanged at 0.175 s/step and are now 36 % of the step. N = 32: 0.261 ->
  0.099 s/step (0.38).
- The ratio to plicRDF's serial step (0.209 s, the first profile) is 2.3
  (5.3 before).

(b) The scaling with N after (the full turns at N = 32 and 64; the N = 128
serial window was not rerun: the first profile's x 3.8 to 4.7 for the swept
parts stands): the clips x 3.82 / 3.88, the pre-test x 3.76, the ring search
x 3.63, the trace's velocity lookup x 3.97 and its cell search x 4.37, the
surface build x 3.88, the reconstruction x 6.5, the mesh-wide parts x 8.6 to
9.7 (the point velocity, the alpha update, setFlow). As before, no swept part
grows faster than the interface-cell count (x 4.00).

(c) Parallel, N = 64 on 8 ranks (697 steps) and N = 128 on 64 ranks (1401
steps), one turn each, the comparison's scotch decompositions:

| | N 64, 8 ranks before | after | N 128, 64 ranks before | after |
|---|---|---|---|---|
| step (CPUTIME per rank / steps; the wall clock) | 0.545 s (0.551) | 0.207 s (0.212) | 1.341 s (1.358) | 0.508 s (0.519) |
| MPI library, share of all CPUTIME | 72.1 % | 64.5 % | 91.2 % | 84.6 % |
| MPI per rank: mean / max / min | 274 / 359 / 72 s | 93 / 123 / 27 s | 1715 / 1841 / 929 s | 603 / 674 / 332 s |
| compute per rank: max / mean (max/mean) | 303 / 104 (2.92) | 113 / 49 (2.31) | 931 / 156 (5.98) | 361 / 101 (3.59) |
| compute per rank, s | 235 59 303 58 43 27 85 20 | 96 34 113 34 28 24 43 20 | 931 758 692 691 483 419 385 359 on eight ranks, 42 ranks at 37 to 80 | 361 324 300 286 228 194 182 181 on the same eight ranks (27, 7, 28, 15, 29, 26, 9, 6), 35 ranks at 37 to 80, none below 37 |
| ranks with compute < 10 % of the max | 2 of 8 | 0 of 8 | 42 of 64 | 0 of 64 |
| the balanced step (the sum of compute / ranks / steps) | 0.149 s: the step is 3.6 x | 0.070 s: the step is 3.0 x | 0.111 s: the step is 12 x | 0.072 s: the step is 7.1 x |
| the busiest rank: pieces + face update + pre-test + trace velocity, s/step | 0.204 + 0.083 + 0 + 0.031 | 0.029 + 0.015 + 0.016 + 0.012 | 0.299 + 0.129 + 0 + 0.046 | 0.041 + 0.022 + 0.026 + 0.017 |
| MPI by call site, share of the MPI time | syncProcPatches 90.8 %, sweptExchange 6.4 %, davof::reconstruct 1.5 % | - | syncProcPatches 90.1 %, sweptExchange 6.2 %, davof::reconstruct 1.6 % | syncProcPatches 81.9 %, sweptExchange 9.8 %, davof::reconstruct 2.7 %, the alpha update 1.9 % |
| the exchange proper (pack, map, unpack) | 0.0049 s/step/rank | 0.0049 | 0.0058 | 0.0057 |

The interface ranks' compute fell by the serial factor (8 ranks: rank 2
303 -> 113 s; 64 ranks: rank 27 931 -> 361 s, 2.6 to 2.7 x), the ranks
without interface not (rank 7 of 8: 20 -> 20 s; the 35 to 42 ranks at 37 to
80 s: the mesh-wide parts), so the imbalance eased from 2.92 to 2.31 and
from 5.98 to 3.59, the MPI wait from 72 to 65 % and from 91 to 85 %; the
step is still 3.0 and 7.1 times the balanced step. The decomposition, not
the communication, stays the parallel bottleneck (the exchange proper 2.4 /
1.1 % of the step; its wait, now 9.8 % of the MPI time at 64 ranks, is the
idle ranks reaching the mapDistribute before the busy ones).

(d) The top exclusive hotspots after (serial, N = 64, one turn; 342 s),
with the kind:

| # | % | s | function (file:line; fork 407f969) | kind |
|---|---|---|---|---|
| 1 | 7.0 | 23.8 | Foam::minusOp (ops.H:209): the vector subtraction x - x0 of planeDistances (davofSweptPieces.H:189) and of the pre-test (isoAdvectionSwept.C:966-983) | algorithmic (vertices x planes) |
| 2 | 4.7 | 16.0 | Foam::Vector::inner (VectorI.H:157): the signed distances of the pre-test and the clip | algorithmic |
| 3 | 4.3 | 14.9 | sweptFaceCellsFlux's own lines (970: the pre-test's loop over the vertices, 3.7 s; 903-905: the cell planes, 2.6 s) | algorithmic |
| 4 | 4.2 | 14.2 | sweptPieceVolume's own lines (405-470: the per-triangle loop, the cap bookkeeping) | algorithmic |
| 5 | 4.0 | 13.6 | volPointInterpolation::interpolateInternalField (libfiniteVolume): the point velocity of the whole mesh | mesh-wide, redundant (the band's points would do) |
| 6 | 3.4 | 11.5 | __mempcpy (libc): GeometricField::operator== 3.6 s and operator= 2.8 s (the solver loop's and setFlow's field assignments), the pointField construction 1.8 s | redundant copies, mesh-wide |
| 7 | 2.7 | 9.2 | davof::clipWith (davofSweptPieces.H:216-229): the clip proper | algorithmic |
| 8 | 2.6 (6.3 with its callees) | 9.0 (21.6) | cellPointWeight::findTetrahedron (libfiniteVolume): the velocity at a traced point; cellTetIndices and faceTetIndices and their lists under it | redundant + allocation (the cell's tetrahedra rebuilt per lookup; once per point now, but still per cell visit) |
| 9 | 1.9 (5.3 with magSqr, mag) | 6.4 (18.2) | plusEqMagSqrOp2, magSqr, mag: the edge lengths and the coplanar tests | algorithmic, in part redundant (L per polygon recomputed) |
| 10 | 1.8 | 6.2 | fvc::surfaceIntegrate (the alpha update) | mesh-wide |
| 11 | 1.8 | 6.2 | Foam::Vector::cross (VectorI.H:169): areaVector, the triangle normals | algorithmic |
| 12 | 1.6 | 5.3 | surfaceInterpolationScheme::dotInterpolate (libfiniteVolume): phi from U (setFlow, the Courant number) | mesh-wide |

Then: inCellTol (isoAdvectionSwept.C:221-222, the cell search) 3.2 s,
planeDistances' own lines 2.7 s, the curvedFlux face loop 3.3 s,
kbIntegrateFace (the surface build) 1.0 s; the allocation (malloc, free
and their callees) 13.7 s = 4.0 % of the step (20.3 s = 2.6 % before: the
trace's tetrahedron lists and the DynamicLists of sweptFaceCellsFlux). The
line table also lists 22.4 s (6.5 %) at a PLT stub of IOstream::fatalCheck
in libOpenFOAM: an attribution artefact of OpenFOAM's libraries built
without -g (the samples sit in their callers' categories above).

The next changes with the largest expected gain (not implemented):
1. The pre-test on a hex mesh (serial, 7.6 % + the ring search 2.6 %): the
   cell planes of a hex are axis-aligned, so the separating-plane test
   reduces to the overlap of the swept volume's vertex box with the cell's
   box (six comparisons instead of 60 dot products and subtractions:
   hotspots 1 to 3 carry most of the pre-test's 0.037 s/step); for a
   general polyhedron the box test first and the planes only for the boxes
   that overlap. Expected: the pre-test to ~ 0.01 s/step, the step from
   0.49 to ~ 0.46 s (6 %).
2. The point velocity of the band only (mesh-wide, 4.5 %): volPointInterpolation
   of U over the whole mesh (hotspot 5, x 8.6 from N = 32 to 64) where only
   the points of the swept faces are read; interpolating the band's points
   (the vertex ring of the active cells, ~ 5 % of the points at N = 64)
   removes ~ 4 % of the step and grows with N like the interface.
3. The cell's tetrahedra once per cell per step (serial, 6.2 % + 3.1 %):
   findTetrahedron rebuilds the decomposition per lookup (hotspot 8); a
   per-cell cache of the tetrahedra (or the cell-centre + face-centre
   decomposition of the cell kept per step for the band's cells) and the
   cell search started from the traced point's own cell would take the
   trace's 0.046 s/step to ~ 0.02 (5 %).
4. In parallel the imbalance remains the loss (the step 3.0 x the balanced
   step at 8 ranks): the static interface weighting and the renewed
   decomposition measured no net gain in the efficiency work (the renewal's
   13.9 s at N 128 against the gain per step); a cheaper renewal
   (redistributePar's cost is the mesh rebuild) or a two-level scheme (the
   band's faces distributed apart from the mesh) is the next parallel item.

What this demonstrates: with no key in the dictionary the fork at 407f969
runs the measured best configuration (the final wisp rule, the consistent
variant, the pyramid, ccs with the anchor skipped, the pre-test, the trace
cache): the plane gate's verdicts and round-off maxima, the spin, the
parallel identity and the transported fields of the rotation and the lines
(E_shape and L1(alphaf) within 1.4e-6, L1(x) within 0.62 %) equal the
opt-in runs', bounded (min alpha -6.6e-11, the excess as the comparison's,
6.2e-10 at most) and conservative to 1e-11 %; the wall clock at equal ranks
is 2.1 to 2.8 times lower than the comparison's in every cell longer than
15 s (1.7 to 2.0 in the 4 to 11 s cells); the serial N = 64 step is 0.49 s
(2.3 times plicRDF's, 5.3 before) with the swept parts at 64 % of it and
the mesh-wide parts shared with plicRDF at 36 %; the 8-rank step 0.21 s and
the 64-rank step 0.51 s with the imbalance eased to 2.3 and 3.6 (the MPI
wait 65 and 85 %). The limitation: the end-state L1(n) of the exact row is
up to 43 % above the comparison's (the planes of the corner slivers under
ccs with the anchor skipped), the ring row's up to 1.5 %.


## Profile round 3: the flux against isoAdvector's, the two structural savings, the flux isolated (2026-11-05)

Tomislav (2026-11-05, via the coordinating session): more speed-up from the
closed forms, more HPCToolkit benchmarking of the swept-flux advection; the
decomposition deferred. Reference: the reference setting 407f969 and its
profile (prof/runs2). Three parts: A, the flux against isoAdvector's
function by function with in-code counters per face class; B, the two
structural savings the profile points at (the donor fast path, the face
update from the volume clip), each verified bit for bit on the eight cases;
C, the flux accuracy isolated (the two loops with the same exact normal).

### C: the flux accuracy isolated (PRE-REGISTERED 2026-11-05, 11:30, before the runs)

The preprint's plane-clipped isoAdvector loop (the rows davof-exact and
davof-ringArea: fluxIntegration prism, clippingSurface none, fitAnchor
faces; the reconstruction at the reference setting, ccs with the anchor
skipped, as the swept loop's) on the comparison's cells: the rotation (one
turn) and the sphere along the x-z diagonal (END 0.45), N 32 / 64 / 128 on
1 / 8 / 64 ranks, OpenFOAM Co 1 and 0.2, three repetitions, from the
comparison's rendered cases with the three keys edited and the
reference-setting keys stripped (workflow/scripts/davof_prof/ref-run.sbatch,
MODE cmpC); scored as the comparison (L1(x), E_shape, L1(alphaf), L1(n),
the bounds, the mass, the wall clock at equal ranks). The swept rows are the
reference rerun's (the exact row one repetition, the ring row three) and
plicRDF's the comparison's.
- F-1 (the flux's share of the exact-normal advantage): at OpenFOAM Co 0.2
  the plane-clipped loop with the exact normal has L1(x) and E_shape within
  10 % of the swept loop's on the rotation and the diagonal at every N, so
  the flux's share of the exact-normal rows' advantage over plicRDF
  ((L1(x)_plane-exact - L1(x)_swept-exact) / (L1(x)_plicRDF -
  L1(x)_swept-exact)) is below 0.15: the normal carries the advantage (the
  fair rows at OpenFOAM Co 0.126 had Ec within 4 % between the two loops
  with the ring normal). FALSIFIED by a share above 0.15 in any cell.
- F-2 (OpenFOAM Co 1): the plane-clipped loop (isoAdvector's flux with the
  exact plane) loses boundedness on the diagonal in parallel at N 64 / 128
  (min alpha below -1e-6) as plicRDF + isoAdvector did (C-4), and its L1(x)
  there exceeds the swept loop's by more than 2 x; on the rotation it stays
  bounded (min alpha >= -1e-9) and within 20 % of the swept loop's L1(x).
- F-3 (the ring-normal rows at Co 1): the two loops with the ring normal
  within 10 % in L1(x) and E_shape on the rotation at Co 0.2 and 1; at Co 1
  on the diagonal the plane-clipped loop unbounded in parallel as in F-2.
- F-4 (the cost): the plane-clipped loop's wall clock at equal ranks 1.2 to
  2 times plicRDF's, i.e. 1.3 to 2.2 times below the swept loop's.

### A: the flux against isoAdvector's, function by function (2026-11-05)

Fork: da04ece (key sweptClassStats, default false: per face class and step
the faces, the touched cells, the pieces, the half-space clips of the volume
clip and of the face update, the triangles entered, and the wall time in
sweptPieceVolume, sweptPieceArea, the trace and the first pass, by
steady_clock around the calls; pieceResult carries nClip and nTri) and
c453c47 (a fourth class, the faces without a flux). Counters only: bit for
bit on the eight cases with the key on (jobs 55964204, 55965249: 11
identical field files per case, max |diff| 0). The HPCToolkit profile of
the reference setting (prof/runs2) and the first profile's plicRDF +
isoAdvector databases (prof/runs/serial/rdf64_full, rdf32_full; the
plicRDF loop reads none of the nine keys) re-attributed alike with
rules_cmp.py (prof/an3/{ring64,rdf64,ring32,rdf32}.*): the reconstruction
split into the plane positioning (surfaceIteratorPLIC / cutCellPLIC, the
same code in both; ours ccs on the same cut) and the normal, the flux into
its sub-parts, the point velocity and the stencil machinery named alike.
The in-code counters: the full turn of the rotation N 64 (job 55965249,
697 steps, the key on: 0.4986 s/step against 0.4847 without it, the timers
3 %) and N 32 (55965250).

The rotation N 64 serial, OpenFOAM Co 1, one turn, s/step (N 32 in
brackets):

| part | swept (ours) | plicRDF + isoAdvector |
|---|---|---|
| reconstruction: the plane positioning (vofCutCell / ccs on cutCellPLIC) | 0.0072 (0.0017) | 0.0260 (0.0063) |
| reconstruction: the normal (ours the ring normal; theirs the RDF, its interpolation, gradient, residual, the initial normals) | 0.0002 | 0.0341 (0.0043) |
| reconstruction: the rest (own lines, field algebra) | 0.0225 (0.0029) | 0.0084 (0.0007) |
| ours: the surface build (fit, anchor, volume rule) | 0.0155 (0.0040) | - |
| reconstruction total | 0.0454 (0.0086) | 0.0684 (0.0113) |
| flux: the trace (the cell search 0.0153, the velocity at the traced points 0.0302, the midpoint trace and the list 0.0107) / theirs the velocity at the face points | 0.0562 (0.0138) | 0.0019 (0.0004) |
| flux: the ring search 0.0127, the cell planes, state and wisp rule 0.0291, the pre-test 0.0369 | 0.0787 (0.0206) | - |
| flux: the clips (the pieces 0.0672, the face update 0.0353) / theirs the plane sweep per face (cutFaceAdvect) | 0.1025 (0.0267) | 0.0019 (0.0004) |
| flux: the fan, sweptFaceCellsFlux other, the band and the face loop of curvedFlux / theirs timeIntegratedFlux other | 0.0364 (0.0070) | 0.0032 (0.0002) |
| the point velocity (volPointInterpolation of U, the whole mesh; both) | 0.0223 (0.0026) | 0.0113 (0.0013) |
| the stencil machinery (zoneDistribute, zoneCPCStencil, mapDistribute; both) | 0.0089 (0.0022) | 0.0171 (0.0042) |
| isoAdvector's defaultAreaFractions (ours: the faces without a plane) | 0.0144 (0.0017) | 0 |
| flux total | 0.3200 (0.0746) | 0.0354 (0.0065) |
| bounding | 0.0113 (0.0004) | 0.0077 (0.0006) |
| shared mesh-wide (the alpha update 0.0291 / 0.0248, setFlow 0.0309 / 0.0306, the Courant number, I/O 0.0218 / 0.0105, main 0.0297 / 0.0293) | 0.1131 (0.0150) | 0.0967 (0.0116) |
| total | 0.4907 (0.0989) | 0.2090 (0.0303) |

Where the extra 0.28 s/step sits: in the flux, +0.285: the two clips
0.103, the pre-test with the ring search and the cell states 0.079, the
trace 0.056, the band, fan and face loop 0.036, defaultAreaFractions
+0.014, the point velocity +0.011, the stencil machinery -0.008 (plicRDF's
RDF stencil costs more than our exchange's). The reconstruction is cheaper
by 0.023 (the ring normal against the RDF, ccs against vofCutCell), the
bounding +0.004, I/O +0.011 (more fields written), the rest +0.005.

Per face class (the in-code counters, N 64, means per step over the turn):

| class | faces/step | touched cells/face | pieces/face | volume clips/face | area clips/face | triangles/face | t volume | t area | t trace | t pre-test (ring, states) | sum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| no-sweep (phi_f = 0: the face update alone) | 5856 (33 %) | 8.17 | 0 | 0 | 165.5 | 0 | 0 | 0.0241 | 0.0173 | 0.0229 | 0.0642 |
| (i) donor-only, sweeping | 0 | - | - | - | - | - | 0 | 0 | 0 | 0 | 0 |
| (ii) multi-cell, sweeping | 11718 (67 %) | 2.04 | 2.04 | 140.3 | 50.2 | 32.7 | 0.0681 | 0.0137 | 0.0528 | 0.0418 | 0.1764 |
| (iii) boundary | 0 | - | - | - | - | - | - | - | - | - | - |

N 32 the same shares (no-sweep 1479 faces, 8.35 cells, 167 area clips,
27 % of the four timed parts; multi-cell 2956 faces, 2.09 cells; donor-only
0). The four timed parts sum to 0.241 s/step of the 0.499; the rest of the
swept parts (the point velocity, the band, the fan, the surface build) and
the shared parts make the step.

What this shows:
- Under the rotation at OpenFOAM Co 1 class (i) is empty: every sweeping
  face's traced vertices leave the donor across a lateral plane (by |u_t|
  dt, the tangential velocity), so a face touches the donor and one
  neighbour, and both hold a piece (2.04 touched = 2.04 pieces: the
  pre-test is exact on the sweeping faces, no clip without a piece). The
  fraction of the flux time on class (i) faces is 0 %.
- One third of the swept faces carry no flux (the faces normal to the
  rotation axis), and they cost 0.064 s/step (13 % of the step, 27 % of
  the four timed parts) for a face fraction without a flux: their
  pre-image lies in the face plane, the pre-test's coplanar exception
  keeps the four cells across that plane which the face update's coplanar
  rule then drops, so each touches 8.2 cells and runs 166 area clips.
- The sweeping faces spend 0.068 in the volume clip (140 clips over 33
  triangles per face), 0.053 in the trace, 0.042 in the pre-test and the
  cell states, 0.014 in the face update.

### B: the structural savings (2026-11-05)

Three changes, each behind a key (default false), each verified on the
eight cases against the reference setting (ref407full: 407f969, no keys;
the rotation N 32 / 64, the plane N 40 exact normal, the diagonal N 32,
the plane gate N 20 / 40 both arms, 20 steps, writePrecision 17), and timed
on the rotation N 64 serial at OpenFOAM Co 1 (200 steps, three
repetitions, the median of the per-step mean; the reference 0.4847 s/step,
0.4805 to 0.4856, the same build with no key: job 55964206).

| change | commit | key | identity on the eight cases | s/step (median of 3) | against the reference |
|---|---|---|---|---|---|
| B.1 the donor fast path | cd9385b, 8933403 | sweptDonorFastPath | bit for bit (55965251: 11 identical field files per case) | 0.4847 (0.4843 to 0.4870) | 1.000 (0 faces take it) |
| B.3 the pre-test's coplanar side | a11ebc1 | sweptPreTestCoplanarSide | bit for bit (55965847) | 0.4622 (0.4622 to 0.4636) | 0.954 |
| B.2 the face update from the volume clip | a11ebc1 | sweptFaceUpdateFromVolume | round-off (55965848; below) | 0.4688 (0.4678 to 0.4712) | 0.967 |
| B.1 + B.2 + B.3 | | the three keys | round-off (B.2's) | 0.4531 (0.4505 to 0.4541) | 0.935 |

B.1, the donor fast path: a face whose swept volume's vertices (x_f, the
face points, the apex, the traced points) have no vertex outside any donor
plane (signed distance <= 0) is clipped by the donor's PLIC plane alone,
with no ring search, no pre-test of the ring and no cell-plane clips. Bit
for bit by construction: a clip changes nothing without a vertex strictly
outside (clipWith), a triangle lying in a donor plane is kept uncut by the
contact rule (the donor's centre inside), a pre-image part lying in the
face's own plane is kept by the coplanar rule on the donor's side, and
every other ring cell is separated from the swept volume by a donor plane
or, the cell across a face parallel to the velocity, adds nothing. Measured
on the eight cases at step 20: 0 faces of 4582 (rotation N 32), 0 of 17528
(N 64), 0 of 4580 (the diagonal), 0 of 367 (the plane N 40, normal arm),
1 of 598 (N 20, oblique), 26 of 5111 (N 40, oblique). Why so few: under
the rotation every sweeping face's traced vertices leave the donor across a
lateral plane (A); under an axis-aligned translation the face points lie on
the donor's lateral planes, where the signed distance (p - x0) . n to the
face centre x0 of that lateral face is +-1e-17 in floating point, so one
vertex is "outside" by round-off on nearly every face. A tolerance would
take those faces but the clips would no longer be no-ops (a vertex at
s = +1e-17 is clipped, the crossing point replaces it), so the result would
be round-off, not identical; the rule "a tolerance that keeps the result
identical" leaves tolerance 0. The path stays in the code (the key), off.

B.2, the face update from the volume clip: the pre-image fan triangles are
every fourth triangle of the closed list, so their wetted parts after the
volume clip's cell-plane and PLIC clips (pieceResult::aPre, covPre) give
the wetted area vector of a sweeping face without the separate
sweptPieceArea clip (50 area clips per sweeping face in A). Not bit for
bit: the triangles enter the clip in the opposite vertex order (outward),
so the crossing points are computed from the other end, and the volume
clip's remnant filter (1e-14 L^2) applies. Measured after 20 steps: max
|dalpha| 1.0e-12 (R32), 2.6e-10 (R64), 1.6e-12 (D), 2.7e-12 (P), 1.5e-17
(P20n), 3.4e-14 (P20o), 2.9e-13 (P40n), 4.7e-13 (P40o); max |dalphaf|
1.2e-11, 1.3e-8 (R64), 3.4e-10, 1.7e-12, 1.5e-17, 1.1e-13, 1.2e-12,
2.9e-12: the round-off of the reordered clip, amplified on the rotation at
N 64 by the rank-fallback faces as the trace cache's first version was
(2026-11-04). The faces without a flux, the dump, the debug face and the
check keep the separate clip, which stays selectable (the key off). Opt-in
by the rules (it alters the numbers).

B.3, the pre-test's coplanar side (from A's class table): a ring cell
separated from the swept volume by one of its planes that holds a
pre-image triangle was kept for the face update's coplanar rule; with the
rule donorSide that part belongs to the cell on the donor's side of the
plane, so a cell across the plane ((n . nRef) <= 0) adds nothing (its
coplanar part is dropped by the rule, a contact is dropped by the volume
clip) and is skipped. Bit for bit on the eight cases. The faces without a
flux touch 4.09 cells instead of 8.17 and run 92 area clips instead of
165; their four timed parts 0.064 -> 0.050 s/step, the sweeping faces'
pre-test 0.042 -> 0.035 (job 55965852, with B.1 and B.2 on: the class
table below). Bit for bit, it can become the default; the key is left
false in this commit so that the reference setting 407f969 stays
reproducible by the builds between.

The class table with the three keys on (N 64, one turn, job 55965852;
0.4623 s/step with the timers, 0.4531 without):

| class | faces/step | touched cells/face | pieces/face | volume clips/face | area clips/face | t volume | t area | t trace | t pre-test | sum |
|---|---|---|---|---|---|---|---|---|---|---|
| no-sweep | 5856 | 4.09 (8.17) | 0 | 0 | 92.0 (165.5) | 0 | 0.0142 (0.0241) | 0.0165 | 0.0195 (0.0229) | 0.0502 (0.0642) |
| (ii) multi-cell, sweeping | 11718 | 2.04 | 2.04 | 140.3 | 0 (50.2) | 0.0681 | 0 (0.0137) | 0.0516 | 0.0345 (0.0418) | 0.1542 (0.1764) |
| the four timed parts | | | | | | | | | | 0.2044 (0.2407) |

### C MEASURED (2026-11-05, jobs 55965344-55; csv/fluxC.csv, csv/fluxC.md)

The plane-clipped loop (the preprint's: the davof reconstruction with the
prism flux and the plane face update) on the comparison's cases with the
exact normal and with the ring normal, three repetitions each; the swept
rows from the reference rerun (ref/cmp: the exact row one repetition, the
ring row three), plicRDF from the comparison. The main clone's build at
da04ece (the reference setting; the counters' key off). A run the abort
guard stopped (the interface cells above the factor of the initial count)
is "diverged" with its step; its norms are not scored.

| cell | OF Co | loop, normal | L1(x) | E_shape | L1(alphaf) | L1(n) | min / excess | mass % | wall median [min, max] s | step s |
|---|---|---|---|---|---|---|---|---|---|---|
| rotation N 32, 1 ranks | 0.200 | swept, exact | 1.9624e-03 | 5.2607e-04 | 2.9807e-02 | 3.1766e-03 | -2.1e-13 / 2.1e-12 | 9.3e-13 | 168.3 [168.3, 168.3] | 0.100 |
| rotation N 32, 1 ranks | 0.200 | swept, ring | 3.9776e-03 | 1.0587e-03 | 5.7035e-02 | 5.7649e-02 | -6.9e-13 / 1.9e-11 | 1.1e-12 | 164.4 [164.2, 164.4] | 0.090 |
| rotation N 32, 1 ranks | 0.200 | plane-clipped, exact | 2.0397e-03 | 5.4593e-04 | 3.0518e-02 | 3.8068e-03 | -9.1e-06 / 2.4e-15 | 1.0e-12 | 40.3 [40.3, 40.7] | 0.020 |
| rotation N 32, 1 ranks | 0.200 | plane-clipped, ring | 4.1106e-03 | 1.0907e-03 | 5.7716e-02 | 5.9394e-02 | -6.8e-13 / 5.1e-15 | 1.5e-12 | 37.5 [37.5, 37.6] | 0.020 |
| rotation N 32, 1 ranks | 0.200 | plicRDF + isoAdvector | 4.3887e-03 | 1.1552e-03 | - | 5.7522e-02 | -5.1e-19 / 2.2e-16 | 1.2e-12 | 45.2 [45.0, 45.4] | 0.020 |
|  | | the flux's share of the exact-normal advantage: L1(x) 0.032, E_shape 0.032; plane-clipped / swept (exact): L1(x) 1.039, E_shape 1.038, wall 0.24 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 1.033, E_shape 1.030, L1(alphaf) 1.012, wall 0.23; plane-clipped-ring / plicRDF wall 0.83 | | | | | | | | |
| rotation N 32, 1 ranks | 1.000 | swept, exact | 1.4072e-03 | 3.7576e-04 | 2.1590e-02 | 3.5835e-03 | -7.1e-13 / 6.1e-13 | 1.3e-12 | 35.6 [35.6, 35.6] | 0.090 |
| rotation N 32, 1 ranks | 1.000 | swept, ring | 2.6449e-03 | 7.0707e-04 | 3.8232e-02 | 4.8100e-02 | -9.2e-13 / 6.0e-13 | 8.9e-13 | 34.8 [34.8, 35.1] | 0.090 |
| rotation N 32, 1 ranks | 1.000 | plane-clipped, exact | - | - | - | - | -1.4e-01 / 9.3e-15 (diverged, the abort guard at step 208 of 208) | 3.8e-13 | 7.3 [7.3, 7.4] | 0.020 |
| rotation N 32, 1 ranks | 1.000 | plane-clipped, ring | - | - | - | - | -1.7e-01 / 7.3e-15 (diverged, the abort guard at step 109 of 109) | 9.2e-01 | 5.0 [5.0, 5.0] | 0.020 |
| rotation N 32, 1 ranks | 1.000 | plicRDF + isoAdvector | 3.3276e-03 | 8.4772e-04 | - | 5.4431e-02 | -5.1e-18 / 2.2e-16 | 1.2e-12 | 11.3 [11.3, 11.3] | 0.030 |
|  | | the flux's share of the exact-normal advantage: L1(x) nan, E_shape nan; plane-clipped / swept (exact): L1(x) nan, E_shape nan, wall 0.20 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.14; plane-clipped-ring / plicRDF wall 0.44 | | | | | | | | |
| rotation N 64, 8 ranks | 0.200 | swept, exact | 6.2757e-04 | 1.7175e-04 | 1.9290e-02 | 8.4615e-04 | -7.7e-13 / 1.0e-10 | 4.5e-12 | 720.4 [720.4, 720.4] | 0.200 |
| rotation N 64, 8 ranks | 0.200 | swept, ring | 1.0577e-03 | 2.9561e-04 | 3.2196e-02 | 2.5174e-02 | -6.2e-14 / 6.2e-10 | 4.2e-12 | 705.1 [704.2, 706.1] | 0.200 |
| rotation N 64, 8 ranks | 0.200 | plane-clipped, exact | 6.4763e-04 | 1.7697e-04 | 1.9783e-02 | 9.7257e-04 | -1.4e-05 / 8.7e-15 | 4.8e-13 | 131.0 [130.9, 131.3] | 0.040 |
| rotation N 64, 8 ranks | 0.200 | plane-clipped, ring | - | - | - | - | -5.1e-02 / 1.1e-14 (diverged, the abort guard at step 3420 of 3420) | 5.5e-03 | 126.2 [126.0, 126.4] | 0.040 |
| rotation N 64, 8 ranks | 0.200 | plicRDF + isoAdvector | 1.1393e-03 | 3.1273e-04 | - | 2.4068e-02 | -2.6e-19 / 2.2e-16 | 1.1e-12 | 192.2 [189.9, 192.3] | 0.050 |
|  | | the flux's share of the exact-normal advantage: L1(x) 0.039, E_shape 0.037; plane-clipped / swept (exact): L1(x) 1.032, E_shape 1.030, wall 0.18 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.18; plane-clipped-ring / plicRDF wall 0.66 | | | | | | | | |
| rotation N 64, 8 ranks | 1.000 | swept, exact | 4.4486e-04 | 1.2085e-04 | 1.3792e-02 | 7.5390e-04 | -5.6e-12 / 5.3e-12 | 4.2e-12 | 147.9 [147.9, 147.9] | 0.200 |
| rotation N 64, 8 ranks | 1.000 | swept, ring | 7.1023e-04 | 1.9738e-04 | 2.1912e-02 | 2.2567e-02 | -5.3e-13 / 6.1e-12 | 4.1e-12 | 145.3 [144.7, 145.5] | 0.200 |
| rotation N 64, 8 ranks | 1.000 | plane-clipped, exact | - | - | - | - | -1.7e-01 / 3.0e-03 (diverged, the abort guard at step 357 of 357) | 7.0e-13 | 18.2 [18.2, 18.3] | 0.040 |
| rotation N 64, 8 ranks | 1.000 | plane-clipped, ring | - | - | - | - | -3.1e-01 / 3.8e-03 (diverged, the abort guard at step 189 of 189) | 3.0e-01 | 12.4 [12.4, 12.5] | 0.040 |
| rotation N 64, 8 ranks | 1.000 | plicRDF + isoAdvector | 8.4673e-04 | 2.2972e-04 | - | 2.2329e-02 | -2.0e-17 / 4.4e-16 | 9.4e-13 | 42.2 [41.3, 44.0] | 0.050 |
|  | | the flux's share of the exact-normal advantage: L1(x) nan, E_shape nan; plane-clipped / swept (exact): L1(x) nan, E_shape nan, wall 0.12 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.09; plane-clipped-ring / plicRDF wall 0.29 | | | | | | | | |
| rotation N 128, 64 ranks | 0.200 | swept, exact | 1.8746e-04 | 5.2188e-05 | 1.1923e-02 | 2.6495e-04 | -7.9e-14 / 5.1e-10 | 8.8e-12 | 3622.8 [3622.8, 3622.8] | 0.470 |
| rotation N 128, 64 ranks | 0.200 | swept, ring | 2.8343e-04 | 8.0418e-05 | 1.8297e-02 | 1.1299e-02 | -1.8e-13 / 7.4e-10 | 7.6e-12 | 3541.6 [3538.6, 3551.4] | 0.470 |
| rotation N 128, 64 ranks | 0.200 | plane-clipped, exact | 1.9343e-04 | 5.3950e-05 | 1.2288e-02 | 2.9719e-04 | -3.1e-03 / 1.2e-14 | 1.4e-12 | 680.2 [678.4, 680.5] | 0.090 |
| rotation N 128, 64 ranks | 0.200 | plane-clipped, ring | - | - | - | - | -3.0e-02 / 1.0e-14 (diverged, the abort guard at step 1205 of 1205) | -2.0e-13 | 117.7 [117.5, 117.8] | 0.080 |
| rotation N 128, 64 ranks | 0.200 | plicRDF + isoAdvector | 2.9492e-04 | 8.3739e-05 | - | 9.5597e-03 | -5.0e-13 / 4.4e-16 | 1.8e-12 | 959.4 [958.6, 970.1] | 0.130 |
|  | | the flux's share of the exact-normal advantage: L1(x) 0.056, E_shape 0.056; plane-clipped / swept (exact): L1(x) 1.032, E_shape 1.034, wall 0.19 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.03; plane-clipped-ring / plicRDF wall 0.12 | | | | | | | | |
| rotation N 128, 64 ranks | 1.000 | swept, exact | 1.3028e-04 | 3.6197e-05 | 8.3619e-03 | 2.3922e-04 | -6.6e-11 / 1.4e-10 | 9.6e-12 | 739.1 [739.1, 739.1] | 0.480 |
| rotation N 128, 64 ranks | 1.000 | swept, ring | 1.8942e-04 | 5.3417e-05 | 1.2274e-02 | 9.9009e-03 | -5.4e-11 / 7.2e-11 | 8.1e-12 | 727.8 [717.6, 851.4] | 0.470 |
| rotation N 128, 64 ranks | 1.000 | plane-clipped, exact | - | - | - | - | -6.1e-01 / 8.0e-02 (diverged, the abort guard at step 803 of 803) | 1.0e-13 | 102.9 [102.9, 103.3] | 0.120 |
| rotation N 128, 64 ranks | 1.000 | plane-clipped, ring | - | - | - | - | -3.1e-01 / 9.3e-02 (diverged, the abort guard at step 331 of 331) | 8.1e-02 | 44.4 [44.4, 44.5] | 0.080 |
| rotation N 128, 64 ranks | 1.000 | plicRDF + isoAdvector | 2.2338e-04 | 6.0959e-05 | - | 9.4322e-03 | -9.8e-13 / 2.9e-13 | 1.9e-12 | 204.3 [201.6, 207.6] | 0.130 |
|  | | the flux's share of the exact-normal advantage: L1(x) nan, E_shape nan; plane-clipped / swept (exact): L1(x) nan, E_shape nan, wall 0.14 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.06; plane-clipped-ring / plicRDF wall 0.22 | | | | | | | | |
| diagonal N 32, 1 ranks | 0.200 | swept, exact | 1.1550e-03 | 2.8940e-04 | 1.6873e-02 | 2.3675e-03 | -3.2e-13 / 2.0e-13 | 7.8e-14 | 18.1 [18.1, 18.1] | 0.090 |
| diagonal N 32, 1 ranks | 0.200 | swept, ring | 2.0126e-03 | 5.7286e-04 | 3.1749e-02 | 5.3562e-02 | -2.0e-14 / 2.8e-13 | 2.0e-13 | 16.1 [16.1, 16.6] | 0.080 |
| diagonal N 32, 1 ranks | 0.200 | plane-clipped, exact | 1.4169e-03 | 3.4274e-04 | 1.9826e-02 | 2.9428e-03 | -1.3e-04 / 1.0e-14 | 3.9e-13 | 5.5 [5.5, 5.6] | 0.020 |
| diagonal N 32, 1 ranks | 0.200 | plane-clipped, ring | 2.1058e-03 | 5.9854e-04 | 3.2753e-02 | 5.5036e-02 | -1.0e-14 / 1.0e-14 | 3.3e-13 | 5.3 [5.3, 5.3] | 0.020 |
| diagonal N 32, 1 ranks | 0.200 | plicRDF + isoAdvector | 1.9296e-03 | 5.1187e-04 | - | 4.2305e-02 | -2.5e-13 / 2.2e-16 | 8.5e-13 | 6.1 [6.0, 6.1] | 0.020 |
|  | | the flux's share of the exact-normal advantage: L1(x) 0.338, E_shape 0.240; plane-clipped / swept (exact): L1(x) 1.227, E_shape 1.184, wall 0.31 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 1.046, E_shape 1.045, L1(alphaf) 1.032, wall 0.33; plane-clipped-ring / plicRDF wall 0.87 | | | | | | | | |
| diagonal N 32, 1 ranks | 1.000 | swept, exact | 4.1866e-04 | 1.1818e-05 | 4.8457e-03 | 2.3924e-03 | -2.3e-14 / 1.4e-13 | 1.6e-13 | 5.4 [5.4, 5.4] | 0.080 |
| diagonal N 32, 1 ranks | 1.000 | swept, ring | 7.8639e-04 | 1.8090e-04 | 1.1178e-02 | 3.8703e-02 | -6.0e-13 / 2.6e-13 | 1.7e-13 | 5.7 [5.3, 6.3] | 0.080 |
| diagonal N 32, 1 ranks | 1.000 | plane-clipped, exact | 2.3603e-02 | 3.9105e-03 | 1.7036e-01 | 3.5782e-02 | -9.3e-02 / 1.6e-02 | 1.4e-13 | 3.2 [3.2, 3.2] | 0.020 |
| diagonal N 32, 1 ranks | 1.000 | plane-clipped, ring | 9.8349e-02 | 3.3791e-03 | 1.8808e-01 | 7.6533e-01 | -1.2e-01 / 1.6e-03 | 8.8e-01 | 3.3 [3.3, 3.3] | 0.020 |
| diagonal N 32, 1 ranks | 1.000 | plicRDF + isoAdvector | 6.4695e-03 | 1.2946e-03 | - | 1.2457e-01 | -1.0e-14 / 9.3e-15 | 1.3e-13 | 3.1 [3.1, 3.2] | 0.030 |
|  | | the flux's share of the exact-normal advantage: L1(x) 3.832, E_shape 3.039; plane-clipped / swept (exact): L1(x) 56.378, E_shape 330.897, wall 0.59 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 125.063, E_shape 18.679, L1(alphaf) 16.826, wall 0.58; plane-clipped-ring / plicRDF wall 1.05 | | | | | | | | |
| diagonal N 64, 8 ranks | 0.200 | swept, exact | 3.2914e-04 | 8.5795e-05 | 9.8911e-03 | 8.4231e-04 | -1.2e-11 / 5.3e-12 | 3.7e-13 | 63.4 [63.4, 63.4] | 0.190 |
| diagonal N 64, 8 ranks | 0.200 | swept, ring | 5.5372e-04 | 1.5684e-04 | 1.7696e-02 | 2.4281e-02 | -2.0e-11 / 5.3e-12 | 3.9e-13 | 60.4 [60.4, 60.5] | 0.180 |
| diagonal N 64, 8 ranks | 0.200 | plane-clipped, exact | 2.6684e-03 | 1.7558e-04 | 3.8586e-02 | 5.5897e-03 | -2.5e-02 / 1.0e-14 | 3.1e-13 | 14.0 [13.9, 14.3] | 0.030 |
| diagonal N 64, 8 ranks | 0.200 | plane-clipped, ring | 6.9686e-02 | 1.9137e-04 | 1.0868e-01 | 5.4833e-01 | -4.5e-02 / 8.4e-03 | 3.9e-02 | 19.9 [19.9, 19.9] | 0.060 |
| diagonal N 64, 8 ranks | 0.200 | plicRDF + isoAdvector | 5.1269e-04 | 1.4356e-04 | - | 1.9224e-02 | -4.5e-10 / 4.7e-13 | 5.3e-13 | 19.9 [19.8, 20.2] | 0.050 |
|  | | the flux's share of the exact-normal advantage: L1(x) 12.744, E_shape 1.554; plane-clipped / swept (exact): L1(x) 8.107, E_shape 2.046, wall 0.22 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 125.852, E_shape 1.220, L1(alphaf) 6.142, wall 0.33; plane-clipped-ring / plicRDF wall 1.00 | | | | | | | | |
| diagonal N 64, 8 ranks | 1.000 | swept, exact | 1.0070e-04 | 3.7843e-06 | 2.3569e-03 | 7.2236e-04 | -3.3e-11 / 1.5e-12 | 2.4e-13 | 16.1 [16.1, 16.1] | 0.190 |
| diagonal N 64, 8 ranks | 1.000 | swept, ring | 1.7639e-04 | 3.8868e-05 | 5.0451e-03 | 1.9075e-02 | -1.2e-13 / 2.0e-12 | 1.9e-13 | 15.9 [15.4, 17.1] | 0.180 |
| diagonal N 64, 8 ranks | 1.000 | plane-clipped, exact | 1.6971e-02 | 2.7172e-03 | 1.9562e-01 | 2.3644e-02 | -1.1e-01 / 1.9e-01 | 1.1e-13 | 6.3 [6.3, 6.5] | 0.040 |
| diagonal N 64, 8 ranks | 1.000 | plane-clipped, ring | 6.9098e-02 | 2.9658e-03 | 1.8403e-01 | 9.1427e-01 | -1.3e-01 / 1.2e-01 | 2.8e-01 | 7.2 [7.2, 7.2] | 0.040 |
| diagonal N 64, 8 ranks | 1.000 | plicRDF + isoAdvector | 5.0766e-03 | 6.3059e-04 | - | 1.3923e-01 | -1.5e-03 / 1.2e-03 | -4.4e-14 | 9.1 [9.0, 9.3] | 0.080 |
|  | | the flux's share of the exact-normal advantage: L1(x) 3.390, E_shape 4.329; plane-clipped / swept (exact): L1(x) 168.522, E_shape 718.017, wall 0.39 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 391.728, E_shape 76.305, L1(alphaf) 36.477, wall 0.45; plane-clipped-ring / plicRDF wall 0.79 | | | | | | | | |
| diagonal N 128, 64 ranks | 0.200 | swept, exact | 9.0167e-05 | 2.4429e-05 | 5.7134e-03 | 2.4190e-04 | -2.4e-11 / 1.2e-11 | 1.5e-12 | 369.6 [369.6, 369.6] | 0.560 |
| diagonal N 128, 64 ranks | 0.200 | swept, ring | 1.2827e-04 | 3.6888e-05 | 8.8200e-03 | 9.7307e-03 | -1.8e-11 / 1.2e-11 | 1.3e-12 | 348.9 [347.6, 350.8] | 0.530 |
| diagonal N 128, 64 ranks | 0.200 | plane-clipped, exact | 1.2379e-04 | 3.3292e-05 | 7.5656e-03 | 3.9157e-04 | -3.9e-02 / 1.5e-03 | 7.2e-13 | 66.5 [66.3, 66.5] | 0.090 |
| diagonal N 128, 64 ranks | 0.200 | plane-clipped, ring | 5.4488e-02 | 7.0303e-05 | 1.1177e-01 | 5.7199e-01 | -4.4e-02 / 6.4e-02 | 5.2e-03 | 90.4 [90.2, 90.5] | 0.130 |
| diagonal N 128, 64 ranks | 0.200 | plicRDF + isoAdvector | 1.2306e-04 | 3.4374e-05 | - | 7.7935e-03 | -1.2e-07 / 1.7e-13 | 8.1e-12 | 99.5 [98.8, 100.0] | 0.140 |
|  | | the flux's share of the exact-normal advantage: L1(x) 1.022, E_shape 0.891; plane-clipped / swept (exact): L1(x) 1.373, E_shape 1.363, wall 0.18 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) 424.775, E_shape 1.906, L1(alphaf) 12.672, wall 0.26; plane-clipped-ring / plicRDF wall 0.91 | | | | | | | | |
| diagonal N 128, 64 ranks | 1.000 | swept, exact | 2.6397e-05 | 8.7709e-07 | 1.2420e-03 | 1.8660e-04 | -1.4e-13 / 1.3e-11 | 1.7e-12 | 82.2 [82.2, 82.2] | 0.580 |
| diagonal N 128, 64 ranks | 1.000 | swept, ring | 3.5490e-05 | 7.1839e-06 | 2.1775e-03 | 8.2923e-03 | -3.6e-11 / 1.3e-11 | 1.6e-12 | 77.5 [77.3, 77.7] | 0.540 |
| diagonal N 128, 64 ranks | 1.000 | plane-clipped, exact | 1.8244e-02 | 2.3967e-03 | 2.3181e-01 | 1.3253e-02 | -2.0e-01 / 2.4e-01 | 2.0e-13 | 27.4 [27.4, 27.5] | 0.130 |
| diagonal N 128, 64 ranks | 1.000 | plane-clipped, ring | - | - | - | - | -2.2e-01 / 1.6e-01 (diverged, the abort guard at step 73 of 73) | 3.0e-13 | 18.0 [17.9, 18.0] | 0.160 |
| diagonal N 128, 64 ranks | 1.000 | plicRDF + isoAdvector | 2.9987e-03 | 3.2206e-04 | - | 1.4836e-01 | -6.6e-03 / 4.1e-03 | 4.3e-13 | 40.5 [40.5, 40.6] | 0.250 |
|  | | the flux's share of the exact-normal advantage: L1(x) 6.129, E_shape 7.459; plane-clipped / swept (exact): L1(x) 691.151, E_shape 2732.546, wall 0.33 | | | | | | | | |
|  | | plane-clipped / swept (ring): L1(x) nan, E_shape nan, L1(alphaf) nan, wall 0.23; plane-clipped-ring / plicRDF wall 0.44 | | | | | | | | |

Verdicts:
- F-1 (the flux's share of the exact-normal advantage at OpenFOAM Co 0.2):
  HELD on the rotation: the plane-clipped loop with the exact normal has
  L1(x) and E_shape 1.03 to 1.04 times the swept loop's at N 32 / 64 / 128,
  the flux's share of the advantage over plicRDF 0.03 to 0.06; the normal
  carries it. FALSIFIED on the diagonal: 1.23 x at N 32 (share 0.34), 8.1 x
  at N 64 (the plane-clipped loop loses its boundedness, min alpha -2.5e-2;
  share 12.7), 1.37 x at N 128 (min alpha -3.9e-2; share 1.02). And the
  plane-clipped loop is not bounded even where its norms are close: min
  alpha -9.1e-6 (rotation N 32), -1.4e-5 (N 64), -3.1e-3 (N 128), -1.3e-4
  (diagonal N 32); the swept loop >= -7.9e-14 there.
- F-2 (OpenFOAM Co 1): FALSIFIED in its rotation half. The plane-clipped
  loop diverges on the rotation too (the abort guard at step 208 of 208 at
  N 32, 357 of 357 at N 64, 803 of 803 at N 128, min alpha -0.14 to -0.61),
  not only on the diagonal (min alpha -0.093 / -0.11 / -0.20 at N 32 / 64 /
  128, L1(x) 56 / 169 / 691 times the swept loop's); the swept loop is
  bounded to 6.6e-11 in every cell. The "prism aborts at N 64 / 128" of
  2026-10-29 (OpenFOAM Co 0.86) returns at Co 1 at every N.
- F-3 (the ring-normal rows): the two loops with the ring normal are within
  3.3 % (L1(x)) and 3.0 % (E_shape) on the rotation at N 32, Co 0.2 (HELD
  there, within 4.6 % on the diagonal N 32); at N 64 / 128 the plane-clipped
  ring loop diverges on the rotation at Co 0.2 (step 3420 of 3420, 1205 of
  1205) and loses its boundedness on the diagonal (L1(x) 126 to 425 times
  the swept loop's, the mass off by 0.5 to 4 %); at Co 1 it diverges in
  every cell. FALSIFIED beyond N 32.
- F-4 (the cost): the plane-clipped loop's wall clock is 0.12 to 1.05 times
  plicRDF's and 0.03 to 0.59 times the swept loop's (the diverged runs
  shorter still): below the band, FALSIFIED; per step 0.02 to 0.16 s against
  the swept loop's 0.08 to 0.58.

What this shows: with the same exact normal, the swept-face advection and
the plane-clipped loop transport alike on the rotation at OpenFOAM Co 0.2
(within 4 %), so the exact-normal rows' advantage over plicRDF there is the
normal's, not the flux's; everywhere else the flux decides: the
plane-clipped loop loses its boundedness on the diagonal at Co 0.2 from
N 64 and diverges at Co 1 on both flows and every N, where the swept loop
stays bounded to 1e-10. The plane-clipped loop costs 2 to 8 times less per
step than the swept loop (its flux is the prism of isoAdvector's kind).

### Round 3b: the re-profile with the three keys (2026-11-05)

The full-turn profiles of the first and the second round repeated for the
fork a11ebc1 with the three keys of B on (sweptDonorFastPath,
sweptPreTestCoplanarSide, sweptFaceUpdateFromVolume; the numbers of B.2's
round-off), the same cases, flags, events and readers (prof/runs3/{serial,
par8, par64}/<run>/hpctoolkit-database, the analyses prof/an4/ with
rules_swept4.py, the line ranges of a11ebc1; csv/hp3_categories.csv,
hp3_ranks.csv, hp3_hotspots.csv, hp3_timing.csv; jobs 55966356-59). The
before column is the reference setting's profile (prof/an2, 407f969).

(a) The serial breakdown per step, N = 64, one turn (697 steps):

| category | 407f969 s/step | % | a11ebc1 + B s/step | % | ratio |
|---|---|---|---|---|---|
| pieces: the clip by the cell planes and the PLIC plane | 0.0660 | 13.4 | 0.0662 | 14.4 | 1.00 |
| pieces: the cap closure and the volume | 0.0012 | 0.2 | 0.0010 | 0.2 | 0.83 |
| face update: the pre-image fan clipped (sweptPieceArea; after: the faces without a flux only) | 0.0353 | 7.2 | 0.0137 | 3.0 | 0.39 |
| swept cells: the ring search | 0.0127 | 2.6 | 0.0089 | 1.9 | 0.70 |
| swept cells: the cell planes, state, wisp rule and the pre-test (one range after) | 0.0291 + 0.0372 | 13.5 | 0.0543 | 11.8 | 0.82 |
| swept cells: the shortcut and the check (off), the fast-path test | 0.0013 | 0.3 | 0.0002 | 0.1 | - |
| swept cells: the pre-image fan and the rest | 0.0122 | 2.5 | 0.0122 | 2.6 | 1.00 |
| trace: the velocity at the traced points (findTetrahedron) | 0.0302 | 6.2 | 0.0293 | 6.4 | 0.97 |
| trace: the cell search (ring, octree) | 0.0153 | 3.1 | 0.0146 | 3.2 | 0.95 |
| trace: the midpoint trace, the triangle list | 0.0107 | 2.2 | 0.0111 | 2.4 | 1.04 |
| the point velocity (volPointInterpolation of U) | 0.0223 | 4.5 | 0.0223 | 4.9 | 1.00 |
| curvedFlux: the band, the wisp surfaces, the exchanges, the face loop | 0.0226 | 4.6 | 0.0234 | 5.1 | 1.04 |
| the remote exchange (serial: the stencil setup) | 0.0031 | 0.6 | 0.0030 | 0.7 | 0.97 |
| the surface build (fit, anchor, volume rule) | 0.0155 | 3.2 | 0.0158 | 3.4 | 1.02 |
| the reconstruction (davof, ring normal, positioning) | 0.0299 | 6.1 | 0.0301 | 6.6 | 1.01 |
| the bounding | 0.0113 | 2.3 | 0.0112 | 2.4 | 0.99 |
| isoAdvector: timeIntegratedFlux | 0.0150 | 3.0 | 0.0149 | 3.2 | 0.99 |
| isoAdvector: advect, the rest (alpha update, counters) | 0.0291 | 5.9 | 0.0290 | 6.3 | 1.00 |
| setFlow | 0.0309 | 6.3 | 0.0305 | 6.6 | 0.99 |
| the Courant number, other function objects | 0.0016 | 0.3 | 0.0016 | 0.3 | 1.00 |
| I/O and start-up | 0.0218 | 4.4 | 0.0218 | 4.8 | 1.00 |
| the solver loop (main) | 0.0355 | 7.2 | 0.0351 | 7.7 | 0.99 |
| total | 0.4907 | 100 | 0.4586 | 100 | 0.93 |

- The step 0.4907 -> 0.4586 s (0.93; the 200-step timing 0.4847 -> 0.4531,
  0.935). The face update's clip fell to the faces without a flux
  (0.0353 -> 0.0137), the pre-test and the ring search by the cells across
  the plane (0.0663 -> 0.0543, 0.0127 -> 0.0089). N 32: 0.0989 -> 0.0921.
- The ratio to plicRDF's serial step (0.2090): 2.19 (2.35 before, 5.3 in
  the first profile). The swept parts are 0.283 s/step (62 % of the step),
  the parts shared with plicRDF 0.175 (38 %): the shared parts alone are
  84 % of plicRDF's whole step.

(c) Parallel, the comparison's scotch decompositions (no balancing):

| | N 64, 8 ranks, 407f969 | after | N 128, 64 ranks, 407f969 | after |
|---|---|---|---|---|
| step (CPUTIME per rank / steps; the wall clock) | 0.207 s (0.212) | 0.190 s (0.194) | 0.508 s (0.519) | 0.471 s (0.481) |
| MPI library, share of all CPUTIME | 64.5 % | 63.2 % | 84.6 % | 83.5 % |
| MPI per rank: mean / max / min | 93 / 123 / 27 s | 84 / 111 / 24 s | 603 / 674 / 332 s | 551 / 621 / 304 s |
| compute per rank: max / mean (max/mean) | 113 / 49 (2.31) | 104 / 46 (2.25) | 361 / 101 (3.59) | 337 / 100 (3.38) |
| the balanced step | 0.070 s | 0.066 s | 0.072 s | 0.071 s |
| the ratio to plicRDF's step at the same ranks (0.0556 / 0.1395 s) | 3.8 | 3.5 | 3.7 | 3.4 |

The busiest rank's compute fell by the serial factor (113 -> 104, 361 ->
337 s), the idle ranks' not, so the step fell by 8 % at 8 ranks and 7 % at
64; the imbalance (2.25 and 3.38) sets the step as before (the
decomposition is deferred).

(d) The top exclusive hotspots after (serial, N = 64, one turn):

| # | % | s | function (file:line; fork a11ebc1) | kind |
|---|---|---|---|---|
| 1 | 5.1 | 16.4 | Foam::minusOp (ops.H:209): the subtractions of planeDistances and of the pre-test | algorithmic |
| 2 | 5.0 | 15.9 | sweptFaceCellsFlux's own lines (the pre-test's loops over the vertices, the cell planes) | algorithmic |
| 3 | 4.5 | 14.2 | sweptPieceVolume's own lines (the per-triangle loop, the cap bookkeeping) | algorithmic |
| 4 | 4.3 | 13.6 | volPointInterpolation::interpolateInternalField: the point velocity of the whole mesh | mesh-wide, redundant (the band's points would do) |
| 5 | 4.1 | 13.0 | Foam::Vector::inner (VectorI.H:157): the signed distances | algorithmic |
| 6 | 3.9 | 12.3 | __mempcpy (libc): the field assignments of the solver loop and setFlow | redundant copies, mesh-wide |
| 7 | 2.7 (6.6 with its callees) | 8.6 (20.9) | cellPointWeight::findTetrahedron: the velocity at a traced point | redundant + allocation (the tetrahedra rebuilt per lookup) |
| 8 | 2.0 | 6.3 | fvc::surfaceIntegrate (the alpha update) | mesh-wide |
| 9 | 1.9 | 5.9 | davof::clipWith (davofSweptPieces.H:216): the clip proper | algorithmic |
| 10 | 1.8 | 5.7 | Foam::Vector::cross: areaVector, the triangle normals | algorithmic |
| 11 | 1.7 | 5.3 | surfaceInterpolationScheme::dotInterpolate (phi from U) | mesh-wide |
| 12 | 1.6 (3.2 with magSqr) | 5.1 (10.2) | plusEqMagSqrOp2, magSqr: the edge lengths | algorithmic, in part redundant |

Total exclusive 319.7 s; the allocation 13.0 s (4.1 %). The clip's
vector subtraction and dot product fell from 23.8 + 16.0 s to 16.4 + 13.0
s with the area clips gone from the sweeping faces and the cells across
the plane gone from the faces without a flux.

The next gains (not implemented; the expected figures from this profile):
1. The pre-test and the cell states of the faces without a flux (one
   third of the faces, 0.050 s/step for a face fraction without a flux):
   a face with phi_f = 0 needs the wetted area of its in-plane shifted
   pre-image in the donor's layer only; its ring search, planes and states
   could be built for that layer alone (the four cells of the donor's
   side), and the trace of its points is shared with the sweeping faces
   through the cache already. Expected ~ 0.02 s/step (4 %).
2. The pre-test on the hex mesh as a box test (the cell planes are
   axis-aligned: six comparisons instead of 60 subtractions and dot
   products per ring cell), ~ 0.02 s/step (4 %); together with 1 the first
   pass (0.054) would halve.
3. The point velocity of the band's points only (0.022 s/step, 4.9 %,
   mesh-wide and growing x 8.6 with N), and the cell's tetrahedra cached
   per cell for the trace (0.029 + 0.015, 10 %): the two largest remaining
   non-clip items, ~ 0.03 s/step together.
4. The volume clip itself (0.066 s/step, 14 %, 140 half-space clips per
   sweeping face over 33 triangles in 2.04 cells, every touched cell with a
   piece): no cell without a piece is clipped any more, so the gain left is
   inside the clip: the sliver neighbour's piece (the second touched cell,
   a thin wedge) costs as much as the donor's; a closed form for a wedge
   between two planes (the lateral plane and the PLIC plane) would halve
   it (~ 7 %), round-off against the clip.
With 1 to 4 the serial N 64 step would reach ~ 0.38 s (1.8 times
plicRDF's); the parts shared with plicRDF (0.175 s) then dominate.

## Round 4: the five remaining gains (2026-11-05)

Tomislav's instruction (via the coordinating session): the five gains named
after round 3b, each verified on the eight cases (the rotation N 32 / 64,
the plane N 40 exact normal, the diagonal N 32, the plane gate N 20 / 40
both arms; 20 steps at writePrecision 17), bit for bit where possible,
round-off stated per case otherwise, opt-in keys for anything not bit for
bit, three timed repetitions. Baseline: a11ebc1 with B.2 and B.3 on
(refB23, job 55972899 for the identity; the 200-step timing of that
configuration 0.4531 s/step, the full-turn profile 0.4586). Fork commits
(tmaric, nothing merged): c884540 (the keys sweptNoFluxLayer,
sweptAxisPlanes, sweptPointVelocityLazy, sweptTetCache; davofVelocityAt.H),
81dcc8a (sweptSliverClosedForm), fbbfaf9 (gain 1 made exact). Jobs:
55973440 (the first form of gain 1), 55973713-18 (the cumulative cycles
g1b, g12, g123, g1234, g12345, g12345stats), 55973719-22 (the round-4
profiles).

### The gains, each with its identity and its step (rotation N 64 serial, OpenFOAM Co 1, 200 steps x 3, the median of the per-step mean)

| gain | key | what | identity on the eight cases | s/step (min .. max) | against the baseline 0.4531 |
|---|---|---|---|---|---|
| 1 | sweptNoFluxLayer | a face without a flux: a ring cell across the face plane skipped before its planes and state when its own face in that plane separates it (the pre-test's test on that one plane, B.3's rule) | bit for bit (fbbfaf9, job 55973713); the first form (c884540, decided from the plane through x_f) moved the plane gate's oblique arm by 3.4e-14 / 6.2e-13 in alpha, 1.3e-13 / 2.1e-12 in alphaf (N 20 / 40) and gave 0.4398 | 0.4534 (0.4518 .. 0.4535) | 1.00: no gain |
| 2 | sweptAxisPlanes | the signed distance to an axis-aligned plane as the one coordinate difference (the clips, the triangle distances, the pre-test, the coplanar exception, the covered-vertex test, the fast path) | bit for bit (55973714): the two zero products and their sums vanish exactly, the product by +-1 is exact | 0.4408 (0.4399 .. 0.4433) | 0.973 |
| 3 | sweptPointVelocityLazy | the point velocity per point on first use: an interior point by the library's inverse-distance weights and sums repeated in order, a point of any patch by the whole-mesh field once; 7063 of 545025 points evaluated per step at N 64 (8182 with gain 4) | bit for bit (55973715), the plane gate's patch points included | 0.4406 (0.4403 .. 0.4421) | 0.972 (alone no gain: interpolationCellPoint still interpolates the whole mesh for its own lookup) |
| 4 | sweptTetCache | the trace's lookup on the cell's tetrahedra cached per cell (davofVelocityAt.H: cellTetIndices once per cell as point triples, cellPointWeight::findTetrahedron's search and interpolationCellPoint::interpolate's weights repeated operation by operation); no interpolationCellPoint on the swept path, so with gain 3 no whole-mesh point interpolation at all | bit for bit (55973716) | 0.3908 (0.3907 .. 0.3923) | 0.862 |
| 5 | sweptSliverClosedForm | the piece of a touched ring cell without a PLIC plane that one plane alone cuts, in closed form (sweptPieceOnePlane: triNegativePolygon per triangle, the volume (1/3) sum (x_T - x0) . a_T about the plane's point x0 where the cap's moment vanishes; the contact rule and the remnant filter kept); 6114 of 23982 pieces per step take it | round-off (55973717): max dalpha 7.4e-13 (R32), 1.1e-9 (R64), 3.3e-12 (D), 4.8e-14 (P20o), 5.0e-13 (P40o), 0 on the normal arms (P, P20n, P40n: no such piece); max dalphaf 7.8e-12, 4.3e-9, 2.7e-10, 1.5e-13, 1.5e-12 | 0.3823 (0.3814 .. 0.3827) | 0.844 |

Gain 1: the cells it spares are exactly those B.3 already skips at their
first plane; the 2.9 % of the first form came from the cells the pre-test
clips for a round-off sliver (their own face centre's distances differ from
the plane through x_f by 1e-17), which no exact rule can skip.
Gain 5's construction: the piece inside the half-space H of the plane is
Omega_f cut by H when every other plane of the cell has no vertex outside
(the pre-test now counts, per touched cell, the planes with a vertex
outside); its boundary is the triangles' parts inside H plus the cap in the
plane, so V = (1/3) sum over the parts of (x_T - x0) . a_T with x0 on the
plane (the divergence theorem on the piece's closed boundary; the cap's
term (x_cap - x0) . a_cap vanishes), and the parts come from the one-plane
closed form of the one-donor clip (triNegativePolygon). No citation beyond
that identity. Its error against the clip: the fields above (the clip
computes the polygon by Sutherland-Hodgman and its area vector by Newell,
the moments about the face centre with the cap solved); on the rotation
N 64 the difference reaches 1.1e-9 in alpha after 20 steps, the
rank-fallback amplification seen for B.2 and the trace cache's first form.

Together (gains 2 to 5 with B.2 and B.3): 0.4531 -> 0.3823 s/step (0.844
of the baseline; 0.789 of the 200-step reference-setting figure 0.4847).

### The round-4 profile (the seven keys on: B.2, B.3, gains 1 to 5)

The full-turn profiles of the rounds before repeated for fbbfaf9 (serial,
job 55973719) and c5c00c5 (8 and 64 ranks: c5c00c5 adds the parallel rule
of gain 3 only, below), the same cases, flags, events and readers
(prof/runs4/{serial,par8,par64}/<run>/hpctoolkit-database, the analyses
prof/an5/ with rules_swept5.py; csv/hp4_*.csv against 407f969, csv/hp4b_*.csv
against a11ebc1 + B). The class statistics with every key on (job
55973718, N 64, one turn, 0.3954 s/step with the timers): the sweeping
faces' volume clip 0.0681 -> 0.0529 (gain 5), their trace 0.0528 -> 0.0301
(gains 3 and 4), the no-flux faces' trace 0.0173 -> 0.0115, the four timed
parts 0.241 -> 0.167 s/step.

(a) The serial breakdown per step, N = 64, one turn (697 steps):

| category | 407f969 s/step | % | a11ebc1 + B | after (fbbfaf9 + B + gains) s/step | % | after / 407f969 |
|---|---|---|---|---|---|---|
| pieces: the clip by the cell planes and the PLIC plane (after: the one-plane closed form included) | 0.0660 | 13.4 | 0.0662 | 0.0459 | 11.3 | 0.70 |
| pieces: the cap closure and the volume | 0.0012 | 0.2 | 0.0010 | 0.0006 | 0.2 | 0.50 |
| face update: the pre-image fan clipped (the faces without a flux only, after B.2) | 0.0353 | 7.2 | 0.0137 | 0.0138 | 3.4 | 0.39 |
| swept cells: the ring search | 0.0127 | 2.6 | 0.0089 | 0.0083 | 2.0 | 0.65 |
| swept cells: the cell planes, state, wisp rule and the pre-test | 0.0291 + 0.0372 | 13.5 | 0.0543 | 0.0616 | 15.2 | 0.93 |
| swept cells: the shortcut, the check, the fast-path test (off) | 0.0013 | 0.3 | 0.0002 | 0.0003 | 0.1 | - |
| swept cells: the pre-image fan and the rest | 0.0122 | 2.5 | 0.0122 | 0.0125 | 3.1 | 1.02 |
| trace: the velocity at the traced points (findTetrahedron / the cached tetrahedra) | 0.0302 | 6.2 | 0.0293 | 0.0100 | 2.5 | 0.33 |
| trace: the cell search (ring, octree) | 0.0153 | 3.1 | 0.0146 | 0.0146 | 3.6 | 0.95 |
| trace: the midpoint trace, the triangle list | 0.0107 | 2.2 | 0.0111 | 0.0107 | 2.6 | 1.00 |
| the point velocity (volPointInterpolation of U / the lazy points) | 0.0223 | 4.5 | 0.0223 | 0.0019 | 0.5 | 0.09 |
| curvedFlux: the band, the wisp surfaces, the exchanges, the face loop | 0.0226 | 4.6 | 0.0234 | 0.0219 | 5.4 | 0.97 |
| the remote exchange (serial: the stencil setup) | 0.0031 | 0.6 | 0.0030 | 0.0031 | 0.8 | 1.00 |
| the surface build (fit, anchor, volume rule) | 0.0155 | 3.2 | 0.0158 | 0.0155 | 3.8 | 1.00 |
| the reconstruction (davof, ring normal, positioning) | 0.0299 | 6.1 | 0.0301 | 0.0299 | 7.4 | 1.00 |
| the bounding | 0.0113 | 2.3 | 0.0112 | 0.0109 | 2.7 | 0.96 |
| isoAdvector: timeIntegratedFlux | 0.0150 | 3.0 | 0.0149 | 0.0157 | 3.9 | 1.05 |
| isoAdvector: advect, the rest (alpha update, counters) | 0.0291 | 5.9 | 0.0290 | 0.0290 | 7.2 | 1.00 |
| setFlow | 0.0309 | 6.3 | 0.0305 | 0.0316 | 7.8 | 1.02 |
| the Courant number, other function objects | 0.0016 | 0.3 | 0.0016 | 0.0016 | 0.4 | 1.00 |
| I/O and start-up | 0.0218 | 4.4 | 0.0218 | 0.0199 | 4.9 | 0.91 |
| the solver loop (main) | 0.0355 | 7.2 | 0.0351 | 0.0350 | 8.7 | 0.99 |
| total | 0.4907 | 100 | 0.4586 | 0.4048 | 100 | 0.82 |

- The step 0.4907 -> 0.4586 -> 0.4048 s (the 200-step timings 0.4847 ->
  0.4531 -> 0.3823); N 32 0.0989 -> 0.0921 -> 0.0764 (28.6 s against 35.2
  for the turn). The ratio to plicRDF's serial step (0.2090): 1.94 (2.19
  after round 3b, 2.35 at 407f969, 5.3 in the first profile).
- What remains specific to the swept flux: 0.229 s/step (57 % of the step):
  the first pass (cell planes, states, the pre-test) 0.062, the volume
  clip 0.046, the trace 0.035, the band and face loop 0.022, the surface
  build 0.016, the face update of the no-flux faces 0.014, the fan and the
  rest 0.013, the ring search 0.008, the point velocity 0.002, the
  exchange setup 0.003, the shared parts of isoAdvector that only the
  swept path reaches (defaultAreaFractions, in timeIntegratedFlux) ~ 0.012.
  What is shared with plicRDF: 0.176 s/step (43 %): the reconstruction
  0.030 (plicRDF's 0.068), the bounding 0.011, the alpha update 0.029,
  setFlow 0.032, the Courant number, I/O 0.020, the solver loop 0.035 =
  84 % of plicRDF's whole step (0.209).
- The pre-test's first pass is now the largest swept item (15 %): the
  states of 2.04 (sweeping) and 4.09 (no-flux) touched cells per face with
  the plane of every interface cell (planeOf) rebuilt per face that
  touches it; a per-step cache of the ring cells' half-spaces and states
  (built once per cell, read per face) is the next structural saving
  (expected ~ 0.03 s/step, 7 %).

(c) Parallel, the comparison's scotch decompositions (no balancing):

| | N 64, 8 ranks, 407f969 | a11ebc1 + B | after (c5c00c5 + B + gains) | N 128, 64 ranks, 407f969 | a11ebc1 + B | after |
|---|---|---|---|---|---|---|
| step (CPUTIME per rank / steps; the wall clock) | 0.207 s (0.212) | 0.190 (0.194) | 0.171 s (0.176) | 0.508 s (0.519) | 0.471 (0.481) | 0.426 s (0.436) |
| MPI library, share of all CPUTIME | 64.5 % | 63.2 % | 62.7 % | 84.6 % | 83.5 % | 83.3 % |
| MPI per rank: mean / max / min | 93 / 123 / 27 s | 84 / 111 / 24 | 75 / 100 / 21 s | 603 / 674 / 332 s | 551 / 621 / 304 | 497 / 562 / 272 s |
| compute per rank: max / mean (max/mean) | 113 / 49 (2.31) | 104 / 46 (2.25) | 94 / 42 (2.23) | 361 / 101 (3.59) | 337 / 100 (3.38) | 307 / 91 (3.35) |
| the balanced step | 0.070 s | 0.066 | 0.061 s | 0.072 s | 0.071 | 0.065 s |
| the ratio to plicRDF's step at the same ranks (0.0556 / 0.1395 s) | 3.8 | 3.5 | 3.2 | 3.7 | 3.4 | 3.1 |

The busiest rank's compute fell by the serial factor again (113 -> 94 s at
8 ranks, 361 -> 307 s at 64), the step by 10 % and 10 % against a11ebc1
+ B; the imbalance (2.23 and 3.35) sets the step as before. In parallel
gain 3 computes the whole point field up front on every rank (c5c00c5: the
patch-point interpolation is collective, and the lazy form deadlocked at 8
ranks, job 55973720), so only gain 4's saving of the interpolationCellPoint
interpolation reaches the parallel step from that pair.

(d) The top exclusive hotspots after (serial, N = 64, one turn):

| # | % | s | function (file:line; fork fbbfaf9) | kind |
|---|---|---|---|---|
| 1 | 5.6 | 15.8 | sweptFaceCellsFlux's own lines (the pre-test's loops, the cell planes and states) | algorithmic |
| 2 | 3.7 | 10.4 | sweptPieceVolume's own lines (the per-triangle loop, the cap bookkeeping) | algorithmic |
| 3 | 3.5 | 9.8 | __mempcpy (libc): the field assignments of the solver loop and setFlow | redundant copies, mesh-wide |
| 4 | 3.4 | 9.5 | Foam::minusOp (ops.H:209): the subtractions of the general-plane distances (the PLIC planes) | algorithmic |
| 5 | 3.2 | 9.0 | davof::halfSpace::dist (davofSweptPieces.H:139): the axis-aligned distances | algorithmic |
| 6 | 2.8 | 7.8 | curvedFlux's own lines (the face loop) | algorithmic |
| 7 | 2.2 | 6.1 | fvc::surfaceIntegrate (the alpha update) | mesh-wide |
| 8 | 1.9 | 5.3 | surfaceInterpolationScheme::dotInterpolate (phi from U, setFlow) | mesh-wide |
| 9 | 1.9 | 5.2 | davof::clipWith: the clip proper | algorithmic |
| 10 | 1.8 | 5.1 | plusEqMagSqrOp2 + magSqr + mag (the edge lengths, the tolerances) | algorithmic |
| 11 | 1.7 | 4.9 | Foam::Vector::cross: areaVector, the triangle normals | algorithmic |
| 12 | 1.6 | 4.4 | setDefaultAreaFraction (isoAdvector: the faces without a plane) | shared, mesh-wide |

Total exclusive 282.1 s; the allocation 8.2 s (2.9 %; 13.0 s before: the
trace's tetrahedron lists are gone). findTetrahedron and
volPointInterpolation have left the list (hotspots 4 and 7 of round 3b);
the line table's largest entry is libOpenFOAM's Field multiply (5.8 %, the
solver loop's field algebra).

### Recommendation on the defaults

Bit for bit and a gain: sweptAxisPlanes (-2.8 %), sweptTetCache with
sweptPointVelocityLazy (-11 % together; the lazy key alone gains nothing
and in parallel computes the whole field up front), sweptPreTestCoplanarSide
(B.3, -4.6 %): these four should join the defaults; they reproduce 407f969
bit for bit on the eight cases and change no number anywhere a plane is
general (the axis test falls back to the dot product). Bit for bit without a
gain: sweptNoFluxLayer and sweptDonorFastPath (B.1): leave off (code paths
that add nothing on these meshes). Round-off: sweptFaceUpdateFromVolume
(B.2, -3.3 %) and sweptSliverClosedForm (-2.2 %): opt-in by the rules; the
differences are 1e-12 after 20 steps except the rotation N 64 (1e-9, the
rank-fallback amplification), the same order as the trace cache's first
form and the processor-boundary round-off of the parallel runs; Tomislav's
call whether a change of that size may enter the defaults.


## The port into leia (2026-11-06)

Tomislav's decision of 2026-11-06: the DAVOF method leaves the TwoPhaseFlow fork and
becomes leia's own library against OpenFOAM v2606's geometricVoF ("I do not want to make
TwoPhaseFlow 'famous' by my work, I want to have it in my own project"). The delta
analysis, the layout, the drop list and the acceptance plan are `docs/davof/PORT-DESIGN.md`
(phase 0, 0ec77adf). The fork at `c5c00c5` (branch `davof-swept`) stays the reference of
every number above; from here on new results are recorded against leia commits (the
version stamp of `libleiaDavof`, printed by every app).

### Phase 1: what was built (leia branch `davof-port`, 4dbd4803 .. 1761c0b9)

- `libleiaDavof` (`src/davof`, 20 500 lines moved from the fork, bit for bit where the
  acceptance below measures it): `reconstruction/` (the davof scheme,
  `davofCcsPositioning.H`), `normalModels/` (cellLocal, youngs, plicRDF, exact, ringArea,
  ringFit), `volumeRules/` (kromerBothe, exactSphere with the vendored MIT overlap library,
  splitQuadrature), `geometry/` (the swept pieces, the remote cells, the point velocity,
  the orbit, the curved geometry), `init/davofInit`, and `advection/`: the base
  `davofAdvectionScheme` (run-time selected by `advectionScheme`, registered as
  `davofAdvectionScheme` for the function objects), `davofSweptAdvection`
  (`advectionScheme davofSwept`: the fork's `isoAdvection` with the plane-clipped loop of
  the preprint v1, the face-clipping models, the one-donor option and the time-quadrature
  fluxes removed (PORT-DESIGN section 3), the keys `surface plane | paraboloid |
  exactSphere` and `flux sweptFaceCells | prism`, the fork's `clippingSurface` and
  `fluxIntegration` as deprecated aliases, the log prefix `davofSwept:`) and
  `isoAdvectorScheme` (`advectionScheme isoAdvector`: OpenFOAM's own `isoAdvection`, the
  comparison arm). The three geometry headers keep their one home in
  `src/leiaLevelSet/davof` (included by path); that library, the level-set side's state
  constructor, is `libleiaDavofState` now.
- `libleiaDavofFunctionObjects` (`src/functionObjects/davof`): areaFractionError,
  volumeFractionError, reconstructionError, linked into the solver.
- The apps `davofAdvectorFoam` (from `advectorVoF`), `davofInitFields`,
  `davofStaticFaceFractions`, `davofPlaneVolumeConsistency`, `davofTestSweptPieces`.
- The build gate, laptop and Lichtenberg, OpenFOAM v2606 (job 56002310, rc 0): no
  TwoPhaseFlow include or library in any Make file; `ldd` of the five apps and the two
  libraries lists OpenFOAM's and leia's libraries only; the DAVOF sources compile without
  a warning. `EIGEN_INC` is set by `etc/leia-env.sh`.
- The workflow: the route is the solver token alone (`solver: davofAdvectorFoam`, or
  `advectorVoF` for the fork through `LEIA_TPF_DIR`); the Snakefile hands `DAVOF_SOLVER`
  to `materialize.py`, which derives the three tokenised lines (`application`,
  `advectionScheme`, the function objects' `libs`); `initDavofFields` maps to
  `davofInitFields` on the leia route; both routes' solve output is `advectorVoF.csv`
  (the converter's file, which reads both logs). `port_render_diff.py` renders a config
  on both routes and found nothing but those three lines differing (`port_rotation_leia`,
  2 cases; `port_plane_leia`, 8 cases).
- Not a second key: the joint branch's `sweptHexClip` is the same device as `c5c00c5`'s
  round-4 `sweptAxisPlanes` (`halfSpace::axis`, the one coordinate difference, local and
  remote cell planes alike), found when extracting it; the port carries `sweptAxisPlanes`.
- The round-4 keys and the two round-off opt-ins are tokens of the three cases
  (`TPF_SWEPT_AXIS_PLANES`, `..._TET_CACHE`, `..._POINT_VELOCITY_LAZY`,
  `..._PRE_TEST_COPLANAR_SIDE`, `..._FACE_UPDATE_FROM_VOLUME`, `..._SLIVER_CLOSED_FORM`,
  the fork's defaults false); the plane case's `writePrecision` is a token.

### Phase 2, the identity set (PORT-DESIGN 4.1): MEASURED 2026-11-06

Workspace `/work/scratch/tm83tomy/davof/agents/port` (the fork at `c5c00c5` through
`LEIA_TPF_DIR`, leia `davof-port` 2033532e for the runs), the studies
`config/davof/port_{rotation,diagonal,plane}_{fork,leia}.yaml`: the reference setting
with the four round-4 keys on, `writePrecision 17`, serial; the rotation N 32 (20 steps)
and N 64 (40 steps), the lineDiagXZ diagonal N 32 (18 steps), the plane N 20 and 40 at
Co 0.5 (8 and 16 steps) for both rows and both arms. Judged by
`workflow/scripts/port_identity_check.py` (the written fields of the last time after the
FoamFile header, the function objects' files and the converter's CSV byte for byte, the
logs' numeric lines after the prefixes and the two renamed keys are normalised):

| study | cases | written fields | function-object files | converter CSV | log numeric lines |
|---|---|---|---|---|---|
| rotation N 32 | 1 | 9 fields identical | identical | identical but the 5 wall-clock columns | 212 of 212 identical |
| rotation N 64 | 1 | 9 fields identical | identical | wall-clock columns and `nAnchorEvals`, `nAnchorEdge` at ONE step (1 vs 0) | 410 of 412; the two lines of that step (below) |
| diagonal N 32 | 1 | 9 fields identical | identical | wall-clock columns | 192 of 192 |
| plane N 20, 40, both rows, both arms | 8 | the t = 0 fields only (the case writes no later time, `TPF_WRITE_STEPS`) | identical | wall-clock columns | 91 or 92 / 171 or 172 of as many, all identical |

The written fields are byte-identical in every case (`alpha.water`, `alphaf.water`,
`areaNormalDavof.water`, `pPlane.water`, the defect and sliver fields, `U`, `phi`), the
per-step errors of the function objects (areaFractionError's L1, L2, Linf, the face and
cell counts; volumeFractionError) identical to the last printed digit in all 11 cases.

The one difference, rotation N 64 at one step: the fork's reconstruction line says
"end-of-interval cells 5 (1 without a plane)", leia's "(0 without a plane)", the
positioning residual mean differs in the 11th digit (2.1213936337790828e-15 against
2.1213936338419322e-15), the anchor evaluated 1 cell (edge 1) in the fork and 0 in leia
(skipped 1702 against 1703), the swept flux clipped 129132 against 129136 triangles;
every field written after that step is identical. This is PORT-DESIGN 1.2 item 1 in its
mildest form: one cell whose positioned plane grazes the cell is classified by
`cutCellPLIC::calcSubCell` as no cut by the fork's threshold (`mag(faceArea) < 10*SMALL`,
TwoPhaseFlow's) and as a cut with a vanishing polygon by OpenFOAM v2606's
(`< ROOTVSMALL`); the cut volume is the same either way (the sliver's volume is below
the residual tolerance), so the plane, the anchor's result and the transported state do
not change; what changes is whether the cell counts as "without a plane", whether
`anchorSkipConverged` skips it, and four clipped triangles of a vanishing polygon.
Tomislav's rule for this case (2026-11-06): v2606's class stays the default; the
fork's threshold comes as an OPT-IN copy (`davofCutCellPLIC`) to prove that the
threshold is the whole difference, and the difference is documented as it is here:
not round-off in the state, which is identical, but a diagnostic count. The opt-in and
its proof run are the next item.

Open after this measurement: the opt-in proof (above); the plane case's last-time
fields were not compared (the template writes at `TPF_WRITE_STEPS` only; the per-step
CSV and log numbers were); the parallel identity (4.2), the ladders (4.3) and the Co 1
comparison (4.4), whose configs are in place (`port_parN32np*_{fork,leia}`,
`port_ladder_*_leia`, `port_rotation_leia_noAxis`).

### Phase 2, the acceptance (PORT-DESIGN 4.1 .. 4.4): MEASURED 2026-11-06

The same workspace; leia `davof-port` 5854a387 (the grazing-cut opt-in) and later for the
runs; every verdict by `workflow/scripts/port_identity_check.py` (bit for bit = the written
fields of the last time after the FoamFile header, every function-object file including the
per-rank `sweptFaceLog` files, the converter's CSV apart from its five wall-clock columns,
and every numeric line of the log after the prefixes, the two renamed keys and the MPI host
lines are normalised). The studies and the scoring scripts are in `config/davof/port_*.yaml`,
`workflow/scripts/port_{identity,parallel,ladder}_check.py`, `davof_cmp/port_cmp_an.py`.

| item | what ran | verdict |
|---|---|---|
| 4.1 the identity set, 11 cases, serial, 20 to 40 steps | `port_{rotation,diagonal,plane}_{fork,leia}` | bit for bit in 10 of 11; the rotation N 64 differs in two anchor counts of one step and in the diagnostic line of that step, every field identical (the grazing cut, PORT-DESIGN 1.2 item 1) |
| 4.1 the plane set with the last step's fields written | `port_planeW_{fork,leia}` (`TPF_WRITE_STEPS 8`) | bit for bit, 8 of 8 (9 fields each at t = 8e-4) |
| 4.1 the grazing-cut proof | `port_rotation_leia_grazingFork` (N 64, `grazingCut fork`) against the fork's N 64 case | bit for bit: the 9 fields, the 3 function-object files, the CSV, all 412 numeric log lines; against the default leia run the two counts (`nAnchorEvals` 1 against 0, `nAnchorEdge` 1 against 0) and that line differ, nothing else |
| 4.1 `sweptAxisPlanes` off against on (leia, N 32 and 64) | `port_rotation_leia_noAxis` | fields, function-object files and CSV identical; the log differs in the startup line that prints the key and in nothing else (the fork's claim for gain 2 holds in leia) |
| 4.2 the parallel identity, np 1, 2, 4, 8, 16, both rows, 20 steps | `port_parN32np*_{fork,leia}` | leia np k against fork np k bit for bit at every np (43 to 323 files per case, the per-rank swept-face logs included); leia np k against np 1 within 1e-12 of L1(alphaf) at every step but np 4 on the exact row (2.4e-12), the fork's own scatter digit for digit |
| 4.3 the ladders, 16 rungs (rotation N 32/64 serial, 128 on 64 ranks; sphere N 32/64, 128 at Co 0.1/0.2/0.5), both rows, a full turn / translation | `port_ladder_*_{fork,leia,leiaGF}` | the default route (v2606's `cutCellPLIC`) is bit for bit in 5 rungs and drifts at the last digit in 11 from the step of the first grazing cut on; the route with `grazingCut fork` is bit for bit in all 16 (the threshold is the whole difference over 160 to 1280 steps) |
| 4.4 the Co 1 comparison, 18 cells, both routes on one node, 3 repetitions of the ring row | `port_cmp*_leia` preparations, `davof_cmp/port-cmp-time.sbatch` | 16 cells scored (table below; the rotation N 64 and N 128 at OF Co 0.2 were still running at the time of writing, the N 128 one at 3255 s per 7000-step run needs one job per repetition as the fork's protocol had) |

The drift of the default route on the ladders, in the scored end norms (the preprint's
scorer, `davof_swept_ladders.py`, fork against leia; `grazingCut fork` reproduces the fork
in every digit): L1(alphaf) and E_shape identical to five digits in all 16 rungs; L1(n)
within 3.2e-3 relative on the rotation ring row at N 64 and within 3.4e-2 on the sphere
exact row at N 128 (the sliver-normal sensitivity of that measure, REPORT 2026-11-05
"L1(n) end +43 % on the exact row (sliver planes, measurement only)"); L1(x) within
9.2e-4 (sphere ring N 64), 4.0e-4 (rotation ring N 128); min alpha and the largest
excess at the 1e-13 .. 1e-11 level in both routes. This is the loop's known
round-off amplification (REPORT 2026-10-24: two serial runs 1e-13 apart differ by 2e-6 at
step 2 and hold 1e-9 .. 2e-7 afterwards) seeded by one classification of one grazing cut;
the state stays bounded and conservative in both routes.

The Co 1 comparison (PORT-DESIGN 4.4; `csv/port_cmp_table.csv`; the values of repetition 1,
the wall clock the median of three for the ring row): the DAVOF arm's end norms (L1(x),
E_shape, L1(n), min alpha, the excess, the mass) agree between the routes to the printed
digits in 12 of the 16 cells and within the ladder's drift in the others (LineX N 64 OF
Co 0.2 exact: L1(x) 1.9338e-4 against 1.9313e-4, L1(n) 1.5876e-3 against 1.5596e-3; LineX
N 128 OF Co 0.2 exact: L1(x) 5.2747e-5 against 5.2740e-5; rotation N 64 OF Co 0.2 ring:
1.0590e-3 against 1.0582e-3); the wall clock of the leia binary is 0.994 .. 1.015 of the
fork's for the DAVOF arm (the same code) and 0.94 .. 1.00 for isoAdvector (OpenFOAM
v2606's against TwoPhaseFlow's). The plicRDF + isoAdvector arm (OpenFOAM v2606's
`plicRDF` and `isoAdvection` against TwoPhaseFlow's, PORT-DESIGN 1.3): L1(x), E_shape and
L1(n) identical to seven digits in 15 of the 16 cells; min alpha and the excess differ at
the 1e-13 .. 1e-14 level (the bounding tolerances `100*SMALL` against `1e-12` / `10*SMALL`),
and in the one unbounded cell (LineDiag N 128, OF Co 1, 64 ranks: both versions lose the
boundedness, REPORT 2026-11-03) the numbers differ: L1(x) 3.0422e-3 against 2.9987e-3,
E_shape 3.2446e-4 against 3.2206e-4, the excess 3.52e-3 against 4.09e-3. The DAVOF arm
against plicRDF on the leia route: wall clock 1.46 .. 3.62 x, E_shape 0.45 .. 0.91 x on
the ring row (the fork's picture of 2026-11-03 unchanged).

The bit-for-bit verdict of the DAVOF arm's cells (`port_identity_check.py` on the first
repetition, the fields of the last written time per processor directory, the
function-object files, the CSV, the log): bit for bit in 15 of the 34 scored runs (every
LineDiag and LineX run at OF Co 1 at N 32, 64, 128, both rows; the LineDiag N 32 and the
LineX N 32 ring runs at OF Co 0.2); identical fields with the two anchor counts
(`nAnchorEvals`, `nAnchorEdge`) and the diagnostic line of a few steps different in 5
(LineDiag N 64 OF Co 0.2 both rows, rotation N 32 OF Co 0.2 both rows and OF Co 1 ring);
the fields different at the round-off level in the remaining 14 (the rotation at N 64 and
N 128, LineX at OF Co 0.2 from N 32 on, LineDiag N 128 OF Co 0.2: e.g. alpha
0.99999999994613487 against 0.99999999994613797, a wisp 3.199e-14 against 3.907e-14, the
LineX N 32 exact run only in `reconstructionError`'s L1(n) 1.9368e-4 against 1.9390e-4 and
L1(x) in the seventh digit). The pattern is the ladders': the default route separates from
the fork at the step of the first grazing cut and the loop amplifies the difference to
the 1e-11 .. 1e-4 relative level in the end norms; the rotation, with the most grazing
cuts, never stays bit for bit beyond N 32 at OF Co 0.2; the straight lines at OF Co 1 do.

Tomislav's rule (2026-11-06) is applied: v2606's `cutCellPLIC` stays the default; the
opt-in `grazingCut fork` (`davofCutCellPLIC`, a parameterised copy of v2606's class,
selected through `davofCellCut`) proves the threshold to be the whole difference; the
drift of the default is documented above with its numbers.

Open after phase 2: the rotation N 64 and N 128 cells at OF Co 0.2 of the comparison (one
job per repetition for N 128); the record in `docs/davof/PORT-DESIGN.md` of the measured
scatter as the acceptance tolerance of long runs on the default route.

### Phase 3, the Enright identity (PORT-DESIGN section 5): MEASURED 2026-11-06

The Enright agent's timed comparison (`enright-time.sbatch`, its branch `davof-swept-enright`
402fcc2e merged into `davof-port` 6a3cdcda, the templates tokenised) run from the leia-route
preparations `port_enright{Swept,RDF}N{32,64,128}_leia` on three routes by
`davof_enright/port-enright-time.sbatch`: `leia` (v2606's `cutCellPLIC`), `leiaGF` (the
opt-in `grazingCut fork`) and `fork` (the three lines rewritten); ring = the swept-face
advection with the ring normal, rdf = plicRDF + isoAdvector; N 32 serial, N 64 on 8 ranks,
N 128 on 64; OF Co 0.5 and 1 (the `co_index` 0 and 1 of the preparations); T = 3, the
abort at min(alpha) < -0.5 as in the Enright record below. Scored by
`davof_enright/port_enright_an.py` (`csv/port_enright_table.csv`):

| N | np | OF Co | ring, fork: steps, status | ring, leia (default): steps, status | ring, leiaGF: identity | rdf, v2606 against TwoPhaseFlow |
|---|---|---|---|---|---|---|
| 32 | 1 | 0.49 | 480, ended (E_shape 1.137e-2) | 443, aborted (min alpha -0.742) | bit for bit, 480 ended | same digits |
| 32 | 1 | 0.98 | 240, ended (1.250e-2) | 225, aborted (-0.797) | bit for bit, 240 ended | same digits |
| 64 | 8 | 0.50 | 628, aborted (-0.545) | 777, aborted (-1.12) | bit for bit, 628 aborted (wall 0.995) | E_shape 3.0181e-3 / 3.0181e-3, L1(x) 5.870e-2 / 5.854e-2, L1(n) 0.5963 / 0.5938 |
| 64 | 8 | 0.99 | 399, aborted (-1.20) | 138, aborted (-0.774) | bit for bit, 399 aborted | E_shape 3.1633e-3 / 3.1633e-3, L1(x) 6.087e-2 / 6.075e-2 |
| 128 | 64 | 0.50 | 125, aborted (-1.25) | 216, aborted (-0.750) | bit for bit, 125 aborted | same digits (E_shape 6.9114e-4, L1(x) 7.4270e-3) |
| 128 | 64 | 0.99 | 106, aborted (-0.508) | 67, aborted (-0.932) | bit for bit, 106 aborted | E_shape 6.908e-4 / 6.906e-4, L1(x) 7.622e-3 / 7.841e-3, excess 8.7e-5 / 9.2e-5 |

Bit for bit = the fields of the last written time in every processor directory (72 at N 64,
576 at N 128), the three function-object files, the CSV but its wall-clock columns, every
numeric log line and the abort line with its numbers; the per-rank `kromerBothe area check
failed` warnings reach the log in the ranks' flush order and are compared as sorted sets
(`port_identity_check.py` 2aa8c65a). The wall clock of the leiaGF binary is 0.994 .. 1.004
of the fork's.

What this says. The Enright case is the loop's sensitivity at its sharpest: the default
route and the fork separate at the first grazing cut (N 32 OF Co 0.5: step 6, "2 clips,
12627 against 12629 quadrature leaves", the mass change of the step 4.9e-16 against
6.1e-16), and from there the two runs are different realisations of the same unstable
fragmentation: the fork's N 32 runs end at T = 3 and the default route's abort at step 443
and 225; at N 64 OF Co 1 the fork aborts at step 399 and the default route at 138; at N 128
OF Co 0.5 at 125 against 216. The abort step of the Enright record below is therefore not a
property of the method but of the realisation, to be read with the round-off sensitivity
in mind (the loop amplifies a 1e-16 difference of one step, REPORT 2026-10-24). With
`grazingCut fork` every cell (6 of 6) is bit for bit, in serial and on 8 and 64 ranks: the
threshold of `cutCellPLIC::calcSubCell` is the whole difference between the port and the
fork on this case too. The plicRDF + isoAdvector arm differs between OpenFOAM v2606's and
TwoPhaseFlow's classes in the third to fourth digit of L1(x) and L1(n) at N 64 (both OF
Co) and N 128 OF Co 1, the unbounded runs (min alpha -4.8e-6 .. -1.1e-3, the bounding
tolerances differ, PORT-DESIGN 1.3), and to the printed digits at N 32 and N 128 OF Co 0.5;
leia's reference row is v2606's by Tomislav's decision, the difference is stated, not
resolved.

For the standing study matrix (Tomislav, 2026-11-07) this fixes two things: the DAVOF
rows run with `grazingCut fork` where a result is to be reproduced against the fork's
record and with the default otherwise, the key written into the stamp of every case; and
the Enright rows' abort step is reported with the realisation's seed (the route and the
threshold) named.

### The standing matrix: the mechanism and the configs (2026-11-07, before the matrix runs)

(1) Recompute only on change: `workflow/scripts/case_stamp.py` (README "Recompute only on
change: the case stamp"; PORT-DESIGN 5b), accepted on the cluster with
`config/davof/port_stamp_test.yaml` (two plane cases of seconds; jobs 56013329, 56013443,
56013459, 56013789, 56014057): a fresh study runs both cases and stores their stamps
(92e360124a45, c884e6ea4e15; the table carries `stampId`, `stampBinary` 547e31ad011d,
`stampGrazingCut`, `stampLeia`); the same launch again schedules nothing ("Nothing to be
done", the guard lists 0 solve jobs); the current stamp files touched schedule nothing
either (Snakemake 9.23 counts an input as updated only when its checksum differs from the
one recorded at the job's completion, `dag.py` `update_needrun`; measured); a changed
dictionary value (TPF_WRITE_PRECISION 16 to 17 through the config) recomputes both cases
with the reason "dictionaries changed" in the dry run and the guard, new ids (d78a09747c4f,
757ad626d9f7) and, after the config is restored, the original ids again (the id is
content-addressed) with the 17-results under `previous/d78a09747c4f/`. Two defects found
and fixed on the way (both measured on the test): the current stamp's record carried the
computation time and moved the file at every parse (rescheduling every finished case; the
file is now rewritten on a changed trigger only, 4f9cdd3d), and an executor job re-parses
the Snakefile, so two jobs raced on the stamp's temporary file (the driver alone writes,
ebf5c70c); a re-materialisation archives the result before it rebuilds the directory
(94b885d4).

(2) The cases. The fork agent's leia track (origin/davof-swept 95952849) merged into
`davof-port` (a71f887e): `cases/davof/isoShear2D` (the sheared disc of Rider and Kothe
1998 on the one-cell-thick mesh, T = 8), `isoShear3D` (the sphere in the shear flow of
Liovic et al. 2006, T = 3), `isoDeformation3D`, `isoRotation3D`, the rows `plicRDF-ccs`
and `davof-ringArea-vofCutCell`, the positioning per row, the tet mesh class, its runner
and tables. The pieces of the fork's `davof-isoAdvector` that the shear cases need, taken
over into leia (face94ab, compiled on the laptop): `davofInit::exactCylinder` (the exact
state of the disc; `davofInitFields`, `areaFractionError davofMethod exactCylinder`,
`volumeFractionError exactState auto` picks it for a cylinder on boxes), the explicit
plane position's Gauss sum closed over the empty faces (davof.C, the consistency app), the
apps `davofInitShearFlow3D` and `davofInitVortex2D` (the fork's `initShearFlow3D` and
`initVortex2D`, phi the exact line integral of the stream function in 2D); the Snakefile
maps the fork's initialiser names to leia's inside a setfields command, the stamp hashes
every app of it. The iso* templates tokenised for the route (the three lines) and extended
with the swept-face advection's reference-setting keys and `grazingCut` (the fork route's
isoAdvection does not read them), so that ONE template per case serves both tracks: on the
fork route the rendering differs from the fork agent's by the added keys and in nothing
else, on the leia route by the three lines (checked for isoShear2D and isoShear3D, both
rows). The repetition token `TPF_REPETITION` (a comment line of controlDict, so every
repetition is its own case with its own stamp) in the three standing cases.

(3) The configs, `config/davof/standing/<family>_N<N>.yaml`, written by
`workflow/scripts/make_standing_configs.py` (the matrix in one table): the families
shear2D (N 64/128/256 on 1/1/8 ranks, the fork agent's decomposition), shear3D (32/64/128
on 1/8/64) and deformation3D (the Enright agent's `sphereDeformation3D`, 32/64/128 on
1/8/64); the rows `davof-planeSurf-ringArea` (the swept-face advection, the reference
setting) and `plicRDF` (v2606's isoAdvector + plicRDF); OF Co 0.2, 0.5, 1; three
repetitions; `grazingCut fork` for the DAVOF rows; the pre-registration S2-1..5, S3-1..5,
D-1..4 in the headers. The tables: `workflow/scripts/make_standing_tables.py <family>`
(one table set per family: the cells, the orders over N, the caption with the stamps,
the ranks, dt and the steps). The smoke before the matrix: `config/davof/port_shear2D_smoke.yaml`
(the swept flux on the 2D mesh with empty patches, 17 steps).

The coordinator's decisions (2026-11-07, after the hand-back): one template per case on
both tracks = the fork agent's iso* templates (the deformation family switched to
`isoDeformation3D`; the Enright agent's `sphereDeformation3D` is its diagnosis case only);
the stamp id without the workflow commit stays; the matrix starts after the 2D smoke
passes, as the BASELINE (the two rows), the variant rows added later as new rows (the
stamps make them additive); the fork agent merges davof-port's template changes.

### The 2D smoke of the swept-face advection: PASSED (2026-11-06, job 56018102)

`config/davof/port_shear2D_smoke.yaml` on the leia route (leia 39cea3d7, built 56014921;
`profiles/slurm-small`, the serial rules at 4 GB: a 16 GB single-core request routes to the
`mem` partition, whose five nodes were draining, and the first submission's mesh jobs sat
pending for an hour): the sheared disc isoShear2D at N 64, OF Co 0.5 (measured 0.4988),
17 steps of T = 8, `grazingCut fork`, the rows `davof-planeSurf-ringArea` and `plicRDF`.
Both rows run through (17 steps, End, no FatalError, no abort). The DAVOF row on the
one-cell-thick mesh with two empty faces per cell: the swept flux "alphaf and dVf set on
532 faces (216 second ring, 0 boundary), 0 skipped, 0 clips", the fit "cells 76, converged
76, stalled 0, singular 0"; the explicit plane position's start defect
max |alpha_plane - alpha| = 0.0126 at t = 0 on the exact cylinder state (the fork agent
measured 0.48 before the empty faces closed the Gauss sum); at the last step min(alpha)
-6.5e-14, max(alpha) = 1 + 7.5e-13, the mass change of the step 5.9e-16 relative, the
converter's mass deviation 1.0e-13 %; plicRDF: min(alpha) -2.7e-19, max 1 + 0, mass
deviation 0. The function objects compare with the disc at rest (U0 = 0), so the numbers
at t = 0.1 are the deformation, not errors; they agree between the rows (E_shape 8.156e-4
against 8.146e-4, L1(n) 0.4091 against 0.4094), and L1(alphaf) of the DAVOF row is finite
(0.475, the face fractions of the moved disc against the initial one). Verdict: the swept
flux, the face update and the function objects work on the 2D mesh with empty patches; the
baseline matrix starts (the coordinator's decision 3).

### The Enright agent's decisions 1-3 in leia: the identity acceptance PASSED (2026-11-06)

The fork's `davof-swept-enright` c5c00c5..25e6ebc (the keys `sweptPieceVolume closure |
closureThenTets | signedTets` with `sweptPieceValidate`, `sweptBoundedFaceFraction none |
scale | recompute`, `interfaceCellRule normalOnly | normalAndAlpha | normalAndAlphaReset`,
the diagnostics `sweptFragmentLog`, `sweptReinitFaceFractions`, `sweptCoverageCheck`,
`sweptPieceDumpTol`; the coplanar tolerance of the tetrahedra clips, 25e6ebc) taken over
into leia (e1d11045: the eight diffs applied to the port's files with the class name and
the log prefix renamed, nothing else changed) and selected by one token
`TPF_ENRIGHT_SETTING` (tpf_rows.ENRIGHT_SETTINGS: reference, closureThenTets, signedTets,
signedTetsValidate, scale, recompute, normalAndAlpha, normalAndAlphaReset, fragmentLog,
recommended = signedTets + scale + normalAndAlpha) that renders the five keys as one block
appended to the `sweptPreTestCoplanarSide` line, EMPTY for reference (002ef32d): every
dictionary rendered before the keys existed is byte for byte the same (checked for four
cases against the templates of e1d11045), so the baseline matrix's stamps stand and the
variant rows are added, not recomputed.

The acceptance (`config/davof/port_enrightKeys_*`, make_port_enright_configs.py; the fork
route against `advectorVoF` built from 25e6ebc in the worktree `TwoPhaseFlow-enright` of
the port workspace, the leia route `davofAdvectorFoam` at 002ef32d; `grazingCut fork`;
jobs 56023550-63 and 56024911/12; `port_identity_check.py`, every file, every numeric log
line):

| study | cases per route | verdict |
|---|---|---|
| rotation N 32 (20 steps) and N 64 (40 steps), the ten settings | 20 | ALL IDENTICAL (with the default `grazingCut openfoam` first: 13 of 20, the 7 at N 64 the known grazing-cut separation, PORT-DESIGN 1.2) |
| diagonal translation N 32, the ten settings | 10 | ALL IDENTICAL |
| plane N 20 and 40, both arms, both rows, the ten settings | 80 | ALL IDENTICAL |
| parallel np 1, 2, 4, 8, the recommended setting, both rows | 2 each | ALL IDENTICAL; leia np k against np 1: max |dL1(alphaf)| 9.5e-14 (ring), 2.6e-12 (exact) at np 4 and 8, the fork's digits |

The validation counter (`sweptPieceValidate true` with signedTets): the largest
|V_closure - V_tets|/V_c over the run 1.8e-13 (rotation N 32, 6448 pieces per step
validated), 6.2e-13 (rotation N 64, 24035 per step), 9.0e-14 (diagonal), 3.6e-14 to
6.1e-13 (the plane, 2417 to 10670 per step): the signed tetrahedra and the closure agree
to round-off on the rotation and the translations (the closure's fourth-cap defect is on
the straining flows, the Enright agent 2026-11-07). The bounded face fractions on the
rotation N 64 (40 steps): `scale` touched 1071 faces, all partial, max |1 - s| 8.6e-7;
`recompute` 1084, max |1 - s| 4.5e-9. `normalAndAlpha`: one cell per step with |m_c|
above the threshold and alpha outside (tol, 1 - tol), no plane for it (the reference
gives it one). Wall clock not measured here (20 to 40 steps); the Enright agent's 2.35 to
2.6 x of signedTets stands.

The variant rows of the standing matrix: `config/davof/standing/<family>_N<N>_variants.yaml`
(the DAVOF row with recommended, signedTets, closureThenTets; `grazingCut fork`; three
repetitions), submitted after this acceptance; `make_standing_tables.py` labels the row
with the setting.

### The baseline matrix, first findings while it runs (2026-11-06)

- The DAVOF reference row on the 3D shear (isoShear3D) DIES with a floating-point
  exception in the bounding (`boundFlux` within `limitFluxes`) at N 32 OF Co 0.2 step 112
  of 1440 (all three repetitions, deterministic) and at N 64 OF Co 0.5 step 108 (rank 4),
  after a step with min(alpha) = -0.053 and "explicit plane defects max |alpha_plane -
  alpha| = 1.0004": the swept flux unbounded on the shear, then the bounding divides by
  zero. S3-1 of the pre-registration is FALSIFIED for the reference setting (it predicted
  completion at N 32 for every Co); the same code runs in the fork at 25e6ebc. Whether the
  recommended setting (signedTets + scale + normalAndAlpha) holds the bound on the shear
  is the variant rows' question. The 2D shear (isoShear2D) completes every cell at N 64,
  128 and the finished N 256 cells (13312 steps at OF Co 0.2 on 8 ranks).
- One plicRDF case (shear3D N 64 OF Co 0.2 rep 1, job 56020056) was CANCELLED by SLURM
  at step 575 (not a result); its files moved to `previous/cancelled-56020056/`, the
  study relaunched (56024903): the stamps run that case alone.
- The shear2D_N256 orchestrator (56018304) died in Snakemake's ILP scheduler
  (pulp/CBC `IndexError` reading an empty solution file) with 8 of 18 cases done;
  relaunched (56025308), the stamps skip the finished eight. The Makefile's launches
  take `--scheduler greedy` from now on (no CBC).
- The N 128 rungs of the 3D families wait for the `mem` partition (16 GB serial steps;
  its nodes draining since the afternoon); every other rung ran on deflt.

## The 3D Enright deformation (2026-11-06)

Tomislav's task of 2026-11-06: the 3D deformation test (LeVeque 1996, Enright et
al. 2002) with the swept-face advection in parallel against isoAdvector +
plicRDF, on the equidistant hexahedral meshes of the level-set case
`cases/3Ddeformation`. Agent clones: fork
`/work/scratch/tm83tomy/davof/agents/enright/TwoPhaseFlow` (branch
`davof-swept-enright` from `davof-swept` c5c00c5, no source change), leia
`.../enright/leia` (branch `davof-swept-enright` from `davof-swept` 930d2328);
the job ledger `.../enright/.my_jobs`.

### The case

`cases/davof/sphereDeformation3D` (new; the dictionaries of
`sphereRotation3D` with the setFlow block replaced): the unit cube with N^3
equidistant hexahedra, one patch `walls`; the sphere R = 0.15 at
(0.35 0.35 0.35); the velocity

    u =  2 sin^2(pi x) sin(2 pi y) sin(2 pi z) cos(pi t / T),
    v = -sin(2 pi x) sin^2(pi y) sin(2 pi z) cos(pi t / T),
    w = -sin(2 pi x) sin(2 pi y) sin^2(pi z) cos(pi t / T),  T = 3,

the 3D field of R. J. LeVeque, High-resolution conservative algorithms for
advection in incompressible flow, SIAM J. Numer. Anal. 33 (1996) 627-665,
doi:10.1137/0733033, with the time reversal of D. Enright, R. Fedkiw,
J. Ferziger, I. Mitchell, A hybrid particle level set method for improved
interface capturing, J. Comput. Phys. 183 (2002) 83-116,
doi:10.1006/jcph.2002.7166. The field is set by OpenFOAM's
`functionObjects::setFlow`, mode `vortex3D` (exactly the three components
above on the cell centres and the face centres, `phi = U(x_f) . S_f`), with
the time factor through its `scale` Function1 `cosine` of period 2T;
advectorVoF evaluates it at the half time of every step (as for the
rotation). The velocity vanishes on the six faces of the cube, so the
boundary flux is zero. The face-centre flux is EXACTLY divergence-free on
the equidistant hexahedra, not only to the midpoint rule: with
A(x) = 2 sin^2(pi x), D(x) = sin(2 pi x), B(y) = sin(2 pi y),
E(y) = sin^2(pi y), C(z) = sin(2 pi z), F(z) = sin^2(pi z) a cell's flux
sum is h^2 [(A(x1) - A(x0)) B(yc) C(zc) - D(xc) (E(y1) - E(y0)) C(zc)
- D(xc) B(yc) (F(z1) - F(z0))], and A(x1) - A(x0) = 2 D(xc) sin(pi h),
E(y1) - E(y0) = B(yc) sin(pi h), F(z1) - F(z0) = C(zc) sin(pi h): the
three terms cancel identically. setFlow's continuity error per step is
1e-13 (measured on the preparation runs). The level-set case
`cases/3Ddeformation` uses the same face-centre flux on hexahedra through
`velocityModel::setVolumetricFlux` (its Helmholtz projection is for
polyhedral meshes); advectorVoF's source `divU alpha` is round-off here.

The step: dt = Co h / 2.5 with 2.5 just above the maximum over the cells of
0.5 sum_f |U(x_f) . n_f| at t = 0 (2.460 / 2.475 / 2.478 at N = 32 / 64 /
128 on the three meshes; the continuum maximum of |u| + |v| + |w| is 2.479,
of |U| 2.000), so that the nominal Co is OpenFOAM's Courant number
(0.5 sum_f |phi_f| dt / V, the maximum over the cells, at the first step)
to within 1.6 % at N = 32 and 1 % above; T = 3 is 3 N 2.5 / Co steps, every
count a multiple of four, so that the field writes at T/4, T/2, 3T/4 and T
land on steps:

| N | ranks | Co 0.5: dt, steps | Co 1: dt, steps |
|---|---|---|---|
| 32 | 1 | 0.00625, 480 | 0.0125, 240 |
| 64 | 8 | 0.003125, 960 | 0.00625, 480 |
| 128 | 64 | 0.0015625, 1920 | 0.003125, 960 |

Rows: ours = the swept-face advection with the ring normal
(`davof-planeSurf-ringArea`) in the reference setting of 407f969
(sweptWispSurface matched / 1e-9 / 0, sweptConsistent, pyramid, ccs
positioning, anchorSkipConverged, sweptPreTest touching, sweptTraceCache)
plus the four bit-for-bit keys of round 4 ON as the study setting
(sweptAxisPlanes, sweptTetCache, sweptPointVelocityLazy,
sweptPreTestCoplanarSide; they reproduce 407f969 bit for bit on the eight
verification cases of round 4); the fork's bounding at its defaults
(nAlphaBounds 3, snapTol 0, clip false). The comparator = isoAdvector +
plicRDF with nAlphaBounds 5, snapTol 0, clip false, nAlphaCorr 1 (the
fork's deformation test `testsuite/advection/test-deformationSphere`; that
test runs Roenby's spiralling field in a 1 x 1 x 2 box and is not this case).
No exact-normal row: a deformed sphere has no exact normal. Same initial
alpha and alphaf (initDavofFields exactSphere), same mesh, dt, scotch
decomposition and ranks (1 / 8 / 64 at N = 32 / 64 / 128, as in the
comparison of 2026-11-03). The abort guard's interface-cell factor is 100
instead of the sphere cases' 10 (the sheet's area alone multiplies the
interface cells; the counts at T/2 and T are reported); its alpha bounds
are unchanged.

Scoring at t = T against the initial sphere (the exact end state; the
function objects' body at rest, U0 = 0, omega = 0): E_shape =
sum |alpha - alpha_exact| V (volumeFractionError, exactState exactSphere),
also as the fraction of the sphere volume 4/3 pi R^3 = 0.0141372 (for sharp
fields E_shape is the volume of the symmetric difference); L1(x) and L1(n)
of the recovered planes (reconstructionError at the write times); for ours
L1(alphaf) against the exact initial state (areaFractionError); min alpha
and the excess over the whole run and the mass deviation from advectorVoF's
per-step line; the interface cells and the wisps at T/2 and T; the swept
counters per step. Timing: the total wall clock per case at equal ranks,
three repetitions on exclusive i01 nodes, the diagnostics at the four write
times only for both methods.

Files: `cases/davof/sphereDeformation3D{,.parameter}`, the velocity arm
`deformation3D` in `workflow/scripts/tpf_rows.py` (U0 = 0), the six
preparation configs `config/davof/sphereDeformation3D_hex{Swept,RDF}N{32,64,128}.yaml`
(short END_TIME; the timed runner sets T), the timed runner
`workflow/scripts/davof_enright/enright-time.sbatch` and the scoring
`workflow/scripts/davof_enright/enright_an.py`.

### Pre-registered before any run

(The same text in `config/davof/sphereDeformation3D_hexSweptN32.yaml`.)

- E-1: both methods complete every case (6 cells each) and conserve the
  volume to 1e-10 % (the log's dev % at every step). Falsified by an abort
  (reported with the counters, not tuned around) or by |dev| > 1e-10 %.
- E-2: ours stays bounded to 1e-9 over the whole run (min alpha >= -1e-9,
  excess <= 1e-9) at both Co and every N. plicRDF's boundedness is
  measured; the expectation from the comparison of 2026-11-03 (unbounded on
  the diagonal at OF Co 1 in parallel): plicRDF leaves [0, 1] by more than
  1e-6 at Co 1 for N >= 64. Falsified for ours by any excursion beyond 1e-9.
- E-3: E_shape at T: the ratio ours / plicRDF in [0.7, 1.5] at every
  N >= 64 and both Co (the ring normal is first order, plicRDF's second; the
  swept flux is bounded where isoAdvector's is not; on the rotation ours was
  7-21 % better in L1(x), on the diagonal at Co 0.2 4-12 % worse).
  Falsified if the ratio exceeds 2 at N = 128 for either Co: the thin sheet
  (thinner than h at T/2 for N <= 64) is where the one-layer ring stencil
  and the wisp rule are stressed, and a ratio above 2 would mean the ring
  normal on the sheet loses more than the swept flux gains.
- E-4: the orders over N of E_shape at T (log2 of the ratio on 32 -> 64 and
  64 -> 128): both methods between 1 and 2.5 on 64 -> 128 at both Co (the
  sphere returns under-resolved at N = 32, so 32 -> 64 may exceed 2).
  Falsified for ours by an order below 1 on 64 -> 128.
- E-5: the wall clock ours / plicRDF at equal ranks between 2 and 4
  (Tomislav's prediction; the reference setting measured 1.9 x serial at
  N 64 and 3.1 x on 64 ranks on the rotation, where the interface is a
  sphere; the sheet at T/2 multiplies the interface cells, and the swept
  flux's cost follows the interface band while isoAdvector's is mostly
  mesh-wide). Falsified if the median ratio exceeds 6 at any N.

- E-6 (the third row, added by Tomislav on 2026-11-06 after E-1 .. E-5
  were registered and before any of its runs: isoAdvector's own plane flux
  with the DAVOF ring normal and the CCS positioning, the fork's
  plane-clipped loop = the preprint's davof-ringArea row plus ccs:
  reconstructionScheme davof, normalModel ringArea, positioning ccs,
  faceClipping plane, clippingSurface none, fitAnchor faces,
  fluxIntegration prism, nAlphaBounds 3, snapTol 0, clip false, the swept
  keys inert; configs `sphereDeformation3D_hexPrismN{32,64,128}.yaml`, the
  method `prism` of the timed runner and the scoring): the expectation from
  the rotation (it diverged there at OF Co 1 and was within 4 % of the
  swept loop at Co 0.2): the plane loop loses boundedness beyond 1e-6 at
  Co 1 for every N and ends by the abort guard or unbounded; at Co 0.5 it
  completes with E_shape within a factor 1.3 of the swept row; its wall
  clock is within 1.5 x plicRDF's (the flux is isoAdvector's). Falsified if
  it stays bounded to 1e-9 at Co 1 for every N.

Also reported, not predicted: L1(x), L1(n), L1(alphaf) at T, the
interface-cell and wisp counts at T/2 and T, the swept counters per step
(rank fallbacks, beyond-layer faces, pyramid delta max, twisted faces, wisp
pieces, bounding passes), setFlow's continuity error, the pictures of the
alpha = 0.5 surface at T/2 and T.

Changes after the registration and before the recorded runs (Tomislav's
notes of 2026-11-06 through the coordinator): writePrecision 17 for the
fields, the function objects and the log (OpenFOAM sets the stream
precision from it), so that the leia re-run after the port can be checked
bit for bit (a first timed pass at precision 12 was cancelled after a few
minutes and kept as `time-p12/`; its observations are repeated below from
the recorded runs); the controlDict application from the token
`TPF_APPLICATION` and the timed runner's `APP` / `ENRIGHT_W` / `STUDIES`
variables, so that only the solver changes for the leia app; the third row.

### Results (2026-11-06)

Jobs (ledger `/work/scratch/tm83tomy/davof/agents/enright/.my_jobs`): build
55999984 (fork c5c00c5); preparation studies at writePrecision 17
56000993-56000998 (Swept, RDF) and 56001210-56001212 (Prism); the timed cells
56001326 (N 32, Co 0.5), 56001327 (N 32, Co 1), 56001328 (N 64, Co 0.5),
56001329 (N 64, Co 1), 56001330 (N 128, Co 0.5), 56001331 (N 128, Co 1), one
exclusive i01 node each (Xeon Platinum 9242), three repetitions of the three
rows, the order rotated; the diagnostics 56000705 (diag1) and 56000828
(diag2); the pictures 56002022 / 56003985 / 56006509. The precision-12 first
pass (56000603-56000608, cancelled after a few minutes) is kept as `time-p12/`.
Every number below is from `csv/enright_runs.csv`, `csv/enright_table.csv`,
`csv/enright_orders.csv` and `csv/enright_counters_N<N>_co<i>_<row>.csv`
(the per-step counters of the davof rows) under the agent workspace; the
cases with their logs, fields and function-object output at writePrecision
17 under `time/N<N>_co<i>/<row>_r<rep>/` (`.source` names the prepared case
they were copied from). The three repetitions of every cell reproduce their
digits (`repsIdentical` true for all 18 cells).

Measured OpenFOAM Courant numbers (the first step): 0.492 / 0.984 at N = 32,
0.495 / 0.990 at N = 64, 0.496 / 0.991 at N = 128. setFlow's continuity
error 4.6e-14 .. 2.6e-13 at every step. Mass: the largest |dev| of
sum(alpha V) over every run of every row 5.0e-12 % (plicRDF N 32),
4.2e-12 % (ours); within E-1's 1e-10 % everywhere.

**Completion and bounds (the whole run; the abort guard's bounds are
[-0.5, 1.5] and 100 x the initial interface cells).**

| N | Co | ours (swept, ring) | plicRDF | prism (plane loop, ring, ccs) |
|---|---|---|---|---|
| 32 | 0.5 | completes; min alpha -2.70e-1, excess 2.95e-4 | completes; -2.3e-17, 4.4e-16 | completes; -2.41e-1, 2.86e-2 |
| 32 | 1 | completes; -1.40e-1, 9.17e-4 | completes; -3.9e-13, 4.4e-16 | completes; -4.68e-1, 4.11e-2 |
| 64 | 0.5 | ABORT step 628, t = 1.9625: min alpha -0.5446 (excess 8.96e-3; 51138 interface cells of 1694) | completes; -4.76e-6 (first below -1e-9 at step 47, t = 0.147), 4.4e-16 | ABORT step 495, t = 1.547: interface cells 169410 > 100 x 1694 (min -2.36e-2 there, over the run -2.31e-1, excess 7.85e-2) |
| 64 | 1 | ABORT step 399, t = 2.4938: min alpha -1.2037 (excess 1.87e-2; 56463 cells) | completes; -5.20e-4 (from step 1), 3.94e-7 | ABORT step 264, t = 1.65: interface cells 169703 (min over the run -4.57e-1, excess 8.59e-2) |
| 128 | 0.5 | ABORT step 125, t = 0.1953: min alpha -1.2535 (excess 5.02e-5; 66741 cells of 6836) | completes; -3.23e-5 (from step 25, t = 0.039), 2.46e-6 | ABORT step 234, t = 0.3656: interface cells 684394 (min over the run -3.72e-1, excess 3.98e-2) |
| 128 | 1 | ABORT step 106, t = 0.3313: min alpha -0.5081 (excess 3.25e-2; 82892 cells) | completes; -1.10e-3 (from step 2), 9.20e-5 | ABORT step 131, t = 0.4094: min alpha -0.5261 (excess 1.56e-1; 390192 cells) |

**Errors at t = T against the initial sphere (the rows that reached T; V_s
= 0.0141372).**

| N | Co | row | E_shape | E_shape / V_s | L1(x) | L1(n) | L1(alphaf) | interface cells T/2, T | partial cells T/2, T | wisps T/2, T |
|---|---|---|---|---|---|---|---|---|---|---|
| 32 | 0.5 | ours | 1.137e-2 | 0.804 | 0.306 | 1.460 | 7.06e-2 | 6019, 5230 | 6638, 5871 | 619, 641 |
| 32 | 0.5 | plicRDF | 7.765e-3 | 0.549 | 0.0845 | 0.692 | - | 1324, 747 | 1324, 747 | 0, 0 |
| 32 | 0.5 | prism | 1.441e-2 | 1.019 | 0.326 | 1.196 | 1.46e-1 | 6328, 5328 | 6866, 5359 | 538, 31 |
| 32 | 1 | ours | 1.250e-2 | 0.884 | 0.313 | 1.479 | 7.88e-2 | 6659, 5305 | 7348, 6136 | 689, 831 |
| 32 | 1 | plicRDF | 8.009e-3 | 0.567 | 0.0952 | 0.749 | - | 1327, 833 | 1327, 833 | 0, 0 |
| 32 | 1 | prism | 1.755e-2 | 1.241 | 0.311 | 1.188 | 1.61e-1 | 5343, 5102 | 5597, 4837 | 254, -265 |
| 64 | 0.5 | plicRDF | 3.018e-3 | 0.2135 | 0.0585 | 0.594 | - | 6245, 2736 | 6245, 2736 | 0, 0 |
| 64 | 1 | plicRDF | 3.163e-3 | 0.2238 | 0.0608 | 0.631 | - | 6291, 2973 | 6291, 2973 | 0, 0 |
| 128 | 0.5 | plicRDF | 6.911e-4 | 0.0489 | 7.43e-3 | 0.187 | - | 27650, 7724 | 27650, 7724 | 0, 0 |
| 128 | 1 | plicRDF | 6.906e-4 | 0.0489 | 7.84e-3 | 0.207 | - | 27645, 8181 | 27645, 8181 | 0, 0 |

The interface cells are reconstructionError's count at the write (for the
davof rows the cells with |m_c| > wispTol V^(2/3), for plicRDF the cells
inside (surfCellTol, 1 - surfCellTol)); the partial cells are the cells of
the written field with 1e-8 < alpha < 1 - 1e-8; the wisps their difference
(negative where face-fraction interface cells carry alpha outside the band).
The davof rows' own per-step count (by the face fractions, alpha
unrestricted) is larger: ours 10813 / 10367 at T at N = 32 (of 32768
cells), the prism row 23198 / 23590. At N = 64 the last writes before the
aborts: ours E_shape 2.105e-2 at T/2 (Co 0.5; 35644 interface cells of
262144) and 2.868e-2 at 3T/4 (Co 1); the prism row 2.225e-2 / 2.247e-2 at
T/2 (43372 / 36310 cells); plicRDF's E_shape at T/2 is the sheet against
the sphere and not an error. Orders of plicRDF (E_shape; L1(x); L1(n)):
32 -> 64: 1.36; 0.53; 0.22 (Co 0.5), 1.34; 0.65; 0.25 (Co 1); 64 -> 128:
2.13; 2.98; 1.66 (Co 0.5), 2.20; 2.95; 1.61 (Co 1). No order for the davof
rows (one rung).

**Wall clock (median of three; the aborted rows to their abort, with the
steps reached).**

| N | ranks | Co | ours s (steps) | plicRDF s (steps) | prism s (steps) | ours / plicRDF | per step ours / plicRDF / prism s |
|---|---|---|---|---|---|---|---|
| 32 | 1 | 0.5 | 819.5 (480) | 26.1 (480) | 72.4 (480) | 31.4 | 1.89 / 0.050 / 0.160 |
| 32 | 1 | 1 | 394.3 (240) | 14.8 (240) | 35.0 (240) | 26.7 | 1.885 / 0.050 / 0.150 |
| 64 | 8 | 0.5 | 1450.8 (628 of 960) | 81.1 (960) | 83.6 (495) | 17.9 (partial) | 2.11 / 0.080 / 0.180 |
| 64 | 8 | 1 | 908.1 (399 of 480) | 46.2 (480) | 46.2 (264) | 19.7 (partial) | 2.17 / 0.080 / 0.180 |
| 128 | 64 | 0.5 | 264.3 (125 of 1920) | 354.5 (1920) | 49.8 (234) | 0.75 (6.5 % of the run) | 1.71 / 0.170 / 0.160 |
| 128 | 64 | 1 | 267.3 (106 of 960) | 210.9 (960) | 32.7 (131) | 1.27 (11 % of the run) | 2.15 / 0.190 / 0.210 |

The per-step ratio of ours to plicRDF is 38 (N 32), 26-27 (N 64) and 10-11
(N 128, the first 6-11 % of the run); the preparation runs on the undeformed
sphere gave 3 / 3.9 / 3.3 (0.09 / 0.29 / 1.48 s against 0.03 / 0.075 /
0.445 s per step at N 32 serial / 64 on 8 / 128 on 64 ranks). The growth is
the fragmentation: the swept flux's cost follows the number of cells with
face-fraction surfaces, which grows 20-fold before T/4 (below).

**The swept counters (rep 1; `enright_counters_*_ring.csv`).** Rank
fallbacks (the one-donor clip of a face whose piece has more than three
caps) from step 1: 463 of 4716 swept faces at N 32 (9.8 %), 994 of 18103
at N 64 (5.5 %), 1699 of 70814 at N 128 (2.4 %), at both Co alike (458 at
Co 0.25 in diag1: the share is set by the oblique 3D motion of the
pre-images, not by the step); they grow with the band to 4.5-4.9e3 (N 32),
2.0-2.5e4 (N 64), 4.9-5.2e4 (N 128) per step, 5 % of the swept faces
throughout. Beyond-layer faces (the uncovered pre-image parts): 0 at every
step except 199 (N 32, Co 0.5, max share 0.00), 12 (N 32, Co 1), 12897 (N
64, Co 0.5) and 2658 (N 64, Co 1) in total, 0 at N 128. Twisted faces 0
everywhere; the pyramid delta at most 0.053 h (N 32, Co 1), 0.013 h (N 64),
0.0034 h (N 128); the anchor never failed; the volume rule fell back
7104 / 6348 (N 32), 114454 / 88718 (N 64), 21787 / 35326 (N 128) times in
total. The second-layer bounding corrects cells from step 1 (25 / 38 at
N 32, 559 / 518 at N 64, 3300 / 3277 at N 128) and then 1.4-1.8e4 (N 32),
0.9-1.4e5 (N 64), 1.9-2.2e5 (N 128) cells per step: the bounded alpha is
held by the bounding, and the excess it cannot remove is the unboundedness
above. The wisp pieces without a neighbour surface grow from 0 at step 1 to
1.4-1.6e4 (N 32), 5.6-5.8e4 (N 64), 2.0-2.5e5 (N 128) per step at their
peak. The face-fraction interface cells grow 410 -> 3625 / 7155 at T/4
(N 32, Co 1 / 0.5), 1694 -> 38713 / 48181 at T/4 (N 64), 6836 -> 66741 /
82892 at the aborts (N 128, t = 0.20 / 0.33).

**Controls (not production rows, N 32 serial, precision 12; diag1 job
56000705, diag2 56000828).** diag1: the four round-4 keys off reproduce
the keys-on run bit for bit (alpha at t = 0.1 identical, the same counters
at every step); Co 0.25 has the same rank-fallback share (458 of 4711 at
step 1) and the same growth of the wisp pieces; sweptCapSet crossing has 0
degenerate pieces, 517 rank fallbacks and min alpha -7.6e-7 by step 3
(final: -2.6e-6 by step 8). diag2 (the full T): the plane-clipped loop with
the ring normal (= the prism row) ends unbounded (-0.278 / -0.135 at Co 1 /
0.5, E_shape 1.755e-2 / 1.441e-2, the same digits as the timed prism row);
the swept flux with plicRDF's normal in the davof scheme (normalModel
plicRDF, everything else ours) aborts at step 26, t = 0.325, at Co 1
(min -0.665) and completes at Co 0.5 with min -6.9e-3, excess 5.3e-3,
E_shape 8.706e-3 (plicRDF 7.765e-3), L1(x) 0.240, L1(alphaf) 5.44e-2 and
1924 interface cells at T (plicRDF 747).

**Pictures** (`/work/scratch/tm83tomy/davof/agents/enright/pics/`,
`N<N>_co<i>_<row>_t<t>.png` for every written time of every row and
`diag_<control>_t<t>.png`; the alpha = 0.5 surface by marching cubes on
the reconstructed fields, `render_alpha_mc.py`; the cluster's ParaView
5.13.3 modules carry no pvpython, `render_sheet.py` is kept for a build
that does): plicRDF at N 128 shows the sheet at T/2 (a thin shell with one
hole, 27650 interface cells) and a sphere with a small lobe at T; ours at
N 32 shows at T/2 a blunt body with detached fragments (7348 partial
cells) and at T a dented blob with two satellites; the rows that aborted
have the writes before their abort only (ours at N 64: T/4 and T/2 at
Co 0.5, T/4 .. 3T/4 at Co 1; nothing at N 128).

### Verdicts

- E-1 FALSIFIED for ours and for the prism row: ours aborts in the four
  cells with N >= 64 (min alpha below -0.5 at t = 0.20, 0.33, 1.96, 2.49),
  the prism row in the same four (three by the interface-cell factor, one
  by min alpha); both complete at N = 32 only. plicRDF completes every
  cell. The mass part holds for every row (5e-12 %).
- E-2 FALSIFIED for ours: min alpha -0.14 .. -1.25, excess 3e-4 .. 3e-2
  (N 32: -0.27 / -0.14 at Co 0.5 / 1). plicRDF: bounded to 4e-13 at N 32;
  below -1e-6 at N >= 64 at BOTH Co (-4.8e-6 / -5.2e-4 at N 64, -3.2e-5 /
  -1.1e-3 at N 128), the Co 1 part as expected, the Co 0.5 part beyond the
  expectation.
- E-3 FALSIFIED by absence: no value of ours at N >= 64. At N 32 the ratio
  is 1.46 (Co 0.5) and 1.56 (Co 1), inside the registered band but for the
  under-resolved rung the band excluded.
- E-4 not evaluable for ours (one rung). plicRDF: E_shape 1.35 then 2.1-2.2.
- E-5 FALSIFIED: the median wall clock ratio is 31 / 27 at N 32, 18-20 for
  the partial runs at N 64, and the per-step ratio 38 / 26 / 10-11 at N
  32 / 64 / 128 against the 2-4 predicted; the cause is the fragmentation
  (the band of face-fraction surfaces grows 20-fold), not the per-face
  cost (3-4 x on the undeformed sphere).
- E-6 (prism row): the loss of boundedness at Co 1 held (every N); the
  completion at Co 0.5 within 1.3 x of the swept row is falsified at
  N >= 64 (aborts) and at N 32 (1.27 x of ours, but 1.86 x plicRDF); the
  wall clock within 1.5 x plicRDF held at N >= 64 (1.0 x, aborted runs) and
  not at N 32 (2.8 / 2.4 x).

### What the counters say about the failure

The swept flux itself does what it did on the rotation and the lines: no
twisted faces, the pyramid correction below 0.053 h, no anchor failures, no
uncovered pre-image parts at N 128. Two things are new under the straining
field. (1) From step 1 a fixed share of the swept faces (2-10 %, falling
with N, independent of Co) has a piece with more than three caps and takes
the one-donor clip, the rule the rotation study showed to be unbounded;
the second-layer bounding has to correct tens to thousands of cells at
step 1 already. (2) The face-fraction state fragments: cells with a
face-fraction surface multiply 20-fold before T/4 while the alpha field's
partial cells do not (N 32 at T/2: 6638 partial cells against plicRDF's
1324; the davof count by the faces alone 7155 at T/4), the wisp pieces
without a neighbour surface rise from none to 1e4-2e5 per step, and the
rank fallbacks and the bounding grow with that band. The prism row, which
shares the davof state and the ring normal but uses isoAdvector's own
flux, fragments faster (23590 face-fraction interface cells of 32768 at
T at N 32; 684394 of 2.1 M at N 128 at t = 0.37) and is as unbounded, and
the swept flux with plicRDF's normal in the davof scheme still aborts at
Co 1: the fragmentation sits in the transported face-fraction state under
straining (the face update by the moved plane / the clipped pieces leaves
faces inconsistent with their cells, m_c = -sum alphaf S_f then marks cells
as interface cells that the alpha field does not have), and the normal
and the flux only change how fast. Not shown here: which of the two (the
face update or the one-donor fallback) starts it; a per-cell history of
|m_c| against alpha on the first fragmenting cells would tell. Nothing was
tuned; the rows ran as registered.

Commits (leia, branch `davof-swept-enright`): df4a8805 (the case, the arm,
the six configs with E-1 .. E-5, the timed runner, the report's setup and
registration), 54638749 (the scoring and the pvpython script), dbecf5fa /
796f4e6d (the exact divergence), 4180c110 / 9bcd55e2 (writePrecision 17,
TPF_APPLICATION, the runner's variables, the marching-cubes render),
8d849ca5 / 72d3f42e / 72210ab4 (the render job), 402fcc2e (the third row,
E-6), 704d709a (the wisps), and this one (the results). The fork branch
`davof-swept-enright` is c5c00c5 without a change (nothing to push). The
templates and configs a re-run needs: `cases/davof/sphereDeformation3D/`,
`cases/davof/sphereDeformation3D.parameter` (every key with its value),
`config/davof/sphereDeformation3D_hex{Swept,RDF,Prism}N{32,64,128}.yaml`;
only `solver` and `TPF_APPLICATION` change for the leia app.

### The diagnosis of the fragmentation (2026-11-06, Tomislav's follow-up)

Short runs only: N = 32, Co 0.5, serial, the full T = 3 (N = 64 only where a
control needs it); fork branch `davof-swept-enright` 5d8e57d (three opt-in
keys, default false, the default path unchanged, built from c5c00c5):

- `sweptFragmentLog`: the per-cell history (1). At the start of every step,
  after the reconstruction, every cell that is a face-fraction interface
  cell (|m_c| > wispTol V^(2/3)) while alpha < 1e-6 or alpha > 1 - 1e-6
  and was not logged before is recorded once: the step, t, the cell and
  its centre, alpha, |m_c|/V^(2/3), the six alphaf in cell-face order,
  whether a face of the cell was computed by the one-donor rank fallback
  (more than three caps or a rank-deficient cap system) in the previous
  step or the one before, whether the cell was in the second-layer
  bounding set then, and whether it was a wisp without a neighbour
  surface then (`postProcessing/sweptFragmentLog/cells_<rank>.csv`, the
  first 200 per rank); the class counts of every step in the log.
- `sweptExactCaps` (control a): where the closure's cap system is not
  solvable (k > 3 caps, or a PLIC plane parallel to a cell face) the piece
  volume is taken exactly instead of the one-donor fallback
  (`davof::sweptPieceExactVolume`, davofSweptPieces.H): the closed
  triangle list of the swept volume is decomposed into the signed
  tetrahedra about the face centre (Zhang & Chen, ICIP 2001,
  doi:10.1109/ICIP.2001.958278), each tetrahedron is clipped as a convex
  polyhedron by the cell's half-spaces and the PLIC half-space face by
  face (Sutherland & Hodgman, Commun. ACM 17 (1974) 32-42,
  doi:10.1145/360767.360802; the cap of each clip from the crossing points
  ordered about their centroid), the volume of each clipped polyhedron by
  the divergence theorem; sum_T vol(tet_T ∩ K) = vol(P ∩ K) for every
  convex K by the winding number, so no cap system is solved. The 2D face
  update (no caps) is unchanged. The ordinary pieces (k <= 3) keep the
  closure path bit for bit.
- `sweptReinitFaceFractions` (controls c, d): at the start of every step
  (from step 2) the face fractions are set from the previous
  reconstruction's planes (cutFacePLIC: the mean of the two clips on a
  face between two interface cells, the one clip by one plane, the mean of
  the two alphas on a face without a plane; the boundary faces by the
  owner), the transported alphaf discarded, the wet-face list rebuilt; the
  reconstruction then forms m_c and the ring normal from those. One
  reconstruction per step; the planes are those positioned on the previous
  alpha, which the next step's planes correct.
- control (e): wispTol 1e-2 in the reconstruction's dictionary (the
  initialisation and the function objects keep 1e-3).
- The runner `workflow/scripts/davof_enright/enright-diag.sbatch` (six
  runs on one node: frag = the production row + the log, exactCaps,
  reinit, reinitPrism = the prism row + reinit, wispTol, exactReinit = a
  + c together as the cross-check of D-2 and D-3), the scoring
  `enright_diag_an.py` (csv/diag_runs_diag3.csv, diag_frag_*.csv,
  diag_history_*.csv).

Pre-registered before the runs:

- D-1 (the seeds): the first fragmenting cells appear by step 3; of the
  first 200, more than 60 % have a face computed by the rank fallback in
  the previous step or the one before, fewer than 20 % belong to none of
  the three classes; more than 80 % are dry cells (alpha < 1e-6) that
  gained a face-fraction surface, with |m_c|/V^(2/3) between 1e-3 and 1e-1.
  Falsified if the fallback-adjacent share is below 40 % or the "none"
  share above 40 %.
- D-2 (a): the rank fallbacks are 0 at every step; the second-layer
  bounding at step 1 drops below 10 cells (from 25); min alpha stays above
  -1e-9 over the first 20 steps; over T the run does NOT stay bounded
  (min alpha below -1e-3 by T/2) and the face-fraction interface cells
  still exceed twice the partial cells at T/2: the fallback contributes,
  it is not the seed. Falsified if the run stays bounded to 1e-9 over T
  (then the fallback IS the seed) or if the first 20 steps already leave
  -1e-9 (then the exact caps do not remove the step-1 inconsistency).
- D-3 (c): the face-fraction interface cells stay within 1.5 x the
  partial cells at every write; min alpha over the run above -1e-2 (ten
  times better than the row's -0.27); E_shape at T within 1.5 x plicRDF's
  7.765e-3 and the sphere returns in the picture. Falsified if the ratio
  exceeds 2 at any write or min alpha falls below -0.1: then the loop's
  geometry (the flux and its fallback), not the transported state,
  fragments.
- D-4 (d): as D-3 for the prism row but less bounded: min alpha above
  -5e-2, the ratio within 1.5 x.
- D-5 (e): the interface count at T/2 falls by 30-50 %, min alpha stays
  below -0.1: a sensitivity, not a cure.
- D-6: no control keeps N = 32 bounded to 1e-9 over T. If (c) or (a + c)
  keeps N = 32 above -1e-2, that control is run at N = 64, Co 0.5 (the
  prepared 8-rank case) and is expected to complete without an abort.

## The isoAdvector package: the fork branch `davof-isoAdvector` (2026-11-06)

Tomislav's decision of 2026-11-06: the method (the swept-face advection, the
paraboloid dual advection, the volume rules, the swept geometry) is ported into
leia (branch `davof-port`); the fork keeps ONLY what the isoAdvector improvement
needs and runs the verification that demonstrates it. The fork branch
`davof-isoAdvector` starts from `davof-swept` at `c5c00c5` and holds one commit,
`fab3805` (pushed to `tmaric`); `davof-swept` is untouched (the history). The
package is defined in the fork's `etc/README-davof.md` ("The package"):

- (a) the DAVOF reconstruction for isoAdvector: the `davof` scheme (the normal
  from the face fractions, the plane positioned by alpha, `positioning` default
  `ccs`), the normal models `ringArea` (the stable model), `cellLocal`, `youngs`,
  `exact` and `plicRDF` (the controls), `ringFit`; the face update of
  isoAdvector's own plane sweep (`storeAreaFractions`, the second-ring sweep,
  `planeSpeaks`/`planeFlux`, the content bound, the second bounding layer, the
  defaults, the processor sync of `alphaf` with `dVf`); `davofInit` with the
  simplex, quadratic-face and sphere geometry (`detrixheAslam`,
  `quadraticFaces`, `exactSphere`); `davofOrbit.H`.
- (b) the CCS positioning for every plane-positioning reconstruction scheme:
  `surfaceIteratorPLIC::vofCutCell` (plicRDF, gradAlpha) reads `positioning
  vofCutCell | ccs`, `positioningTol` (1e-14 of V_c), `positioningMaxIter` (30)
  from the scheme's dictionary; the default `vofCutCell` leaves the code and the
  logs unchanged; `ccs` calls `davof::ccsPositionPlane` on the same `cutCellPLIC`
  cut from the linear estimate between the vertex extrema and prints one Info
  line per reconstruction. `isoAlpha` positions an iso-value of the
  point-interpolated alpha (`surfaceIteratorIso`), not a plane: `cutCellIso`
  gives no slope of the cut volume, so `isoAlpha` is unchanged (the third
  reference row).
- (c) the function objects `areaFractionError`, `volumeFractionError`
  (`exactState exactSphere`), `reconstructionError`; the apps `advectorVoF`
  (`setFlow` at the half step, the abort guard), `initDavofFields`,
  `planeVolumeConsistency`, `staticFaceFractions`, `testCutFaceAdvectEndArea`;
  the benchmark dictionaries `run/benchmark/{advection,reconstruction}`.

Removed (history on `davof-swept`): `isoAdvectionSwept.C` and every `swept*`
key; `isoAdvectionCurved.C` (the paraboloid dual advection, the quadrature
fluxes, the volume anchor; keys `clippingSurface`, `fluxIntegration`, `fit*`,
`anchor*`, `curved*`, `tensorCentroid`, `planeSource`, `defectLevels`,
`writeCurvedFields`); the volume rules with the vendored overlap library and the
Eigen dependency (`volumeRule`); the swept geometry (`davofSweptPieces.H`,
`davofSweptRemote.H`, `davofVelocityAt.H`, `davofCurvedGeometry.H`); the
face-clipping models with the curvature tensor of the `davof` scheme
(`faceClipping`, `tensor*`, `writeCurvatureTensor`); the vector face state of
`davof-warped` (`faceState vector`, `wettedAreaDavof`, in `davof`,
`isoAdvection`, `davofInit`, `initDavofFields`, `areaFractionError`); the chord
of `cutFaceAdvect` (`keepChord`, `chordEnd`, `edgeCrossings`);
`apps/test/testSweptPieces`. The diff of `fab3805` against upstream
`origin/master` (`ba10702`): 225 files, +21006/-83, i.e. the package plus the
`inline namespace twoPhaseFlow` fix of `feature/parasitic-currents` (the
3-line diffs of the cellCuts, surfaceIterators, reconstructionSchemes,
thermoDynamics files and the testsuite controlDicts). `isoAdvection.H` has
upstream's include list again; `isoAdvectionTemplates.C` keeps only the second
bounding layer and the default face fractions. Build: `Allwmake` clean on the
laptop (libVoF, postProcessing, every app and solver) and on Lichtenberg (job
56001915: `iso/TwoPhaseFlow` at `fab3805` in 384 s, `iso/TwoPhaseFlow-swept` at
`c5c00c5` in 389 s, both with the default flags, worktrees of the main clone).

### The identity check (job 56001915, mpsc0039, 2026-11-06)

Keys off = the same code path: 20-step copies (writePrecision 17) of the
plane-clipped loop's rows run with both binaries, the written fields compared
with `cmp` and by max |difference|, the per-step log lines diffed
(`prof/iso-ident.sbatch`; the cases in `iso/ident/{A,B}`). The davof cases are
the comparison's swept cases with `clippingSurface none`, `fluxIntegration
prism`, `fitAnchor faces` (as `ref/cmpC`); the plicRDF and isoAlpha cases lose
the keys `positioning*` (the study templates write `positioning ccs` for every
row, which the new binary's plicRDF would read; see the note below).

| case | row | source study, N, OF Co | ranks | field files identical | alpha / alphaf differing values | per-step log lines differing |
|---|---|---|---|---|---|---|
| R32 | davof-ringArea | hexCmpRotSweptN32, 32, 1 | 1 | 11 of 11 | 0 / 0 | 0 |
| R32x | davof-exact | hexCmpRotSweptN32, 32, 1 | 1 | 11 of 11 | 0 / 0 | 0 |
| D32 | davof-ringArea | hexCmpLineDiagSweptN32, 32, 1 | 1 | 11 of 11 | 0 / 0 | 0 |
| P20n | davof-ringArea | hexSweptPyr, 20, 0.5, normal arm | 1 | 11 of 11 | 0 / 0 | 0 |
| P20o | davof-ringArea | hexSweptPyr, 20, 0.5, oblique arm | 1 | 11 of 11 | 0 / 0 | 0 |
| RDF32 | plicRDF | hexCmpRotRDFN32, 32, 1 | 1 | 6 of 6 | 0 / - | 0 |
| ISO32 | isoAlpha | hexComparators, 32, 0.2 | 1 | 5 of 5 | 0 / - | 0 |
| R32p8 | davof-ringArea | as R32, simple (2 2 2) | 8 | 11 of 11 | 0 / 0 | 0 |

Every row reproduces the `davof-swept` binary bit for bit (min(alpha) after the
bounding equal to the last digit in every log, e.g. R32 -8.0394031111856562e-15,
D32 -0.099564557019122008: the diagonal line at OF Co 1 loses the bound on both
binaries alike, as the comparison C of 2026-11-05 found).

Smoke test of (b), plicRDF with `positioning ccs` (RDF32ccs, the new binary):
runs, 1248 positionings per reconstruction, 5.6 cuts per cell (max 9), the
positioning residual max 8.7e-15, mean 1.1e-16 of V. Against the vofCutCell row
after 20 steps alpha differs in 2519 cells, max |dalpha| 2.8e-4: NOT round-off.
vofCutCell positions to `isoFaceTol` 1e-8 in alpha (the cubic through four cuts,
Newton, then the secant), ccs to 1e-14 of V, and plicRDF's normal iteration (the
RDF from the plane centres) carries the plane differences into the normals. The
pre-registered expectation "ccs cheaper than vofCutCell with identical errors" is
therefore measured in step 2 as the difference of the two rows' errors against
the exact state (L1(n), L1(x), E_shape, with vofCutCell's own plane residual
reported by the same diagnostic), not as a field identity.

Note for the verification's templates: `positioning` is read from the alpha
dictionary by whichever scheme is selected (davof since 2026-11-04, plicRDF and
gradAlpha since `fab3805`); the token `TPF_POSITIONING` must be set per row
(plicRDF: `vofCutCell` for the classical setting, `ccs` for the ccs row),
not study-wide.

## The isoAdvector improvement: the verification in the fork (2026-11-06)

The fork branch `davof-isoAdvector` (fab3805 the package, c699e9e `initShearFlow3D` and the
jet case, 65ad6ae the jet README and comparison scripts, 690a023 / 4da5799 / b09b800 /
ae08dc0 the inflow face fraction, fe5ff45 the `initShearFlow3D` fix, b467390 the 2D sheared
disc: `exactCylinder` and `initVortex2D`) against upstream isoAdvector + plicRDF, on
Lichtenberg (exclusive i01 nodes, 96 cores), every row with isoAdvector's own flux and
bounding (`nAlphaBounds 3`, `snapTol 0`, `clip false`), the same initial state
(`initDavofFields exactSphere`), mesh, dt, decomposition (scotch) and ranks. Rows: RDF =
plicRDF with `positioning vofCutCell` (the classical setting), RDFccs = plicRDF with `ccs`,
RING = `davof` with the ring-area normal and `ccs` (the package's row), RINGvcc = the same
with `vofCutCell`, ISO = isoAlpha. Three repetitions per cell (the order rotated), the
diagnostics at the writes only (quarter turns / T/4), one profiled run per cell (OpenFOAM's
profiling, rank 0) outside the timing. The pre-registration F-1 .. F-5 is in leia
`config/davof/iso*3D_hexN32.yaml` (committed 1bf91bcd before any run); the tables are
`iso/csv/{iso_table,iso_orders,iso_prof,iso_runs}.csv` on the cluster
(`/work/scratch/tm83tomy/davof/agents/swept/iso/`, the cases under `time/`), the scripts
`workflow/scripts/davof_iso/{iso-time.sbatch,iso_an.py,iso_prof.py}`. The hexahedral
rotation, shear and deformation ran with the binary of c699e9e (`iso/TwoPhaseFlow`; the inflow
rule that came later changes nothing on these cases: `zeroGradient` inlets). The inflow face
fraction (Tomislav's rule of 2026-11-06: at a liquid inlet the face fraction is the alpha
condition's value, as the flux is) went through four versions, each checked against
davof-swept c5c00c5 on the eight verification cases (R32, R32x, D32, P20n, P20o, RDF32,
ISO32, R32p8 on 8 ranks; `iso/ident2..5`): the condition's value on every swept inflow face
(690a023) moved the oblique plane gate P20o by 1.8e-2 in alpha within 20 steps; the liquid
inlet's faces (alpha_b = 1) in the sweep, the second ring and the defaults (4da5799, and
b09b800 with the full cells' faces excluded) by 5.7e-6 (312 cells; the plane's cut of an inlet
face the condition calls liquid is below 1 where the plane crosses the inlet); the value in
the defaults alone (ae08dc0: a face no plane sets) leaves all eight cases unchanged bit for
bit (`iso/ident5`, job 56013017: every field file identical, no per-step log line differing).
The other seven cases were unchanged bit for bit under every version. The three families
below do not depend on the rule.

### The solid-body rotation, hexahedral (one turn; interface Co 0.046 / 0.116 / 0.23 = OpenFOAM Co 0.2 / 0.5 / 1)

Wall clock in seconds (the median of three repetitions; the spread max-min over the median
is 0.1-2.5 % on every cell), the per-step time from the run's wall clock; E_shape, L1(n),
L1(x) at the end of the turn against the exact sphere; min(alpha) over the run.

| N (ranks) | OF Co | row | status | steps | E_shape | E/RDF | L1(n) | L1(n)/RDF | L1(x) | min alpha | mass dev. % | wall s | ms/step | step/RDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 (1) | 0.2 | RDF | completed | 1728 | 1.155e-3 | 1 | 0.0575 | 1 | 4.39e-3 | -7e-20 | 6e-8 | 47.7 | 27.6 | 1 |
| | | RDFccs | completed | 1728 | 1.155e-3 | 1.00 | 0.0575 | 1.00 | 4.39e-3 | -7e-20 | 9e-8 | 47.0 | 27.2 | 0.99 |
| | | RING | completed | 1728 | 1.091e-3 | 0.94 | 0.0594 | 1.03 | 4.11e-3 | -6.9e-13 | 4.0e-4 | 39.9 | 23.1 | 0.84 |
| | | RINGvcc | completed | 1728 | 1.090e-3 | 0.94 | 0.0587 | 1.02 | 4.10e-3 | -9.4e-13 | 2.7e-4 | 42.5 | 24.6 | 0.89 |
| | | ISO | completed | 1728 | 1.172e-3 | 1.01 | 0.0737 | 1.28 | 4.37e-3 | -2e-20 | 5e-8 | 33.4 | 19.3 | 0.70 |
| 32 (1) | 0.5 | RDF | completed | 691 | 1.039e-3 | 1 | 0.0541 | 1 | 3.92e-3 | -3e-18 | 6e-8 | 21.6 | 31.3 | 1 |
| | | RDFccs | completed | 691 | 1.038e-3 | 1.00 | 0.0542 | 1.00 | 3.92e-3 | -6e-19 | 9e-8 | 21.5 | 31.2 | 1.00 |
| | | RING | ABORT step 497 (4270 interface cells of 416) | 497 | - | - | - | - | - | -0.064 | 0.083 | 14.6 | 29.3 | 0.94 |
| | | RINGvcc | completed, unbounded | 691 | 1.065e-3 | 1.02 | 0.361 | 6.7 | 0.041 | -0.071 | 2.1e-3 | 20.2 | 29.2 | 0.93 |
| | | ISO | completed | 691 | 1.052e-3 | 1.01 | 0.0713 | 1.32 | 3.99e-3 | -3e-18 | 5e-8 | 15.9 | 23.0 | 0.73 |
| 32 (1) | 1 | RDF | completed | 346 | 8.478e-4 | 1 | 0.0544 | 1 | 3.33e-3 | -1e-17 | 5e-8 | 13.1 | 37.9 | 1 |
| | | RDFccs | completed | 346 | 8.467e-4 | 1.00 | 0.0545 | 1.00 | 3.32e-3 | -1e-17 | 9e-8 | 13.1 | 37.8 | 1.00 |
| | | RING | ABORT step 111 | 111 | - | - | - | - | - | -0.21 | 0.33 | 6.0 | 54.1 | 1.43 |
| | | RINGvcc | ABORT step 231 | 231 | - | - | - | - | - | -0.22 | 76 | 10.6 | 45.7 | 1.21 |
| | | ISO | completed | 346 | 8.746e-4 | 1.03 | 0.0720 | 1.32 | 3.38e-3 | -5e-18 | 5e-8 | 10.1 | 29.1 | 0.77 |
| 64 (8) | 0.2 | RDF | completed | 3456 | 3.127e-4 | 1 | 0.0241 | 1 | 1.139e-3 | -5e-19 | 2e-8 | 191.2 | 55.3 | 1 |
| | | RDFccs | completed | 3456 | 3.127e-4 | 1.00 | 0.0240 | 1.00 | 1.138e-3 | -5e-19 | 2e-8 | 189.8 | 54.9 | 0.99 |
| | | RING | completed | 3456 | 3.016e-4 | 0.96 | 0.0256 | 1.06 | 1.087e-3 | -3.3e-6 | 1.7e-4 | 127.3 | 36.8 | 0.67 |
| | | RINGvcc | completed | 3456 | 3.016e-4 | 0.96 | 0.0256 | 1.06 | 1.085e-3 | -9.9e-13 | 1.7e-4 | 139.8 | 40.5 | 0.73 |
| | | ISO | completed | 3456 | 3.154e-4 | 1.01 | 0.0375 | 1.56 | 1.137e-3 | -1e-19 | 2e-8 | 103.4 | 29.9 | 0.54 |
| 64 (8) | 0.5 | RDF | completed | 1383 | 2.789e-4 | 1 | 0.0224 | 1 | 1.007e-3 | -5e-18 | 2e-8 | 80.4 | 58.1 | 1 |
| | | RDFccs | completed | 1383 | 2.789e-4 | 1.00 | 0.0224 | 1.00 | 1.007e-3 | -5e-18 | 2e-8 | 80.1 | 57.9 | 1.00 |
| | | RING | ABORT step 588 (18159 of 1760) | 588 | - | - | - | - | - | -0.086 | 1.1e-3 | 27.3 | 46.4 | 0.80 |
| | | RINGvcc | ABORT step 1364 | 1364 | - | - | - | - | - | -0.15 | 0.022 | 61.4 | 45.0 | 0.77 |
| | | ISO | completed | 1383 | 2.827e-4 | 1.01 | 0.0361 | 1.61 | 1.016e-3 | -3e-18 | 1e-8 | 44.6 | 32.3 | 0.55 |
| 64 (8) | 1 | RDF | completed | 691 | 2.287e-4 | 1 | 0.0218 | 1 | 8.39e-4 | -1e-17 | 1e-8 | 43.8 | 63.4 | 1 |
| | | RDFccs | completed | 691 | 2.285e-4 | 1.00 | 0.0220 | 1.01 | 8.42e-4 | -2e-17 | 2e-8 | 43.9 | 63.5 | 1.00 |
| | | RING | ABORT step 156 | 156 | - | - | - | - | - | -0.19 | 0.033 | 11.6 | 74.4 | 1.17 |
| | | RINGvcc | ABORT step 292 | 292 | - | - | - | - | - | -0.12 | 29 | 19.7 | 67.5 | 1.06 |
| | | ISO | completed | 691 | 2.336e-4 | 1.02 | 0.0355 | 1.63 | 8.60e-4 | -1e-17 | 1e-8 | 25.5 | 36.8 | 0.58 |
| 128 (64) | 0.2 | RDF | completed | 6911 | 8.364e-5 | 1 | 9.59e-3 | 1 | 2.937e-4 | -6e-13 | 5e-9 | 947.1 | 137 | 1 |
| | | RDFccs | completed | 6911 | 8.364e-5 | 1.00 | 9.60e-3 | 1.00 | 2.937e-4 | -1e-18 | 7e-9 | 947.1 | 137 | 1.00 |
| | | RING | ABORT step 670 (67049 of 6680) | 670 | - | - | - | - | - | -0.035 | 1.6e-3 | 69.6 | 104 | 0.76 |
| | | RINGvcc | ABORT step 2131 | 2131 | - | - | - | - | - | -0.039 | 7.5e-3 | 223.2 | 105 | 0.76 |
| | | ISO | completed | 6911 | 8.304e-5 | 0.99 | 0.0202 | 2.11 | 2.935e-4 | -6e-13 | 5e-9 | 488.1 | 70.6 | 0.52 |
| 128 (64) | 0.5 | RDF | completed | 2766 | 7.447e-5 | 1 | 9.28e-3 | 1 | 2.614e-4 | -1e-17 | 4e-9 | 389.1 | 141 | 1 |
| | | RDFccs | completed | 2766 | 7.446e-5 | 1.00 | 9.28e-3 | 1.00 | 2.614e-4 | -1e-17 | 6e-9 | 389.3 | 141 | 1.00 |
| | | RING | ABORT step 835 | 835 | - | - | - | - | - | -0.11 | 1.4e-3 | 85.8 | 103 | 0.73 |
| | | RINGvcc | completed, min -1.2e-3 | 2766 | 7.367e-5 | 0.99 | 0.0109 | 1.17 | 2.590e-4 | -1.2e-3 | 2.9e-4 | 283.4 | 102 | 0.73 |
| | | ISO | completed | 2766 | 7.410e-5 | 1.00 | 0.0203 | 2.19 | 2.630e-4 | -3e-13 | 4e-9 | 203.1 | 73.4 | 0.52 |
| 128 (64) | 1 | RDF | completed, min -1.5e-4 | 1382 | 6.065e-5 | 1 | 9.34e-3 | 1 | 2.221e-4 | -1.5e-4 | 7e-8 | 201.4 | 146 | 1 |
| | | RDFccs | completed | 1382 | 6.066e-5 | 1.00 | 9.43e-3 | 1.01 | 2.227e-4 | -1.1e-4 | 5e-8 | 203.0 | 147 | 1.01 |
| | | RING | ABORT step 334 | 334 | - | - | - | - | - | -0.16 | 0.026 | 46.7 | 140 | 0.96 |
| | | RINGvcc | ABORT step 508 | 508 | - | - | - | - | - | -0.15 | 33 | 68.6 | 135 | 0.93 |
| | | ISO | completed | 1382 | 6.093e-5 | 1.00 | 0.0201 | 2.15 | 2.254e-4 | -4e-5 | 1e-8 | 108.0 | 78.2 | 0.54 |

Orders over N (E_shape; L1(n)) at OF Co 0.2: RDF 1.89 / 1.90 (32-64 / 64-128); 1.26 / 1.33;
RING 1.85; 1.21 (32-64, the only completed pair); ISO 1.89 / 1.93; 0.97 / 0.89. ISO's normal
is first order and 1.3-2.2 x plicRDF's, its shape error equal to plicRDF's: the shape error
of isoAdvector on this case is set by the flux, not by the normal.

The profiled runs (rank 0; the time of the step and its parts per step, ms; the share of
`reconstruct()` in `time.run()`):

| N (ranks) | OF Co | row | step | reconstruct | share | advect | flux sweep | ring normal | 2nd ring | content bound | defaults |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 (1) | 0.2 | RDF | 25.7 | 13.3 | 52 % | 18.7 | 2.5 | | | | |
| | | RDFccs | 25.3 | 13.0 | 51 % | 18.3 | 2.5 | | | | |
| | | RING | 21.2 | 5.1 | 24 % | 13.8 | 3.4 | 0.43 | 0.63 | 0.05 | 2.2 |
| | | RINGvcc | 22.9 | 6.7 | 30 % | 15.4 | 3.4 | 0.43 | 0.62 | 0.06 | 2.2 |
| | | ISO | 17.9 | 4.8 | 27 % | 10.7 | 2.4 | | | | |
| 64 (8) | 0.2 | RDF | 54.3 | 35.7 | 66 % | 45.4 | 5.2 | | | | |
| | | RDFccs | 53.8 | 35.2 | 65 % | 44.9 | 5.1 | | | | |
| | | RING | 35.6 | 11.2 | 32 % | 26.0 | 7.5 | 3.6 | 0.89 | 0.10 | 2.3 |
| | | RINGvcc | 39.0 | 14.7 | 38 % | 29.4 | 7.5 | 3.5 | 0.89 | 0.10 | 2.2 |
| | | ISO | 28.8 | 7.0 | 24 % | 19.9 | 7.9 | | | | |
| 128 (64) | 0.2 | RDF | 136 | 86.0 | 63 % | 113 | 13.1 | | | | |
| | | RDFccs | 136 | 86.5 | 64 % | 113 | 12.8 | | | | |
| | | RING (670 steps) | 93.5 | 28.0 | 30 % | 66.0 | 19.4 | 0.003 | 0.008 | 0.68 |
| | | RINGvcc (2131 steps) | 101 | 35.4 | 35 % | 73.0 | 20.4 | 0.002 | 0.008 | 0.69 |
| | | ISO | 69.6 | 5.2 | 8 % | 46.4 | 26.1 | | | | |

(The cells of Co 0.5 and 1 are in `iso_prof.csv`; the N 128 rows are the resubmitted
profiled runs of 2026-11-06, 12:30.) The DAVOF reconstruction costs 0.38 x (N 32), 0.31 x
(N 64) and 0.33 x (N 128, up to the abort) of plicRDF's: 2.6 to 3.2 x cheaper, not the
registered fifth; the ring normal is a third of it at N 64 (3.6 of 11.2 ms, the vertex-ring
exchange), the ccs positioning the rest (6.7 - 5.1 = 1.6 ms less than vofCutCell at N 32).
The face update adds the plane sweep's extra faces (the flux sweep 3.4 vs 2.5 ms at N 32,
7.5 vs 5.2 at N 64, 19.4 vs 13.1 at N 128), the second ring (0.6-0.9 ms at N <= 64, 0.003
ms at N 128 on 64 ranks), the content bound (0.1 ms) and the defaults (2.2 ms at N <= 64,
0.7 ms at N 128: a pass over the faces whose value can differ): 3.6-6.8 ms per step. The
step of the RING row is 0.84 x (N 32, serial), 0.67 x (N 64, 8 ranks) and 0.69 x (N 128,
64 ranks, up to the abort) of plicRDF's wall clock; the whole step's saving is the
reconstruction's. isoAlpha is cheaper still (0.70 / 0.54 / 0.51 x).

### The 3D deformation (LeVeque / Enright), hexahedral (T = 3; OF Co 0.2 / 0.5 / 1)

| N (ranks) | OF Co | row | status | steps | E_shape | E/RDF | L1(n) | L1(x) | min alpha | interface cells (first / max) | wall s | ms/step | step/RDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 (1) | 0.2 | RDF | completed | 1200 | 7.93e-3 | 1 | 0.692 | 0.079 | -1.2e-8 | | 62.0 | 51.6 | 1 |
| | | RDFccs | completed | 1200 | 7.93e-3 | 1.00 | 0.690 | 0.080 | -1.2e-7 | | 60.5 | 50.4 | 0.98 |
| | | RING | completed, unbounded | 1200 | 1.21e-2 | 1.5 | 1.20 | 0.355 | -0.12 | 410 / 26628 | 197.7 | 165 | 3.2 |
| | | RINGvcc | completed, unbounded | 1200 | 1.25e-2 | 1.6 | 1.20 | 0.363 | -0.11 | 410 / 27359 | 190.3 | 159 | 3.1 |
| | | ISO | completed | 1200 | 9.07e-3 | 1.14 | 0.742 | 0.077 | -3e-17 | | 43.2 | 36.0 | 0.70 |
| 32 (1) | 0.5 | RDF | completed | 480 | 7.78e-3 | 1 | 0.719 | 0.084 | -1.6e-4 | | 26.6 | 55.3 | 1 |
| | | RDFccs | completed | 480 | 7.77e-3 | 1.00 | 0.695 | 0.083 | -1.6e-4 | | 26.6 | 55.3 | 1.00 |
| | | RING | completed, unbounded | 480 | 1.44e-2 | 1.9 | 1.20 | 0.326 | -0.24 | 410 / 24350 | 71.8 | 149 | 2.7 |
| | | RINGvcc | completed, unbounded | 480 | 1.49e-2 | 1.9 | 1.19 | 0.340 | -0.22 | 410 / 23692 | 64.8 | 135 | 2.4 |
| | | ISO | completed | 480 | 9.04e-3 | 1.16 | 0.753 | 0.079 | -1e-9 | | 18.9 | 39.4 | 0.71 |
| 32 (1) | 1 | RDF | completed | 240 | 8.06e-3 | 1 | 0.783 | 0.101 | -5.5e-4 | | 15.0 | 62.5 | 1 |
| | | RDFccs | completed | 240 | 7.95e-3 | 0.99 | 0.766 | 0.098 | -5.0e-4 | | 14.9 | 61.9 | 0.99 |
| | | RING | completed, unbounded | 240 | 1.76e-2 | 2.2 | 1.19 | 0.311 | -0.47 | 410 / 23590 | 34.8 | 145 | 2.3 |
| | | RINGvcc | completed, unbounded | 240 | 1.77e-2 | 2.2 | 1.19 | 0.275 | -0.47 | 410 / 22259 | 28.6 | 119 | 1.9 |
| | | ISO | completed | 240 | 8.56e-3 | 1.06 | 0.802 | 0.104 | -3.8e-4 | | 11.0 | 45.6 | 0.73 |
| 64 (8) | 0.2 | RDF | completed | 2400 | 3.07e-3 | 1 | 0.618 | 0.061 | -5.9e-5 | | 194.1 | 80.9 | 1 |
| | | RDFccs | completed | 2400 | 3.05e-3 | 0.99 | 0.608 | 0.060 | -5.9e-5 | | 190.1 | 79.2 | 0.98 |
| | | RING | ABORT step 1222 (t 1.53) | 1222 | - | - | - | - | -0.094 | 1694 / 169630 | 203.1 | 166 | 2.1 |
| | | RINGvcc | ABORT step 1211 (t 1.51) | 1211 | - | - | - | - | -0.11 | 1694 / 170240 | 201.8 | 167 | 2.1 |
| | | ISO | completed | 2400 | 3.27e-3 | 1.06 | 0.551 | 0.038 | -4.9e-6 | | 113.2 | 47.2 | 0.58 |
| 64 (8) | 0.5 | RDF | completed | 960 | 2.98e-3 | 1 | 0.597 | 0.059 | -2.1e-4 | | 82.2 | 85.6 | 1 |
| | | RDFccs | completed | 960 | 3.10e-3 | 1.04 | 0.609 | 0.060 | -2.7e-4 | | 79.7 | 83.0 | 0.97 |
| | | RING | ABORT step 495 (t 1.55) | 495 | - | - | - | - | -0.23 | 1694 / 169410 | 82.2 | 166 | 1.9 |
| | | RINGvcc | ABORT step 493 | 493 | - | - | - | - | -0.22 | 1694 / 170555 | 74.4 | 151 | 1.8 |
| | | ISO | completed | 960 | 3.32e-3 | 1.11 | 0.546 | 0.036 | -5.3e-4 | | 48.6 | 50.6 | 0.59 |
| 64 (8) | 1 | RDF | completed | 480 | 3.13e-3 | 1 | 0.633 | 0.061 | -8.3e-3 | | 46.5 | 96.9 | 1 |
| | | RDFccs | completed | 480 | 3.16e-3 | 1.01 | 0.629 | 0.062 | -7.1e-3 | | 46.9 | 97.6 | 1.00 |
| | | RING | ABORT step 264 (t 1.65) | 264 | - | - | - | - | -0.46 | 1694 / 169703 | 45.5 | 172 | 1.8 |
| | | RINGvcc | ABORT step 295 | 295 | - | - | - | - | -0.34 | 1694 / 169622 | 42.8 | 145 | 1.5 |
| | | ISO | completed | 480 | 3.36e-3 | 1.08 | 0.624 | 0.062 | -0.05 | | 27.0 | 56.1 | 0.58 |
| 128 (64) | 0.2 | RDF | completed | 4800 | 6.89e-4 | 1 | 0.196 | 7.9e-3 | -2.4e-5 | | 861.0 | 179 | 1 |
| | | RDFccs | completed | 4800 | 6.85e-4 | 0.99 | 0.179 | 7.1e-3 | -2.5e-5 | | 859.8 | 179 | 1.00 |
| | | RING | ABORT step 646 (t 0.40) | 646 | - | - | - | - | -0.15 | 6836 / 685612 | 90.0 | 139 | 0.78 |
| | | RINGvcc | ABORT step 862 | 862 | - | - | - | - | -0.099 | 6836 / 684118 | 132.0 | 153 | 0.85 |
| | | ISO | completed | 4800 | 5.39e-4 | 0.78 | 0.131 | 6.4e-3 | -2.7e-5 | | 406.9 | 84.8 | 0.47 |
| 128 (64) | 0.5 | RDF | completed | 1920 | 6.90e-4 | 1 | 0.190 | 7.7e-3 | -2.3e-4 | | 357.6 | 186 | 1 |
| | | RDFccs | completed | 1920 | 6.87e-4 | 1.00 | 0.187 | 7.5e-3 | -2.2e-4 | | 355.4 | 185 | 0.99 |
| | | RING | ABORT step 234 (t 0.37) | 234 | - | - | - | - | -0.37 | 6836 / 684394 | 50.0 | 214 | 1.1 |
| | | RINGvcc | ABORT step 441 | 441 | - | - | - | - | -0.28 | 6836 / 684785 | 89.3 | 203 | 1.1 |
| | | ISO | completed | 1920 | 5.43e-4 | 0.79 | 0.146 | 7.1e-3 | -5.2e-4 | | 172.4 | 89.8 | 0.48 |
| 128 (64) | 1 | RDF | completed | 960 | 6.88e-4 | 1 | 0.196 | 7.6e-3 | -5.4e-3 | | 207.4 | 216 | 1 |
| | | RDFccs | completed | 960 | 6.92e-4 | 1.01 | 0.207 | 8.1e-3 | -6.3e-3 | | 211.5 | 220 | 1.02 |
| | | RING | ABORT step 131 (t 0.41, min alpha) | 131 | - | - | - | - | -0.53 | 6836 / 390192 | 31.2 | 238 | 1.1 |
| | | RINGvcc | ABORT step 219 | 219 | - | - | - | - | -0.50 | 6836 / 480647 | 46.5 | 213 | 0.98 |
| | | ISO | completed | 960 | 5.50e-4 | 0.80 | 0.173 | 8.4e-3 | -2.9e-3 | | 92.2 | 96.0 | 0.44 |

Orders of E_shape (32-64 / 64-128) at OF Co 0.2: RDF 1.37 / 2.16, RDFccs 1.38 / 2.15, ISO
1.47 / 2.60; the same within 0.05 at Co 0.5 and 1. The profiled RING runs at N 32: the
reconstruction 59 ms of a 144 ms step (41 %; plicRDF 26 of 51 ms) because the state
fragments into 24,000-27,000 interface cells by |m_c| against 410 at the start (the
transported face fractions, the wisps and the second-layer bounding: the Enright agent's
finding of 2026-11-06 on the swept loop, the same here), the flux sweep 34 ms against 2.8,
the second ring 15 ms.


### The 3D shear flow (Liovic et al. 2006), hexahedral (T = 3; OF Co 0.2 / 0.5 / 1)

The sphere of radius 0.15 at (0.5, 0.75, 0.25) in the 1 x 1 x 2 box, the field of
`initShearFlow3D` (phi = U(x_f) . S_f at t = 0; the first batch of 2026-11-06 ran with the
app's face loop indexing the neighbour list beyond its size, four corrupted fluxes at N 128
and a segfault at N 32/64, fixed in fe5ff45 and re-prepared) scaled by cos(pi t / 3), scored
at t = 3 against the initial sphere; the timed runs of 12:14 were lost to a SLURM step
refusal (the submitter job's TRES inherited by the timed jobs; the submitter now clears
its SLURM variables) and re-run at 12:21 from the login node.

| N (ranks) | OF Co | row | status | steps | E_shape | E/RDF | L1(n) | L1(x) | min alpha | mass dev. % | interface cells (first / max) | wall s | ms/step | step/RDF |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 32 (1) | 0.2 | RDF | completed | 720 | 4.71e-03 | 1.00 | 0.293 | 0.0177 | -1.8e-06 | 1.7e-08 |  | 28.03 | 38.9 | 1.00 |
|  |  | RDFccs | completed | 720 | 4.71e-03 | 1.00 | 0.293 | 0.0177 | -1.8e-06 | 1.5e-08 |  | 27.32 | 37.9 | 0.97 |
|  |  | RING | ABORT step 382 | 382 | - | - | - | - | -0.081 | 0.18 | 416 / 41614 | 65.24 | 171 | 4.39 |
|  |  | RINGvcc | ABORT step 384 | 384 | - | - | - | - | -0.073 | 0.36 | 416 / 41807 | 57.41 | 150 | 3.84 |
|  |  | ISO | completed | 720 | 4.67e-03 | 0.99 | 0.29 | 0.0172 | -1.7e-06 | 4.7e-09 |  | 16.93 | 23.5 | 0.60 |
| 32 (1) | 0.5 | RDF | completed, min -1.4e-04 | 288 | 4.60e-03 | 1.00 | 0.291 | 0.0174 | -1.4e-04 | 3.6e-08 |  | 14.42 | 50.1 | 1.00 |
|  |  | RDFccs | completed, min -1.4e-04 | 288 | 4.60e-03 | 1.00 | 0.29 | 0.0174 | -1.4e-04 | 5.4e-08 |  | 14.18 | 49.2 | 0.98 |
|  |  | RING | ABORT step 169 | 169 | - | - | - | - | -0.18 | 0.18 | 416 / 41692 | 30.78 | 182 | 3.64 |
|  |  | RINGvcc | ABORT step 211 | 211 | - | - | - | - | -0.14 | 0.54 | 416 / 41603 | 32.89 | 156 | 3.11 |
|  |  | ISO | completed, min -1.3e-04 | 288 | 4.58e-03 | 0.99 | 0.285 | 0.0171 | -1.3e-04 | 1.3e-07 |  | 9.57 | 33.2 | 0.66 |
| 32 (1) | 1 | RDF | completed, min -3.1e-03 | 144 | 4.45e-03 | 1.00 | 0.303 | 0.0217 | -3.1e-03 | 8.9e-05 |  | 9.91 | 68.8 | 1.00 |
|  |  | RDFccs | completed, min -5.9e-03 | 144 | 4.46e-03 | 1.00 | 0.302 | 0.0205 | -5.9e-03 | 8.9e-05 |  | 9.76 | 67.8 | 0.98 |
|  |  | RING | ABORT step 111 | 111 | - | - | - | - | -0.38 | 0.26 | 416 / 41755 | 21.82 | 197 | 2.86 |
|  |  | RINGvcc | completed, unbounded | 144 | 0.0187 | 4.20 | 1.24 | 0.436 | -0.38 | 13 | 416 / 32849 | 20.83 | 145 | 2.10 |
|  |  | ISO | completed, min -1.3e-03 | 144 | 4.43e-03 | 0.99 | 0.314 | 0.0273 | -1.3e-03 | 1.5e-06 |  | 7.03 | 48.8 | 0.71 |
| 64 (8) | 0.2 | RDF | completed, min -3.0e-05 | 1440 | 1.34e-03 | 1.00 | 0.141 | 4.82e-03 | -3.0e-05 | 3.7e-10 |  | 128 | 89 | 1.00 |
|  |  | RDFccs | completed, min -3.0e-05 | 1440 | 1.34e-03 | 1.00 | 0.141 | 4.82e-03 | -3.0e-05 | 5.4e-10 |  | 122 | 84.7 | 0.95 |
|  |  | RING | ABORT step 265 | 265 | - | - | - | - | -0.097 | 4.2e-03 | 1760 / 177187 | 30.45 | 115 | 1.29 |
|  |  | RINGvcc | ABORT step 367 | 367 | - | - | - | - | -0.093 | 0.15 | 1760 / 176601 | 49.59 | 135 | 1.52 |
|  |  | ISO | completed, min -1.0e-04 | 1440 | 1.42e-03 | 1.06 | 0.144 | 5.11e-03 | -1.0e-04 | 1.1e-09 |  | 54.94 | 38.2 | 0.43 |
| 64 (8) | 0.5 | RDF | completed, min -3.4e-04 | 576 | 1.31e-03 | 1.00 | 0.138 | 4.80e-03 | -3.4e-04 | 1.2e-07 |  | 58.5 | 102 | 1.00 |
|  |  | RDFccs | completed, min -3.4e-04 | 576 | 1.31e-03 | 1.00 | 0.138 | 4.79e-03 | -3.4e-04 | 1.2e-07 |  | 56.65 | 98.4 | 0.97 |
|  |  | RING | ABORT step 130 | 130 | - | - | - | - | -0.23 | 8.8e-11 | 1760 / 178058 | 23.51 | 181 | 1.78 |
|  |  | RINGvcc | ABORT step 216 | 216 | - | - | - | - | -0.23 | 0.03 | 1760 / 176442 | 35.09 | 162 | 1.60 |
|  |  | ISO | completed, min -2.6e-04 | 576 | 1.38e-03 | 1.05 | 0.145 | 5.05e-03 | -2.6e-04 | 7.0e-08 |  | 26.1 | 45.3 | 0.45 |
| 64 (8) | 1 | RDF | completed, unbounded | 288 | 1.27e-03 | 1.00 | 0.265 | 0.0149 | -0.011 | 1.0e-05 |  | 38.26 | 133 | 1.00 |
|  |  | RDFccs | completed, min -9.2e-03 | 288 | 1.27e-03 | 1.00 | 0.226 | 0.0139 | -9.2e-03 | 1.6e-05 |  | 37.45 | 130 | 0.98 |
|  |  | RING | ABORT step 106 | 106 | - | - | - | - | -0.42 | 1.6e-12 | 1760 / 176436 | 25.57 | 241 | 1.82 |
|  |  | RINGvcc | ABORT step 30 (max alpha) | 30 | - | - | - | - | -0.23 | 5.3e-13 | 1760 / 9216 | 5.59 | 186 | 1.40 |
|  |  | ISO | completed, unbounded | 288 | 1.31e-03 | 1.03 | 0.228 | 0.0118 | -0.02 | 1.2e-05 |  | 16.54 | 57.4 | 0.43 |
| 128 (64) | 0.2 | RDF | completed, min -4.6e-05 | 2880 | 3.31e-04 | 1.00 | 0.053 | 1.17e-03 | -4.6e-05 | 1.3e-09 |  | 462 | 160 | 1.00 |
|  |  | RDFccs | completed, min -4.7e-05 | 2880 | 3.31e-04 | 1.00 | 0.053 | 1.17e-03 | -4.7e-05 | 1.3e-09 |  | 452 | 157 | 0.98 |
|  |  | RING | ABORT step 328 | 328 | - | - | - | - | -0.11 | 1.2e-12 | 6680 / 675893 | 69.55 | 212 | 1.32 |
|  |  | RINGvcc | ABORT step 726 | 726 | - | - | - | - | -0.11 | 8.4e-07 | 6680 / 670097 | 117 | 162 | 1.01 |
|  |  | ISO | completed, min -1.9e-05 | 2880 | 3.78e-04 | 1.14 | 0.0608 | 1.33e-03 | -1.9e-05 | 1.3e-10 |  | 215 | 74.8 | 0.47 |
| 128 (64) | 0.5 | RDF | completed, min -9.5e-04 | 1152 | 3.29e-04 | 1.00 | 0.0564 | 1.23e-03 | -9.5e-04 | 1.8e-07 |  | 214 | 185 | 1.00 |
|  |  | RDFccs | completed, min -1.1e-03 | 1152 | 3.29e-04 | 1.00 | 0.0574 | 1.24e-03 | -1.1e-03 | 1.9e-07 |  | 212 | 184 | 0.99 |
|  |  | RING | ABORT step 181 | 181 | - | - | - | - | -0.37 | 1.2e-12 | 6680 / 669254 | 55.26 | 305 | 1.65 |
|  |  | RINGvcc | ABORT step 297 | 297 | - | - | - | - | -0.27 | 1.7e-12 | 6680 / 668695 | 73.72 | 248 | 1.34 |
|  |  | ISO | completed, min -2.7e-04 | 1152 | 3.67e-04 | 1.11 | 0.0617 | 1.33e-03 | -2.7e-04 | 8.3e-08 |  | 96.57 | 83.8 | 0.45 |
| 128 (64) | 1 | RDF | completed, unbounded | 576 | 3.38e-04 | 1.00 | 0.292 | 0.0116 | -0.015 | 3.6e-06 |  | 153 | 265 | 1.00 |
|  |  | RDFccs | completed, min -9.6e-03 | 576 | 3.38e-04 | 1.00 | 0.289 | 0.0114 | -9.6e-03 | 3.4e-06 |  | 159 | 276 | 1.04 |
|  |  | RING | ABORT step 24 (max alpha) | 24 | - | - | - | - | -0.32 | 3.4e-13 | 6680 / 20364 | 9.92 | 413 | 1.56 |
|  |  | RINGvcc | ABORT step 39 (max alpha) | 39 | - | - | - | - | -0.46 | 2.2e-13 | 6680 / 34119 | 12.12 | 311 | 1.17 |
|  |  | ISO | completed, unbounded | 576 | 3.65e-04 | 1.08 | 0.185 | 6.74e-03 | -0.015 | 1.2e-06 |  | 55.84 | 96.9 | 0.37 |

Three repetitions on every completed cell (the spread of the wall clock 0.1-3.5 %); the
aborted DAVOF cells ran once (Tomislav's rule of 2026-11-06), the profiled run repeats their
abort step exactly. Orders of E_shape over N (32-64 / 64-128): RDF 1.81 / 2.02 (Co 0.2),
1.81 / 2.00 (Co 0.5), 1.81 / 1.91 (Co 1); ISO 1.72 / 1.91; of L1(n): RDF 1.05 / 1.42 at Co
0.2, 1.08 / 1.29 at Co 0.5, 0.19 / -0.14 at Co 1. At OF Co 1 plicRDF's normal error does
not converge on this flow (L1(n) 0.30 / 0.27 / 0.29 at N 32 / 64 / 128 against 0.29 / 0.14 /
0.053 at Co 0.2) and plicRDF, plicRDF-ccs and isoAlpha leave [0, 1] by 0.011-0.015 at N 64
and 128 (RDFccs 9.2e-3 / 9.6e-3), while E_shape still converges at 1.9: the shape error is
the flux's, the normal's error at Co 1 is a wider band of surface cells.

The DAVOF rows abort at every cell by the interface-cell guard (100 x the initial count: 416
to 41,600 at N 32, 1760 to 177,000 at N 64, 6680 to 676,000 at N 128) between t = 0.34
(N 128, Co 0.2) and t = 2.31 (N 32, Co 1), or, at Co 1 on N >= 64, by the max-alpha guard
within 24-39 steps; the one exception is RINGvcc at N 32 Co 1, which completes unbounded
(min alpha -0.38, E_shape 4.2 x RDF, the mass deviation 13 %). The fragmentation is the same
mechanism as on the deformation (the transported face fractions, the wisps and the second
bounding layer under straining; the Enright agent's finding on the swept loop) and sets in
earlier the finer the mesh. The profiled runs: at N 32 Co 0.2 the RING step is 160 ms
against 34.9 (4.6 x) with the reconstruction at 68.7 ms over 41,000 interface cells (43 %),
the flux sweep 41.7 against 3.1 ms and the second ring 19 ms; at N 128 Co 0.2 on 64 ranks
189 against 156 ms (1.2 x; the reconstruction 77.6 ms, 41 %) up to the abort at step 328.

### The verdicts on the hexahedral rotation, shear and deformation

- F-1 (completion and bounds): FALSIFIED for the DAVOF rows beyond OF Co 0.2 on the
  rotation and at N 128 even there: RING completes bounded (min alpha -6.9e-13 / -3.3e-6)
  at N 32 and 64 at OF Co 0.2 only; at Co 0.5 it aborts at every N (the interface-cell
  guard, 4270 / 18159 / 67945 cells of 416 / 1760 / 6680), at Co 1 as well, and at N 128
  at Co 0.2 at step 670 of 6911. RINGvcc completes two more cells (N 32 Co 0.5 with
  min alpha -0.07 and L1(n) 0.36, N 128 Co 0.5 with -1.2e-3). On the deformation the
  DAVOF rows complete at N 32 unbounded (-0.12 to -0.47) and abort at N 64 (t 1.5-1.65)
  and N 128 (t 0.37-0.41), as registered. On the shear they abort at every cell (the
  registration expected completion at N 32 and at N 64 for Co 0.2: FALSIFIED in the
  pessimistic direction), RINGvcc completing N 32 Co 1 unbounded. RDF, RDFccs and ISO
  complete every cell; plicRDF leaves [0, 1] by 1.5e-4 at N 128 Co 1 on the rotation, by
  8.3e-3 / 5.4e-3 at Co 1 on the deformation (ISO 0.05 at N 64 Co 1) and by 0.011 / 0.015
  at N 64 / 128 Co 1 on the shear (ISO 0.02 / 0.015).
- F-2 (errors where RING completes): HELD. E_shape(RING)/E_shape(RDF) 0.94 (N 32) and 0.96
  (N 64), L1(n) 1.03 / 1.06, L1(x) 0.94 / 0.95; the order 32-64 of E_shape 1.85 (RDF 1.89),
  of L1(n) 1.21 (RDF 1.26).
- F-3 (ccs against vofCutCell): HELD for the errors and the cost. plicRDF: E_shape, L1(x)
  within 0.1 % of each other at every rotation cell, L1(n) within 1 %; on the deformation
  E_shape within 1 % (4 % at N 64 Co 0.5), L1(n) within 3 % (8 % at N 128 Co 0.2); the wall
  clock within 1.5 % (RDFccs 0.98-1.02 x RDF). The DAVOF rows: RING and RINGvcc equal to
  three digits where both complete (N 32 / 64 Co 0.2), ccs 6-8 % faster. The ccs
  positioning: 5.6-6.0 cuts per cell, the residual max 1e-14 .. 3e-14 of V. The residual
  of vofCutCell's own planes, measured with `planeVolumeConsistency` on the written end
  states of the DAVOF rows (the only rows that write `alphaf`, which the app needs for its
  interface-cell set; job 56011148, `pvc_{vofCutCell,ccs}.dat` under each case's
  `postProcessing/planeVolumeConsistency/`), with the ring normal the run used (`field`)
  and the area normal, at N 32 (435 cells) and N 64 (1758 cells) on the rotation and on the
  N 32 deformation (21,880 cells by |m_c|, 5,800-6,200 inside [tol, 1 - tol]): vofCutCell
  at its exit tolerance 1e-8 leaves |V_cut - alpha V|/V up to 9.0-10.0e-9 (mean 7-11e-10,
  two thirds of the cells above 1e-14, signed sum 1e-12 (rotation) to 1.9e-10 (deformation)
  against N eps V 1.3e-17), 5-6 cuts per cell, 3.7-4.2 us per positioning; ccs on the same
  cells and normals 6e-15 .. 1.6e-14 (mean 1e-15, 0-5 cells above 1e-14 at N 64, signed sum
  1e-18 .. 1e-17), 2.85-2.95 cuts per cell (2.0 / 4.8 on the deformation), 1.8-2.1 us per
  positioning. So plicRDF's planes (vofCutCell, the same iteration) carry a volume residual
  of the order of its tolerance 1e-8, six decades above ccs's, at twice ccs's cost per
  positioning, and the transport errors of the two plicRDF rows differ by less than 1 % all
  the same: a plane residual of 1e-8 is below what the flux sees.
- F-4 (wall clock): FALSIFIED in the registered direction, the measurement below the
  registered range: RING / RDF per step 0.84 (N 32, serial) and 0.67 (N 64, 8 ranks)
  against the registered [1.0, 2.5]; the reconstruction's share 24 / 32 % (RING) against
  52 / 66 % (RDF), the DAVOF reconstruction 2.6 x (N 32) and 3.2 x (N 64) cheaper than
  plicRDF's (the registered factor 3 held at N 64 and N 128, missed at N 32); the face
  update's extra cost 3.6-6.8 ms per step as itemised above. On the deformation RING / RDF
  2.3-3.2 x (N 32) with 24,000 interface cells, 1.8-2.1 x (N 64) up to the abort; on the
  shear 2.9-4.4 x (N 32), 1.3-1.8 x (N 64), 1.0-1.7 x (N 128), every DAVOF cell up to its
  abort with 40,000-680,000 interface cells.
- F-5 (parallel): HELD on the completed cells: the per-step ratio RING / RDF falls from
  0.84 (serial) to 0.67 (8 ranks); at 64 ranks 0.69-0.76 up to the abort.

Also measured: the mass deviation of the DAVOF rows over the run 1.7e-4 .. 4e-4 % (RDF,
RDFccs, ISO 1e-8 .. 1e-7 %), present on every completed DAVOF cell: a conservation defect
of the plane-flux rule or of the second bounding layer to be pinned (not registered, not
tuned). The wall-clock spread (max - min over the median of three repetitions) is 0.1-3.5 %
on 55 of the 63 completed hexahedral cells; eight cells carry one slow repetition (up to
14 % on the deformation RDF N 64 Co 1, 41 % on the aborted RING N 128 Co 0.5: 120 s against
85 and 86 s), and the tables report the median.

### The liquid jet in cross flow (`run/rhoVoF/jetCrossFlow`, the coarse mesh, 32 ranks)

(a) The meshes of the data archive (doi:10.48328/tudatalib-1277; cfMesh): coarse 524,288 and
fine 4,194,304 cells, 100 % hexahedra, `checkMesh` OK, max non-orthogonality 4.77 degrees,
skewness 0.065, aspect ratio 1.26; faces warped up to 8.7e-3 h (2.0 % of the coarse faces
above 1e-6, 0.54 % of the fine), all within 0.45-6.5 cells of the box boundary, none near the
jet; the cells are not exactly axis-aligned (`exactSphere` refuses them; the static test
used `quadraticFaces`). The static normal test (the sphere of radius 4 mm at (12, 10, 0) mm):
the ring normal L1(n) 1.2027e-2 (coarse) and 5.977e-3 (fine), identical to a blockMesh
control at the same h; plicRDF 9.012e-3 / 4.328e-3 (0.75 / 0.72 of the ring's); the ccs and
vofCutCell positionings agree to 5-6 digits. The ring normal is undegraded by the cfMesh
hexahedra.

(b) The timing window: 2000 steps of dt 5e-7 s (t = 1 ms) on 32 ranks (simple 8 2 2), the
published settings (`interFlow`, the coded `jetInlet`, `nAlphaBounds 3`, the density-ratio
branch's p-U coupling) except the `visRDF` function object (a compile at every start) and the
interface surface written as `reconstructedInterface` of `libpostProcess` (the published
`interface` key loads OpenFOAM's `geometricVoF` next to `libVoF` and crashed at the first
write). plicRDF (vofCutCell): 2045.5 s, 1.02 s per step, min alpha -4.1e-6; plicRDF-ccs:
2078.2 s, 1.04 s per step, -1.9e-6 (the first repetition each; two more and the 8 ms runs are
in the queue). The DAVOF row (ring + ccs, the binary b09b800; the initial `alphaf` is the
initial alpha's state: 0 on every internal face, 1 on the nozzle faces of the `jetInlet`
patch, whose alpha condition is 1) ran 623 steps (t = 0.31 ms, 650 s, 1.04 s per step) and
died with a floating-point exception in the pressure solve (rc 143): the plane state
drifted from the step the liquid entered, the explicit plane defect max |alpha_plane -
alpha| growing from 0 to 0.03 (step 50), 0.16 (200), 0.45 (450) and 0.77 (600), the
interface cells from 32 to 150 and then, within the last 20 steps, to 480,000 with alpha at
1e86; min alpha stayed above -4.6e-4 until step 500. Reported as measured; the row is
re-run once with ae08dc0 (the inflow value in the defaults alone, identical to the swept
branch on the eight cases), and its full 8 ms run with it (jobs 56013018 / 56013019). The windward-trajectory extraction against
Desclaux et al. (2020) at t = 7.1 ms follows the 8 ms runs of the plicRDF rows
(`run/rhoVoF/jetCrossFlow/comparison/`).

### The standing study set (Tomislav, 2026-11-07) and the stamps

The three families above (hex) are the first mesh class of the standing set: 2D shear
(the sheared disc of Rider & Kothe 1998), 3D shear, 3D deformation; three resolutions;
OpenFOAM Co 0.2, 0.5 and 1; the rows plicRDF (vofCutCell), plicRDF-ccs, DAVOF ring + ccs
(+ vofCutCell), isoAlpha; perturbed and polyhedral (and tet) classes next. Added for it on
2026-11-07: the fork's `davofInit::exactCylinder` (the exact DAVOF state of a cylinder on
axis-aligned hexahedra: the 2D disc on a one-cell-thick mesh, selectable in
`initDavofFields`, `areaFractionError` and `volumeFractionError`) and `initVortex2D` (U and
phi of the single vortex, phi the exact line integral of the stream function around every
face; b467390); leia's `cases/davof/isoShear2D` (the isoShear3D template with the N x 1 x N
mesh of cubic cells, front and back empty, R 0.15 at (0.5, 0.75), T = 8), the configs
`isoShear2D_hexN{64,128,256}` (1 / 1 / 8 ranks) with the pre-registration S-1..S-5
(4dfb4f3c, written before any run), Co 1 on the mesh-class configs (0c0ddb84, 0228ea15),
the runner's FLOW Shear2D and the submitter. The first preparation of `isoShear2D_hexN64`
(15 cases with the fork at b467390 in `iso/TwoPhaseFlow-c`): the mesh 64 x 1 x 64, the
exact disc by `exactCylinder` (the exact state of `volumeFractionError` taken
automatically), OpenFOAM Co 0.2163 at the target 0.2 with the provisional reference speed
1.2, so REFERENCE_SPEED is 1.3 (3299c7cf); the DAVOF rows' explicit plane position
(the Gauss sum over the cell's faces) missed the two empty faces of every cell, a start
defect of 0.48 in alpha_plane for the ccs iteration (which converges all the same), and
the sum now takes alpha_c on an empty face (the cross-section is the cell's; a closed cell
is unchanged bit for bit) before the three 2D studies are prepared. The 2D and 3D shear templates
(`cases/davof/isoShear2D`, `cases/davof/isoShear3D`, their `.parameter` files) are the ones
to share with the port track: one template per case, the rows in `tpf_rows.py`. The stamps
of every timed cell (the fork commit and the sha256 of its `advectorVoF` and `libVoF.so`,
the leia commit that last touched the case template, the study configuration and the runner,
the sha256 over the rendered dictionaries of the prepared cases) are in
`iso/csv/iso_stamps.csv`; the recompute-on-change rule itself is the port agent's stamp
mechanism in leia's workflow, which this track will use when it lands.

## The iterative face-fraction corrector (2026-11-09, leia branch `davof-iter`)

Tomislav's instruction (2026-11-09): "Start it." The DAVOF state is (alpha_k, alpha_f).
The transport step sets alpha_k^{n+1} and alpha_f^{n+1} (the swept pieces' clips of the
faces, the predictor); the reconstruction then forms m_k = -sum_f alpha_f S_f, the ring
area normal d_k = sum_{j in R(k)} m_j and positions the plane by CCS from alpha_k. The
predictor's face fractions and the planes they produce are not consistent with each
other: the plane of a cell clips its faces to fractions that differ from the transported
ones. The corrector iterates on this inconsistency inside the reconstruction, the
advection untouched:

1. reconstruct from the current alpha_f (ring area normal, CCS from alpha_k);
2. the face fractions the planes induce: a face between two interface cells takes
   1/2 (alpha_f[Pi_owner] + alpha_f[Pi_neighbour]), alpha_f[Pi] the face polygon
   clipped by the plane Pi (wet area over face area); a face shared with a cell that has
   no plane takes that cell's alpha (0 or 1 within the tolerance) in place of its clip;
   a face touching no interface cell and every non-coupled boundary face are unchanged;
   a processor face receives the identical value on both sides (the half-sum);
3. alpha_f := the induced values; back to 1, up to `correctorIterations` times or until
   the largest change of alpha_f falls below `correctorTolerance` (then no replacement).

The interface cells of the iteration are the cells with eps < alpha_k < 1 - eps
(eps = surfCellTol 1e-8; `interfaceCellRule alphaOnly`, forced while the corrector is
on), not the cells flagged by their face fractions. The final alpha_f is the state
carried to the next step and written; the final planes are the reconstruction the next
step's flux uses (the first reconstruct() of the next time index returns without
recomputing). Keys of the alpha dictionary: `faceFractionCorrector none | mean`
(default none: bit for bit the behaviour before), `correctorIterations` (1),
`correctorTolerance` (1e-12), `correctorLog` (false; the per-step line: iterations
used, max and L1 change of alpha_f per iteration, interface cells, faces, cells without
a plane, and the identity |sum_k m_k - sum_b alpha_f S_f| / sum_k |m_k| over the
non-coupled boundary faces, which holds to round-off exactly when every face carries one
value). Code: `src/davof/reconstruction/davof.{H,C}` (`correctFaceFractions`,
`inducedFaceFractions`, the rule `alphaOnly`, the skip of the next step's first
reconstruction), the hook at the end of `davofSweptAdvection::advect` (the flux
untouched), and `isoAdvectorScheme` (the DAVOF reconstruction under OpenFOAM's
isoAdvection: the predictor alpha_f from the UPWIND cell's plane translated by U_k dt,
the fork's plane-clipped face update in the wrapper, then the corrector). The study
tokens: `TPF_FACE_FRACTION_CORRECTOR none | mean | mean<k>` renders the keys through
`TPF_CORRECTOR_KEYS` (empty for none: every dictionary rendered before stays byte for
byte, the stamps recompute nothing); `TPF_ADVECTION_SCHEME auto | isoAdvector`.

Prior art (cited in the code): B. Swartz, The second-order sharpening of blurred smooth
borders, Math. Comp. 52 (1989) 675-714, doi:10.1090/S0025-5718-1989-0983313-8, and
S. J. Mosso, B. K. Swartz, D. B. Kothe, R. C. Ferrell, A parallel, volume-tracking
algorithm for unstructured meshes, Parallel CFD '96 (Capri, 1996), Elsevier 1997, pp.
368-375, LA-UR-96-2420 (no DOI; the author list, title and conference verified from the
reference list of Maric, Marschall, Bothe, arXiv:1305.3417, 2026-11-09; the proceedings'
pages from the standard citation, not re-verified) iterate the NORMAL of a cell from
the positions of the neighbours' planes (the line through the neighbouring interface
points). Ours iterates the FACE FRACTIONS the neighbouring planes induce; the normal is
never touched directly and follows from alpha_f by the area identity and the ring sum,
so that every iterate is a DAVOF state (alpha_k, alpha_f).

### Pre-registered (2026-11-09, before the first run)

The cases: the plane translation `planeTranslation3D` (both arms: normal U = n and
oblique U = (1 0 0); Co 0.2 / 0.5 / 1; N 20 / 40 / 80), the sphere rotation
`sphereRotation3D` (one turn; OF Co 0.2 = interface Co 0.0463, OF Co 1 = 0.2315; N 32 /
64 / 128 on 1 / 8 / 64 ranks), the 3D deformation `isoDeformation3D` (the standing
matrix's template, OF Co 0.5, N 32 / 64 on 1 / 8 ranks, `grazingCut fork` so that the
none row reproduces the standing record), every DAVOF row the port's reference setting
(`davof-planeSurf-ringArea`: the swept flux with the round-4 keys, the ring area normal,
CCS); the corrector none | mean with 1, 2 and 4 iterations (tolerance 1e-12); the
rotation at N 32 / 64, OF Co 0.2, also under `isoAdvector` with the DAVOF reconstruction.
Configs `config/davof/iter_*.yaml` (`workflow/scripts/iter_make_configs.py`), tables
`workflow/scripts/iter_tables.py`, cluster directory
`/work/scratch/tm83tomy/davof/agents/iter`.

- **P1, the plane gate (MUST hold).** On the plane translation, both arms, every Co and
  N, the corrector changes nothing: 0 or 1 iteration per step with max |dalpha_f| <=
  1e-12 at the first evaluation of every step; the run maxima of L1(alphaf) and L1(n)
  <= 1e-12 as without the corrector; alpha, alphaf at the end of the np 8 run equal the
  np 1 run's to round-off (max |difference| <= 1e-12, 20 steps, plane N 40 oblique Co
  0.5 and rotation N 32 OF Co 0.2). FALSIFIED by a step with max |dalpha_f| > 1e-12 on
  the plane (the planes of a planar state must reproduce the transported fractions;
  a larger change is a defect of the clip, the sync or the alphaOnly set), by a gate
  value above 1e-12, or by a parallel difference above 1e-12.
- **P2, convergence on the sphere.** Rotation N 32 / 64 at OF Co 0.2 and 1, deformation
  N 32 / 64 at OF Co 0.5, 4 iterations: the max change of alpha_f falls by at least a
  factor 2 from one iteration to the next over the first 4 iterations at every step
  (the median and the worst factor per iteration are reported), and the iterations to
  reach max |dalpha_f| 1e-6 and 1e-10 are reported (">4" when not reached). FALSIFIED
  by a step whose max change grows from one iteration to the next (factor < 1) or whose
  factors stay below 2 at every iteration after the first; a factor below 2 at a few
  steps is reported with the located cells, not a falsification.
- **P3, accuracy.** Rotation N 32 / 64 / 128 at OF Co 0.2 and deformation N 32 / 64 at
  OF Co 0.5, mean with 1, 2 and 4 iterations against none: E_shape and L1(n) at the end
  not worse than none by more than 5 % with 1 iteration; L1(n) improves with the
  iterations (the ring normal is first order, L1(n) about 1e-2 at N 128 on the sphere
  as the ladders of 2026-11-05 measured; the expected size of the improvement is a few
  per cent to a factor 2, not an order change); boundedness (min alpha >= -1e-9, excess
  <= 1e-9 on the rotation) and the mass deviation as without the corrector (round-off);
  the interface-cell count of every corrected run equals the alpha-based count (no
  face-fraction-flagged cells: the log's count equals the number of cells with
  1e-8 < alpha < 1 - 1e-8). FALSIFIED by E_shape or L1(n) worse than none by more than
  5 % at 1 iteration on any rung, by L1(n) growing with the iterations, by an
  unbounded or non-conservative run that none completes bounded, or by an
  interface-cell count differing from the alpha-based one.
- **P4, cost.** Wall clock per step at equal ranks (rotation N 64, OF Co 0.2, 8 ranks,
  one exclusive node, two repetitions, the diagnostics at the end only): one corrector
  iteration adds at most one reconstruction's cost to the step; the reconstruction's
  share of the step is measured first (OpenFOAM's profiling, `reconstruct` against
  `advect`; the fork measured the DAVOF reconstruction at 0.49-0.61 x plicRDF's at N
  128). FALSIFIED if (t_mean1 - t_none) per step exceeds 1.2 x the measured
  reconstruction time per step, or if t_mean4 - t_none exceeds 4.8 x of it.
- **P5, the first correction as a diagnostic.** On the rotation N 32 / 64 / 128 at OF
  Co 0.2 the max and the L1 change of alpha_f at iteration 1 (the median over the
  steps) scale as O(h): the least-squares order over N between 0.7 and 1.3 for the L1
  change (the advection's face-fraction inconsistency is first order); the order of
  the max change is reported. Under isoAdvector (N 32 / 64) the same quantity is
  reported next to the swept flux's. FALSIFIED by an order of the L1 change below 0.5
  or above 1.5.

Also reported, not predicted: the iterations used per step (median, max), the identity
defect, the cells without a plane, L1(alphaf) at the end and its order, L1(x), the
interface-cell counts, the total wall clock of every case (shared nodes, one run; the
exclusive-node timing of P4 is the cost measurement), the job ids and the stamps.

### MEASURED (2026-11-09 to 2026-11-11; three runs, the third the one reported)

Branch `davof-iter` (leia): bb14d0c5 the corrector and the pre-registration, 967b2657 the
unit app `davofTestCorrector`, 9c0d983d the unit-normal clip, 44da32c8 the parallel check's
CSV fallback, 914505cd the rule `correctorNoPlaneRule`. Cluster
`/work/scratch/tm83tomy/davof/agents/iter` (the clone `leia/`, the tables `csv/run3/`
with `iter_cases.csv`, `iter_orders.csv`, `timing_N64.csv`, `profiling_{none,mean1}.txt`,
`par_checks.txt`; the earlier runs kept in `csv/run2_full/` and `csv/run2_keep/`). Jobs:
run 1 build 56128315, orchestrators 56130420-31; run 2 build 56168740, orchestrators
56168741-54, unit 56168755, timing 56168756; run 3 build 56170303, orchestrators
56170304-17, unit 56170318, timing 56170319; the run-2 tables 56170339.

**Three runs, two located defects, one rule changed after the data.** Run 1 (2026-11-09,
bb14d0c5) failed the plane gate P1 with changes of alpha_f of exactly 0.5 at N 40 and 1
at N 80 and three divergences; the unit check on the exact planar state (N 40, 840
interface cells) reproduced it with the EXACT normal (pass 1: max change 0.5, L1 4.0e-3
over 3203 faces, i.e. about 25 faces; passes 2-3 round-off), identical to seven digits
with the ring normal: a defect of the clip, not of the normals or the positioning.
`cutFacePLIC::calcSubFace` snaps a point's signed value (x - base).n to zero below SMALL
(1e-15, absolute) and calls a face without a positive point a gas face; `cutCellPLIC`
passes the unit normal, so that snap acts at 1e-15 m, but the corrector passed the cut
polygon's AREA vector `normal_`, down to 1e-12 for a sliver cut, so every point of a
fully wet face of such a cell snapped and the face was induced as 0 (mean(1, 0) = 0.5
with its full neighbour); the number of sliver cuts grows with N (8.8e-4 at N 20, 0.5 at
N 40, 1 at N 80). Fixed in 9c0d983d (the clip takes the unit normal; the same defect was
in the port's `reinitFaceFractions` diagnostic of 2026-11-06 and is fixed on this branch).
Run 2 (2026-11-10, 9c0d983d) then showed changes of 1e-10 to 1.7e-5 per step at every N
and a divergence at N 80, Co 1 (normal arm, step 14-16), every worst face being a face
shared by an interface cell and a cell the plane GRAZES (alpha within surfCellTol 1e-8 of
1, e.g. 1 - 7.6e-9, no plane) whose shared face is still cut by the plane (the
transported value and the owner's clip agree at 0.99998958): the pre-registered rule
"that cell's alpha enters the mean in place of its clip" is wrong there by
O(alpha^(2/3)) ~ 1e-5. The plane that cuts the face is the only information about its
wet area; the neighbour's alpha says nothing about the face. Changed in 914505cd, AFTER
the data: `correctorNoPlaneRule clip` (the default: a face touching one plane takes that
plane's clip, a face touching no plane is unchanged), the pre-registered rule kept as
`alpha`. Run 3 (2026-11-11, 914505cd) is reported below; run 2's tables (the alpha rule)
are kept on the cluster and differ from run 3 on the sphere by less than 1 % in every
error (the rule touches grazed faces only).

**The default is unchanged (demonstrated).** The `none` rows of the deformation
(grazingCut fork, OF Co 0.5) reproduce the port's standing record bit for bit in every
run: N 32 every one of the 480 per-step bounding lines and E_shape 1.13689512056824021e-02
identical, N 64 the abort at step 628 with the identical min alpha -0.54457534646676198.

**The unit check (davofTestCorrector, run 3).** Exact planar state N 40: pass 1 changes
alpha_f by 4.0e-13 (max) and 2.9e-15 (L1), both arms, ring and exact normal: the planes
of a planar state reproduce the transported fractions to round-off. Static exact sphere
N 32, exact normal: pass 1 the realizability defect (max 0.0236, L1 6.6e-3), then
2.2e-4, 1.0e-5, 6.9e-7, 5.8e-8, 5.2e-9 (a contraction of 10-20 per pass: the exact model
takes the normal at the plane centroid of the first pass, which moves with m_c, so the
fixed point is approached, not hit); ring normal: 0.0270, 1.4e-3, 1.6e-4, 3.7e-5, 5.9e-6,
8.3e-7, 2.2e-7, 3.8e-8, 9.9e-9, 2.1e-9 (contraction 4-19 per pass). The identity
|sum_k m_k - sum_b alpha_f S_f| / sum_k |m_k| is 5e-15 to 5e-17 in every pass.

**P1, the plane gate: FAILED in its strict form, HELD in its gate values except one cell;
the mechanisms located.** Table "The plane gate" below (run 3; it = iterations used per
step median/max, d1 = the change of alpha_f at iteration 1, median over the steps; the
none rows are the gate as before, L1(alphaf) and L1(n) 1e-14 to 1.6e-11). (i) The
corrector does not change nothing: the first evaluation changes alpha_f by 8e-13 (N 20)
to 1.2e-10 (N 80) at Co 0.2 / 0.5 and 1.9e-9 (N 40) / 5e-8 (N 80) at Co 1, so 1 to 4
iterations are used and the strict clause (max change <= 1e-12 every step) is falsified
at N >= 40. Located: at N 80 the worst face is shared with an interface cell of alpha
1 - 6.7e-7 whose plane is a sliver (|normal_| 6e-13 against 3e-9 for a full cut); CCS
positions it to 1e-14 of V, which for a corner of 7e-7 V is a plane shift of 1e-10 h,
and the clip of the shared face inherits the 1e-10. (ii) The gate values hold: L1(alphaf)
and L1(n) at the end are at the none row's level (1e-14 to 3.3e-12) in 17 of the 18
mean4 cells, the np 8 run's per-step function-object numbers equal the np 1 run's to
3e-14 (L1), 4e-14 (L1(n)), 3e-24 (E_shape); Linf(alphaf) differs by up to 5e-12 in 7 of
20 steps. (iii) The exception: normal arm, Co 1, N 80 ends at L1(alphaf) 9.7e-3, L1(n)
0.133, min alpha -9.2e-4. Located: at step 13 a cell with alpha 1.49e-8, just inside the
tolerance band (surfCellTol 1e-8), is an interface cell under `alphaOnly` and gets a
sliver plane (|normal_| 4e-14) whose clip of its face (1.0e-5) replaces the transported
1.5e-8 next to a cell at alpha -9.2e-9; the none row itself holds cells at -1.2e-8 on this
cell (the swept flux's round-off wisps at Co 1). The alpha-based interface set admits
round-off wisps at the band's edge and gives them meaningless sliver planes: a
limitation of the alpha-only rule (the normal-based rule excludes them by |m_c|), not of
the clip. Parallel: the rotation pair np 1 / np 8 (20 steps, mean4) differs by 4.2e-10
per value in alpha (54 cells above 1e-12) and 3.4e-10 in alpha_f (178 faces), L1 2.8e-14;
the none pair differs by 1.1e-9 (181 cells) and 2.4e-9 (487 faces), L1 1.9e-13 / 4.1e-13:
the per-value scatter is the swept flux's own (2.5e-9 on the ring row, 2026-11-02), the
corrector adds none, and the pre-registered per-value 1e-12 was below the baseline's
scatter (the port's criterion was the L1, which holds).

**P2, convergence on the sphere: HELD on the rotation, FALSIFIED on the deformation; the
mechanism located.** Rotation (N 32 / 64 / 128 at OF Co 0.2, N 32 / 64 at Co 1, mean4):
the contraction factors of the max change, median / worst over the steps, are 18.9 /
4.3, 8.3 / 2.6, 6.0 / 2.8 (N 32, Co 0.2), 41 / 8.8, 5.9 / 1.8, 6.2 / 2.5 (N 64), 51 /
8.4, 4.8 / 1.6, 6.0 / 2.1 (N 128), 17 / 4.8, 8.7 / 2.6, 5.7 / 2.8 (N 32, Co 1), 44 / 11,
7.2 / 2.0, 6.3 / 2.5 (N 64, Co 1); under isoAdvector 22 / 4.2, 9.4 / 2.9, 5.9 / 2.7 (N
32) and 50 / 13, 7.7 / 2.0, 6.0 / 2.1 (N 64). No step grows; the worst factor falls
below 2 at the second iteration on a few steps at N 64 and 128 (1.6 to 2.0). The first
correction is 0.15 to 0.20 (max) and the fourth about 1e-4, so 1e-6 needs about 7 and
1e-10 about 12 iterations (not reached within 4; ">4" in the table). Deformation (N 32
serial, N 64 on 8 ranks, OF Co 0.5): the max change is 1 at EVERY iteration of every
step after the sheet forms (factor 1 / 1), while the L1 change contracts (step 300: L1
0.053, 0.043, 0.038, 0.036; step 1: 5.5e-3, 3.5e-4, 5.5e-5, 1.0e-5). Located: the faces
that change by 1 are (a) faces between two interface cells of the sub-cell sheet (alpha
0.09 and 0.08) whose planes both clip the shared face to 0 while the transported value
is 1, i.e. a sheet thinner than a cell that one plane per cell cannot carry across the
face, and (b) faces of wisp cells (alpha 1e-4 to 1e-8) whose ring normal is axis-aligned
so that the CCS plane lies ON a cell face (|normal_| = h^2 exactly) and the clip of that
face is 0 or 1 by the SMALL snap, flipping between iterations; next to them cells at
alpha -7e-6 (unbounded by the swept flux, no plane). These are the realizability limit of
one plane per cell on the sheet and the wisps of the fragmentation, not a defect of the
iteration.

**P3, accuracy: HELD on the rotation for the 5 % clause, FALSIFIED for "L1(n) improves
with iterations", FALSIFIED on the deformation (E_shape), with a positive finding on the
deformation's completion.** Rotation, OF Co 0.2, E_shape none / mean1 / mean4: N 32
1.059e-3 / 1.097e-3 (+3.6 %) / 1.102e-3 (+4.1 %); N 64 2.95e-4 / 3.02e-4 (+2.4 %) /
3.02e-4; N 128 8.05e-5 / 8.14e-5 (+1.1 %) / 8.14e-5. L1(n): N 32 0.0575 / 0.0596 (+3.7 %)
/ 0.0593; N 64 0.0252 / 0.0265 (+4.8 %) / 0.0265; N 128 0.0113 / 0.0116 (+2.5 %) /
0.0116. L1(x) +3 to +4 %. The corrector makes every error slightly WORSE, within the 5 %
clause at one iteration, and the iterations do not move it (mean4 = mean1 to three
digits): the fixed point of (reconstruct, induce) on the sphere is a state whose normals
are no better than the transported state's; the ring-area normal is first order and the
corrector cannot make it second order. The orders over N are unchanged: E_shape 1.88
(none 1.86), L1(n) 1.18 (1.17), L1(x) 1.93 (1.90), L1(alphaf) 0.84 (0.82); at Co 1 the
same picture (+5 % in E_shape, +3 % in L1(n)). Boundedness as without: min alpha -1e-12
to -1e-14, the excess over 1 up to 8.6e-10 at N 128 (none 6.5e-10); mass to 8e-12 %;
the interface-cell count of every corrected run is the alpha-based count (424 against
423 flagged by |m_c| at N 32, 6970 against 6870 at N 128). Under isoAdvector (the DAVOF
reconstruction, the predictor from the translated planes) the same: E_shape +3.7 %
(N 32), +2 % (N 64), L1(n) +2 to +4 %. Deformation N 32: E_shape none 1.137e-2, mean1
1.407e-2 (+24 %), mean2 1.486e-2, mean4 1.537e-2: falsified; L1(n) 1.46 / 1.30 / 1.29 /
1.40 and the run minimum of alpha -0.27 (none) against -0.05 / -0.06 / -0.06: the
corrector improves the bound and the normals on the sheet and worsens the shape. N 64
(8 ranks): none aborts at step 628 (the standing record), mean1 at step 889, mean2 and
mean4 COMPLETE T = 3 with E_shape 7.84e-3 and 7.72e-3 (orders 32 -> 64 of 0.92 and
0.99), L1(n) 1.30 and 1.31, min alpha -0.065 (mean2) and -0.315 (mean4), mass to 1e-12 %:
one realisation each, not predicted, reported as a lead for the fragmentation work and
nothing more. The identity defect on the deformation N 64 runs is 7.6e-5 to 1.3e-4
(relative) against 1e-15 serial and 1e-12 to 1e-11 on the rotation in parallel: the
state of these parallel deformation runs carries double-valued processor faces among
the faces no plane touches (the swept flux's own sync on the sheet's fragments); the
faces the corrector sets are single-valued by construction. To be located in the swept
flux, not here.

**P4, cost: HELD.** Rotation N 64, OF Co 0.2, 8 ranks, one exclusive i01 node
(mpsc0203), 864 steps, the diagnostics at the end, two repetitions within 0.6 %: the
median step 0.190 s (none), 0.200 s (mean1), 0.210 s (mean2), 0.2375 s (mean4), i.e.
+10, +20, +48 ms per step. The profile (OpenFOAM's profiling, one run each): none spends
9.83 s of 156.3 s of `advect` in `reconstruct` over 864 steps, 11.4 ms per step (6.3 % of
the step; `curvedFlux` 137.5 s); mean1 spends 20.87 s in `correctFaceFractions` (two
reconstructions 17.9 s, the induced fractions 2.6 s, the exchange of the ring normal
inside) and 0.013 s in the skipped start-of-step reconstruction, so one iteration adds
12.8 ms net = 1.12 reconstructions (the clause: at most 1.2), and four add 48 ms = 4.2
(the clause: at most 4.8). The DAVOF reconstruction's share of the swept step is small;
the corrector's cost is the reconstructions it repeats.

**P5, the first correction as a diagnostic: HELD for the L1, the max is O(1).** The L1
change of alpha_f at iteration 1 (median over the steps) is 5.41e-3, 2.55e-3, 1.28e-3 at
N 32, 64, 128 (swept flux, OF Co 0.2): order 1.04; at Co 1 5.90e-3, 2.83e-3: order 1.06;
under isoAdvector 5.23e-3, 2.49e-3: order 1.07. The two fluxes' inconsistency is the same
size within 3 % (the swept pieces' clips and isoAdvector's translated planes leave the
same first-order inconsistency with the planes reconstructed from alpha^{n+1}). The max
change is 0.175, 0.161, 0.200 (order -0.09): it is set by isolated faces (grazing cuts,
wisps) at every N, not by h.

**What is demonstrated, what is a limitation, what is open.** Demonstrated: the
corrector as implemented leaves a planar state unchanged to round-off at the first
evaluation (the unit check) and to 1e-12 .. 1e-10 under transport at every N and Co
except the one cell where the swept flux's own round-off wisps enter the band; the
default path is bit for bit the old one; the iteration contracts on the sphere by a
factor 5 to 50 per iteration (median) with the first correction first order in h and the
same for the swept and the isoAdvector flux; the cost is one reconstruction per
iteration; the identity holds to round-off wherever every face is single-valued, and it
locates double-valued processor faces where it is not. Limitations: on the sphere the
fixed point is no more accurate than the transported state (every error +2 to +5 %, no
gain with iterations: the ring-area normal's first order is not improved by making the
face fractions consistent with it), and on the deformation the one-plane-per-cell
representation cannot carry a sub-cell sheet across a face, so the max change never
contracts and the shape error grows by 24 % at N 32, while the bound improves and N 64
completes with 2 and 4 iterations. Open: the interface set at the tolerance band's edge
(round-off wisps of alpha ~1e-8 get sliver planes; the normal-based rule excludes them),
the sliver positioning tolerance (1e-14 of V is 1e-10 h for a corner of 1e-7 V), the
double-valued processor faces of the parallel deformation runs (the swept flux's), and
whether a different normal model at the fixed point (the ring fit, second order) gains
what the ring area normal cannot.

Side findings for the other tracks: the port's `reinitFaceFractions` diagnostic
(2026-11-06) clipped with the area vector and returned 0 on fully wet faces of sliver
cuts (fixed here in 9c0d983d, the same line); the none np 1 / np 8 pair of the rotation
differs by up to 2.4e-9 per face after 20 steps (L1 4e-13), the swept flux's own
parallel scatter, which the port's L1 criterion does not see.

#### The plane gate (run 3, 914505cd; mean4 = the corrector with 4 iterations, tolerance 1e-12; it = iterations used per step, median/max; d1 = the change of alpha_f at iteration 1, median over the steps, max over the faces / L1 over the set faces; ident = |sum_k m_k - sum_b alpha_f S_f| / sum_k |m_k|, max over the steps; serial)

| arm | Co | N | corrector | it med/max | d1 max | d1 L1 | L1(alphaf) end | L1(alphaf) max | E_shape | L1(n) | min alpha | max alpha - 1 | interface cells | no plane | ident | steps | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| normal | 0.2 | 20 | mean4 | 0/1 | 8.37e-13 | 3.72e-14 | 3.75e-14 | 7.09e-14 | 8.73e-24 | 6.63e-14 | -1.62e-14 | 1.85e-13 | 19 | 0 | 9.06e-16 | 20 | completed |
| normal | 0.2 | 40 | mean4 | 1/4 | 1.84e-11 | 1.89e-13 | 1.25e-13 | 3.18e-13 | 1.32e-23 | 1.43e-13 | -1.16e-13 | 8.84e-13 | 71 | 0 | 4.60e-15 | 40 | completed |
| normal | 0.2 | 80 | mean4 | 4/4 | 1.16e-10 | 7.17e-13 | 1.16e-12 | 1.74e-12 | 7.50e-23 | 5.43e-13 | -1.04e-14 | 3.27e-12 | 287 | 0 | 2.41e-14 | 80 | completed |
| normal | 0.2 | 20 | none | -/- | - | - | 4.33e-14 | 4.33e-14 | 9.42e-24 | 6.05e-14 | -2.75e-14 | 1.84e-13 | 24 | - | - | 20 | completed |
| normal | 0.2 | 40 | none | -/- | - | - | 7.05e-14 | 9.57e-14 | 1.25e-23 | 8.97e-14 | -9.54e-14 | 8.84e-13 | 80 | - | - | 40 | completed |
| normal | 0.2 | 80 | none | -/- | - | - | 1.06e-12 | 1.06e-12 | 7.66e-23 | 4.62e-13 | -1.12e-14 | 3.27e-12 | 302 | - | - | 80 | completed |
| normal | 0.5 | 20 | mean4 | 0.5/1 | 1.05e-12 | 4.09e-14 | 1.05e-14 | 2.99e-14 | 6.20e-24 | 2.41e-14 | -2.65e-14 | 4.75e-14 | 19 | 0 | 6.68e-16 | 8 | completed |
| normal | 0.5 | 40 | mean4 | 2/4 | 5.76e-11 | 4.81e-13 | 3.30e-13 | 1.58e-12 | 2.53e-23 | 5.69e-13 | -6.51e-13 | 3.27e-13 | 71 | 0 | 4.63e-15 | 16 | completed |
| normal | 0.5 | 80 | mean4 | 4/4 | 9.51e-10 | 1.95e-12 | 1.60e-12 | 1.60e-12 | 7.26e-23 | 3.23e-12 | -5.83e-13 | 1.74e-12 | 287 | 0 | 1.10e-14 | 32 | completed |
| normal | 0.5 | 20 | none | -/- | - | - | 1.06e-14 | 2.87e-14 | 6.44e-24 | 1.94e-14 | -3.59e-14 | 4.75e-14 | 32 | - | - | 8 | completed |
| normal | 0.5 | 40 | none | -/- | - | - | 3.28e-13 | 3.82e-13 | 1.76e-23 | 7.06e-13 | -1.64e-13 | 3.28e-13 | 94 | - | - | 16 | completed |
| normal | 0.5 | 80 | none | -/- | - | - | 1.20e-12 | 1.31e-12 | 5.07e-23 | 2.43e-12 | -6.95e-13 | 1.74e-12 | 330 | - | - | 32 | completed |
| normal | 1 | 20 | mean4 | 1/1 | 5.50e-12 | 2.23e-13 | 4.60e-13 | 4.60e-13 | 6.29e-23 | 6.97e-13 | -2.66e-13 | 3.06e-14 | 19 | 0 | 6.97e-16 | 4 | completed |
| normal | 1 | 40 | mean4 | 3.5/4 | 1.94e-09 | 1.80e-11 | 4.91e-11 | 6.14e-11 | 6.94e-21 | 5.45e-11 | -3.43e-10 | 2.78e-13 | 71 | 0 | 3.92e-15 | 8 | completed |
| normal | 1 | 80 | mean4 | 4/4 | 5.05e-08 | 2.55e-10 | 0.009704 | 0.009704 | 1.13e-15 | 0.133 | -9.23e-04 | 6.65e-12 | 325 | 0 | 7.42e-15 | 16 | completed |
| normal | 1 | 20 | none | -/- | - | - | 1.10e-13 | 1.10e-13 | 5.38e-23 | 1.31e-13 | -1.76e-13 | 3.06e-14 | 48 | - | - | 4 | completed |
| normal | 1 | 40 | none | -/- | - | - | 1.57e-12 | 1.57e-12 | 6.07e-22 | 2.33e-12 | -1.52e-11 | 2.78e-13 | 121 | - | - | 8 | completed |
| normal | 1 | 80 | none | -/- | - | - | 1.19e-11 | 1.19e-11 | 2.25e-19 | 1.59e-11 | -1.22e-08 | 6.65e-12 | 378 | - | - | 16 | completed |
| oblique | 0.2 | 20 | mean4 | 1.5/4 | 4.71e-12 | 7.90e-14 | 7.59e-14 | 2.13e-13 | 6.92e-23 | 1.01e-13 | -1.11e-14 | 2.31e-13 | 132 | 0 | 7.06e-16 | 20 | completed |
| oblique | 0.2 | 40 | mean4 | 4/4 | 1.11e-10 | 6.27e-13 | 4.36e-13 | 7.72e-13 | 2.12e-22 | 5.25e-13 | -1.81e-14 | 2.66e-12 | 521 | 0 | 4.82e-15 | 40 | completed |
| oblique | 0.2 | 80 | mean4 | 4/4 | 7.86e-10 | 2.17e-12 | 1.71e-12 | 3.17e-12 | 5.12e-22 | 1.99e-12 | -3.71e-14 | 7.21e-12 | 2.09e+03 | 0 | 1.24e-14 | 80 | completed |
| oblique | 0.2 | 20 | none | -/- | - | - | 8.41e-14 | 8.41e-14 | 7.06e-23 | 1.29e-13 | -9.87e-15 | 2.31e-13 | 131 | - | - | 20 | completed |
| oblique | 0.2 | 40 | none | -/- | - | - | 3.04e-13 | 4.20e-13 | 2.03e-22 | 3.33e-13 | -3.06e-14 | 2.66e-12 | 525 | - | - | 40 | completed |
| oblique | 0.2 | 80 | none | -/- | - | - | 1.25e-12 | 1.43e-12 | 4.92e-22 | 1.23e-12 | -3.91e-14 | 7.02e-12 | 2.06e+03 | - | - | 80 | completed |
| oblique | 0.5 | 20 | mean4 | 1/1 | 1.54e-12 | 4.26e-14 | 3.08e-14 | 3.54e-14 | 3.33e-23 | 3.74e-14 | -1.10e-14 | 1.87e-13 | 132 | 0 | 5.37e-16 | 8 | completed |
| oblique | 0.5 | 40 | mean4 | 3.5/4 | 7.91e-11 | 3.65e-13 | 2.85e-13 | 2.94e-13 | 9.31e-23 | 4.14e-13 | -3.94e-14 | 7.37e-13 | 521 | 0 | 4.00e-15 | 16 | completed |
| oblique | 0.5 | 80 | mean4 | 4/4 | 3.11e-10 | 1.09e-12 | 8.14e-13 | 9.11e-13 | 2.34e-22 | 8.85e-13 | -1.15e-13 | 4.18e-12 | 2.09e+03 | 0 | 9.91e-15 | 32 | completed |
| oblique | 0.5 | 20 | none | -/- | - | - | 2.89e-14 | 3.30e-14 | 3.43e-23 | 3.07e-14 | -1.15e-14 | 1.87e-13 | 136 | - | - | 8 | completed |
| oblique | 0.5 | 40 | none | -/- | - | - | 1.43e-13 | 1.71e-13 | 8.85e-23 | 1.40e-13 | -4.06e-14 | 7.37e-13 | 534 | - | - | 16 | completed |
| oblique | 0.5 | 80 | none | -/- | - | - | 7.06e-13 | 7.12e-13 | 2.29e-22 | 7.31e-13 | -1.15e-13 | 4.18e-12 | 2.10e+03 | - | - | 32 | completed |
| oblique | 1 | 20 | mean4 | 1/1 | 1.30e-12 | 3.59e-14 | 4.20e-14 | 4.20e-14 | 2.52e-23 | 5.62e-14 | -1.40e-14 | 6.37e-14 | 132 | 0 | 2.77e-16 | 4 | completed |
| oblique | 1 | 40 | mean4 | 2/4 | 4.28e-11 | 2.34e-13 | 2.45e-13 | 2.47e-13 | 7.46e-23 | 2.91e-13 | -1.26e-13 | 4.58e-13 | 521 | 0 | 4.64e-15 | 8 | completed |
| oblique | 1 | 80 | mean4 | 4/4 | 2.80e-10 | 1.14e-12 | 1.57e-12 | 1.57e-12 | 1.83e-22 | 1.93e-12 | -2.08e-13 | 2.95e-12 | 2.09e+03 | 0 | 8.82e-15 | 16 | completed |
| oblique | 1 | 20 | none | -/- | - | - | 4.68e-14 | 4.68e-14 | 2.68e-23 | 6.12e-14 | -1.66e-14 | 6.33e-14 | 147 | - | - | 4 | completed |
| oblique | 1 | 40 | none | -/- | - | - | 2.26e-13 | 2.33e-13 | 7.57e-23 | 2.76e-13 | -8.43e-14 | 4.58e-13 | 553 | - | - | 8 | completed |
| oblique | 1 | 80 | none | -/- | - | - | 1.24e-12 | 1.24e-12 | 1.83e-22 | 1.32e-12 | -2.15e-13 | 2.95e-12 | 2.14e+03 | - | - | 16 | completed |
| normal | 0.2 | mean4 | -1.552 | -1.517 | -1.513 | -3.56 ||||||||||||
| normal | 0.2 | none | -1.512 | -1.467 | -0.8392 | - ||||||||||||
| normal | 0.5 | mean4 | -1.775 | -3.533 | -1.219 | -4.913 ||||||||||||
| normal | 0.5 | none | -1.487 | -3.484 | -0.7899 | - ||||||||||||
| normal | 1 | mean4 | -12.05 | -18.74 | -17.21 | -6.582 ||||||||||||
| normal | 1 | none | -6.016 | -3.462 | -1.985 | - ||||||||||||
| oblique | 0.2 | mean4 | -1.443 | -2.154 | -2.176 | -3.691 ||||||||||||
| oblique | 0.2 | none | -1.399 | -1.623 | -1.5 | - ||||||||||||
| oblique | 0.5 | mean4 | -1.404 | -2.281 | -1.914 | -3.829 ||||||||||||
| oblique | 0.5 | none | -1.368 | -2.287 | -1.661 | - ||||||||||||
| oblique | 1 | mean4 | -1.43 | -2.551 | -2.003 | -3.872 ||||||||||||
| oblique | 1 | none | -1.385 | -2.217 | -1.579 | - ||||||||||||

#### The sphere rotation (one turn; OF Co 0.2 = interface Co 0.0463, OF Co 1 = 0.2315; N 32 / 64 / 128 on 1 / 8 / 64 ranks), the deformation (isoDeformation3D, OF Co 0.5, grazingCut fork, N 32 / 64 on 1 / 8 ranks) and the np 1 / np 8 pairs (20 steps); run 3. fac_i = the contraction factor max|dalpha_f|(i) / max|dalpha_f|(i+1), median / worst over the steps; wall = ClockTime at the end on shared nodes (the timed comparison is the table below)

| study | scheme | Co | N | corrector | it med/max | d1 max | d1 L1 | fac1 med/worst | fac2 | fac3 | L1(alphaf) end | E_shape | L1(x) | L1(n) | min alpha | max alpha - 1 | mass dev % | interface cells | no plane | ident | steps | wall s | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| iter_def_N32 | davofSwept | 0.5 | 32 | mean1 | 1/1 | 1 | 0.04082 | -/- | -/- | -/- | 0.03764 | 0.01407 | 0.275 | 1.299 | -0.05041 | 0.00129 | -6.88e-13 | 4.49e+03 | 51 | 1.41e-15 | 480 | 301 | completed |
| iter_def_N32 | davofSwept | 0.5 | 32 | mean2 | 2/2 | 1 | 0.05077 | 1/1 | -/- | -/- | 0.04155 | 0.01486 | 0.296 | 1.289 | -0.06278 | 0.001915 | -8.44e-13 | 4.48e+03 | 54 | 1.15e-15 | 480 | 267 | completed |
| iter_def_N32 | davofSwept | 0.5 | 32 | mean4 | 4/4 | 1 | 0.05206 | 1/1 | 1/0.9999 | 1/0.407 | 0.05065 | 0.01537 | 0.284 | 1.396 | -0.06074 | 0.004057 | -6.66e-14 | 4.73e+03 | 101 | 1.24e-15 | 480 | 285 | completed |
| iter_def_N32 | davofSwept | 0.5 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.07063 | 0.01137 | 0.3063 | 1.46 | -0.2701 | 2.95e-04 | -8.44e-13 | 1.08e+04 | - | - | 480 | 466 | completed |
| iter_def_N64 | davofSwept | 0.5 | 64 | mean1 | 1/1 | 1 | 0.03768 | -/- | -/- | -/- | - | - | - | - | -0.7715 | 0.001959 | 3.77e-13 | 2.64e+04 | 133 | 7.91e-05 | 889 | 446 | diverged |
| iter_def_N64 | davofSwept | 0.5 | 64 | mean2 | 2/2 | 1 | 0.05633 | 1/1 | -/- | -/- | 0.04164 | 0.007841 | 0.2792 | 1.297 | -0.06511 | 0.0118 | 1.13e-12 | 2.87e+04 | 223 | 7.64e-05 | 960 | 559 | completed |
| iter_def_N64 | davofSwept | 0.5 | 64 | mean4 | 4/4 | 1 | 0.05221 | 1/1 | 1/0.8029 | 1/0.5992 | 0.0483 | 0.007724 | 0.2865 | 1.314 | -0.315 | 0.008982 | 3.44e-13 | 2.96e+04 | 391 | 1.32e-04 | 960 | 668 | completed |
| iter_def_N64 | davofSwept | 0.5 | 64 | none | -/- | - | - | -/- | -/- | -/- | - | - | - | - | -0.5446 | 0.00896 | -2.00e-13 | 5.11e+04 | - | - | 628 | 828 | diverged |
| iter_par_plane_np1 | davofSwept | 0.5 | 40 | mean4 | 3/4 | 8.71e-11 | 3.70e-13 | 39.2/6.135 | 1.751/0.5046 | 1.481/0.5672 | 1.97e-13 | 7.90e-23 | 8.22e-16 | 2.25e-13 | -3.94e-14 | 7.37e-13 | 60.39 | 452 | 0 | 4.00e-15 | 20 | 5 | completed |
| iter_par_plane_np8 | davofSwept | 0.5 | 40 | mean4 | 4/4 | 8.64e-11 | 3.80e-13 | 41.44/5.867 | 1.469/0.7695 | 1.008/0.782 | 1.96e-13 | 7.82e-23 | 8.08e-16 | 2.26e-13 | -2.86e-14 | 7.37e-13 | 60.39 | 452 | 0 | 3.13e-15 | 20 | 2 | completed |
| iter_par_rot_np1 | davofSwept | 0.0462996198086 | 32 | mean4 | 4/4 | 0.04227 | 0.003093 | 13.42/4.309 | 7.014/5.585 | 7.377/5.3 | 0.007105 | 4.96e-05 | 5.62e-04 | 0.0385 | -1.10e-14 | 5.73e-13 | 7.77e-14 | 444 | 0 | 4.76e-16 | 20 | 10 | completed |
| iter_par_rot_np1_none | davofSwept | 0.0462996198086 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.007048 | 5.01e-05 | 5.58e-04 | 0.03834 | -7.62e-15 | 5.73e-13 | 1.22e-13 | 433 | - | - | 20 | 11 | completed |
| iter_par_rot_np8 | davofSwept | 0.0462996198086 | 32 | mean4 | 4/4 | 0.04227 | 0.003093 | 13.42/4.309 | 7.014/5.585 | 7.377/5.3 | 0.007105 | 4.96e-05 | 5.62e-04 | 0.0385 | -1.06e-14 | 5.73e-13 | -4.44e-14 | 444 | 0 | 8.37e-16 | 20 | 5 | completed |
| iter_par_rot_np8_none | davofSwept | 0.0462996198086 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.007048 | 5.01e-05 | 5.58e-04 | 0.03834 | -6.76e-15 | 5.73e-13 | -8.88e-14 | 433 | - | - | 20 | 5 | completed |
| iter_rotIso_N32 | isoAdvector | 0.0462996198086 | 32 | mean1 | 1/1 | 0.1744 | 0.005227 | -/- | -/- | -/- | 0.06019 | 0.001131 | 0.004237 | 0.06012 | -5.12e-19 | 2.22e-16 | 1.12e-12 | 422 | 0 | 8.79e-16 | 1728 | 239 | completed |
| iter_rotIso_N32 | isoAdvector | 0.0462996198086 | 32 | mean4 | 4/4 | 0.1745 | 0.005206 | 21.72/4.183 | 9.356/2.923 | 5.852/2.687 | 0.06042 | 0.001135 | 0.004251 | 0.0603 | -1.29e-19 | 2.22e-16 | 1.32e-12 | 422 | 0 | 5.97e-16 | 1728 | 257 | completed |
| iter_rotIso_N32 | isoAdvector | 0.0462996198086 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.05831 | 0.001091 | 0.004119 | 0.05893 | -6.79e-14 | 2.22e-16 | 8.44e-13 | 423 | - | - | 1728 | 235 | completed |
| iter_rotIso_N64 | isoAdvector | 0.0462996198086 | 64 | mean1 | 1/1 | 0.1601 | 0.002491 | -/- | -/- | -/- | 0.03343 | 3.08e-04 | 0.001109 | 0.02658 | -5.14e-19 | 2.22e-16 | 7.44e-13 | 1.77e+03 | 0 | 8.41e-16 | 3456 | 952 | completed |
| iter_rotIso_N64 | isoAdvector | 0.0462996198086 | 64 | mean4 | 4/4 | 0.1594 | 0.002488 | 50.36/13.03 | 7.675/1.995 | 5.993/2.13 | 0.03346 | 3.08e-04 | 0.001111 | 0.02666 | -1.29e-19 | 2.22e-16 | 9.99e-13 | 1.77e+03 | 0 | 9.36e-16 | 3456 | 642 | completed |
| iter_rotIso_N64 | isoAdvector | 0.0462996198086 | 64 | none | -/- | - | - | -/- | -/- | -/- | 0.03272 | 3.02e-04 | 0.001088 | 0.02565 | -9.51e-14 | 2.22e-16 | 1.44e-13 | 1.75e+03 | - | - | 3456 | 898 | completed |
| iter_rot_N128 | davofSwept | 0.0462996198086 | 128 | mean1 | 1/1 | 0.1981 | 0.00128 | -/- | -/- | -/- | 0.01851 | 8.14e-05 | 2.87e-04 | 0.01159 | -3.22e-14 | 6.64e-10 | 7.58e-12 | 6.97e+03 | 0 | 1.40e-12 | 6912 | 3.24e+03 | completed |
| iter_rot_N128 | davofSwept | 0.0462996198086 | 128 | mean2 | 2/2 | 0.1997 | 0.00128 | 50.94/8.459 | -/- | -/- | 0.01851 | 8.14e-05 | 2.87e-04 | 0.0116 | -1.11e-12 | 8.61e-10 | 7.67e-12 | 6.97e+03 | 0 | 1.41e-12 | 6912 | 3.36e+03 | completed |
| iter_rot_N128 | davofSwept | 0.0462996198086 | 128 | mean4 | 4/4 | 0.2002 | 0.00128 | 50.86/8.396 | 4.752/1.601 | 6.04/2.084 | 0.01851 | 8.14e-05 | 2.87e-04 | 0.01161 | -1.13e-12 | 5.20e-10 | 7.88e-12 | 6.97e+03 | 0 | 1.44e-12 | 6912 | 4.60e+03 | completed |
| iter_rot_N128 | davofSwept | 0.0462996198086 | 128 | none | -/- | - | - | -/- | -/- | -/- | 0.01829 | 8.05e-05 | 2.84e-04 | 0.01131 | -2.07e-14 | 6.45e-10 | 7.45e-12 | 6.87e+03 | - | - | 6912 | 3.14e+03 | completed |
| iter_rot_N32 | davofSwept | 0.0462996198086 | 32 | mean1 | 1/1 | 0.1749 | 0.005411 | -/- | -/- | -/- | 0.05933 | 0.001097 | 0.004127 | 0.05961 | -8.74e-13 | 1.57e-12 | 2.40e-12 | 424 | 0 | 3.22e-15 | 1728 | 306 | completed |
| iter_rot_N32 | davofSwept | 0.0462996198086 | 32 | mean2 | 2/2 | 0.1758 | 0.005397 | 18.76/4.321 | -/- | -/- | 0.05933 | 0.001101 | 0.004173 | 0.0593 | -9.76e-13 | 5.13e-11 | 2.43e-12 | 422 | 0 | 3.38e-15 | 1728 | 294 | completed |
| iter_rot_N32 | davofSwept | 0.0462996198086 | 32 | mean4 | 4/4 | 0.1761 | 0.005398 | 18.85/4.309 | 8.329/2.622 | 5.959/2.773 | 0.0596 | 0.001102 | 0.004175 | 0.05933 | -6.89e-13 | 4.88e-12 | 2.10e-12 | 422 | 0 | 3.27e-15 | 1728 | 326 | completed |
| iter_rot_N32 | davofSwept | 0.0462996198086 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.05704 | 0.001059 | 0.003972 | 0.05749 | -8.12e-13 | 1.99e-11 | 2.14e-12 | 423 | - | - | 1728 | 270 | completed |
| iter_rot_N32 | davofSwept | 0.231498099043 | 32 | mean1 | 1/1 | 0.1535 | 0.005898 | -/- | -/- | -/- | 0.04032 | 7.40e-04 | 0.002748 | 0.04976 | -9.99e-13 | 9.09e-13 | 2.03e-12 | 420 | 0 | 3.29e-15 | 346 | 61 | completed |
| iter_rot_N32 | davofSwept | 0.231498099043 | 32 | mean2 | 2/2 | 0.1536 | 0.005889 | 17.09/4.782 | -/- | -/- | 0.04051 | 7.43e-04 | 0.002759 | 0.04992 | -9.74e-13 | 9.82e-13 | 2.25e-12 | 420 | 0 | 3.57e-15 | 346 | 68 | completed |
| iter_rot_N32 | davofSwept | 0.231498099043 | 32 | mean4 | 4/4 | 0.1536 | 0.005892 | 17.1/4.766 | 8.714/2.645 | 5.727/2.771 | 0.04053 | 7.44e-04 | 0.00276 | 0.04995 | -9.76e-13 | 9.89e-13 | 2.33e-12 | 420 | 0 | 3.64e-15 | 346 | 72 | completed |
| iter_rot_N32 | davofSwept | 0.231498099043 | 32 | none | -/- | - | - | -/- | -/- | -/- | 0.03824 | 7.07e-04 | 0.002653 | 0.04828 | -8.67e-13 | 8.47e-13 | 2.32e-12 | 415 | - | - | 346 | 59 | completed |
| iter_rot_N64 | davofSwept | 0.0462996198086 | 64 | mean1 | 1/1 | 0.1606 | 0.002549 | -/- | -/- | -/- | 0.03278 | 3.02e-04 | 0.001087 | 0.02645 | -6.61e-14 | 5.51e-10 | 4.47e-12 | 1.76e+03 | 0 | 2.45e-12 | 3456 | 984 | completed |
| iter_rot_N64 | davofSwept | 0.0462996198086 | 64 | mean2 | 2/2 | 0.1608 | 0.002548 | 40.95/8.845 | -/- | -/- | 0.03281 | 3.02e-04 | 0.001087 | 0.02648 | -2.19e-14 | 1.70e-10 | 4.39e-12 | 1.76e+03 | 0 | 2.92e-12 | 3456 | 938 | completed |
| iter_rot_N64 | davofSwept | 0.0462996198086 | 64 | mean4 | 4/4 | 0.1609 | 0.002548 | 40.97/8.808 | 5.893/1.837 | 6.219/2.525 | 0.03281 | 3.02e-04 | 0.001087 | 0.02648 | -1.64e-14 | 1.17e-10 | 4.57e-12 | 1.76e+03 | 0 | 2.81e-12 | 3456 | 1.01e+03 | completed |
| iter_rot_N64 | davofSwept | 0.0462996198086 | 64 | none | -/- | - | - | -/- | -/- | -/- | 0.03211 | 2.95e-04 | 0.00106 | 0.02523 | -2.02e-14 | 1.72e-10 | 4.17e-12 | 1.76e+03 | - | - | 3456 | 981 | completed |
| iter_rot_N64 | davofSwept | 0.231498099043 | 64 | mean1 | 1/1 | 0.134 | 0.002825 | -/- | -/- | -/- | 0.02235 | 2.02e-04 | 7.28e-04 | 0.02351 | -4.25e-13 | 3.22e-11 | 2.52e-12 | 1.77e+03 | 0 | 1.10e-11 | 691 | 191 | completed |
| iter_rot_N64 | davofSwept | 0.231498099043 | 64 | mean2 | 2/2 | 0.1338 | 0.002824 | 44.05/11.33 | -/- | -/- | 0.02242 | 2.02e-04 | 7.29e-04 | 0.02357 | -3.88e-13 | 5.22e-11 | 2.95e-12 | 1.77e+03 | 0 | 1.09e-11 | 691 | 215 | completed |
| iter_rot_N64 | davofSwept | 0.231498099043 | 64 | mean4 | 4/4 | 0.1338 | 0.002825 | 44.04/11.33 | 7.193/1.973 | 6.338/2.52 | 0.02242 | 2.02e-04 | 7.29e-04 | 0.02357 | -4.02e-13 | 4.87e-11 | 2.95e-12 | 1.77e+03 | 0 | 1.11e-11 | 691 | 196 | completed |
| iter_rot_N64 | davofSwept | 0.231498099043 | 64 | none | -/- | - | - | -/- | -/- | -/- | 0.02174 | 1.97e-04 | 7.07e-04 | 0.02257 | -6.07e-13 | 1.89e-10 | 3.05e-12 | 1.75e+03 | - | - | 691 | 213 | completed |
| iter_def | davofSwept | 0.5 | mean1 | - | - | - | 0 | 0.1155 ||||||||||||||||
| iter_def | davofSwept | 0.5 | mean2 | 0.9219 | -0.009266 | 0.08423 | 0 | -0.15 ||||||||||||||||
| iter_def | davofSwept | 0.5 | mean4 | 0.9927 | 0.08771 | -0.0129 | 0 | -0.004169 ||||||||||||||||
| iter_def | davofSwept | 0.5 | none | - | - | - | - | - ||||||||||||||||
| iter_rot | davofSwept | 0.0462996198086 | mean1 | 1.877 | 1.182 | 1.923 | -0.08987 | 1.04 ||||||||||||||||
| iter_rot | davofSwept | 0.0462996198086 | mean2 | 1.879 | 1.177 | 1.931 | -0.09193 | 1.038 ||||||||||||||||
| iter_rot | davofSwept | 0.0462996198086 | mean4 | 1.879 | 1.177 | 1.931 | -0.09226 | 1.038 ||||||||||||||||
| iter_rot | davofSwept | 0.0462996198086 | none | 1.859 | 1.173 | 1.903 | - | - ||||||||||||||||
| iter_rot | davofSwept | 0.231498099043 | mean1 | 1.874 | 1.082 | 1.916 | 0.1959 | 1.062 ||||||||||||||||
| iter_rot | davofSwept | 0.231498099043 | mean2 | 1.877 | 1.083 | 1.92 | 0.1984 | 1.06 ||||||||||||||||
| iter_rot | davofSwept | 0.231498099043 | mean4 | 1.878 | 1.083 | 1.921 | 0.1988 | 1.061 ||||||||||||||||
| iter_rot | davofSwept | 0.231498099043 | none | 1.848 | 1.097 | 1.908 | - | - ||||||||||||||||
| iter_rotIso | isoAdvector | 0.0462996198086 | mean1 | 1.878 | 1.177 | 1.933 | 0.1233 | 1.069 ||||||||||||||||
| iter_rotIso | isoAdvector | 0.0462996198086 | mean4 | 1.881 | 1.178 | 1.936 | 0.131 | 1.065 ||||||||||||||||
| iter_rotIso | isoAdvector | 0.0462996198086 | none | 1.855 | 1.2 | 1.92 | - | - ||||||||||||||||

#### The orders over N (least squares; d1 = the first correction)

| family | scheme | Co | corrector | N | order L1(alphaf) | order E_shape | order L1(n) | order L1(x) | order d1 max | order d1 L1 |
|---|---|---|---|---|---|---|---|---|---|---|
| iter_def | davofSwept | 0.5 | mean1 | 32 64 | - | - | - | - | 0 | 0.1155 |
| iter_def | davofSwept | 0.5 | mean2 | 32 64 | -0.003361 | 0.9219 | -0.009266 | 0.08423 | 0 | -0.15 |
| iter_def | davofSwept | 0.5 | mean4 | 32 64 | 0.06852 | 0.9927 | 0.08771 | -0.0129 | 0 | -0.004169 |
| iter_def | davofSwept | 0.5 | none | 32 64 | - | - | - | - | - | - |
| iter_rot | davofSwept | 0.0462996198086 | mean1 | 32 64 128 | 0.8404 | 1.877 | 1.182 | 1.923 | -0.08987 | 1.04 |
| iter_rot | davofSwept | 0.0462996198086 | mean2 | 32 64 128 | 0.8402 | 1.879 | 1.177 | 1.931 | -0.09193 | 1.038 |
| iter_rot | davofSwept | 0.0462996198086 | mean4 | 32 64 128 | 0.8436 | 1.879 | 1.177 | 1.931 | -0.09226 | 1.038 |
| iter_rot | davofSwept | 0.0462996198086 | none | 32 64 128 | 0.8204 | 1.859 | 1.173 | 1.903 | - | - |
| iter_rot | davofSwept | 0.231498099043 | mean1 | 32 64 | 0.8511 | 1.874 | 1.082 | 1.916 | 0.1959 | 1.062 |
| iter_rot | davofSwept | 0.231498099043 | mean2 | 32 64 | 0.8537 | 1.877 | 1.083 | 1.92 | 0.1984 | 1.06 |
| iter_rot | davofSwept | 0.231498099043 | mean4 | 32 64 | 0.8543 | 1.878 | 1.083 | 1.921 | 0.1988 | 1.061 |
| iter_rot | davofSwept | 0.231498099043 | none | 32 64 | 0.8145 | 1.848 | 1.097 | 1.908 | - | - |
| iter_rotIso | isoAdvector | 0.0462996198086 | mean1 | 32 64 | 0.8482 | 1.878 | 1.177 | 1.933 | 0.1233 | 1.069 |
| iter_rotIso | isoAdvector | 0.0462996198086 | mean4 | 32 64 | 0.8527 | 1.881 | 1.178 | 1.936 | 0.131 | 1.065 |
| iter_rotIso | isoAdvector | 0.0462996198086 | none | 32 64 | 0.8337 | 1.855 | 1.2 | 1.92 | - | - |

#### P4: rotation N 64, OF Co 0.2, 8 ranks, one exclusive i01 node (mpsc0203), 864 steps, the function objects at the end; the median of the per-step ExecutionTime differences

| setting | rep | steps | wall s | median step s |
|---|---|---|---|---|
| none | 1 | 864 | 168.61 | 0.1900 |
| mean1 | 1 | 864 | 178.12 | 0.2000 |
| mean2 | 1 | 864 | 189.15 | 0.2100 |
| mean4 | 1 | 864 | 211.29 | 0.2400 |
| mean1 | 2 | 864 | 179.13 | 0.2000 |
| mean2 | 2 | 864 | 189.70 | 0.2100 |
| mean4 | 2 | 864 | 210.63 | 0.2350 |
| none | 2 | 864 | 168.13 | 0.1900 |
| none (profiled) | 1 | 864 | 167.26 | 0.1800 |
| mean1 (profiled) | 1 | 864 | 178.20 | 0.2000 |


## Normal order versus normal bias (2026-11-12, leia branch `davof-normal-order`)

Tomislav's study of 2026-11-12: separate the ORDER of the interface normal from its BIAS,
and report the error |1 - n_k . n_e| with its convergence orders wherever a normal error
appears. Branch `davof-normal-order`, based on `davof-iter` 76db7b64, which carries the
unit-normal clip fix 9c0d983d (the port's ee6fddc8 carries the same fix in the
`reinitFaceFractions` diagnostic only; the standing pieces path does not call the clip,
so both bases run the same flux; davof-iter adds the face-fraction corrector, default
`none`, bit for bit the old path, verified on the standing deformation record).

### The argument under test (pre-registered before the first run)

CCS positions the plane inside the cell from alpha_k for a given normal, so the plane
pivots about a point of the cell: a normal whose angle error is theta_k = O(h^p) places
the plane with an in-cell position error O(h^{p+1}) (the lever arm is O(h)). The volume
flux of a face is the liquid volume in the slab the face sweeps in one step; the error a
tilted plane makes there is the sliver between the true and the tilted plane inside that
slab, O(Co h^{p+1}) relative to the cell volume, per interface cell and per step. Over the
O(1/(Co h)) steps of a fixed travel, the GLOBAL shape error has an order between p and
p + 1: close to p when the per-step errors keep their sign from step to step (a BIASED
normal estimator, whose error accumulates), close to p + 1 when they are unbiased and the
next reconstruction, from the transported alpha_k, absorbs them. A first-order normal
therefore does not imply a first-order scheme: the shape order is set by the bias of the
normal, not by its order. This is a heuristic error budget, not a proof; the study tests
its four consequences. The data so far straddle it: on the rotation (the corrector study,
2026-11-11) the ring-area normal converges at 1.18 in L1(n) while E_shape converges at
1.88; on the sphere translation ladder of 2026-10-25 E_shape converged at 1.35 with the
ring normal and 2.84 with the exact normal.

### The norms (every row, every N, every Co; least-squares and pairwise orders)

All over the interface cells of the reconstruction at the end of the run (the
`reconstructionError` function object), with n_k the unit normal of the positioned plane
and n_e the exact unit normal at the plane's centroid x_k (the gradient of the exact
implicit function there):

(a) the angle theta_k = acos(n_k . n_e), computed as 2 asin(|n_k - n_e|/2) (acos loses
half the digits at small angles): L1 (the mean) and Linf (the maximum);
(b) 1 - n_k . n_e: L1 and Linf. For unit vectors 1 - n_k . n_e = |n_k - n_e|^2/2 exactly,
and 1 - cos theta = theta^2/2 (1 + O(theta^2)): when the angle converges at order p this
norm reads 2p, so a first-order angle appears as second order in it. Consequently
L1(1 - n.n_e) = L2(|n - n_e|)^2/2 and Linf(1 - n.n_e) = Linf(|n - n_e|)^2/2 exactly, their
orders are exactly twice the L2 and Linf orders of |n - n_e|, and every finished case that
stored the L2 and Linf of |n - n_e| gives this norm without a re-run;
(c) the existing column L1(n) = `LNormalDiff1` of `reconstructionError`: the mean over the
interface cells of |n_k - n_e|, the length of the difference of the two unit vectors (not
the angle; |n_k - n_e| = 2 sin(theta_k/2) = theta_k (1 - theta_k^2/24 + ...), so it has
the angle's order); `LNormalDiff2` its root mean square, `LNormalDiffInf` its maximum;
(d) the in-cell position error: the distance of the plane centroid x_k from the exact
sphere, |R - |x_k - o||, L1 over the interface cells (`LCentreDiff1`);
(e) the per-step flux error on the translation: for every face of the interface cells, the
computed volume flux dVf of the step (after the swept flux, before the bounding) against
the exact one, the liquid volume of the prism between f - U dt and f under the exactly
moved sphere at t^n (the exact sphere-tetrahedron overlap of Strobl, Formella and
Poeschel 2016, the `exactSphere` volume rule, exact for a uniform translation),
|dVf - dVf_exact|/V of the upwind cell, the mean over the faces per step (`fluxErrL1`), and
the max over the steps of that mean (the median over the steps reported next to it);
(f) E_shape (`volumeFractionError`, against the exact end state), E_mass (the volume
deviation in %), min and max alpha over the run;
plus the normal error of the plicRDF benchmark, defined by H. Scheufler and J. Roenby,
Accurate and efficient surface reconstruction from volume fraction data on general meshes,
J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009, arXiv:1801.05382,
Section 4.1, Eq. (21), verbatim: "For the interface position, x_S, and normal, n̂_S,
inside a cell we define the two norms, L∞^x = max_i(|x_S,i − x^exact_S,i|) and
L∞^n = max_i(1 − n̂_S,i · n̂^exact_S,i), (21) where the maximum is taken over all interface
cells, x^exact_S,i is the exact interface position, and n̂^exact_S,i is the exact interface
unit normal." It is the Linf of (b) (their reference point of the exact normal is not
stated; ours is the plane centroid), reported as its own column so that the comparison
with their figures is like for like; they report it, not the angle, and so read a
first-order angle as second order. Wall clock: the ClockTime at the end of every case, at
equal ranks (8; the serial N 32 check on 1), on shared nodes, one run per case.

### The rows and the tests

Three normals on identical initial states (initDavofFields `exactSphere`), the same
swept-face scheme for all three, the pieces path in the standing Co 1 setting (the
round-4 keys; the Enright setting `recommended`: `sweptPieceVolume signedTets`,
`sweptBoundedFaceFraction scale`, `interfaceCellRule normalAndAlpha`; CCS positioning to
1e-14 of V; `grazingCut openfoam`): (1) `davof-planeSurf-exact`, the exact normal (the
control); (2) `davof-planeSurf-ringArea`, the ring area normal d_k = sum_{j in R(k)} m_j;
(3) `davof-planeSurf-youngs`, a deliberately biased first-order normal, Youngs' normal
-grad(alpha_k)/|grad(alpha_k)| (the code's sign convention points into the liquid) from
OpenFOAM's Gauss linear gradient of alpha_k (the gradient of OpenFOAM's gradAlpha), the
existing `youngs` model with the new key `youngsGradScheme "Gauss linear"` (D. L. Youngs,
Time-dependent multi-material flow with large fluid distortion, in: K. W. Morton,
M. J. Baines (eds.), Numerical Methods for Fluid Dynamics, Academic Press, 1982,
pp. 273-285; no DOI exists). Tests: the sphere translation (`sphereTranslation3D`, R 0.15,
U = (1 0 0), one unit of travel; OpenFOAM Co = the interface Co for the uniform velocity)
and the solid-body rotation of the standing family (`sphereRotation3D`, one turn); N 32,
64, 128 on 8 ranks; OpenFOAM Co 0.2 and 1 (the rotation: interface Co 0.0463 and 0.2315);
N 32 also serial. Configs `config/davof/nor_{trans,rot}_N{32,64,128}[_np1].yaml`
(`workflow/scripts/nor_make_configs.py`), tables `workflow/scripts/nor_tables.py`, cluster
directory `/work/scratch/tm83tomy/davof/agents/normal`.

### Predictions per row (p = the least-squares order of L1(theta) over N 32-128, per row, test and Co)

- **N-1 position:** the order of L1(d) is p + 1 within 0.3, for all three rows.
- **N-2 flux (translation):** the order over N of the max over the steps of the per-step
  L1 flux error is p + 1 within 0.3, for all three rows.
- **N-3 shape:** with q the order of E_shape, q - p >= 0.6 ("near p + 1") for the exact and
  the ring normal and q - p <= 0.4 ("near p") for the Youngs normal; 0.4 < q - p < 0.6 is
  undecided and reported as such.
- **N-4 the 1 - n.n_e norm:** the order of L1(1 - n.n_e) is 2p within 0.3, for all three rows
  (it is exactly twice the L2 order of |n - n_e|, so this tests whether the L2 and the L1
  of the angle share an order).
- **N-5 serial/parallel:** at N 32 the end values of E_shape, L1(theta), L1(d), L1(alphaf)
  of np 1 and np 8 agree to four significant digits (relative 1e-4; the swept flux's
  parallel scatter is 1e-9 per value per 20 steps, its growth over a whole run unmeasured).
- Expected sizes, not predictions: p about 2 for the exact normal (the exact normal at a
  centroid that is itself O(h^2) off), p about 1 for the ring and the Youngs normal.

Falsifiers stated with the predictions: the ring row's E_shape near p on BOTH tests means
the ring normal is biased (the fix is then a normal, not a flux); the ring row near p + 1 on
the translation but near p on the rotation means the flux path carries a rotation error
the plane gate does not see. N-1 or N-2 failing for a row falsifies the error budget's
first two steps for that row and voids its N-3 verdict.

The one-plane path (davof-kiss, KISS round 3: one plane per swept face) enters as rows 4-6
in a second pass if its plane gate passes; it is not waited for.

#### Addendum (2026-11-12, 22:50, after the first N 32 results and before any oblique run)

The axial translation at OpenFOAM Co 1 is an exact transport: U = (1 0 0) and dt = h move
the sphere by exactly one cell per step, the swept volume of every face is its upwind cell,
and the N 32 runs end with E_shape 2.5e-15 to 1.2e-13, the per-step flux error 7.6e-14 and
L1(alphaf) 7.7e-14 for all three normals (the normal errors at the end are those of the
initial state, shifted). The axial Co 1 rows therefore carry no order information about the
normal; N-2 and N-3 are not evaluated on them. Added, with the same predictions N-1 .. N-4:
the oblique translation `nor_transOb_N{32,64,128}` (8 ranks), U = (0.3 0.5 0.8124), |U| = 1,
over 0.3 (the sphere stays in the box), OpenFOAM Co = (|Ux| + |Uy| + |Uz|) dt/h = 1.6124
Co_target, so the targets 0.124039 and 0.620194 give OpenFOAM Co 0.2 and 1.

### MEASURED (2026-11-12/13; 54 cases on 8 ranks, 12 serial; every job completed)

Jobs: build 56255681; orchestrators 56255682-89 (axial translation N 32/64/128 and N 32
serial, rotation N 32/64/128 and N 32 serial), 56258347-49 (the oblique translation, post
hoc). Tables: `/work/scratch/tm83tomy/davof/agents/normal/csv/run1/{nor_cases.csv,
nor_orders.csv,nor_tables.md}` (every norm (a)-(f), least-squares and pairwise orders, wall
clock); the 1 - n.n_e columns of the finished studies in `csv/existing/{standing,corrector,
ladders}/`. Two rotation cases ended early: the Youngs row at N 128 aborted by the solver's
guard at both Courant numbers (min alpha -0.85 at step 4919 of 6912, -0.72 at step 658 of
1382); they are not resubmitted, the abort is the row's result.

The four orders, least squares over N 32, 64, 128 (p = L1 of the angle; d = L1 of the
in-cell position error; flux = max over the steps of the per-step L1 flux error relative to
the cell volume; q = E_shape; the last column L1(1 - n.n_e), exactly twice the L2 order of
|n - n_e|). OF Co is the OpenFOAM Courant number.

| test | OF Co | normal | p | d | flux | q | 1 - n.n_e | wall s at N 128 (8 ranks) |
|---|---|---|---|---|---|---|---|---|
| axial translation | 0.2 | exact | 1.63 | 1.92 | 0.97 | 1.95 | 2.28 | 314 |
| | | ring | 0.96 | 1.88 | 0.89 | 1.86 | 1.76 | 305 |
| | | Youngs | 0.02 | 1.02 | 0.10 | 1.08 | 0.00 | 469 |
| oblique translation | 0.2 | exact | 2.53 | 1.94 | 0.91 | 1.92 | 6.06 | 1040 |
| | | ring | 1.29 | 2.16 | 1.11 | 2.14 | 2.38 | 1060 |
| | | Youngs | 0.15 | 0.73 | -0.35 | 0.66 | 0.20 | 1360 |
| oblique translation | 1 | exact | 2.33 | 1.99 | 0.93 | 1.93 | 4.77 | 266 |
| | | ring | 1.22 | 2.37 | 1.34 | 2.34 | 2.27 | 282 |
| | | Youngs | 0.13 | 0.77 | -0.28 | 0.75 | 0.22 | 234 |
| rotation | 0.2 | exact | 1.83 | 1.69 | - | 1.67 | 3.61 | 9660 |
| | | ring | 1.17 | 1.90 | - | 1.86 | 2.14 | 9330 |
| | | Youngs (N 32, 64; N 128 aborted) | 0.18 | 1.22 | - | 1.16 (pair) | 0.33 | aborted |
| rotation | 1 | exact | 1.73 | 1.72 | - | 1.69 | 3.36 | 1760 |
| | | ring | 1.13 | 1.90 | - | 1.87 | 2.06 | 2030 |
| | | Youngs (N 32, 64; N 128 aborted) | 0.14 | 1.08 | - | 1.00 (pair) | 0.27 | aborted |

The axial translation at OF Co 1 is an exact one-cell shift per step (E_shape 1e-15, flux
error 1e-13 to 5e-12 for all three normals) and carries no order information (addendum).

Readings needed before the verdicts. (i) The Youngs normal from the Gauss linear gradient is
not a first-order normal under transport: its angle error is 0.12 to 0.22 at every N (order
0.02 to 0.18), a zero-order normal, as the static tables of the preprint find for gradAlpha.
The row tests the argument at p = 0, not at p = 1. (ii) The exact-normal row's angle column
is not a normal error of the scheme: the model takes the exact normal at the first-pass
centroid, and the score compares it with the exact normal at the final centroid, so the
column measures the centroid's displacement. Its mechanical verdicts against "p" are void;
the row is the control with no normal error.

**N-1 (position order p + 1): HELD for the ring normal on every test** (d 1.88, 2.16, 2.37,
1.90, 1.90 against p + 1 = 1.96, 2.29, 2.22, 2.17, 2.13; within 0.3) and for the Youngs
normal on the axial translation and the rotation (1.02 against 1.02; 1.22 and 1.08 against
1.18 and 1.14), FAILED for Youngs on the oblique translation (0.73 and 0.77 against 1.15 and
1.13). The position is one order above the angle: the pivot of the CCS-positioned plane.

**N-2 (per-step flux order p + 1): FALSIFIED for all three rows.** Measured relative to the
cell volume the per-step flux error has the order of the angle, not one more: ring 0.89,
1.11, 1.34 (p 0.96, 1.29, 1.22), Youngs 0.10, -0.35, -0.28 (p 0.02, 0.15, 0.13). The control
with the exact normal has order 0.91 to 0.97: a plane through a sphere leaves the sagitta
sliver in every swept slab, O(Co h/R) of the cell volume per step, whatever the normal. The
pre-registered budget lost one power of h in the normalisation: a tilt theta changes a face's
wetted fraction by O(theta), and the flux through it by O(Co theta) of the cell volume.

**N-3 (shape order near p + 1 for exact and ring, near p for Youngs): HELD for the ring
normal on the axial translation, on both oblique translations and on both rotations**
(q - p = 0.90, 0.85, 1.13, 0.69, 0.74); **FALSIFIED for the Youngs normal**, whose shape
error converges one order above its angle as well (q - p = 1.06 on the axial translation,
about 1 on the rotation pairs; 0.52 and 0.62 on the oblique translation, undecided and
p + 1). In every row and test the order of E_shape equals the order of the in-cell position
error within 0.05 (1.95/1.92, 1.86/1.88, 1.08/1.02, 2.14/2.16, 2.34/2.37, 0.66/0.73,
0.75/0.77, 1.67/1.69, 1.86/1.90, 1.87/1.90). Neither stated falsifier occurred: the ring
row is near p + 1 on the translation and on the rotation.

**N-4 (1 - n.n_e of order 2p): HELD for the ring normal on every test** (1.76, 2.38, 2.27,
2.14, 2.06 against 2p = 1.92, 2.58, 2.44, 2.34, 2.26) and trivially for Youngs; the measure
is exactly half the square of the L2 of |n - n_e| (verified on every case to 5e-13
relative), so its order is twice the L2 order, which lags the L1 order by 0.05 to 0.1 for
the ring normal. A first-order ring normal reads as second order in this measure: at N 128
on the rotation L1(theta) = 1.13e-2 and L1(1 - n.n_e) = 1.10e-4, orders 1.17 and 2.14.

**N-5 (serial against 8 ranks at N 32): HELD for the exact and the ring normal** (E_shape,
L1(theta), L1(d), L1(alphaf) agree to 1e-9 relative or better on the translation and the
rotation), **FAILED for the Youngs normal on the rotation at OF Co 0.2** (L1(alphaf) differs
by 2e-3 relative, min alpha -1.6e-6 on 8 ranks against -4e-11 serial; the other norms to
6e-7).

Bounds and mass: the exact and the ring rows are bounded to 1e-10 and conservative to
1e-11 % on every case; the Youngs row loses the bound on the rotation (-1.6e-6 at N 32,
-7.3e-5 at N 64, abort at N 128).

**What is demonstrated, what is a limitation, what is open.** Demonstrated on the sphere
translation (axial and oblique) and rotation, N 32 to 128, OF Co 0.2 and 1, pieces path: a
first-order normal does not give a first-order scheme; with the plane positioned by CCS from
alpha_k the in-cell position and the shape error converge one order above the angle of the
normal, for the ring normal (angle 1.0 to 1.3, shape 1.9 to 2.3) and for a zero-order normal
(Youngs: angle 0 to 0.2, shape 0.7 to 1.1). Limitation: the "biased first-order normal" row
did not exist as specified (Youngs is zero order), so the bias hypothesis itself, shape
order near p for a biased normal of order p >= 1, is not tested; no row showed a shape order
near p. The per-step flux budget of the pre-registration is wrong by one power of h. Open: a
normal that is first order AND biased (for example the ring normal with a prescribed O(h)
tilt of fixed sign) to test the bias clause; the mechanism by which per-step flux errors of
order p sum to a shape error of order p + 1 (the next reconstruction re-positions the plane
from alpha_k; not measured step by step here); the one-plane path's rows 4-6 (the KISS
round 3 gate jobs 56250533-34 completed, their verdict is the KISS agent's to report; not
run here).

The 1 - n.n_e columns of the finished studies (derived exactly from the stored L2 and Linf
of |n - n_e|; no case needed a re-run): the standing matrix (405 cases; 90 diverged or
crashed DAVOF cases report no normal error at all), the corrector tables (34 cases), the
sphere ladders (38 cases); in the preprint tables on branch `davof-normal-order-article`
(040aad92, 1d423d24).

## The bias test (2026-11-13, leia branch `davof-normal-order`)

Tomislav's instruction of 2026-11-13: "Do the bias test." The normal-order study (above)
measured that the plane positioned by CCS from alpha_k has a position and a shape error one
order above the angle of the normal, for the ring normal (angle order 1.0 to 1.3, shape 1.9
to 2.3), and that the Youngs normal, which was meant as a biased first-order normal, is a
zero-order normal. The hypothesis that the shape order is set by the BIAS of the normal,
not by its order, is therefore untested. This test gives the bias a known size and sign.

### The model (src/davof/normalModels/ringTiltNormal.{H,C}, key normalModel ringTilt)

The ring area normal n_k rotated by theta_b = c h / R, h = V_k^(1/3), R = 0.15 the radius
of the exact sphere, about the axis a = e_z x n_k / |e_z x n_k| (where |e_z x n_k| < 1e-6,
e_x x n_k): n' = n_k cos(s theta_b) + (a x n_k) sin(s theta_b). The rotation turns n_k
away from e_z in the plane of e_z and n_k, and keeps |n_k|. biased: s = +1 in every cell:
one sense relative to the fixed direction e_z everywhere, so the tilt is first order and
has a fixed sign. alternating: s = +1 or -1 by the parity of floor(x/h + 1/2) +
floor(y/h + 1/2) + floor(z/h + 1/2) of the cell centre, the same magnitude and order, no
bias (the parity is geometric, independent of the cell numbering and the decomposition).
c is chosen so that theta_b at N 32 is comparable to the ring normal's own L1 angle error
of 0.041 to 0.058 (the oblique translation and the rotation, measured above):
c = 0.125, 0.25, 0.5 give theta_b = 0.026, 0.052, 0.104 at N 32 (c = 0.5, 1, 2 would give
0.10, 0.21, 0.42, two to eight times the ring error; the instruction's two clauses
conflict and the comparability clause is followed).

### The tests

The oblique translation (U = (0.3 0.5 0.8124), over 0.3) and the rotation (one turn), N 32,
64, 128, OpenFOAM Co 0.2 and 1, 8 ranks, the pieces path of the normal-order study, with
the position diagnostic. Configs `config/davof/nor_bias{Ob,Rot}_N*.yaml` and
`nor_balance_N64.yaml` (`workflow/scripts/nor_make_bias_configs.py`). Same measures as the
normal-order study: the angle, 1 - n.n_e, the position error, the per-step flux error (the
translation), E_shape, the orders over N, the wall clock at 8 ranks.

### Predictions (stated before any run; p = the L1 order of the angle error)

- **B-1 (the tilt):** the L1 angle error of the biased and the alternating rows is
  theta_b + (the ring error): order about 1 (0.7 to 1.3), and the L1 angle of the tilted
  rows exceeds the plain ring normal's by at least the factor 1.5 at N 32 for c = 0.25 and
  0.5.
- **B-2 (the biased row):** E_shape of the biased tilt converges near p = 1: the order lies
  in [0.7, 1.4] for c = 0.25 and 0.5.
- **B-3 (the controls):** E_shape of the alternating tilt and of the plain ring normal
  converges near p + 1: the order lies in [1.6, 2.5].
- **B-4:** the order of L1(1 - n.n_e) is 2p within 0.3 for all rows.

Falsifier: the biased row also converges near 2 (order above 1.6 for c = 0.25 and 0.5). The
mechanism is then that CCS re-positioning removes a fixed-sign first-order tilt too, and
the shape order is p + 1 regardless of bias; the bias clause of the argument is then
falsified, and the tilt sets the constant, not the order.

### The step-resolved balance (oblique translation, N 64, OpenFOAM Co 0.2, c = 0.25)

Per step, for the exact, the ring, the biased and the alternating normal: (i) the per-step
flux error (L1 over the faces of the interface cells, relative to V_k of the upwind cell,
and its signed mean), (ii) the position error of the plane at the CCS start (the closed-form
position along the used normal, clamped into the vertex range) and of the final plane, both
the distance of the plane centroid from the exact sphere, L1 over the interface cells, and
the mean CCS shift |p_final - p_start|/h (key positionDiag). It shows where the order is
gained: the position error before CCS carries the tilt (order p), and after CCS the tilt
only through the pivot (order p + 1) if the position order of the normal-order study is the
mechanism. Not predicted numerically; the ratio before/after is reported.

### MEASURED, the bias test (2026-11-13/14; 78 cases on 8 ranks)

Jobs: build 56305790; orchestrators 56305791 (balance), 56305792-94 (oblique translation N 32/64/128), 56305795-97 (rotation N 32/64/128); all completed. Tables: `/work/scratch/tm83tomy/davof/agents/normal/csv/bias/{nor_cases.csv,nor_orders.csv,nor_tables.md}`; preprint table body `table_bias.tex` on development.

Orders over N = 32, 64, 128 (p: L1 angle; d: position; q: E_shape; omd: L1(1 - n.n_e)); the plain ring normal for reference (oblique 0.2: p 1.29, d 2.16, q 2.14, omd 2.38; oblique 1: 1.22, 2.37, 2.34, 2.27; rotation 0.2: 1.17, 1.90, 1.86, 2.14; rotation 1: 1.13, 1.90, 1.87, 2.06).

| test, OF Co | tilt | c | p | d | q | omd | L1(theta) N 32 / ring | E_shape N 128 / ring |
|---|---|---|---|---|---|---|---|---|
| oblique 0.2 | biased | 0.125 / 0.25 / 0.5 | 1.18 / 1.13 / 1.12 | 2.17 / 2.13 / 2.05 | 2.09 / 2.06 / 2.02 | 2.35 / 2.29 / 2.25 | 1.13 / 1.52 / 2.68 | 0.95 / 1.16 / 1.95 |
| | alternating | | 1.18 / 1.09 / 1.03 | 2.14 / 2.11 / 2.08 | 2.12 / 2.10 / 2.06 | 2.32 / 2.21 / 2.11 | 1.11 / 1.37 / 2.19 | 1.03 / 1.08 / 1.26 |
| oblique 1 | biased | | 1.11 / 1.08 / 1.08 | 2.27 / 2.17 / 2.07 | 2.22 / 2.12 / 2.04 | 2.21 / 2.19 / 2.17 | 1.17 / 1.63 / 2.85 | 1.29 / 1.83 / 3.12 |
| | alternating | | 1.12 / 1.06 / 1.02 | 2.34 / 2.29 / 2.18 | 2.31 / 2.26 / 2.14 | 2.22 / 2.14 / 2.06 | 1.15 / 1.51 / 2.51 | 1.05 / 1.15 / 1.46 |
| rotation 0.2 | biased | | 1.21 / 1.18 / 1.15 | 1.92 / 1.94 / 1.99 | 1.88 / 1.90 / 1.94 | 2.30 / 2.30 / 2.27 | 1.27 / 1.64 / 2.58 | 1.01 / 1.02 / 1.05 |
| | alternating | | 1.12 / 1.06 / 1.03 | 1.90 / 1.91 / 1.92 | 1.86 / 1.86 / 1.87 | 2.14 / 2.11 / 2.08 | 1.08 / 1.29 / 2.00 | 1.00 / 1.00 / 1.02 |
| rotation 1 | biased | | 1.15 / 1.13 / 1.10 | 1.93 / 1.97 / 1.99 | 1.89 / 1.91 / 1.96 | 2.21 / 2.22 / 2.19 | 1.29 / 1.71 / 2.80 | 1.01 / 1.02 / 1.07 |
| | alternating | | 1.07 / 1.03 / 1.01 | 1.90 / 1.91 / 1.93 | 1.87 / 1.87 / 1.89 | 2.06 / 2.05 / 2.04 | 1.11 / 1.41 / 2.33 | 1.00 / 1.01 / 1.03 |

Wall clock at N 128, 8 ranks, shared nodes: oblique 0.2 1180-1490 s, oblique 1 244-306 s, rotation 0.2 8700-10700 s, rotation 1 1680-1880 s (the position diagnostic on in all of these rows; every row within 20 % of the plain ring normal's 1059 / 282 / 9327 / 2034 s except the biased c = 0.5 rows at Co 0.2 on the translation, 1492 s).

**B-1 (the tilt): held in order, held in the factor for the biased tilt.** The angle order of the tilted rows is 1.01 to 1.21. The L1 angle at N 32 exceeds the ring normal's by 1.5 to 1.7 (biased c = 0.25) and 2.6 to 2.9 (biased c = 0.5), as registered; the alternating rows reach 1.29 to 1.51 at c = 0.25 (below the registered 1.5 on three of four tests) and 2.0 to 2.5 at c = 0.5.
**B-2 (biased shape order near 1): FALSIFIED.** The biased shape order is 2.02 to 2.22 on the translation and 1.88 to 1.96 on the rotation for every c; the falsifier (above 1.6 for c = 0.25 and 0.5) occurred in all eight cells.
**B-3 (controls near 2): held.** Alternating 2.06 to 2.31 (translation) and 1.86 to 1.89 (rotation); ring 2.14, 2.34, 1.86, 1.87; all in [1.6, 2.5].
**B-4 (1 - n.n_e of order 2p): held** within 0.3 for every row (omd 2.04 to 2.35 against 2p of 2.02 to 2.42).

The mechanism implied by the falsifier: CCS from alpha_k does not remove the fixed-sign tilt either, and neither is the shape order p + 1 produced by a cancellation of per-step errors of alternating sign. The step-resolved balance (nor_balance_N64: oblique translation, N 64, OF Co 0.2, c = 0.25, 155 steps, medians over the steps): per-step flux error relative to V_k 2.5e-4 (exact), 4.0e-4 (ring), 4.1e-4 (biased), 4.6e-4 (alternating); position error L1(d) of the plane at the CCS start (the closed-form position along the used normal, clamped into the vertex range) 1.33e-4, 2.08e-4, 2.52e-4, 2.98e-4 and after CCS 1.48e-4, 1.74e-4, 1.95e-4, 2.04e-4; mean CCS shift |p_final - p_start|/h 0.008, 0.009, 0.014, 0.012. The ratio before/after CCS is 0.9 (exact; CCS makes it slightly worse), 1.2, 1.3 and 1.5: CCS reduces the position error of the tilted planes by at most a third and does not change its order (over N: before 2.02 to 2.14, after 2.03 to 2.29 on the translation; 1.88 to 1.95 before and 1.90 to 2.01 after on the rotation). The order p + 1 of the position error is already in the plane at the CCS start; it is the lever arm of the normal (the position error O(h theta_k)), not the re-positioning, that carries it. The signed mean of the per-step flux error is 1.3e-4 (exact), 1.4e-4 (ring), 1.4e-4 (alternating) and -1.2e-6 (biased): a near-zero signed mean for the biased tilt, a change of the per-step flux bias by two orders of magnitude, which the shape order does not follow; its mechanism was not located.

**Demonstrated** (sphere, oblique translation and rotation, N 32 to 128, OF Co 0.2 and 1, pieces path): a first-order normal error of fixed sign gives the same shape order as an alternating one and as the plain ring normal, and a larger constant (up to 3.1 times at c = 0.5); the position error has order p + 1 before CCS. **Limitation:** one form of bias (a tilt toward a fixed direction, not correlated with the curvature of the interface), one geometry, the pieces path only. **Open:** a tilt correlated with the interface curvature (toward the centre of the sphere); the one-plane path; why the signed mean of the per-step flux error vanishes for the biased tilt while the shape error does not follow.


## Signed flux error of the biased tilt (2026-11-14, leia branch `davof-normal-order`)

The bias test measured the signed mean of the per-step flux error (sweptExactFluxError, oblique translation, N 64, OpenFOAM Co 0.2, c = 0.25, 155 steps) as 1.3e-4 (exact normal), 1.4e-4 (ring), 1.4e-4 (alternating tilt) and -1.2e-6 (biased tilt), while E_shape of the biased row is not smaller than the others (order about 2, larger constant). This section registers three hypotheses for that observation before any new run.

### The definition found in the code (located before any run; no run needed)

`davofSweptAdvection::exactFluxError` (src/davof/advection/davofSweptAdvectionSurfaces.C, from the line `void Foam::davofSweptAdvection::exactFluxError`; called from davofSweptAdvectionTemplates.C after `curvedFlux`, before the bounding). Faces: every face of every interface cell (`surf_->interfaceLabels()`), internal and processor faces, the processor faces flagged on both sides and weighted 1/2 per side; empty patch faces skipped. For each face: pv = U dt . S_f; the upwind cell `up` is owner for pv >= 0 and neighbour for pv < 0 (owner on boundary faces); vEx = +/- the exact liquid volume of the prism f - U dt .. f (sphere-tetrahedron overlap, exactSphere), sign + for pv >= 0, - for pv < 0, zero when |pv| < 1e-10 |U dt||S_f|; diff = (dVf_f - vEx)/V_up, with dVf_f the flux volume through the face along S_f (the owner-to-neighbour orientation), after the flux and before the bounding. L1 = sum w |diff| / sum w, signed mean = sum w diff / sum w, Linf = max |diff|; one line per step. The sign of diff is therefore the owner-to-neighbour orientation of S_f. For a flow with all components of U positive on the blockMesh (the oblique translation) every face has pv > 0, so this orientation is the flow direction.

### Hypotheses (pre-registered)

H1 (coincidence of two terms): the signed mean is s(c) = s_0 + s_1 c, with s_0 ~ 1.3e-4 the curvature (sagitta) term common to all normals and s_1 c the tilt's contribution, which for the biased tilt has a fixed sign and happens to cancel s_0 near c = 0.25. Prediction: s is linear in c over c in {-0.5, -0.25, -0.125, 0, 0.125, 0.25, 0.5} (negative c = opposite tilt sense), zero crossing near c = 0.25, and s(-0.25) ~ 2 s_0. For the alternating tilt s is independent of c.

H2 (the face mean is not the quantity that sets the shape error): the shape error is driven by the per-cell NET flux error e_k = (1/V_k) sum over the cell's faces of the outward-signed (dVf - dVf_exact), not by a face mean whose sign depends on the face orientation convention. Prediction: L1 over interface cells of |e_k| per step is of the same size for biased and ring rows (ratio within the ratio of their E_shape, 1.0 to 1.3 at c = 0.25) and does NOT vanish for the biased row; the sum of e_k V_k over all cells is zero to round-off (conservation).

H3 (convention artefact): the signed mean changes when the sign convention changes (owner-to-neighbour orientation vs flux-direction sign vs outward from the wetter cell).

Falsifiers: s not linear in c or no zero crossing near 0.25 (H1 false); L1|e_k| of the biased row much smaller than the ring's (H2 false; then the shape error has another source and it is reported).

### The runs and the measures

The balance case (oblique translation U = (0.3 0.5 0.8124), OpenFOAM Co 0.2, 8 ranks, 155 steps): exact, ring, biased c in {-0.5, -0.25, -0.125, 0.125, 0.25, 0.5}, alternating c in {0.125, 0.25, 0.5} at N 64; exact, ring, biased c = +-0.25, alternating c = 0.25 at N 32 and N 128. Print-only diagnostics (key sweptFluxErrorDiag, default off, the transport unchanged): per step the signed mean under the three conventions (as implemented; by the sign of the face flux; outward from the wetter cell), the per-cell net error e_k (L1, L2, Linf and signed mean over the interface cells, sum of e_k V_k over all cells touched), and the signed mean per face orientation group (x-, y-, z-normal faces). Statement of an identity found while writing this registration: the sign by the face flux and the sign outward from the upwind cell are the same quantity (the flux leaves the upwind cell), so the three conventions of H3 are implemented as: orientation of S_f (as implemented), flow direction (equal to outward from the upwind cell), outward from the wetter cell (the larger alpha of the two cells at t^n). The sum of e_k V_k over cells is a telescoping sum over the internal faces; it tests the boundary and processor faces and the conservation of the computed dVf, not the exact flux.


### MEASURED, the signed flux error (2026-11-14; 21 cases on 8 ranks, all completed)

Build 56325885 (24b4a3bc; the first build 56314737 failed on a `FixedList(int)` initialisation, fixed in 24b4a3bc); orchestrators 56326064-66 (N 64 base, bia, alt), 56326067-69 (N 32), 56326071-73 (N 128). Tables: `/work/scratch/tm83tomy/davof/agents/normal/csv/sgn/{sgn_cases.csv,sgn_tables.md}` (`workflow/scripts/sgn_tables.py`). The dictionaries of every existing config render byte-identical (checked on nor_balance_N64: all files equal except case_params.json: the commit stamp, the run date and the new default token TPF_FLUX_ERROR_DIAG "false"). Means over the steps (77, 155, 310 steps at N 32, 64, 128). s: the signed mean of (dVf - dVf_exact)/V_up over the faces of the interface cells, orientation convention as implemented.

H1, N 64, s(c) of the biased tilt (c = 0: the ring row): c = -0.5, -0.25, -0.125, 0, 0.125, 0.25, 0.5 gives s = 4.183e-4, 2.804e-4, 2.109e-4, 1.408e-4, 7.05e-5, -7.5e-8, -1.408e-4. Least-squares line s = s_0 + s_1 c: s_0 = 1.400e-4, s_1 = -5.595e-4, R^2 = 1.0000, max residual 1.5e-6, zero crossing c* = 0.250; s(-0.25) = 2.804e-4 against 2 s_0 = 2.800e-4. Alternating, independent of c: 1.405e-4, 1.400e-4, 1.388e-4 (c = 0.125, 0.25, 0.5); exact 1.276e-4. The same at N 32 and 128 (two points c = +-0.25): s(-0.25), s(0.25) = 5.556e-4, 9.4e-6 (N 32); 1.428e-4, 1.6e-6 (N 128); the zero crossing of the two-point line is c* = 0.259 (N 32) and 0.256 (N 128); s_1 = -1.09e-3, -5.6e-4, -2.8e-4 at N 32, 64, 128 (order 1), s_0 follows (2.4e-4 exact / 2.9e-4 ring, 1.4e-4, 7.2e-5 ring).

H3, signed mean under three conventions, N 64 (orientation / flow / wetter): exact 1.276e-4 / 1.305e-4 / 1.293e-4; ring 1.408e-4 / 1.445e-4 / 1.225e-4; alternating c = 0.25 1.400e-4 / 1.433e-4 / 1.222e-4; biased c = 0.25 -7.5e-8 / 3.6e-6 / -4.06e-5; biased c = -0.25 2.804e-4 / 2.841e-4 / 2.830e-4. N 32 biased c = 0.25: 9.4e-6 / 6.2e-6 / -7.2e-5; N 128: 1.6e-6 / 2.0e-6 / -1.63e-5. Orientation and flow differ by at most 3.7e-6.

H2, per-cell net error e_k (L1 over the interface cells, mean over the steps) and E_shape at the end, ratios to the ring row:

| N | row | c | e_k L1 | e_k L1 / ring | E_shape / ring | e_k L1 / face L1 | wall s |
|---|---|---|---|---|---|---|---|
| 64 | exact | 0 | 9.32e-4 | 0.68 | 0.57 | 3.46 | 71 |
| 64 | ring | 0 | 1.372e-3 | 1.00 | 1.00 | 3.31 | 71 |
| 64 | biased | -0.5 | 3.423e-3 | 2.50 | 2.60 | 3.08 | 74 |
| 64 | biased | -0.25 | 2.308e-3 | 1.68 | 1.70 | 3.20 | 80 |
| 64 | biased | -0.125 | 1.771e-3 | 1.29 | 1.29 | 3.26 | 88 |
| 64 | biased | 0.125 | 1.212e-3 | 0.88 | 0.91 | 3.29 | 93 |
| 64 | biased | 0.25 | 1.310e-3 | 0.95 | 1.06 | 3.03 | 70 |
| 64 | biased | 0.5 | 1.979e-3 | 1.44 | 1.72 | 2.69 | 83 |
| 64 | alternating | 0.125 / 0.25 / 0.5 | 1.459e-3 / 1.628e-3 / 2.083e-3 | 1.06 / 1.19 / 1.52 | 1.02 / 1.06 / 1.19 | 3.39 / 3.44 / 3.36 | 81 / 90 / 86 |
| 32 | ring / biased 0.25 / biased -0.25 | | 3.131e-3 / 2.974e-3 / 4.840e-3 | 1.00 / 0.95 / 1.55 | 1.00 / 1.05 / 1.65 | | 17 / 17 / 18 |
| 128 | ring / biased 0.25 / biased -0.25 | | 5.993e-4 / 6.022e-4 / 1.095e-3 | 1.00 / 1.00 / 1.83 | 1.00 / 1.16 / 1.89 | | 1029 / 1387 / 1392 |

Position error L1 (median over the steps, before / after CCS) at N 64: exact 1.33e-4 / 1.48e-4, ring 2.08e-4 / 1.74e-4, biased c = 0.25 2.52e-4 / 1.95e-4, biased c = -0.25 3.46e-4 / 2.86e-4, alternating c = 0.25 2.98e-4 / 2.04e-4. Sum of e_k V_k over all cells touched: at most 3.4e-16 (N 64), 7.2e-16 (N 128), against sum of |dVf - dVf_exact| over the faces of 7e-6 to 3e-5 (N 64): round-off.

Split by face normal at N 64 (orientation convention): x, y, z for exact 9.3e-5, 1.29e-4, 1.61e-4; ring 9.7e-5, 1.42e-4, 1.83e-4; biased c = 0.25 9.3e-5, 1.35e-4, -2.28e-4; biased c = -0.25 1.02e-4, 1.47e-4, 5.93e-4; biased c = 0.5 8.7e-5, 1.28e-4, -6.38e-4; alternating c = 0.25 9.7e-5, 1.41e-4, 1.82e-4. The x and y means change by at most 15 % over c in [-0.5, 0.5]; the z mean (the faces whose normal is e_z, the tilt axis direction) carries the change.

Orders over N 32/64/128 (least squares of log y against log N): exact: s 0.93, e_k L1 0.99, face L1 0.93, E_shape 1.92; ring 1.00, 1.19, 1.09, 2.14; biased c = -0.25 0.98, 1.07, 1.00, 2.04; biased c = 0.25 1.30 (s changes sign: +, -, +; |s| from 9.4e-6 to 7.5e-8 to 1.6e-6), 1.15, 1.08, 2.06; alternating c = 0.25 0.99, 1.09, 1.05, 2.10.

**H1: held.** s is linear in c at N 64 (R^2 = 1.0000, residual 1.5e-6), the zero crossing is c* = 0.250, s(-0.25) = 2.804e-4 against 2 s_0 = 2.800e-4, the alternating tilt is independent of c within 1.2 %. The zero crossing stays at 0.256 to 0.259 from the two-point lines at N 32 and 128: s_0 and s_1 both have order 1 and their ratio is N-independent.
**H2: held.** L1 of e_k of the biased c = 0.25 row is 0.95, 0.95, 1.00 of the ring's (N 32, 64, 128), E_shape 1.05, 1.06, 1.16; the ratio of e_k L1 follows the ratio of E_shape over c (2.50 / 2.60, 1.68 / 1.70, 1.29 / 1.29, 0.88 / 0.91 for c = -0.5, -0.25, -0.125, 0.125), and it does not vanish for the biased row. The ratio at c = 0.25 is 0.95 to 1.00, below the registered 1.0 to 1.3 by 0.05 at most. The sum of e_k V_k is zero to round-off; it is an identity of the internal faces (each counts once with each sign), so it tests the boundary and processor faces only.
**H3: held as stated, not an explanation.** The wetter-cell convention changes the signed mean of the biased c = 0.25 row from -7.5e-8 to -4.06e-5 (N 64) and of the ring row from 1.41e-4 to 1.23e-4; orientation and flow differ by at most 3.7e-6. In no convention does the biased c = 0.25 row reach the 1.2e-4 to 1.4e-4 of the other rows (|s| at most 7.2e-5 at N 32, 4.1e-5 at N 64, 1.6e-5 at N 128).

Unexpected: (i) the orientation and the flow conventions differ although every face of the oblique translation has U dt . S_f > 0 on the blockMesh: the processor faces are seen with S_f pointing outward from the local cell, so the sign of the orientation convention on the processor faces (weight 1/2 per side) depends on the decomposition. (ii) The N 64 biased c = 0.25 value is -7.5e-8 here (mean over 155 steps) against the -1.2e-6 of the bias test; the earlier value is not reproduced by this mean (the earlier statistic of the step values was not recorded in the table).

**Demonstrated** (sphere, oblique translation, N 32 to 128, OF Co 0.2, pieces path): the signed per-step flux error of the face mean is the sum of a tilt-independent term s_0 (the sagitta term, 1.4e-4 at N 64) and a term s_1 c linear in the tilt, carried by the faces normal to the tilt axis; the biased tilt of c = 0.25 cancels s_0 to 5e-4 of its value at N 64 and to 0.02 to 0.03 of it at N 32 and 128; the per-cell net error follows the shape error over all c and does not vanish at the zero crossing. **Limitation:** one geometry, one tilt axis (e_z), one Co; the cancellation at c = 0.25 is a coincidence of the constants s_0 and s_1 (their ratio), not a property of the method; the means over the steps hide the step-to-step variation (L1 over the interface cells of e_k, not its sign, relates to E_shape). **Open:** why s_1/s_0 = 0.25 here (the tilt c is registered in units of h/R); the processor-face sign of the orientation convention; the same decomposition of the signed mean on the rotation.
