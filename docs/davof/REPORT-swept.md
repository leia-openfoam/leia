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
