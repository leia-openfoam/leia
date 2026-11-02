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
