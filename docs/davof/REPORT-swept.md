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
N = 32 and 64 at
the end is 1.98; the rotation's Eshape grows by 26 % / 47 % from T/4 to T.
