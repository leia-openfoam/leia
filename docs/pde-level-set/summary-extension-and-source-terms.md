# Velocity extension and gradient source terms in leia: what we tried, and the roadmap

Summary of 2026-11-13, for D. Bothe and K. Soga and for the supervisor. Numbers and links:
`STATUS.md` (sections 9 to 11.30), `docs/plan-velocity-extension-verification.md`, the deck
`docs/pde-level-set/pde-level-set-presentation/`, and the knowledge base
(`docs/knowledge-base/hubs/gradient-control.md`).

## What we tried

### July 2026: the velocity-extension line (Eulerian FVM, `leiaLevelSetFoam`)

- Eight extension models in one runtime-selectable family: `anisotropicDiffusion`, `pseudoTime`,
  `steadyUpwind`, `steadyUpwindLinear`, `closestPoint`, `meshWave`, and later `haloLimited` and
  `projectionBlend`.
- A static test at t = 0 measured only whether the extended velocity is constant along the normals.
  `closestPoint` was named the static winner, at about second order.
- On the reversed 2D vortex the extensions were judged "dominated": the best accuracy needed 12 to 27
  times the cost. That comparison was against the semi-Lagrangian method, read at the end of a
  reversed flow, and ran in parallel before two seam defects were fixed (2026-08-26).

### August 2026: the SDPLS strain source (Eulerian FVM)

- The source psi * a with the normal strain a = n.grad(u).n keeps |grad psi| near one in kinematic
  flows: band gradient order +0.74 against -0.26 without it.
- Under two-phase coupling it fails: it reads the strain of the spurious current and amplifies the
  mesh-locked current 260 times; with exact curvature the current stays zero. The seed is the
  curvature error.
- Two variants were withdrawn: `beta` (its target is structurally wrong) and the material form `Rdiv`
  (diverges kinematically).

### September 2026: the dossiers in the semi-Lagrangian line (first gradient-control campaign)

- Six candidates from the two dossiers: the soft-wall source alone (S1), the halo-limited extension
  alone (HL0), the extension with a linear law in q or in z = q^2 (HL1q, HL1z), the extension with
  the soft wall (HL2), and the closest-point extension (FP0).
- Every candidate failed the 2D method gate. On the shear arm the shape error grew 14 to 410 times;
  HL1q and HL1z destroyed the stationary droplet.
- The technical report found the main cause in the source's own discretisation: a centred gradient of
  q lets the grid-scale mode grow. The halo-limited extension at R = h moves the normal strain into a
  shell instead of removing it.
- FP0 depended on the decomposition, through an unfixed halo surrogate in its implementation.

### October 2026: the combined source F = a + lambda (1 - q) (semi-Lagrangian source step)

- The centred q fails even at zero velocity; a one-sided, upwind-weighted q (`upwindFit`) passes.
- The combined source and both single terms then lose the shape order on the kinematic arms
  (0.88 / 1.24 / 0.26 against 3.25). The cause, measured: the strain term depends on psi through the
  normal n = grad(psi)/|grad psi|; with a fixed normal the shape stays within 1.2 times the baseline.
- Not run: the next form, a strain term with the normal taken from the zero set.

### October 2026: the static test of every Eulerian source law (zero velocity)

- The first run was void: the psi solver did zero iterations below its tolerance in 68 % of the
  solves. Fixed with `minIter 1` and re-run.
- Result: every law (Soga's linear law, the z laws, the cubic laws, the soft wall) lets grid-scale
  extrema of psi grow at the law's own rate. The defect is the estimate of q, not the law.

### Late October 2026: D. Bothe's algebraic rescaling and blended projection extension (Eulerian FVM)

- The second-order distance d-hat = 2 psi / (|grad psi| + sqrt(|grad psi|^2 - 2 psi d_nn psi)),
  applied once, is second order with the zero set and the volume unchanged (exact for a quadratic
  psi). Applied every step it diverges with the centred gradient.
- The blended projection extension `projectionBlend` halved the band gradient error at T/2 on the
  reversed vortex and the 3D shear, but made the shape and the volume at T worse. On the unresolved 3D
  deformation it improved every metric.

### November 12-13, 2026: re-verification from the bottom (this direction)

- **An audit of every earlier verdict.** No extension had ever been compared with the exact extension
  u(Pi(x)). Most losses were read at T of a reversed flow, where a transport with a prescribed
  velocity cancels its own error and an extension that depends on psi breaks the cancellation. Read
  at T/2, the extensions won in every data set found. The advected extension verdicts are re-opened.
- **Two code defects found by the audit.** (1) Every extension except `none` and `haloLimited` rebuilt
  the far-field flux too (`interpolate(Uext) & Sf` on every face). Repaired with a correction form,
  phi + interpolate(Uext - U) & Sf: an inert extension now equals no extension bit for bit. Measured
  size of the defect: shape 0.1 to 0.2 %, volume 1 to 4 % at T. Real, but small. (2) `closestPoint`
  and `meshWave` write nothing beyond their band, so the velocity jumps at the band edge.
- **The extension band assumed a distance.** The models' own band and fade used |psi| <= (n+1) h.
  The token is now rendered, and the topological `layer` band is the setting from step B on. The
  cases' narrow band is `signChange` throughout.
- **The Eulerian transport passes where it can be checked exactly** (one-way translation, where psi
  must stay a signed distance): band gradient error order 2.3 to 3.9, shape about 2.2, volume 2 to 4.
  The serial/np 4 difference is solver convergence noise (it falls by 2e4 with a tighter tolerance),
  not a seam defect.
- **Nine extensions against the exact u(Pi(x)) at one instant** (circle in a shear flow, N 32 to
  256): `closestPoint` is the only second-order extension (360 times better than no extension);
  `meshWave` is first order; `pseudoTime`, `anisotropicDiffusion` and `haloLimited` are barely better
  than no extension in the band; `steadyUpwind` does not converge; `steadyUpwindLinear` is worse than
  no extension. `projectionBlend` is second order inside 1.5 h, where its blend weight is zero.
- **The face flux, which the transport actually uses,** is second order inside the band and first
  order only at the band edge: the jump, not the interpolation.
- **Two of my readings were corrected by the supervisor.** The parallel-surface curvature inverse does
  not read psi as a distance (its offset is d-hat from the cell's quadratic fit). interFoam is no
  reference; the goal is an accurate, resolved, adaptive solver, judged against exact solutions and
  observed orders.

## What we conclude so far

- The extension was never given a fair test before November 12; the source terms were, and their
  defect is the estimate of q.
- The Eulerian transport itself is sound.
- `closestPoint` is the extension to build on, once its band edge is repaired.
- In a two-phase run only the scalar face flux and the cell velocity come from Navier-Stokes. An
  extended face flux must therefore be phi_f + [u(Pi(x_f)) - u(x_f)] . S_f, with both samples from
  one reconstruction of the cell velocity; it keeps the Navier-Stokes flux exactly at the interface.

## Roadmap

Each step has an exact reference, a pre-registered prediction, at least three resolutions with the
observed order, a serial/np 4 check, and runs on Lichtenberg.

1. **Step A0, finish the repairs.** The band edge of `closestPoint` and `meshWave` (extend to the end
   of the fade, or fade inside the band); the `closestPoint` far-field switch and a true Newton
   descent; the fixable settings of `pseudoTime`, `anisotropicDiffusion` and `steadyUpwind`. Each
   gated by bit identity where it claims to be inert.
2. **Step A, the extension operator against u(Pi(x)), completed.** Inputs with |grad psi| != 1 (d g,
   r^2 - R^2, d-hat), an off-lattice ellipse, 3D, perturbed hexahedra and polyhedra. New candidates:
   the oracle extension, the Newton point of Bothe and Soga sampled by cell location, Bothe's capped
   and blended form, the face-flux form, and a normal-component-only extension.
3. **Step 0.5, the q range (waits for the supervisor).** A static gate that measures what the
   production geometry operators need, against the same operators on psi = d. Prediction: the level of
   q does not matter, its variation along the interface does.
4. **Step B, transport with the exact extension (oracle).** Flows with a closed-form interface at
   every instant: translation, rotation, and linear flows that map a circle to an ellipse (strain,
   no reversal). Reference: the signed distance to the moved interface. This isolates the transport
   of a non-solenoidal velocity.
5. **Step C, transport with the computed extension.** The same cases with the best two or three
   extensions from step A; the cost of the loop psi -> foot point -> velocity -> psi; the time step;
   the decomposition. A reversed flow enters only here, read at T/2 against a marker reference.
6. **Step D, the source term.** Alone at zero velocity against the exact logistic law z' = mu z (1 - z),
   with the one-sided q ported to the Eulerian source; then with the oracle extension, where it only
   has to remove the initial error of q and the discretisation drift.
7. **Step E, assembly and coupling.** The best extension with the best source on the kinematic cases;
   then the stationary, translating and oscillating droplets against their exact solutions (Laplace
   jump, zero velocity, the Lamb frequency and damping), first with exact curvature; then on statically
   refined meshes, because the goal is an adaptive resolved solver. 3D and polyhedra only on a pass.
8. **Open questions for D. Bothe and K. Soga.** The extension variants a) to d) of the notes of
   2026-10-29; what sets beta; whether the velocity-extension theorem holds with the blend and the cap;
   which range of |grad psi| the analysis needs; and, in the blend, which term carries delta.
