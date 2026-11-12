# Plan: verified building blocks of a PDE level-set method on the Eulerian FVM line

Velocity extension first, the gradient source second, then the assembly. The model is the one of
D. Bothe and K. Soga. Written 2026-11-12. Style: ASD-STE100. Status: the decisions of section 0 are
taken; nothing of sections 4 to 10 is implemented yet. Record: STATUS.md 11.29.

---

## 0. Decisions taken (supervisor, 2026-11-12)

| # | question | decision |
|---|---|---|
| 1 | the goal | assemble a PDE level-set method in OpenFOAM that is usable for multiphase flow; "better than interFoam" is enough, Basilisk accuracy is not the target |
| 2 | the line | the Eulerian finite-volume line only (`leiaLevelSetFoam`, `leiaLevelSetTwoPhaseFoam`); the semi-Lagrangian line is out of scope for this plan |
| 3 | the order | split the model into sub-problems that each have an exact reference; the velocity extension first (its verification test, then the existing models, then the best one), the source terms after it, then the assembly |
| 4 | the attitude to old verdicts | every earlier verdict is re-read critically: the work moved fast, and a bug, a setting or a bottleneck may have decided a verdict that the method did not deserve |
| 5 | the record for the supervisor | a reveal.js deck in the style of a PhD student's progress report, updated with each sub-problem (section 11); the numbers live in STATUS.md and the curated CSVs |
| 6 | the reading of D. Bothe's notes (deck, stack 2) | right (supervisor, 2026-11-12, second round) |
| 7 | the advected velocity-extension verdicts | re-opened and analysed (section 12) |
| 8 | the target range of q | not yet decided; the supervisor asked for more material: section 4b and the deck's q stack |
| 9 | start | step 0 and step A0 start on Lichtenberg |

---

## 1. The target model

From the two dossiers (`agent-input/`, not in git; cite by title) and from D. Bothe's meeting notes
("Dieter-Besprechung-LevelSet", 2026-10-29):

$$\partial_t\phi + \mathbf v^{\rm ext}\cdot\nabla\phi = \beta\,\phi\,g(\lvert\nabla\phi\rvert),\qquad
g(q) = 1 - q \ \ \text{(Soga)}\quad\text{or}\quad \tfrac12(1 - q^2)\ \ \text{(the dossier's $z$ law)}.$$

1. **Velocity extension** (Bothe and Soga, J. Differential Equations 468, 2026): the velocity at a
   point is the physical velocity at its foot point on the interface, so the extended field is
   constant along the normals, $\mathbf n\cdot\nabla\mathbf v^{\rm ext} = 0$. Then pure transport no
   longer stretches the level set: $D_t q = -q\,\mathbf n\cdot\nabla\mathbf v\,\mathbf n$ becomes zero.
   With an exact extension, $\phi(t)$ stays the signed distance to $\Gamma(t)$ in a tube, for all time.
2. **The foot point.** Analytic: $\Pi_\Sigma(x)$, the closest point. Algebraic (Newton point):
   $P_x = x - \phi\nabla\phi/\lvert\nabla\phi\rvert^2$. Capped (Bothe's notes):
   $x - \delta(x)\nabla\phi/\lvert\nabla\phi\rvert_{\max}$. Geometric (joint dossier): the closest point
   on the reconstructed zero iso-surface, which decouples the sample point from the noisy cell
   gradient.
3. **The blend** (Bothe's notes): $\mathbf v^{\rm ext} = \delta\,\mathbf v(x) + (1-\delta)\,\mathbf v(P_x)$, with
   $\delta = 0$ near $\Sigma$ and $\delta = 1$ far away.
4. **The source** (Bothe, Fricke and Soga, Math. Ann. 390, 2024; the joint dossier): it drives $q$ to 1
   and multiplies $\phi$, so it does not move the zero set in the continuum. With an exact extension
   the interface strain is already zero, so the source only has to remove the initial error of $q$
   and the discretisation drift. In the ideal extension, $z = q^2$ obeys $z' = \mu z(1 - z)$ exactly.
5. **Algebraic improvements** (Bothe's notes, item 0): $\phi \to \phi/\lvert\nabla\phi\rvert$ (first order)
   and $\phi \to \hat d = 2\phi/(\lvert\nabla\phi\rvert + \sqrt{(\lvert\nabla\phi\rvert^2 - 2\phi\,\partial_{nn}\phi)^+})$
   (second order); the curvature from $\hat d$.

The model has five parts that can fail separately: the foot point, the velocity sample at the foot
point, the transport of $\phi$ with a non-solenoidal velocity, the estimate of $q$ in the source, and
the coupling of all of these with a flow solver. This plan verifies them one at a time.

---

## 2. What we know, re-read critically

Source: two read-only audits of 2026-11-12 (the code of the nine extension models and of
`eulerianAdvection`; the evidence chain of every negative verdict). Details in section 12. The
summary:

1. **No extension was ever measured against the exact extension $\mathbf u(\Pi(x))$.** The static test
   `leiaTestVelocityExtension` measures only $\lvert\hat{\mathbf n}\cdot\nabla\mathbf u_e\rvert$, which a
   constant field passes. Every verdict comes from a transport or a coupled study.
2. **The flows were reversed and the shape was read at $T$.** The vortex and the shear reverse, and
   a linear transport with a prescribed velocity cancels much of its error on the way back. An
   extension depends on $\phi$, so it breaks that cancellation. The `projectionBlend` result of
   2026-10-31 (STATUS 11.28 item 9: band gradient error halved at $T/2$, shape at $T$ 4.5 to 5 times
   worse in 2D) fits this reading, and on the unresolved 3D deformation the same extension improves
   every metric. HYPOTHESIS, not tested: the reversed-flow read-out is biased against every extension.
3. **The semi-Lagrangian verdicts do not transfer.** FP0 and HL0 ran in the semi-Lagrangian solver
   with a source step and the `projectedFlux` trace; the "legacy models are dominated" verdict
   compares their accuracy with the semi-Lagrangian transport. None of these is a verdict on an
   extension in the Eulerian FVM.
4. **The Eulerian baseline itself is not verified.** Its band gradient error does not converge
   (orders -0.3 to 0.0, STATUS 11.28 item 9), and serial and np 4 of the 2D vortex baseline differ
   by 1.6e-6 relative at step 1 (STATUS 11.28 item 5, open). A transport that is not verified cannot
   judge an extension.
5. **Solver traps have voided data before.** The static source study of 2026-10-14 was void because
   the psi solver did zero iterations below its tolerance (STATUS 11.27 part 2). Every gate of this
   plan counts the solver iterations.
6. **The defect of the sources so far is the estimate of $q$, not the law** (STATUS 11.27 part 3,
   MEASURED): the centred gradient lets grid-scale extrema grow at $G(0)$ under every law. The
   one-sided `upwindFit` q passed at zero velocity in the SL step (STATUS 11.26). The Eulerian source
   has not been tested with a one-sided q.
7. **The rescaling is a one-time map, not a per-step operation** (STATUS 11.28 item 8): applied once
   it gives a second-order distance with the zero set unchanged; fed back every step it diverges. It
   has never been used as a read-only geometry evaluator (section 9, idea 3).

8. **A flux confound in every extension model except `none` and `haloLimited`** (code audit,
   CONFIRMED, `interfaceExtension.C` 435-468, `projectionBlend.C` 172): the transport flux is
   rebuilt as `interpolate(Uext) & Sf` on every face, so even far from the interface, where
   `Uext = U`, the exact prescribed flux is replaced by an interpolated one. Every comparison against
   `none` therefore also compared two far-field fluxes. The `projectionBlend` smoke check used a
   uniform translation, where the two fluxes are equal, so it could not see this.
9. **A hard jump at the band edge** in `closestPoint` and `meshWaveExt` (CONFIRMED): the fade blends
   `Uext` with `U` between L and 2L, but these models write values only inside L, so the velocity
   jumps at |psi| = L. The code's own comment records a 10x growth of the band gradient error from
   this jump.
10. **At T/2 the extensions win in every data set the audit found** (the bench at N = 512, the
    deck's steady vortex, `projectionBlend`): the losses are read at T of reversed flows.

Consequence: the earlier verdicts on velocity extension are WEAK or VOID for this plan's question.
They are kept in the knowledge base as what they are, verdicts on a configuration, and the models
are measured again from the bottom.

---

## 3. Rules of this campaign

1. **Every sub-problem has an exact reference.** No verdict comes from a comparison with another
   method's numerical result.
2. **Non-reversing flows first.** The kinematic references are flows whose interface $\Gamma(t)$ is
   known in closed form at every $t$: translation, rigid rotation, and the linear flows
   $\mathbf u = \mathbf A(\mathbf x - \mathbf x_0)$, which map a circle to an ellipse with known axes and
   orientation. A reversed flow is used only with a forward reference at $T/2$.
3. **Oracles separate the parts.** An extension computed from the exact interface (the oracle)
   isolates the transport; a transport with the exact $\phi$ isolates the extension.
4. **Pre-registration in the config header**: the prediction, the threshold, the falsifying result,
   and "how this gate could pass for the wrong reason" with the control that excludes it.
5. **Every gate counts the solver iterations** (no psi solve with a nonzero residual and zero
   iterations), runs serial against np 4 for the seam, reads the library stamps, and checks the
   patches of the mesh.
6. **At least three resolutions, the order next to the error, L1 and L2 only.** Uniform hexahedra
   first, then perturbed hexahedra, then cfMesh polyhedra.
7. **The whole vector**: the field error of $\phi$ in a band of fixed physical width, the band error
   of $q$, the zero-set position, the volume, the curvature, the cost. For the coupled cases add the
   spurious current and the pressure jump.
8. **Everything runs on Lichtenberg**; the laptop builds and runs unit tests only if the supervisor
   allows it (current rule: no runs on the laptop).
9. **A bug found is a result.** It is recorded in STATUS.md and in the knowledge base the same day,
   and every number it touched is re-measured.

---

## 4. Sub-problem 0: the Eulerian transport without an extension

**Question.** Does `leiaLevelSetFoam` with `ADVECTION eulerian`, `VELOCITY_EXTENSION none` transport a
signed distance correctly where the exact solution stays a signed distance?

**Why first.** Every later sub-problem measures a change against this transport. Item 4 of section 2
says it is not verified.

**Cases and references.**

| case | velocity | exact $\phi(t)$ | what it tests |
|---|---|---|---|
| S0.1 | translation, $\mathbf u = (1, 0.5)$ | $d(\mathbf x - \mathbf u t)$, a signed distance for all $t$ | the order of the transport alone; $q$ must stay 1 to the discretisation error |
| S0.2 | rigid rotation about a point off the circle centre | the rotated signed distance | a curved trajectory, no strain |
| S0.3 | linear strain $\mathbf u = a(x, -y)$ | $\phi_0(e^{-\mathbf A t}\mathbf x)$: not a signed distance; $q$ drifts by the exact law | that the transport reproduces the exact drift of $q$ (the reason the extension exists) |
| S0.4 | manufactured: $\phi_{\rm ex}(\mathbf x, t)$ chosen, the source $\partial_t\phi_{\rm ex} + \mathbf u\cdot\nabla\phi_{\rm ex}$ added | $\phi_{\rm ex}$ | the order of the whole discretisation with a non-solenoidal $\mathbf u$ (the material form) |

**Checks.** The order of the field error in the band and of the $q$ error in the band; the
`nDefCorr` loop (does a fourth correction change the result?); the div scheme and its limiter; the
material term $-\mathrm{Sp}(\nabla\cdot\Phi)$ (exactly zero for a solenoidal prescribed flux?); serial
against np 4, with the open 1.6e-6 dependence of STATUS 11.28 item 5 located.

**Prediction (to pre-register).** S0.1 and S0.2: the field error converges at order 2 or close to it
(`linearUpwind` with deferred correction) and the $q$ band error at order 1 or better. S0.3: the
computed $q$ follows the exact $q$ to the same order. If the $q$ band error does not converge on
S0.1, the Eulerian baseline has a defect, and that defect explains the non-convergence of
STATUS 11.28 item 9.

---

## 4b. Step 0.5 (proposed): what range of q is enough? A static gate

**Why.** The supervisor decides the target of the extension and the source. DERIVED: the zero
set, the normal and the curvature are invariant under $\phi \to f(\phi)$, $f$ increasing; the
operators that read $\phi$ as a function (the Detrixhe-Aslam phase indicator, the normal, the cell
curvature) are discretely invariant under a uniform $\phi \to c\phi$ (a linear gradient, a linear
fit), and they suffer only from a profile that bends across a cell, $\varepsilon = h\lvert\partial_n q\rvert/q$
(the zero of the linear fit moves by about $\tfrac{h^2}{8}\lvert\partial_n q\rvert/q$). The operators
that read $\phi$ as a distance (`cellCentreInverse`, a `distance` narrow band, the mollifier, the
$\phi$ width of a smeared indicator, the Newton point) need $q$ near 1. So the decision may need
two numbers: a robustness bound $q_{\min} \le q \le q_{\max}$ and an accuracy bound on $\varepsilon$.

**The gate.** No flow; a circle and an ellipse, three resolutions; $\phi = c\,d$ with
$c \in \{0.25, 0.5, 1, 2, 4\}$; $\phi = d(1 + s\,d/h)$ with $s \in \{0, 0.01, 0.03, 0.1, 0.3\}$; $\phi = d\,g(\theta)$.
Read-out with the production operators: the volume and shape error of $\alpha$, the normal error,
the curvature error with and without `cellCentreInverse`; interFoam's curvature error on the same
mesh as the bar. Prediction: $c$ changes nothing except through `cellCentreInverse`; the errors
grow with $\varepsilon$. Waits for the supervisor's go.

## 4a. Step A0: repair what the audit found, before any measurement

Each repair is an inert-by-default change, gated by bit identity where it claims to be inert.

1. **The flux.** The correction form `phiExt = phi + interpolate(Uext - U) & Sf` for every model of the
   `interfaceExtension` base and for `projectionBlend`. Gate: each model with a zero extension
   (`Uext = U`) reproduces `none` bit for bit on a non-uniform flow (the 2D vortex), serial and np 4.
2. **The band edge.** `closestPoint` and `meshWaveExt` compute the extension out to the end of the
   fade (2L), or fade inside L. Gate: the velocity jump at the band edge is zero to round-off in the
   static test.
3. **The `closestPoint` failure switch.** The far-field steady solve runs always or never, not when
   any cell on any rank fails. The descent becomes a true Newton step on distance.
4. **The fixable settings.** `pseudoTime` iterates until it reaches the band edge; `anisotropicDiffusion`
   iterates its lagged cross term and gets a dimensionally correct `eps`; `steadyUpwind` stops
   repeating the same solve.
5. **The static test** gets the exact reference of section 5 and drops L_inf.

The Eulerian transport (`eulerianAdvection`, sub-problem 0) is NOT changed in A0; `Sp` against
`SuSp` for the material term is an arm of sub-problem B.

## 5. Sub-problem A: the extension operator, static

**Question.** Given $\phi$ and $\mathbf u$ at one instant, how close is the computed $\mathbf u_e$ to the
exact extension $\mathbf u(\Pi(\mathbf x))$, and at what order and cost?

**The test application.** `leiaTestVelocityExtension` gains:

1. an exact reference: the analytic surface of `levelSetImplicitSurfaces` (circle and sphere in
   closed form; `signedDistanceEllipse` and `signedDistanceEllipsoid`, which already carry a verified
   closest point) gives $\Pi(\mathbf x) = \mathbf x - d\nabla d$ and $\mathbf u_e^* = \mathbf u(\Pi(\mathbf x))$ from the
   prescribed velocity model evaluated at that point;
2. the metrics $E_{\rm ext} = \lVert\mathbf u_e - \mathbf u_e^*\rVert$ in L1 and L2, over a band of FIXED
   physical width $\lvert d\rvert \le W$ and over $\lvert d\rvert \le k h$; the error of the face flux
   $\lVert\Phi_e - \mathbf u_e^*(\mathbf x_f)\cdot\mathbf S_f\rVert$, because the transport sees the face flux,
   not the cell velocity; the foot-point error $\lVert\Pi_h - \Pi\rVert$ where a model exposes it;
   $\hat{\mathbf n}\cdot\nabla\mathbf u_e$ as the old secondary metric; the wall time per call;
3. four input fields with the same zero set: $d$; $d\,g(\mathbf x)$ with a smooth $g$ between 0.8 and
   1.25 ($q \ne 1$, smooth); $r^2 - R^2$; and $\hat d$ of $r^2 - R^2$. Input 2 separates the Newton
   point from the geometric foot point, which is the joint dossier's argument;
4. a CSV that the gate tooling reads, one row per case.

**Interfaces and velocities.** A circle and an ellipse, centre off the cell lattice; velocities: the
linear flows (rotation, strain, simple shear; a new prescribed model `linearFlow`, $\mathbf u = \mathbf A(\mathbf x
- \mathbf x_0) + \mathbf b$), and the vortex `shear2D` (nonlinear). 3D: sphere and ellipsoid with
`deformation3D`.

**Candidates.** The nine models of `libleiaVelocityExtension` as they are, each with the settings
that the audit (section 12) finds most favourable and with its default; and three new ones:

| new model | what it computes | why |
|---|---|---|
| `exactFootPoint` | $\mathbf u(\Pi(\mathbf x))$ from the analytic surface and the prescribed velocity | the oracle for sub-problem B; also the 0 of the metric |
| `newtonPoint` | $\mathbf u_h(P_x)$, $P_x = \mathbf x - \phi\nabla\phi/\lvert\nabla\phi\rvert^2$, with $\mathbf u_h$ sampled by locating the cell of $P_x$ and interpolating (not by a Taylor step from the cell centre) | Bothe and Soga's continuum extension, discretised directly |
| `newtonPointCapped` | as `newtonPoint` with Bothe's cap $\delta(\mathbf x)$ and blend | his notes, page 1 |

`projectionBlend` stays as it is (its sample is a first-order Taylor step from the cell centre);
the comparison with `newtonPoint` measures what that Taylor step costs.

**Prediction (to pre-register per model).** `exactFootPoint` at round-off. `closestPoint` and
`newtonPoint` at order 2 on input 1; `newtonPoint` loses order on input 2 where $q$ is far from 1,
the geometric foot point does not. `projectionBlend` at order 1 at best (its Taylor step).
`steadyUpwind` and `meshWave` at order 1. `anisotropicDiffusion` and `pseudoTime`: undetermined; if
they do not converge, the audit tells whether a setting causes it.

**Exit.** A ranking of the models by error, order and cost per call, on hexahedra and on polyhedra,
serial and np 4. The two or three best go to sub-problem C. A model that fails here is not tested
in transport until the cause is known.

---

## 6. Sub-problem B: transport with the exact extension (oracle)

**Question.** With the exact extension, does the Eulerian transport keep $\phi$ a signed distance,
and does it put the interface where the flow puts it?

**Reference.** With an exact extension, $\phi(\mathbf x, t) = d(\mathbf x, \Gamma(t))$, with $\Gamma(t)$ the
interface moved by $\mathbf u$. Closed form for translation, rotation and every linear flow (a circle
becomes an ellipse with axes and orientation from $e^{\mathbf A t}$; its signed distance from
`signedDistanceEllipse` in the rotated frame). For the vortex, a marker polygon integrated with RK4
on the analytic velocity gives $\Gamma(t)$ to round-off; the reference is read FORWARD, at $t \le T/2$.

**Arms.** The material form (the code: $\mathrm{ddt} + \mathrm{div}(\Phi_e, \phi) - \mathrm{Sp}(\nabla\cdot\Phi_e,
\phi)$) against the conservative form; the div scheme; the time step law $\Delta t \sim h$.

**Prediction.** The field error of $\phi$ in the band converges at the transport's order (sub-problem
0), and the $q$ band error converges to zero, which the baseline without an extension cannot do on
S0.3. If the oracle does not keep $q$ near 1, the transport of a non-solenoidal velocity is the
defect, and no discrete extension can be judged before it is repaired.

---

## 7. Sub-problem C: transport with the discrete extension

**Question.** How much of B survives when the extension is computed from $\phi_h$?

The same cases as B with the best models from A. The difference C minus B is the cost of the loop
$\phi_h \to \Pi_h \to \mathbf u_e \to \phi_h$ that the joint dossier names. Additional checks: the time
step (the loop is explicit in the extension), the band width, the seam, the polyhedral mesh. A
reversed flow enters here for the first time, read at $T/2$ against the marker reference and at $T$
only for the shape, and the baseline at $T$ is reported with the warning of section 2 item 2.

---

## 8. Sub-problem D: the source term

**D1, the source alone.** $\mathbf u = 0$, $\phi_0 = c\,d$ (q constant $= c$): along every ray $z = q^2$
follows the logistic law exactly, so the whole field has a closed form at every $\tau$. The harness
of STATUS 11.27 (`cases/2DstaticCircle`, `config/staticSource/`, `static_source_metrics.py`, minIter
1) is reused. Arms: the laws `linearQ` and `linearZ` with the centred $q$ (the known defect, as the
control), with the one-sided q of STATUS 11.26 ported to the Eulerian source, and with $q$ from
$\hat d$. Prediction: the one-sided q follows the logistic law at order 1 or better and does not
grow the grid-scale mode; the centred q reproduces the growth at $G(0)$ of STATUS 11.27.

**D2, the source with the oracle extension**, on the cases of B with $\phi_0 = c\,d$ or $d\,g$: the
source must bring $q$ to 1 at the rate of the law without moving the zero set beyond the transport
error.

---

## 9. Out of the box: hypotheses that are cheap to test inside the sub-problems

1. **Only the normal component of the extension matters.** $\mathbf u_e\cdot\nabla\phi = (\mathbf u_e\cdot\mathbf
   n)\lvert\nabla\phi\rvert$: the tangential part of $\mathbf u_e$ does nothing in the continuum, but in the
   FVM it is a non-solenoidal flux that the transport must carry. Arm: $\mathbf u_e = (\mathbf u(\Pi)\cdot\mathbf
   n)\,\mathbf n$ in A and C. If it gives the same $\phi$ with less divergence, the extension is simpler.
2. **Extend the face flux, not the cell velocity.** The transport uses $\Phi_e$ at the faces. A model
   that evaluates $\mathbf u(\Pi(\mathbf x_f))\cdot\mathbf S_f$ at each face removes one interpolation and the
   kink of the cell field at the medial axis. Arm in A (the face-flux metric) and C.
3. **Read-only geometry.** The rescaling $\hat d$ diverged when fed back every step and works once
   (STATUS 11.28). For a two-phase solver the level set is needed for the normal, the curvature and
   the phase indicator, not for its own sake. Arm: transport $\phi$ without any control, and compute
   the normal, the curvature and $\alpha$ from $\hat d(\phi)$ each step without changing $\phi$. If
   $q$ stays within, for example, $[0.5, 2]$ in the band over the run, this may be enough for the
   multiphase target and is the simplest method. It is the control arm of sub-problem E.
4. **The target is a bound, not $q = 1$.** The dossier's target is $0.9 \le q \le 1.1$; the multiphase
   need may be weaker. Sub-problem E measures which $q$ range keeps the curvature and the volume at
   or below interFoam's error, and that range becomes the criterion.
5. **The extension only where it is needed.** With Bothe's blend, $\delta = 1$ outside a few cells:
   the extension acts where the phase indicator and the curvature read $\phi$.

---

## 10. Sub-problem E: assembly and coupling

1. **E1, kinematic**: the best extension of C plus the best source of D on the cases of B and C,
   against the oracle; the control arm of section 9 idea 3.
2. **E2, coupled**: `leiaLevelSetTwoPhaseFoam` on the stationary, the translating and the oscillating
   droplet, three resolutions, np 4, against interFoam (`interFoamDroplet2D`) on the same mesh, time
   step and fluids, and against the Eulerian baseline. The whole vector: shape, volume, $q$ band,
   curvature, spurious current, pressure jump, cost. Escalation to 3D and polyhedra only on a pass.
3. Before E2: the well-balanced gate with exact curvature, which isolates the force balance from the
   level-set method (the CLAUDE.md research ladder).

---

## 11. The deck

`docs/pde-level-set/pde-level-set-presentation/pde-level-set.template.html` (a theme folder of its own
for this direction, `docs/pde-level-set/README.md`), built with `workflow/scripts/export_html.py`.
It also carries the transcription of D. Bothe's notes of 2026-10-29 and their reading. One vertical stack per sub-problem, each with: the question, the
test and why it is exact, the prediction, the result (a figure, few numbers), the reading (MEASURED /
DERIVED / HYPOTHESIS), and a "questions for the supervisor" box. The first stacks: the model, the
critical review of the old verdicts, the plan, sub-problem 0. Every later sub-problem adds its stack
when its gate lands.

---

## 12. The audits of 2026-11-12

Two read-only audits of the up-to-date tree (456cae4f): the evidence chain of every negative verdict
on velocity extension, and the code of the nine models, of `eulerianAdvection` and of the static
test. Labels: STRONG = clean, isolated, several resolutions; WEAK = confounded or one rung; VOID = a
setup later found wrong.

### 12.1 The verdicts

| verdict | where | main confound | strength |
|---|---|---|---|
| the legacy extensions give the worst accuracy at 12 to 27 times the cost | `benchVortexVET2/T8`, only `closestPoint`, reversed vortex, N 32 to 512, np 4, 2026-07-30 | compared with the SL method, not with Eulerian `none`; against `none` the volume error at T is better (1.11e-2 against 2.11e-2) and the band gradient error at T/2 is better (0.87 against 1.67 at N = 512); np 4 before the seam fixes of 2026-08-26 (808d384b); non-monotone errors at N = 512 | WEAK, close to VOID |
| `closestPoint` 2185 times worse on the band gradient | the same row, read at T | the time-reversal cancellation of `none` (1.67 at T/2 to 0.0013 at T); a serial re-measure after the fixes is better | VOID as stated |
| `anisotropicDiffusion` cannot converge (static) | `config/staticExtension.yaml`, serial | the lagged tensor cross term is never iterated (code audit) | STRONG for the code, WEAK for the method |
| `pseudoTime` cannot converge (static) | the same | its reach (5 steps of 0.5 h) is shorter than the band | WEAK: a setting |
| `steadyUpwind` degrades, `steadyUpwindLinear` diverges | the same | the deferred correction on an undamped steady equation; `steadyUpwindLinear` is linear along the normal by design, not `u(Pi)` | MEDIUM |
| `meshWave` stays first order | the same | seed quantisation, intrinsic | STRONG (static) |
| translating-droplet survival times | SL two-phase, N = 32, 2026-07-26 | the closed-box setup voided on 2026-09-02 | VOID |
| FP0 fails the 2D method gate (shear 34x, order 1.38; seams) | `methodGate2D`, SL | reversed shear read at T; the extension AND the trace changed together; the seam failure is the halo Taylor surrogate that 808d384b names as not fixed | WEAK on accuracy; STRONG that this implementation depends on the decomposition |
| FP0 stationary and oscillating time out | the same | a cost bottleneck, not accuracy | not a verdict |
| HL0 fails (shear 14x, order 1.39) | the same | reversed shear at T; trace changed with the extension; R <= h hard-coded, the R ladder never ran; the sampler extrapolates up to 1.5 h | MEDIUM for "R = h relocates the strain", WEAK otherwise |
| `projectionBlend` B1 fails | STATUS 11.28 item 9, Eulerian | reversed flows read at T (the gradient at T/2 halves, the 2D volume at T/2 is 16x lower); the flux confound (12.2); U(P) is a first-order Taylor step up to 4 h; on the unresolved 3D deformation it improves every metric | WEAK |

Also found: earlier "all three models abort" records were configuration errors (the missing `alpha`
key, the type name `meshWaveExt`, the missing `UextFinal` solver entry). The skeleton article
`docs/velocity-extension/velocity-extension-article/velocityExtension.tex` makes no claim; the deck of
record (July 2026, Eulerian, before the seam fixes) names `closestPoint` the static winner at about
order 2, and on its steady (non-reversing) vortex the extensions keep the interface alive longer than
`none`.

### 12.2 The code

| model | what it computes | main defects (CONFIRMED unless marked) | parallel |
|---|---|---|---|
| base `interfaceExtension` | `phiExt = interpolate(Uext) & Sf` everywhere; seeds = cells with 0 < alpha < 1; foot point by `nDescent` steps `x -= psi nHat` (not a Newton step, no 1/q); fade from L to 2L | the far-field flux replacement; seeds stop at a rank boundary; alpha is stale between write times when `reportAtWriteTimesOnly` is true | seeds depend on the decomposition |
| `anisotropicDiffusion` | one solve of `laplacian(nn + eps I, Uext) = 0` over the domain | the tangential tensor part is explicit and never iterated; `eps` has units of length in `cellSize` mode | safe |
| `pseudoTime` | 5 implicit pseudo-steps of 0.5 h, upwind | reach 2.5 h shorter than the band | safe |
| `steadyUpwind` | steady `(w.grad) Uext = 0`, 5 passes | with upwind the passes repeat the same solve; first order | safe |
| `steadyUpwindLinear` | `u(Pi) + d dn u`, linear along the normal | not `u(Pi)` by design | safe |
| `closestPoint` | descent with a 0.9 h step limit, `cellPoint` sample at the foot; a global far-field solve when any cell fails | no fade (band-edge jump); the global `nFailed` switch changes the far field between steps and decompositions; not a Newton step; the acceptance test reads `cellPoint` psi in every mode | collectives correct; decomposition-dependent through the one-layer halo and the switch |
| `meshWaveExt` | the foot point of the nearest sign-change face, by FaceCellWave | internal faces only; constant per Voronoi region; no fade | decomposition-dependent |
| `haloLimited` | `(1-w) u(x) + w u_h(Y)`, |Y - x| < R <= h, flux in correction form | `u(Pi)` only for |d| <= h by design; first-order distance | safe |
| `projectionBlend` | `delta U + (1-delta)(U + dx.grad U)`, `dx = -dist n` | far-field flux replacement; first-order Taylor sample, error up to about 8 h^2 |grad grad U|; `h` from the bounding box | safe |

`eulerianAdvection`: `correct()` once per step from psi^n and U^{n+1}; the material form
`ddt + div(phiExt, psi) - Sp(div phiExt, psi)` with 3 deferred-correction passes (consistent for a
non-solenoidal flux; `SuSp` is the safer choice when `div phiExt > 0`). The static test
`leiaTestVelocityExtension` measures only `|nHat . grad Uext|` inside the band: it cannot see a wrong
interface value, the flux replacement, the band-edge jump, time behaviour, or the decomposition; and
it reports L_inf.

### 12.3 What this means

The supervisor's suspicion holds for most of the record: of the negative verdicts, only the static
verdicts on `meshWave` and part of those on `steadyUpwind` and `anisotropicDiffusion` are clean. The
advected verdicts carry at least one of: the reversed-flow read-out, two changes at once, the flux
confound, the seam defects of before 2026-08-26, or a voided setup. They are re-measured after A0.

---

## 13. Infrastructure to build, in order

1. `linearFlow`, a prescribed velocity model $\mathbf u = \mathbf A(\mathbf x - \mathbf x_0) + \mathbf b$, inert by
   default, with the exact $\Gamma(t)$ for a circle start (a header with the ellipse axes and angle).
2. A header of exact references: $\Pi(\mathbf x)$, $d(\mathbf x, \Gamma(t))$, $\mathbf u_e^*(\mathbf x)$ for the analytic
   surfaces and the prescribed velocities.
3. A function object that writes, per step, the field error of $\phi$ against the exact $d(\mathbf x,
   \Gamma(t))$ in a fixed physical band, the $q$ band error, the zero-set position error and the volume
   error; the solver iteration counts.
4. `leiaTestVelocityExtension` upgraded (section 5); the new models `exactFootPoint`, `newtonPoint`,
   `newtonPointCapped`, each a new class in the extension family (CLAUDE.md "Extension without
   modification").
5. A marker reference for the vortex (`workflow/scripts/marker_reference.py`, RK4, a polygon, its
   signed distance on the cell centres).
6. Study configs under `config/veVerification/`, one per sub-problem, with the pre-registered
   read-out; a report script that writes the tables and the deck figures.

Each item gets its unit gate (exact solution, exits nonzero) before any study uses it.

---

## 14. Decisions the supervisor takes during the campaign

1. After sub-problem 0: is the baseline transport good enough, or is its repair the first result?
2. After A: which two or three extensions go on.
3. After B and C: whether the extension needs the source at all for the multiphase target, and the
   $q$ range that is enough (section 9 idea 4).
4. Before E2: the interFoam comparison set (cases, resolutions, fluids).
