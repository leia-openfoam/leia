# Normal-Projected Semi-Lagrangian Level Set with Per-Step Geometric Renormalization

Method description, current understanding. Status: design stage — derived
and internally checked, not implemented, not measured. Companion documents:
the baseline solver configuration note (METHOD.md) and the stabilized 3D
foot-point algorithm (stable_foot_point_3d.md, needed only for the optional
refinement in §7.2).

Working name: **normal-projected SL** (nSL). The scheme replaces the
baseline semi-Lagrangian value transport
$\psi^{n+1}(\mathbf x_c)=\psi^n(\mathbf x_d)$ by a purely geometric update
in the interface-normal direction, with a built-in, zero-set-preserving
renormalization of the distance property in every cell at every step.

---

## 1. Summary of the scheme

Per time step, per cell $c$ with centre $\mathbf x_c$:

1. Perform the quadratic weighted-least-squares reconstruction of $\psi^n$
   over the cell's stencil (identical to the baseline scheme). This yields
   the value $\psi_c=\psi^n(\mathbf x_c)$, the gradient
   $\mathbf g=\nabla R_c(\mathbf x_c)$ and the Hessian
   $\mathbf H=\nabla\nabla R_c$ of the local fit.
2. Build the **normal-projected velocity**
   $\mathbf u_n=(\mathbf u\cdot\mathbf n)\,\mathbf n$ with
   $\mathbf n=\mathbf g/|\mathbf g|$, at the current and previous time
   levels.
3. Compute the departure point $\mathbf x_d$ with the same second-order
   backward Taylor trace as the baseline scheme, but driven by
   $\mathbf u_n$ instead of $\mathbf u$.
4. Convert the stored value $\psi_c$ into a **geometric offset** $d_c$
   (signed distance from $\mathbf x_c$ to the zero set of the local
   quadratic model), using the stable quadratic root formula with a
   first-order fallback (§5).
5. Measure the **normal displacement** $\delta$ of the trace (§6).
6. Write

$$\boxed{\;\psi^{n+1}(\mathbf x_c)=d_c+\delta\;}$$

No value of $\psi$ is ever interpolated at $\mathbf x_d$. The reconstruction
is consulted only for directions, speeds and geometry — never to resample
the transported field.

---

## 2. Why the normal-only trace is exact, not an approximation

Decompose the velocity at any point into components normal and tangential
to the local level set of $\psi$:

$$\mathbf u=v_n\,\mathbf n+\mathbf u_t,\qquad v_n=\mathbf u\cdot\mathbf n,
\qquad \mathbf n=\frac{\nabla\psi}{|\nabla\psi|} .$$

Because $\nabla\psi$ is normal to its own level sets,
$\mathbf u_t\cdot\nabla\psi=0$ **identically**, so

$$\mathbf u\cdot\nabla\psi=v_n\,|\nabla\psi|=\mathbf u_n\cdot\nabla\psi,
\qquad \mathbf u_n:=v_n\,\mathbf n .$$

Hence the level-set equation is equivalent to its Hamilton–Jacobi form
$\partial_t\psi+v_n|\nabla\psi|=0$, and $\psi$ is **exactly constant along
trajectories of $\mathbf u_n$**, for any field $\psi$, drifted or not:

$$\frac{\partial\psi}{\partial t}+\mathbf u_n\cdot\nabla\psi=0 .$$

Tracing backward along $\mathbf u_n$ therefore solves the same continuous
transport problem as tracing along $\mathbf u$. The tangential velocity
component carries no information about $\psi$ and is discarded without
loss.

Consequences:

- **Pure tangential motion is an exact no-op.** If $\mathbf u\perp\mathbf n$
  everywhere, then $\mathbf u_n=0$, $\mathbf x_d=\mathbf x_c$, $\delta=0$,
  and $\psi^{n+1}=d_c$ — the field is renormalized but not moved. In the
  baseline scheme, the same situation still resamples $\psi$ through the
  reconstruct-and-evaluate operator every step.
- **The zero-set amplification operator leaves the value path.** The
  grid-scale corrugation growth of the baseline scheme is a property of
  repeatedly resampling $\psi$ through the local fit. In nSL, the written
  value is composed of the cell's own stored value (converted to geometric
  units) plus a geometric length; the fit's noise enters only through
  directions and speeds, multiplied by $\Delta t$ — a qualitatively weaker
  coupling.

---

## 3. Required fields and their sources

All from data already computed by the baseline solver — no new
reconstructions:

| quantity | definition | source |
|---|---|---|
| $\psi_c$ | stored cell value at $t^n$ | field storage |
| $\mathbf g$, $\mathbf H$ | gradient, Hessian of local fit at $\mathbf x_c$ | quadratic WLS fit (baseline §2.2) |
| $\mathbf n$ | $\mathbf g/|\mathbf g|$ | fit |
| $h_{nn}$ | $\mathbf n^{\mathsf T}\mathbf H\,\mathbf n$ | fit |
| $\mathbf u^n$, $\mathbf u^{n-1}$ | cell-centred velocity, two levels | flow solver |
| $\nabla\mathbf u$ | velocity gradient | flow solver |
| $\mathbf u_n^{n}$, $\mathbf u_n^{n-1}$ | normal-projected velocities | computed, §4; the $n{-}1$ level must be stored from the previous step (projected with the *previous* step's normals) |

---

## 4. The normal-projected trace

### 4.1 Projected velocity and its gradient

$$\mathbf u_n=(\mathbf u\cdot\mathbf n)\,\mathbf n=v_n\,\mathbf n$$

For the second-order convective term of the trace, the gradient of
$\mathbf u_n$ is needed. With $\otimes$ the outer product:

$$\nabla\mathbf u_n=\nabla v_n\otimes\mathbf n+v_n\,\nabla\mathbf n$$

$$\nabla v_n=(\nabla\mathbf u)^{\mathsf T}\mathbf n+(\nabla\mathbf n)^{\mathsf T}\mathbf u,
\qquad
\nabla\mathbf n=\frac{(\mathbf I-\mathbf n\mathbf n^{\mathsf T})\,\mathbf H}{|\mathbf g|}$$

Every factor is available at the cell: $\nabla\mathbf u$ from the momentum
solver, $\mathbf g,\mathbf H$ from the fit. (Index convention: whichever
convention is used for $\nabla\mathbf u$ in the code, the same must be used
consistently here; the formulas above are written for
$(\nabla\mathbf u)_{ij}=\partial u_i/\partial x_j$, i.e.
$\mathrm d\mathbf u=(\nabla\mathbf u)\,\mathrm d\mathbf x$.)

### 4.2 Departure point

Identical in structure to the baseline second-order backward trace, with
$\mathbf u\to\mathbf u_n$:

$$\mathbf x_d=\mathbf x_c
-\frac{\Delta t}{2}\left(3\,\mathbf u_n^{\,n}-\mathbf u_n^{\,n-1}\right)
+\frac{\Delta t^{2}}{2}\left(\mathbf u_n^{\,n}\cdot\nabla\right)\mathbf u_n^{\,n}$$

The Adams–Bashforth combination $\tfrac12(3\mathbf u_n^n-\mathbf u_n^{n-1})$
provides the time-centred (midpoint) speed to second order; the convective
term corrects the midpoint sample to the trajectory's spatial midpoint.
Under a varying time step, the same $\Delta t/\Delta t_0$ ratio as in the
baseline applies to the difference quotient.

**Storage note.** $\mathbf u_n^{\,n-1}$ is the *projected* velocity of the
previous step, projected with the previous step's normals. It must be
stored at the end of each step; re-projecting the old velocity with the new
normals is not equivalent and mixes time levels of the geometry.

### 4.3 CFL restriction

Unchanged from the baseline: the trace must remain within the cell's
stencil support (in practice, within one halo layer), which the capillary
time-step restriction already guarantees by a wide margin. Since
$|\mathbf u_n|\le|\mathbf u|$ pointwise, the nSL trace is never longer than
the baseline trace.

---

## 5. Geometric conversion of the stored value

### 5.1 Why a conversion is required

The update $d_c+\delta$ adds a geometric length $\delta$ (metres) to the
cell's offset. If the raw stored value $\psi_c$ were used in place of a
geometric offset, the scheme would mix units whenever $|\nabla\psi|\ne1$:
a 1D field $\psi=a\,x$ with $a\ne1$ under uniform normal speed $v$ would
propagate its zero crossing at $v/a$ instead of $v$ — an error largest at
the interface and growing with drift. The stored value must therefore be
converted to the **signed geometric distance from $\mathbf x_c$ to the zero
set of the local model** before the geometric increment is added.

### 5.2 Second-order conversion (stable quadratic root)

The pipeline is second order throughout (quadratic fit, second-order
trace); the conversion must be too. Restrict the local quadratic model to
the normal ray through $\mathbf x_c$:

$$R(s)=\psi_c+s\,|\mathbf g|+\tfrac12\,s^{2}\,h_{nn},
\qquad h_{nn}=\mathbf n^{\mathsf T}\mathbf H\,\mathbf n .$$

The signed distance from $\mathbf x_c$ to the model's zero level along the
ray is $d_c=-s^\ast$ where $R(s^\ast)=0$, taking the root that reduces to
$-\psi_c/|\mathbf g|$ as $h_{nn}\to0$. In the numerically stable,
cancellation-free form:

$$D=|\mathbf g|^{2}-2\,\psi_c\,h_{nn},
\qquad
d_c=\frac{2\,\psi_c}{\,|\mathbf g|+\sqrt{D}\,}\quad(\text{if accepted, §5.3}).$$

Expansion: $d_c=\dfrac{\psi_c}{|\mathbf g|}
\Bigl(1+\dfrac{\psi_c\,h_{nn}}{2|\mathbf g|^{2}}+\dots\Bigr)$ — the
first-order (linear) conversion plus the field-curvature correction that
the linear form drops. Cost: one square root per cell; both inputs
($|\mathbf g|$, $h_{nn}$) are already computed for the curvature
evaluation.

Properties:

- Exactly zero-preserving: $\psi_c=0\Rightarrow d_c=0$ for any
  $\mathbf g,\mathbf H$. The interface position is never altered by the
  conversion.
- On a clean signed-distance field ($|\mathbf g|=1$ and the fit exact),
  $d_c=\psi_c$ and the conversion is the identity.
- On a drifted field, the conversion is a per-step multiplicative
  reparametrization of $\psi$ — the legitimate class of corrections that
  changes the representative of the level-set family without touching the
  interface — and it drives $|\nabla\psi^{n+1}|$ back toward 1 to the
  accuracy of the local quadratic distance estimate.

### 5.3 Guard and first-order fallback

The root is real only if $D\ge0$; near $D=0$ the cell sits near the local
model's turning point, the conversion becomes ill-conditioned
($\partial d_c/\partial\psi_c\to\infty$), and the formula must be abandoned
*before* the singularity, not at it. Acceptance test with threshold
$\beta$ (suggested $\beta=\tfrac14$):

$$d_c=\begin{cases}
\dfrac{2\,\psi_c}{\,|\mathbf g|+\sqrt{D}\,} & D\ \ge\ \beta\,|\mathbf g|^{2}
\\[2.2ex]
\dfrac{\psi_c}{|\mathbf g|} & \text{otherwise (first-order fallback)}
\end{cases}$$

**Fallback rationale.** The two available failure behaviours are not
symmetric. The first-order conversion $\psi_c/|\mathbf g|$ remains
well-defined and degrades gracefully: worst case, that cell's
renormalization is first-order accurate for one step. Writing the *plain
SL value with no conversion* instead would place an algebraic value next
to neighbours holding geometric values — a units mismatch between adjacent
cells, i.e. an artificial tangential jump in the effective
$|\nabla\psi|$, which is exactly the quantity that drives spurious
capillary flow. Every cell therefore receives *some* zero-preserving
conversion every step; only its order varies.

With $\beta=\tfrac14$ the fallback engages once $\psi_c$ passes roughly
$3/8$ of the distance to the model's turning point — comfortably inside
the analogous half-local-radius guard of the curvature offset correction.

### 5.4 Consistency with the curvature offset correction

The curvature path uses an offset $d$ of the same meaning (distance from
the cell centre to the interface) in its parallel-curve correction
$\kappa=\kappa_d/(1-d\kappa_d)$. Whatever branch of §5.3 produced $d_c$
for the transport must feed the *same* $d_c$ to the curvature correction
in that cell and step, so transport and curvature never disagree about the
cell's offset. (Independently of nSL, replacing the curvature path's
first-order $d=\psi_c/|\mathbf g|$ with the §5.2 root is a candidate
improvement for the baseline scheme as well.)

---

## 6. The normal displacement δ

$\delta$ is the signed geometric displacement of the trace along the
normal:

$$\delta=(\mathbf x_d-\mathbf x_c)\cdot\mathbf n
\;=\;-v_n^{\,\text{mid}}\,\Delta t+O(\Delta t^{3})$$

Sign convention: with $\mathbf n$ pointing toward increasing $\psi$, a
flow moving the interface in the $+\mathbf n$ direction ($v_n>0$) has the
departure point on the $-\mathbf n$ side, $\delta<0$, and the cell's new
value is smaller — the interface has moved toward/past it in the
$+\mathbf n$ direction. Limiting checks: $v_n=0\Rightarrow\delta=0$
(renormalization only); on a clean SDF with uniform $v_n$,
$d_c+\delta=\psi_c-v_n\Delta t$, the exact translation.

**Projection vs. arc length.** The second-order convective term of the
trace is not exactly parallel to $\mathbf n$ (the $v_n\nabla\mathbf n$
part of $\nabla\mathbf u_n$ rotates the direction), so $\mathbf x_d$
acquires a small tangential component of size $O(\Delta t^{2})$. The
definition above *projects* it out. Open point (§9.4): whether the
projection or the full arc length along the traced path is the correct
second-order-consistent measure; the difference is $O(\Delta t^{3})$ per
step under the CFL restriction and therefore below the scheme's global
order, but the choice should be fixed and documented in the
implementation.

**Order discipline.** The scheme has no first-order component left, so
every geometric quantity in the update — $d_c$ and $\delta$ separately —
must carry local error $O(\cdot^{3})$ or better ($\cdot$ = offset for
$d_c$, $\Delta t$ for $\delta$). A first-order shortcut in either
re-caps the whole scheme.

---

## 7. Relationship to previously considered alternatives

### 7.1 Strain-rescaling ODE

The exact transport of the gradient magnitude,
$\mathrm D|\nabla\psi|/\mathrm Dt=-|\nabla\psi|\,\varepsilon_{nn}$ with
$\varepsilon_{nn}=\mathbf n\cdot\mathbf D\cdot\mathbf n$, motivated a
candidate update $d^{n+1}=d^n(\mathbf x_d)(1+\varepsilon_{nn}\Delta t)$ —
a first-order Lagrangian integrator of the stretching. In nSL this
correction is **subsumed**: the drift that the ODE would integrate in time
is instead *measured geometrically and removed* each step by the §5
conversion. The identity behind it: on a signed-distance field, the normal
speed varies along the normal as
$\partial v_n/\partial n=\varepsilon_{nn}$, so transporting each cell with
its own local $v_n$ (as nSL does) accrues exactly the strain drift per
step — which the conversion then renormalizes away. No separate rescaling
stage, and no time-integration-order limit on the renormalization, since
it is geometric rather than integrated.

### 7.2 Foot-point search against the local iso-level

An earlier variant proposed measuring the offset drift by running the
stabilized foot-point algorithm from $\mathbf x_d$ onto the iso-surface
$\{\psi^n=\psi^n(\mathbf x_c)\}$. In nSL this search is unnecessary in the
common case: with the trace already along the normal ray, the §5.2 ray
root *is* the foot-point distance up to a tangential deviation of relative
size $O((d\,\kappa)^2)$ — third order in the offset, below the scheme's
order for band-scale offsets. The full stabilized foot-point iteration
(documented separately) remains the systematic upgrade path if deep-band
or strongly drifted cells ever need it; one cycle of it on the local
quadratic model is the natural next refinement.

### 7.3 Global redistancing

Explicitly avoided, for the standing reason: an Eikonal-type redistancer
rebuilds distance to whatever zero set it is handed, including grid-scale
advection noise. nSL never interrogates the remote zero set: every
quantity in the update is local to the cell's own stencil, and the
zero-crossing of the local model is anchored by the cell's own stored
value.

---

## 8. Properties (expected, to be verified)

- **Interface kinematics.** Exact-in-principle: the normal-projected trace
  solves the same continuous transport as the full-velocity trace (§2).
  Zero-preservation of the conversion guarantees the renormalization never
  moves the front.
- **Distance maintenance.** $|\nabla\psi|$ is driven toward 1 in every cell
  at every step, to second order in the offset, without any global
  operator, and without ever having been allowed to drift far (the
  correction acts before drift accumulates).
- **Corrugation coupling.** The reconstruct-and-evaluate resampling of
  $\psi$ — the operator with measured amplification $>1$ at
  $\lambda\approx2h$ — is removed from the value path. Fit noise enters
  only via $\mathbf n$, $v_n$, $h_{nn}$, i.e. multiplied by $\Delta t$ (in
  the trace) or by the offset (in the conversion). Whether the resulting
  closed loop (noise → normals → trace/conversion → new field → next
  fit) has amplification $\le1$ for the $\lambda\approx2h$ mode is **the**
  question the first numerical test must answer; the acceptance criterion
  is the same as for the baseline (amplification factor $\le1$, not "mode
  is small").
- **Cost.** Per cell per step, relative to baseline: saves the polynomial
  evaluation at $\mathbf x_d$; adds the projection $\mathbf u_n$, the
  $\nabla\mathbf u_n$ assembly for the convective term, one square root,
  and storage of one extra vector field ($\mathbf u_n^{\,n-1}$). Same
  single halo exchange; embarrassingly parallel as before.

---

## 9. Open questions

1. **Closed-loop stability.** The trace direction now depends on the
   reconstructed $\mathbf n$ of the field being transported — the discrete
   scheme is quasi-linear where the baseline was linear. The corrugation
   test (rigid translation, tangential and normal, on hex and perturbed
   meshes) must be rerun; acceptance: amplification $\le1$ at
   $\lambda\approx2h$.
2. **Behaviour at $|\mathbf g|\to$ small.** Far from the interface, or at
   level-set skeletons/ridges, $\mathbf n$ is ill-defined. Cells outside
   the narrow band of interest can be frozen or handled by any crude
   update, since only the force band (≈3 cells) feeds the physics — but
   the band/far-field interface needs a defined treatment.
3. **Fallback clustering.** Isolated fallback cells are benign (§5.3);
   whether *contiguous patches* of fallback cells (e.g. around high-κ
   features) reintroduce a tangential $|\nabla\psi|$ signal at the patch
   boundary should be checked in the droplet cases.
4. **δ as projection vs. arc length** (§6): fix one definition, verify the
   $O(\Delta t^{3})$ claim numerically, document.
5. **Two-level normal consistency.** The AB2 combination mixes
   $\mathbf u_n^{\,n}$ and $\mathbf u_n^{\,n-1}$ projected with *different*
   normals. For smoothly evolving interfaces the difference is
   $O(\Delta t)$ inside an $O(\Delta t)$ term (harmless); for rapidly
   rotating interface elements this should be checked against the
   oscillating-droplet case.
6. **3D validation order.** All reasoning above is
   dimension-independent, but the curvature offset machinery it must stay
   consistent with (§5.4) is currently 2D-only in the baseline; nSL tests
   should follow the same 2D-first path.

---

## 10. Minimal test sequence

1. **1D drifted translation** (analytic): $\psi=a x$, uniform $v$ — checks
   the unit-conversion argument end to end; expected exact front speed and
   $|\nabla\psi^{n+1}|=1$ after one step.
2. **Rigid translation of a circle**, tangential and normal separately —
   the tangential case must be an exact no-op modulo renormalization; the
   normal case measures the closed-loop corrugation amplification (§9.1).
3. **2D vortex / shear transport** — geometric shape error vs. baseline at
   matched CFL; expected: comparable or better order, with band
   $\min/\max|\nabla\psi|$ pinned near 1 throughout instead of spreading.
4. **Stationary droplet** — the anti-convergence check: settled
   $\max|\mathbf U|$ at $N=64,128$; the hypothesis to falsify is that
   pinning band $|\nabla\psi|$ restores convergence under refinement.
5. **Oscillating droplet** — the case that currently fails; also exercises
   §9.5.

---

## 11. First measured results (2026-08-03, sigma = 0 rigid translation, N = 64)

Implementation: `slScheme normalProjected` with `renormalization
geometric|strain|none` (normalProjectedScheme.{H,C}); the Sec. 5 conversion is
`slReconstruction::signedOffset`, shared with the curvature path's new
`offsetCorrection quadraticRoot` mode. All numbers below from the sigma = 0
rigid-translation gate (transISTDroplet2D, N = 64, Co = 0.01, 1600 steps,
16 h total displacement; velocity exact free stream to ~3e-15 m/s in every
run, so all effects are pure transport).

**The Sec. 5 geometric write-back is unstable, for both conversion orders.**
Band error after 100 steps: quadratic root 1.11 h, first-order fallback
1.11 h, raw transport (no renormalization) 0.012 h. Writing ANY fit-derived
offset back into psi feeds neighbour noise into the next step's fit at order
one per step -- d(d_c) ~ psi_c d|g|/|g|^2 is not multiplied by dt -- and the
loop compounds at roughly x1.7 per 10 steps. Two further amplification
channels of the same conversion were measured and gated before that
conclusion: at large offsets the root's psi_c h_nn term amplifies
fitted-Hessian noise linearly in the offset (~30 h one-step errors at 45 h
offset, domain-corner cells), and cells whose stencils reach a level-set
skeleton (a circle's cone tip) seed the runaway even at 3.8 h offset -- hence
the shipped gates (root only where |psi_c|/|g| <= 1 stencil radius) and the
raw far field. The `geometric` mode is retained selectable for the record.

**The trace itself is clean.** One full step of the normal-projected update
against the exact translated distance field: 0.017 h in the band, 0.027 h
everywhere. Over the full 1600-step gate the front radius, phase volume,
zero-set error and band |grad psi| match the pointValue baseline
(volume 1.17e-2 vs 1.02e-2; band |grad psi| [0.82, 1.38] vs [0.84, 1.37];
zero-set L2 8.7e-5 vs 7.5e-5 m).

**The Sec. 8 corrugation hypothesis is falsified.** The m > 4 zero-set
corrugation at 16 h displacement: pointValue 0.209 h, normalProjected
0.223 h -- the SAME exponential growth. Removing the value resampling does
not remove the amplifier, because the corrugation re-enters through the
fitted normals, which steer both the trace and the written increment
delta . n with the same per-displacement gain. Open question 9.1 is answered:
the closed-loop amplification of the lambda ~ 2h mode is the same as the
baseline's, not <= 1.

**Status.** Per the pre-registered acceptance rule (amplification <= 1 at
lambda ~ 2h), the scheme in its present form does not pass the go/no-go gate;
the advection ladders and droplet tests were not run. What survives: the
normal-projected trace (exact free-stream preservation, baseline-equal
transport), the strain-renormalization mode (noise gain O(dt) by
construction; inert on this uniform-translation gate since grad(u) = 0 --
untested where it acts), and the shared second-order offset conversion, which
remains valid WHERE THE FIT SAMPLES THE ZERO SET and is available to the
curvature path as `offsetCorrection quadraticRoot`.

## 12. Variant (a) measured (2026-08-03): the instability is engine-independent

Variant (a) is implemented: `offsetEngine footPoint` runs the stabilized
foot-point algorithm (stable-foot-point-3d.md) on the cell's own quadratic
model for BOTH legs of the geometric update -- d_c as the foot-point distance
of x_c to {R_c = 0}, and the normal displacement as the foot-point distance of
x_d to the cell's own iso-surface {R_c = psi_c}, replacing the flat projection
(x_d - x_c).n. All controls are runtime dictionary entries; per-cell fallback
to the ray root when a search fails its guards (counted; measured 0-2 cells
per step out of ~112 band cells, i.e. the engine genuinely ran).

Band error after 100 steps of the sigma = 0 rigid translation (N = 64,
velocity exact free stream in all runs):

| in-band update                      | band error |
|---|---|
| geometric, quadratic-root engine    | 1.107 h    |
| geometric, first-order conversion   | 1.11 h     |
| geometric, stabilized foot-point    | 1.113 h    |
| raw transport (no renormalization)  | 0.012 h    |

Three offset engines of very different accuracy, one divergence rate. The
write-back of ANY fit-derived offset into psi is the amplifier; the engine's
accuracy is irrelevant to the loop gain. On the drifted-profile translation
(psi scaled by 1.3) the foot-point variant ends with front error -0.43 h at
t = 0.0031 -- the renormalization moves the front through the same feedback.

Reversed 2D vortex, N = 64, CFL 0.5, T = 2 (the first basic advection test):

| scheme | shape error | volume error | band min per-cell gradient magnitude at T |
|---|---|---|---|
| pointValue baseline        | 4.29e-4 | 1.56e-2 | 0.98 |
| nSL, geometric + footPoint | 2.28e-2 | 2.12    | 0.00 |
| nSL, strain (full ladder)  | anti-convergent, order -1.18 | up to 4.7 | 0.00 |

Status: the value path (full quadratic resampling at the foot) remains the
most accurate transport measured in this framework; every geometric
substitute tested so far -- including Variant (a) -- is worse on deforming
flows, and every per-step renormalizing write-back diverges regardless of the
offset engine. Next per plan: the third configuration -- the UNCHANGED
baseline value transport with the Sec. 7.1 strain factor applied in the band
as a runtime switch -- which changes exactly one thing relative to the
verified scheme and whose correction carries the O(dt) noise gain.
