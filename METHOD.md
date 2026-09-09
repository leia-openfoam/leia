# Semi-Lagrangian level set for two-phase flow — current best configuration

**THIS FILE IS THE RECORD OF THE BEST CONFIGURATION.** It says what the current best
configuration is, why each choice was made, and which gate measured it. Any gate that
changes a setting updates this file IN THE SAME COMMIT, and every entry names the config
that decided it and the number it was decided on. A setting nobody can trace to a
measurement is folklore, and this file exists to keep it out.

**Its EXECUTABLE form is the `.parameter` layering, not a separate file.** A document
drifts from the code; the token files cannot, because the workflow renders every case
from them. Three layers, later overriding earlier:

1. `cases/default.parameter` — the global default, flat (`SL_FIT normalEquations;`).
   A token no case and no study varies takes its value from here.
2. `cases/<case>.parameter` — the per-case default, inside a `values { ... }` block
   (`SL_RECONSTRUCTION ( uncachedQuadraticWeightedLeastSquares );`). This layer exists
   because some settings are genuinely case-dependent, and it must be used rather than
   silently picking one global winner — see `CURVATURE_EXTENSION` in Sec. 10.
3. `axes_override` in `config/<study>.yaml` — the per-study sweep.

So changing the best configuration means editing a `.parameter` file, and every study
that does not override that axis inherits the change. That also makes it TESTABLE: the
change must be gated exactly like any other, because it moves every inheriting study.

**A `config/best.yaml` base file does NOT work and must not be added.** Snakemake merges
multiple `--configfile` arguments SHALLOWLY, so a study's `axes_override` replaces the
base's entirely instead of merging into it. MEASURED 2026-09-09: with a base setting
`SL_FIT` and a study setting only `N_CELLS`, the rendered case silently lost `SL_FIT`.
A base config would drop axes without warning, which is worse than having none.

Discussion note, 2026-07-31, kept current. Describes the configuration that currently
gives the best measured result, the equations it actually solves, and what is still open.
Numbers are measured, not estimated; the provenance of each is given.

Solver `leiaSemiLagrangianLevelSetTwoPhaseFoam`, library `libleiaLevelSet`,
OpenFOAM-v2512. Collocated cell-centred finite volume, arbitrary polyhedra.

---

## 1. Continuous problem

One-field incompressible two-phase flow, phase indicator $\alpha\in[0,1]$:

$$\nabla\cdot\mathbf u = 0,\qquad
\frac{\partial(\rho\mathbf u)}{\partial t} + \nabla\cdot(\rho\,\mathbf u\mathbf u)
= -\nabla p_{\mathrm d} - (\mathbf g\cdot\mathbf x)\nabla\rho
+ \nabla\cdot\!\left[\mu(\nabla\mathbf u + \nabla\mathbf u^{\mathsf T})\right] + \mathbf f_\sigma$$

$$\rho = \alpha\rho_1 + (1-\alpha)\rho_2,\qquad \mu = \alpha\mu_1 + (1-\alpha)\mu_2,
\qquad \mathbf f_\sigma = \sigma\kappa\,\nabla\alpha$$

Interface transport is the vanishing material derivative of the level set,

$$\frac{\mathrm D\psi}{\mathrm Dt} = \frac{\partial\psi}{\partial t} + \mathbf u\cdot\nabla\psi = 0 .$$

No reinitialisation is performed at any point. $\psi$ is initialised as a signed
distance and is **not** maintained as one.

---

## 2. Semi-Lagrangian transport

Because $\psi$ is constant along characteristics, the update is an assignment:

$$\boxed{\;\psi^{n+1}(\mathbf x_c) = \psi^{n}(\mathbf x_d)\;}$$

evaluated at every cell centre. No flux, no divergence, no linear system.

### 2.1 Departure foot — explicit in the current and previous velocity

The coupled solver advances the interface **before** its momentum solve, so only
$\mathbf u^{n}$ and $\mathbf u^{n-1}$ are available. The backward trace is
therefore expanded about the **departure** time, which is second order from those
two levels alone:

$$\boxed{\;\mathbf x_d = \mathbf x_c
- \frac{\Delta t}{2}\left(3\mathbf u^{n} - \mathbf u^{n-1}\right)
+ \frac{\Delta t^{2}}{2}\left(\mathbf u^{n}\cdot\nabla\right)\mathbf u^{n}\;}$$

The linear displacement uses the time-centred velocity
$\tfrac12(3\mathbf u^n-\mathbf u^{n-1}) = \mathbf u^{n+1/2}+O(\Delta t^2)$,
obtained by Adams–Bashforth extrapolation. Relative to the arrival-time form,
**only the temporal-derivative term changes sign**; the convective term keeps its
plus sign under both centrings.

Implemented at the call site: the kernel is handed
$\mathbf u^{\ast} = \mathbf u^{n} + (\Delta t/\Delta t_0)(\mathbf u^{n}-\mathbf u^{n-1})$
in the slot it documents as $\mathbf u^{n+1}$, with $\mathbf u^{n}$ as the old
level, and reduces algebraically to the boxed expression. The $\Delta t/\Delta t_0$
ratio keeps the difference quotient correct under a varying step.

*Why it matters:* feeding $(\mathbf u^n,\mathbf u^{n-1})$ into the arrival-time
form leaves a $+\Delta t^{2}\partial_t\mathbf u$ foot error and first-order-in-time
transport. Invisible in every steady or uniform test ($\partial_t\mathbf u=0$ makes
the two forms identical), but 2–4% of the per-step displacement early in an
oscillating droplet and 35–47% by $t=0.02$.

### 2.2 Reconstruction at the foot — quadratic weighted least squares

Per-cell fit over the arrival cell's neighbour stencil $\mathcal S_c$, with
$\mathbf d_i = \mathbf x_i - \mathbf x_c$ and the **constant-free** quadratic basis
($m=5$ in 2D, $9$ in 3D):

$$\mathbf b(\mathbf d) = \bigl(d_a,\ \tfrac12 d_a^{2},\ d_a d_b\bigr)_{a<b},
\qquad R_c(\mathbf x_c+\mathbf d) = \psi_c + \mathbf a\cdot\mathbf b(\mathbf d)$$

$$\mathbf a = \arg\min_{\mathbf a}\sum_{i\in\mathcal S_c} w_i^{2}
\bigl(R_c(\mathbf x_i)-\psi_i\bigr)^{2},\qquad w_i = 1/|\mathbf d_i|$$

solved as the SPD normal system $M\mathbf a = \mathbf r$,
$M = A^{\mathsf T}W^{2}A$, $\mathbf r = A^{\mathsf T}W^{2}\Delta\boldsymbol\psi$,
by in-place Cholesky. Dropping the constant makes $R_c(\mathbf x_c)=\psi_c$ exactly.

The coefficient vector splits as $\mathbf a = (\mathbf g,\mathbf h^{\mathrm d},\mathbf h^{\mathrm o})$
with $\mathbf g = \nabla R_c(\mathbf x_c)$ and, because of the $\tfrac12$ on the
diagonal terms, $\mathbf H = \nabla\nabla R_c$ read directly off the fit.

Two properties are load-bearing and both are measured: the fit must interpolate
stencil **values** (gradient-injecting Taylor variants amplify grid-scale error and
diverge), and the degree must be $\ge 2$ (a linear value fit drives
$\||\nabla\psi|-1\|$ to $10^{21}$).

Stencil: cell–point–cell on hexahedra (6 face neighbours are too few for the
9-term 3D fit), cell–face–cell on polyhedra, where face valence already
over-determines it.

---

## 3. Phase indicator

The indicator is recovered **geometrically**, without assuming $|\nabla\psi|=1$.
A local signed-distance plane is fitted over the cell and its neighbours,

$$(\mathbf n_c,d_c) = \arg\min_{\mathbf n,d}\!\!\sum_{k\in\{c\}\cup N_c}\!\!
\bigl(\psi_k - (\mathbf n\cdot\mathbf x_k + d)\bigr)^{2},
\qquad \phi(\mathbf x) = \mathbf n\cdot\mathbf x + d,\ \ |\mathbf n|=1$$

then $\Omega_c$ is fanned into tetrahedra and each tet's $\{\phi<0\}$ fraction is
filled by the closed form of Detrixhe & Aslam, giving
$\alpha_c = \sum_t f_t V_t/\sum_t V_t$. Second order, tolerance-free.

---

## 4. Curvature — cell-centred symbolic fit with the parallel-curve correction

**This is the part that changed most recently, and it is where the current gain comes from.**

### 4.1 Symbolic curvature from the same fit

No re-differentiation of any field; $\mathbf g$ and $\mathbf H$ come from the fit
coefficients:

$$\kappa_d = \nabla\cdot\!\left(\frac{\nabla\psi}{|\nabla\psi|}\right)
= \frac{\operatorname{tr}(\mathbf H)\,|\mathbf g|^{2}
- \mathbf g^{\mathsf T}\mathbf H\,\mathbf g}{|\mathbf g|^{3}}$$

evaluated **at the cell centre**. There is no normal extension: the foot-point
Newton projection is off (`curvatureExtension none`), having been measured to
amplify band curvature 40–85× on a drifted $\psi$.

### 4.2 The offset (parallel-curve) correction

$\kappa_d$ is the curvature of the level contour *through the cell centre*, not of
the interface. For parallel curves at signed normal offset $d$,
$\kappa_d = \kappa/(1+d\kappa)$, so

$$\boxed{\;\kappa = \frac{\kappa_d}{1 - d\,\kappa_d},
\qquad d = \frac{\psi_c}{|\mathbf g|}\;}$$

guarded: fall back to $\kappa_d$ when $|1-d\kappa_d|\le\tfrac12$, i.e. once the
cell centre passes half the local radius of curvature.

$d = \psi_c/|\mathbf g|$, **not** $d=\psi_c$: it is the first-order distance to the
fit's zero contour and needs no $|\nabla\psi|=1$ assumption, which matters because
the band $|\nabla\psi|$ was measured spreading over $0.84$–$1.37$ under transport.
$|\mathbf g|$ is already computed one line above for $\kappa_d$, so the correction
is one divide.

**Why it dominates.** The force band is about three cells thick, so the uncorrected
band norm carries $\overline{|d|}\,\kappa \approx \tfrac34 h/R$ — a *first-order*
error. Predicting that from the band thickness alone, with no free parameters,
reproduces the measured clean-circle band $L_1$ error at every resolution:

| $N$ | cells/$R$ | predicted | measured |
|---|---|---|---|
| 64 | 6.4 | 11.66% | 11.48% |
| 128 | 12.8 | 6.02% | 6.03% |
| 256 | 25.6 | 3.01% | 2.98% |
| 512 | 51.2 | 0.88% | 0.89% |

and it explains why the full estimator and its Laplacian simplification
$\operatorname{tr}(\mathbf H)$ agree to 0.1%: the dominant term depends on *where*
curvature is evaluated, not on which formula is used.

Corrected, on a clean signed-distance circle (`leiaTestMeanCurvature`, band $L_1$
relative to $\kappa=1000$):

| $N$ | cells/$R$ | uncorrected | $d=\psi_c$ | $d=\psi_c/|\mathbf g|$ |
|---|---|---|---|---|
| 32 | 3.2 | 27.56% | 2.31% | 2.11% |
| 64 | 6.4 | 11.48% | 0.457% | **0.477%** |
| 128 | 12.8 | 6.03% | 0.110% | **0.111%** |
| 256 | 25.6 | 2.98% | 0.0275% | 0.0276% |
| 512 | 51.2 | 0.890% | 0.00703% | 0.00703% |

Ratios 4.15, 4.00, 3.91 per mesh halving: **second order**. The previously reported
$O(h^{1.2})$ rate was the offset error masking the fit's own Hessian accuracy, so
the explanation in the article — that $O(h)$ curvature is intrinsic to
differentiating a quadratic value fit twice — is not what the data shows.

The two forms coincide on a clean signed distance, as they must where
$|\nabla\psi|=1$; they separate only once $\psi$ drifts.

### 4.3 Dimensional caveat (open)

The boxed correction is the **2D** parallel-curve inverse. In 3D, with
$\kappa=\kappa_1+\kappa_2$ and Gaussian curvature $K=\kappa_1\kappa_2$, the exact
parallel-surface relation is

$$\kappa_d = \frac{\kappa + 2dK}{1 + d\kappa + d^{2}K}$$

which reduces to $\kappa = \kappa_d/(1-\tfrac12 d\kappa_d)$ under local sphericity.
$K$ is also available from the same fit, so the exact inversion is reachable, but
**only the 2D form is implemented and only 2D cases have been run.**

---

## 5. Surface tension delivery — balanced force

The force is kept as an **integrated oriented scalar face flux**, never a cell
vector:

$$G_{\sigma,f} = \sigma\,\kappa_f\,\mathrm{snGrad}(\alpha)\,|\mathbf S_f|$$

with arithmetic cell→face interpolation of $\kappa$. It enters the pressure
equation on exactly the faces where the pressure gradient acts,

$$\nabla\cdot\bigl(r_{\!A f}\nabla p_{\mathrm d}\bigr)
= \nabla\cdot\bigl(\phi_{HbyA} + r_{\!A f}\,G_{\sigma,f}\bigr),
\qquad
\mathbf U = \frac{\mathbf H}{A} + r_{\!A}\,
\mathrm{reconstruct}\!\left(G_{\sigma,f} - \frac{\phi_{p,f}}{r_{\!A f}}\right)$$

so a spatially **constant** $\kappa$ is a discrete gradient that the pressure
absorbs exactly, leaving zero spurious velocity. Verified: exact constant
curvature on a frozen circle gives $\max|\mathbf U| \sim 3\times10^{-11}$ m/s.

What remains is governed entirely by the **tangential variation** of the delivered
curvature, because

$$\nabla\times\mathbf f_\sigma = \sigma\,\nabla\kappa\times\mathbf n
= \sigma\,\nabla_{\!t}\kappa\times\mathbf n$$

and $\nabla\times\nabla p\equiv 0$, so no pressure field can balance it. The
surviving part is divergence-free and appears as a vortex at the interface.
Discretely: around a closed loop of cells a pressure gradient's face values
telescope to zero, and a tangentially varying capillary force's do not.

---

## 6. Mass–momentum consistency (rhoLENT)

An auxiliary conservative density equation is solved on **every** outer iteration
with the *same* face flux that appears in the momentum convection,

$$\frac{\rho_c^{n+1}-\rho_c^{n}}{\Delta t}
+ \frac{1}{V_c}\sum_f \rho_f^{n+1}F_f^{o} = 0,
\qquad \rho_f = \alpha_f\rho_1 + (1-\alpha_f)\rho_2$$

with $\alpha_f$ from the same reconstructed plane, and after the outer loop
$\rho \leftarrow \alpha\rho_1+(1-\alpha)\rho_2$ so that $\rho^{\text{old}}$ at the
next step is geometrically consistent. The auxiliary $\rho$ is never clipped.
Measured relative mass residual $\sim10^{-13}$, including right up to every crash —
mass transport is not implicated in any failure observed so far.

Interface advanced **once** per step, on the first outer iteration, then held fixed
across the pressure–velocity correctors.

---

## 7. Diagnostics

Pressure jump over the **pure phases**, not either side of $\alpha=\tfrac12$:

$$\Delta p = \langle p\rangle_{\alpha>1-\varepsilon} - \langle p\rangle_{\alpha<\varepsilon},
\qquad \varepsilon = 10^{-2}$$

Partitioning at $\tfrac12$ puts the whole transition into the drop average and
biases the jump low by ~3%; with exact constant curvature and
$\max|\mathbf U|\sim3\times10^{-11}$ m/s the old probe read $70.4456$ Pa ($-3.15\%$)
while the true discrete jump is $72.7394$ Pa ($-0.00\%$). The old metric was also
non-monotone in solver quality, so any previously reported Laplace-jump convergence
over a mesh ladder (e.g. 4.2% → 1.5%) is that artefact, not a solver result.

Also recorded per step: $\max|\mathbf U|$, band $\|\kappa-\kappa_{\text{exact}}\|$,
band $\min/\max|\nabla\psi|$, phase volume, zero-set radial error, mode-2 amplitude,
rhoLENT residual.

---

## 8. Current measured status

**Transport (prescribed velocity, settled).** Geometric shape error converges at
2nd–3rd order to CFL 1, at the *same* rate on hexahedra and cfMesh polyhedra:
2.97 / 2.59 (2D vortex, CFL ½ / 1), 2.95 / 3.28 (3D shear, hex / poly),
1.36 / 1.46 (3D deformation, filament-limited). $20\times$ more accurate than the
best Eulerian line at half the wall clock at $512^2$, $T=8$.

**Curvature (clean circle).** 0.477% at 6.4 cells per radius, second order.

**Stationary droplet, uniform hex, water/air, $R=1$ mm.**

| $N$ | cells/$R$ | settled $\max|\mathbf U|$ | trend | $\Delta p$ error | band $\min|\nabla\psi|$ |
|---|---|---|---|---|---|
| 64 | 6.4 | $4.71\times10^{-6}$ m/s | decaying | exact | 0.989 |
| 128 | 12.8 | $2.90\times10^{-4}$ m/s | decaying | exact | 0.9935 |

$N=64$ corresponds to $\mathrm{Ca}=\mu|\mathbf U|/\sigma\approx6.5\times10^{-8}$,
about $380\times$ below the best published translating-drop figure
($\mathrm{Ca}\approx2.4\times10^{-5}$).

---

### 8.1 Settled decisions, with the gate that decided each (kept current)

Every row is executable: the value is what the named `.parameter` layer holds, so a study
that does not override the axis runs it. "gate" is the config whose measurement decided it.

| token | value | layer | gate that decided it | the number |
|---|---|---|---|---|
| `MOMENTUM_DDT_SCHEME` | `backward` | default | BDF2-vs-Euler matched windows | gain moved +11.1/+2.9/-3.0 % (sign-flipping = noise); volume and shape within 1.2 %. BDF2 costs nothing and is formally right |
| `RHO_DDT_SCHEME` | `backward` | default | as above | — |
| `SL_RECONSTRUCTION` | `uncachedQuadraticWeightedLeastSquares` | per-case | transport ladders | 2nd–3rd order shape error to CFL 1 on hex AND cfMesh poly (Sec. 8) |
| `SL_CORRECTION` | `direct` | default | — | `deferredCorrection` is a research path; no study selects it |
| `SL_TRACE_VELOCITY` | `projectedFlux` | default | `stationaryDropletFootEval*` | the win is the RECONSTRUCT OPERATOR, not solenoidality (STATUS 2026-08-31) |
| `SL_FOOT_INTEGRATOR` | `taylor` | default | — | — |
| `SL_FIT` | `normalEquations` | default | `popinet3D_La12000_poly_dump4_qr` | `householderQR` blows up IDENTICALLY (step-3 phase volume 0.017512193 vs 0.017512208, ~6 significant figures). The amplifying cells are WELL conditioned (min pivot 0.757 at Λ = 1.2608), so there is nothing for better arithmetic to repair |
| `SL_CLIP` | `false` | default | `popinet3D_poly_sigma0_clipGate` (G4) | **the clip is NOT in the best configuration.** With the extremum exemption it fails at step 506 against the unclipped control's 527; without it, it removes the polyhedral failure but costs +30.4 % volume error on the 2D hex translating droplet |
| `SL_CLIP_REGION` | `all` | default | `popinet2D_clipRegionGate` | the band exclusion removes ZERO firings on hex and fires zero cells on the poly stationary rung: no measured benefit anywhere |
| `SL_CLIP_KEEP_EXTREMA` | `false` | default | G4 | see `SL_CLIP`; 59.2 % of the cells the bound must act on are themselves stencil extrema |
| `PSI_FILTER` | `none` | default | filter-off scoring rule | a filter is a research instrument; production must be stable without it |
| `VOLUME_CORRECTION` | `noVolumeCorrection` | default | — | — |
| `MASS_FLUX` | `rhoLENT` | default | consistency studies | — |
| `PHASE_INDICATOR` | `detrixheAslam` | default | — | — |
| `SURFACE_TENSION_FORCE` | `reconstructedCurvature` | per-case | balanced-force gates | — |
| `FACE_CURVATURE_SOURCE` | `model` | per-case | — | — |
| `CURVATURE_EXTENSION` | **case-dependent** | per-case | — | `none` for the Popinet translating family, `cellCentreInverse` for the stationary droplet family. Do NOT collapse these to one global value |
| mesh family | **hexahedral** | — | `sl_fit_amplification.csv` | Λ_max 1.0527 blockMesh / 1.0549 cartesianMesh / 1.0566 snappy+layers / 1.0538 snapped round channel / 1.0265 and 1.0222 at 2:1 refinement transitions, against **1.2608** on pMesh. cfMesh's own HEX mesher on the SAME meshDict and surface is clean, so the defect is polyhedral cells, not cfMesh and not boundary layers |

### 8.2 The transport operator amplifies on EVERY mesh (2026-09-09)

Sec. 9 item 3 asked for "amplification factor ≤ 1 for the λ ≈ 2h mode" as the acceptance
criterion. That number is now measured directly, by
`applications/test/leiaTestTransportSpectrum` (frozen velocity, clip off, so the transport
stage is exactly linear) driven by `workflow/scripts/transport_spectrum_probe.sh`.

| mesh | cells | h [m] | ρ(B) | ρ−1 | ln ρ/Δt [1/s] |
|---|---|---|---|---|---|
| blockMesh hex | 65 536 | 7.81e-05 | 1.00439 | 4.39e-03 | 236.9 |
| blockMesh hex | 221 184 | 5.21e-05 | 1.00441 | 4.41e-03 | 356.6 |
| pMesh | 95 969 | 7.74e-05 | 1.01108 | 1.11e-02 | 595.9 |
| pMesh | 312 975 | 5.16e-05 | 1.01133 | 1.13e-02 | 922.2 |

All at matched Courant 0.0164, 2000 power iterations, increment converged to 0.

Three things follow, and they change how Sec. 9 item 3 must be read.

1. **ρ > 1 on hexahedra too.** Hexahedral meshes are not stable in the ρ ≤ 1 sense; they
   carry a 2.5× smaller exponent. The difference between a run that completes and one that
   fails is the exponent times the horizon: over 1563 steps, exp(1563 × 0.00439) ≈ 9.6e2 on
   hex takes an O(h²) initial error of ~1e-8 to ~1e-5 (invisible), while exp(1563 × 0.011)
   ≈ 3e7 takes it to ~0.3 on pMesh — the observed failure scale.
2. **ρ−1 is INVARIANT under refinement at fixed Courant** (1.004× and 1.023× for a 1.5×
   refinement), but Δt ∝ h, so the growth rate per unit PHYSICAL TIME scales as 1/h:
   measured ratios 1.505 and 1.548 against the exact 1.5. Refining does not make each step
   worse; it buys more steps per second. That is the mechanism behind the recorded 3D wide
   ladder destabilising at R/h = 15.8 while R/h 10.0 and 12.7 are stable.
3. **Λ is a ranking proxy, not a verdict.** Λ > 1 does not imply growth: the anchored
   quadratic at x_i − Ch on a uniform 1D grid IS Lax–Wendroff, with Λ = 1 + C − C² > 1 and
   |G|² = 1 − 4C²(1−C²)sin⁴(θ/2) ≤ 1. The ordering of ρ−1 does track Λ−1 (2.5× against
   4.0×), so Λ remains a cheap one-mesh-pass screen — but only ρ decides.

CAVEAT ON SCOPE: these arms drive the transport with the CELL velocity. Production uses
`projectedFlux`; the projectedFlux arms are the next measurement.

---

## 9. Open

1. **Anti-convergence.** $N=64\to128$ makes the settled current $62\times$ worse.
   The offset correction flipped the $N=128$ trend from growing to decaying and
   improved both resolutions, but did **not** restore convergence under refinement.
2. **$|\nabla\psi|$ drift.** Reinitialisation-free transport does not maintain the
   distance property; band $|\nabla\psi|$ spreads to $0.84$–$1.37$ under transport.
   Through the offset term, a *tangentially varying* $|\nabla\psi|$ produces a
   tangentially varying effective offset and hence tangentially varying $\kappa$ —
   the quantity that drives flow. Across a 10-run matrix, band $\min|\nabla\psi|$
   ordered every outcome and its collapse *preceded* the current growth by ~8 ms.
   Candidate that is **not** redistancing: transport the distance and rescale by the
   normal strain, $d^{n+1} = d^{n}(\mathbf x_d)\,(1+\varepsilon_{nn}\Delta t)$ with
   $\varepsilon_{nn} = \mathbf n\cdot\mathbf D\cdot\mathbf n$, which never rewrites
   the band from geometry and so cannot rebuild distance to advection noise.
3. **Grid-scale corrugation of the transported zero set.** MEASURED as of 2026-09-09,
   Sec. 8.2: the amplification factor this item asks for is ρ(B) = 1.0044 on hexahedra and
   1.0111 on pMesh, so it is ABOVE ONE ON EVERY MESH TESTED and the acceptance criterion
   is not met anywhere. Original text: in an essentially exact
   rigid translation the $\psi=0$ contour develops RMS deviation $0.475h$
   (max $0.985h$) by $t=0.05$, growing exponentially at a rate fixed per unit
   *time* ($\approx0.28$ per cell of displacement), only weakly dependent on CFL
   over $0.007$–$0.5$. This is an amplification property of the reconstruct-and-evaluate
   operator; the acceptance criterion for any fix is amplification factor $\le1$ for
   the $\lambda\approx2h$ mode, not "make the mode small".
4. **No non-oscillatory guarantee, and it cannot be bought with a clip.** MEASURED
   2026-09-09 (G4): `clipToStencilBounds` removes the polyhedral failure but costs +30.4 %
   volume error on the 2D hex translating droplet, because a quasi-monotone bound cannot
   represent the level set's own extrema (the distance-cone apex, the far box corner);
   exempting those extrema restores the interface metrics to +1.03 % but returns the
   failure at step 506. And there is an IMPOSSIBILITY result: applying a nonnegative,
   partition-of-unity, quadratic-exact rule to p(x) = |x − x_d|² forces
   0 = Σ_j a_j |x_j − x_d|² with every term positive. So quadratic exactness and Λ = 1 are
   mutually exclusive at an off-node point, and negative weights are unavoidable.
   `clipToStencilBounds` is off; the quadratic
   value fit has no maximum principle, and hard limiters collapse the convergence
   order (Barth–Jespersen $3.0\to0.1$, Venkatakrishnan $3.0\to0.9$). This is also
   the precondition that makes Eikonal redistancing unsafe here: a redistancer
   rebuilds distance to whatever zero set it is handed, including advection noise.
5. **Moving cases not yet rerun** with the offset correction. Translating and
   oscillating droplets currently fail; those results predate the correction.
6. **3D offset correction** not implemented (§4.3).
7. **Skew meshes.** With exact constant curvature on 10%-perturbed meshes a
   refinement-growing residual remains ($\approx3.5\times10^{-5}$ at $N=64$,
   $7.7\times10^{-5}$ at $N=128$) that is insensitive to the force form, geometry,
   corrector count, $r_{\!Af}$, momentum predictor, solver and tolerance. The one
   operator common to every case and never removed is the `fvc::reconstruct` in the
   velocity correction.

---

## 10. Reproducing the best configuration

`fvSolution`, `levelSet` sub-dictionary:

```
semiLagrangian
{
    reconstruction      uncachedQuadraticWeightedLeastSquares;
    correction          direct;
    offsetCorrection    psiOverGradPsi;   // parallel-curve correction, Sec. 4.2
    trajectoryVelocity  input;
}
surfaceTensionForce { type reconstructedCurvature; faceInterpolation arithmetic;
                      faceCurvatureSource model; forceWeight alpha; }
curvatureExtension  { type none; }        // no foot-point projection
massFlux            { type rhoLENT; faceInterpolation central; }
redistancer         { type noRedistancing; }
```

Capillary step $\Delta t \lesssim \tfrac14\Delta t_\sigma$ with
$\Delta t_\sigma = \sqrt{(\rho_1+\rho_2)h^{3}/(2\pi\sigma)}$, `adjustTimeStep no`.
