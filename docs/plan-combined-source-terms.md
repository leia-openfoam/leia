# plan-combined-source-terms.md

Implementation and gating plan for the zero-level-set-preserving **combined
source term**

$$\frac{\mathrm D\psi}{\mathrm Dt} \;=\; \psi\,F,\qquad
F \;=\; \underbrace{\mathbf n^{\mathsf T}\mathbf D(\mathbf u)\,\mathbf n}_{a}
\;+\;\lambda\bigl(1-|\nabla\psi|\bigr),\qquad
\mathbf n=\frac{\nabla\psi}{|\nabla\psi|},$$

in leia. Working document for a coding agent — **v0.1**.

Source of the model: `docs/combined-source-terms/levelset_combined_source_note.tex`
(technical note, 6 Aug 2026). Companion article:
`docs/combined-source-terms/cst-article/combinedSourceTerms.tex`. Deck:
`docs/combined-source-terms/cst-presentation/combined-source-terms.template.html`.

Scope: the level-set **profile-maintenance** path — how $|\nabla\psi|$ is kept
near 1 in the cells that feed the curvature stencil, without redistancing and
without moving the interface. Out of scope here (separate plans): the curvature
estimator and delivery (`docs/plan-curvature-stabilization.md` WP1–WP3, WP5),
redistancing (retired, §3), velocity extension (§3).

---

## 0. Ground truth — measured facts the agent must not re-litigate

### 0.1 Both halves of the combined source have already been measured here — separately — and each failed

This is the single most important fact in this plan. leia already contains a
runtime-selectable SDPLS source hierarchy
(`src/leiaLevelSet/sdplsSource/`, Reitzel 2022) with exactly the two
constituents of $F$ implemented as **separate** models, and both have been run
on the method-decision benchmark.

`sdplsR`      → $f_{nl}=R=\mathbf n\!\cdot\!\nabla\mathbf u\!\cdot\!\mathbf n = a$
(`sdplsR.C:47`, `sdplsSource.C:137`)
`sdplsBeta`   → $f_{nl}=\beta-|\nabla\psi|$, $\beta$ from the dictionary,
default `SDPLS_BETA 1.0` (`sdplsBeta.C:48`, `cases/default.parameter:52`)

**Eulerian reversed vortex, $N=512$**
(`docs/method-comparison/method-comparison-article/data/tables/benchVortex_decision.tex`):

| run | $T$ | $E_{geom}$ | $E_{vol}$ | **band grad** | clock [s] |
|---|---|---|---|---|---|
| euler (noSource) | 2 | 1.68e-05 | 3.24e-04 | **1.46e-04** | 3.23e+02 |
| euler + SDPLS:R | 2 | 2.56e-05 | 8.09e-04 | **7.87e-01** | 3.27e+02 |
| euler + SDPLS:beta | 2 | 2.40e-05 | 5.93e-04 | **5.16e-01** | 3.29e+02 |
| euler (noSource) | 8 | 3.04e-04 | 2.11e-02 | **7.69e-02** | 1.17e+03 |
| euler + SDPLS:R | 8 | 5.83e-04 | 6.79e-02 | **9.62e-01** | 1.17e+03 |
| euler + SDPLS:beta | 8 | 9.75e-04 | 1.36e-01 | **9.94e-01** | 1.18e+03 |

Both sources made the band gradient error **three to four orders of magnitude
worse** than no source at all, at ~1% extra wall clock. Neither improved the
shape error.

**Further SDPLS negative results** (`sdpls-level-set-negative-results.template.html`,
`sdplsLevelSet.tex` §Negative results, `config/sdplsStability.yaml`):

- `beta` + `explicit` discretization **diverges at every CFL** (0.1, 0.5, 1.0 at
  $N=128$), and the divergence is $\Delta t$-independent: an explicitly lagged
  relaxation source turns its stable fixed point into an anti-damped
  oscillator. Every implicit variant of the same source is stable. **The
  relaxation term must be carried by the implicit diagonal (or by a
  positivity-preserving exponential update).**
- `R` **over-flattens near the resolution limit**: $T=8$, $N=128$,
  $E_{vol}=0.488$ (R) vs $0.238$ (noSource); final band $\max|\nabla\psi|=0.19$.
- `R` **controls steepening but not flattening**: band $\max|\nabla\psi|$
  $5.5\text{–}7.4\to\le1.36$ under R; band $\min|\nabla\psi|\to0$ by $t\approx1$
  **unchanged** under R.

**Semi-Lagrangian form** (`normalProjectedScheme.C:291`, `nslConv2Dvortex`,
`docs/normal-projected-semi-lagrangian/npsl-article/`): the `strain` mode
$\psi^{n+1}=(\psi_c+\delta|\mathbf g|)(1+\varepsilon_{nn}\Delta t)$ — which is
$F=a$ integrated to first order along the characteristic — is **anti-convergent**
on the 2D reversed vortex, $T=2$, Co $=\tfrac12$:

| $N$ | `pointValue` | nSL `strain` |
|---|---|---|
| 32 | 3.05e-04 | 2.26e-03 |
| 64 | 4.39e-05 | 5.14e-03 |
| 128 | 4.61e-06 | 1.29e-02 |
| 256 | 7.85e-07 | 2.53e-02 |
| **order** | **≈2.9** | **−1.24** |

with band $\min|\nabla\psi|\to0$ at every resolution and $E_{vol}$ growing to
$O(1)$. An earlier **unscoped** (whole-domain) application of the same factor
compounded $\exp\!\int\!\varepsilon\,\mathrm dt$ in the far field: band gradient
errors $10^{6}$–$10^{34}$ across the ladder. **The band gate is not optional.**

### 0.2 Why the combination is not "more of the same" — the diagnosis

The measured failures above are exactly the two failure modes the combined
source is constructed to remove, and each half removes the other's:

| source | interfacial law | equilibrium | failure it produces | leia measurement |
|---|---|---|---|---|
| $F=\hat a$ | $\mathrm D^\Sigma_t g = g\,\delta$ | none — pure **integrator** of the residual $\delta=\hat a-a$ | $g=g_0e^{\int\delta}$ drifts without bound, either sign | band grad 7.87e-01 (Eulerian); order −1.24, $\min g\to0$ (SL) |
| $F=\beta-g$ | $\mathrm D^\Sigma_t g=g(\beta-g-a)$ | $g_\ast=\beta-a$ — **velocity-dependent** | with $\beta=1$ and $a\sim O(\pi)$ in the reversed vortex, $g_\ast$ is $O(1)$ wrong wherever the interface strains | band grad 5.16e-01 |
| $F=\hat a+\lambda(1-g)$ | $\mathrm D^\Sigma_t g=g[\lambda(1-g)+\delta]$ | $g_\ast=1+\delta/\lambda$ — **residual-dependent only** | — | **not yet measured — this plan** |

The combination replaces an unbounded integrator by a first-order lag whose
offset is set by the *cancellation residual* $\delta$, not by the strain $a$
itself. For $|\delta|\le\eta<\lambda$ the interval
$[1-\eta/\lambda,\,1+\eta/\lambda]$ is invariant (note, Prop. 6.1). This is a
structural repair of the exact defect leia measured, not a re-parameterisation
of it — and it is the only reason the gates below are worth running.

**It is not a prediction of success.** §1 states what could still kill it.

### 0.3 Noise gain — the acceptance mode, derived

The repo's standing acceptance criterion for any profile-maintenance operator
(`METHOD.md` §9.3) is: *closed-loop amplification $\le1$ for the
$\lambda_{\text{mode}}\approx2h$ mode*. Linearise the frozen-coefficient
exponential source step $\psi^+=e^{\Delta t\,r}\psi$ about a clean local signed
distance in the normal coordinate $\zeta$, $\psi=\zeta+\epsilon\sin k\zeta$,
with $\hat a=0$ and $\chi=1$:

$$g = 1+\epsilon k\cos k\zeta,\qquad r=-\lambda\epsilon k\cos k\zeta,$$
$$\psi^{+} = \zeta+\epsilon\bigl[\sin k\zeta-(\lambda\Delta t)\,k\zeta\cos k\zeta\bigr]+O(\epsilon^{2}),$$
$$g^{+} = 1+\epsilon k\cos(k\zeta)\,(1-\lambda\Delta t)+\lambda\Delta t\,\epsilon k^{2}\zeta\sin k\zeta .$$

Three consequences, all quantitative and all testable offline (WP0):

1. **At the interface** ($\zeta=0$) the mode in $g$ is **damped** by
   $(1-\lambda\Delta t)$. That is the intended relaxation, and it acts exactly
   where the curvature stencil lives.
2. **At the band edge** ($|\zeta|=\varepsilon_b=Wh$, $k=\pi/h$) a quadrature
   component of relative amplitude $\lambda\Delta t\,k\varepsilon_b
   =\lambda\Delta t\,\pi W$ is **created**. The per-step sup-norm gain over the
   band is bounded by
   $$\boxed{\;G_{\text{source}}\;\le\;1+\lambda\Delta t\,\pi W\;}$$
   $k\varepsilon_b=\pi W$ is **mesh-independent**, so this is a rate in time,
   not a per-cell defect. It is also the derived reason for the flat-core
   narrow band and for the measured far-field blow-up of the unscoped strain
   factor (§0.1).
3. **Tangential modes** ($k\perp\mathbf n$) give $\mathbf n\cdot\nabla\varphi=0$,
   hence $r=O(\epsilon^{2})$ and amplification exactly 1. The source is
   **neutral** on the zero-set corrugation mode. *Prediction to be gated:* the
   $m>4$ corrugation figure (`pointValue` $0.209h$ at $16h$ displacement) comes
   out **unchanged**. If it changes, the implementation has a bug or the band
   gate is leaking.

**Design window.** (2) binds against (1):

$$\lambda\Delta t\;\lesssim\;\frac{1}{\pi W},\qquad\text{with }W=3:\quad
\lambda\Delta t\lesssim0.106,\quad \tau_r=\lambda^{-1}\gtrsim10\,\Delta t .$$

The note's illustrative $\lambda\Delta t=0.2$ is already marginal at $W=3$.

**Feasibility condition of the whole method.** Attraction to a useful
equilibrium needs $\lambda>|\delta|_\infty$; the band gain needs
$\lambda\Delta t\lesssim1/(\pi W)$. Together:

$$\boxed{\;|\delta|_\infty\,\Delta t\;\lesssim\;\frac{1}{\pi W}\;}$$

**Measuring $|\delta|_\infty$ therefore decides the method before any coupled
run.** That is WP0.

### 0.4 The multiplicative corrector is not a complementary route — it is the $\lambda\Delta t\to1$ limit

The note (§9) offers $\psi^{n+1}=M\psi^{A}$, $M=g^{\text{tar}}/g^{A}$, as an
alternative that avoids differentiating $\mathbf u$. With the immediate target
$g^{\text{tar}}=1$,

$$\delta\!\left(\frac{\psi}{g}\right)=-\frac{\psi\,\delta g}{g^{2}},$$

which is *algebraically identical* to the amplifier the repo already measured to
diverge — the nSL geometric write-back, whose failure note reads
"$\delta d_c\sim\psi_c\,\delta|\mathbf g|/|\mathbf g|^{2}$ **carries no factor of
$\Delta t$**": divergence $\times1.7$ per 10 steps, independent of $\Delta t$,
and **identical for three offset engines of very different accuracy** (band error
after 100 steps: quadratic root $1.107h$, first-order conversion $1.11h$,
stabilized foot point $1.113h$, raw transport $0.012h$). Setting $\lambda\Delta t=1$
in §0.3 gives $G\approx1+\pi W\approx10.4$ at $W=3$ — the same $O(1)$/step class.

**Consequence for this plan:** the corrector is not a separate work package. It
is the $\lambda\Delta t\to1$ arm of the same one-parameter family (WP6), costs
nothing extra to run, and its predicted divergence is a *falsification test of
§0.3*, which is why it is run at all.

### 0.5 Non-regression gates (inherited verbatim — none of these may regress)

From `docs/plan-curvature-stabilization.md` §0 and `METHOD.md` §8:

| case | mesh | shape order | volume order |
|---|---|---|---|
| 2D reversed vortex, $T=2$, Co $\tfrac12$ | hex 32→256 | 2.97 | 2.98 |
| 2D reversed vortex, Co 1 | hex | 2.59 | 2.94 |
| 3D shear, $T=3$ | hex 32→128 | 2.50 | 3.25 |
| 3D deformation (LeVeque) | hex 32→160 | 1.36 | 1.86 |
| 3D deformation | poly (cfMesh) | 1.46 | 2.51 |

- Curvature statics: 2D circle re-referenced delivery $\ge O(h^{2.04})$;
  3D sphere $\approx O(h^{1.95})$ **with** the Gaussian term (the K-less scalar
  inverse is $O(h^{1.02})$ — never reintroduce it).
- Serial ↔ np=4 to ~1e-12; **np4 ↔ np8** coupled agreement.
- `redistancer { type noRedistancing; }` stays. No band rewrite of any kind.
- Clean-circle curvature 0.477% at 6.4 cells/radius, second order.

### 0.6 What the coupled failure actually looks like (the target)

`docs/plan-curvature-stabilization.md` §6 (WP0 retrodiction, $N=128$,
arithmetic delivery, np 4):

| observable | onset $t$ |
|---|---:|
| `minGradPsiBand < 0.95` (slow drift) | **0.028** |
| $\max|\mathbf U|$ | 0.0655 |
| A2h band $L_2$ | 0.0763 |
| `kErrL2Band` | 0.0784 |
| blow-up (FPE) | 0.0803 |

The **slow $|\nabla\psi|$ drift is the earliest available signal**, leading the
velocity onset by ~37 ms and the blow-up by ~52 ms. The grid-scale mode
explosion is the endgame, not the driver. A source term that holds the band
$|\nabla\psi|$ inside a tube is aimed exactly at the earliest link in the chain
— and, through the parallel-curve offset correction $d=\psi_c/|\mathbf g|$
(`METHOD.md` §4.2), at the *tangential variation* of the delivered $\kappa$,
which is the only part of the CSF force no pressure field can balance
($\nabla\times\mathbf f_\sigma=\sigma\nabla_{\!t}\kappa\times\mathbf n$).

Under transport without redistancing, band $|\nabla\psi|$ is measured spreading
over $[0.84,1.37]$ (`METHOD.md` §9.2).

### 0.7 Machinery that already exists and must be reused, not rebuilt

| what | where | note |
|---|---|---|
| RTS source hierarchy (Eulerian) | `src/leiaLevelSet/sdplsSource/sdplsSource.{H,C}` | `TypeName("noSource")`, `declareRunTimeSelectionTable(autoPtr, sdplsSource, Dictionary, (const dictionary& dict, const fvMesh& mesh), (dict, mesh))`; derived override `nonLinearPart(R, psi, U)` only |
| $a=\mathbf n^{\mathsf T}\nabla\mathbf u\,\mathbf n$ | `sdplsSource.C:137-149` | already computed as `R_` for **every** model incl. `noSource`, registered and written |
| $F=a$ | `sdplsSource/sdplsR.C` | 77 lines total |
| $F=\beta-g$ | `sdplsSource/sdplsBeta.C` | 82 lines total |
| narrow-band cutoff $\chi$ | `sdplsSource/mollifier/mollifier1.C:50` | **flat core**: $\chi=1$ for $|x|<w_1$, Gaussian taper to $10^{-3}$ at $w_2$ — exactly the note's $\chi$ |
| source linearisations | `sdplsSource/discretization/` | `explicit`, `simpleLinearImplicit`, `strictNegativeSpLinearImplicit` |
| dedicated scheme keywords | `gradPsiSdpls`, `gradUSdpls` | never silently coupled to the advection scheme's limited $\nabla\psi$ |
| $\nabla\mathbf u$ of the **trajectory** velocity | `pointValueScheme.C:353` | `const volTensorField gradU(fvc::grad(*trajectoryNew, "gradU"));` — already built every step; $a$ is one contraction away |
| fit $\mathbf g,\mathbf H$ at the cell centre | `slReconstruction.H:281` `fitDerivatives(c, g, H)` | returns fit order (2/1/0) — the degradation ladder for free |
| the fit is rebuilt on the current $\psi$ | `slAdvection.C:80` (`meanCurvature`), `directCorrector.C:60` | so a post-advection $g$ costs at most one extra `update()` |
| $|\nabla\psi|$ diagnostics | `src/functionObjects/gradPsiErrorCSV` | columns incl. `NARROW_MIN_MAG_GRAD_PSI`, `NARROW_MAX_MAG_GRAD_PSI`, `E_NARROW_L2_GRAD_PSI` |
| study axes | `cases/2Dvortex/system/fvSolution.template:180-186`, `config/sdplsStability.yaml` | `SDPLS_SOURCE`, `SDPLS_BETA`, `SOURCE_SCHEME` already exist |
| $a$ in SL form | `normalProjectedScheme.C:311` | `const scalar epsNN = n & (symm(gU) & n);` — sign convention already validated |

**Nothing in this plan requires a new numerical kernel.** The combined source is
one new derived class in an existing hierarchy (Eulerian) plus one new small
hierarchy composed at one call site (semi-Lagrangian).

---

## 1. Problem decomposition

- **A — continuum correctness.** Does the implementation reproduce the logistic
  law $\mathrm D^\Sigma_t g=\lambda g(1-g)$ and the invariant tube? Cheap,
  offline + affine unit tests. Low risk; a sign or transpose error here
  invalidates everything downstream.
- **B — the residual $\delta=\hat a-a^{\text{eff}}$.** DECISIVE. §0.3 shows the
  method exists only if $|\delta|_\infty\Delta t\lesssim1/(\pi W)$. $\delta$
  collects: truncation error in $\nabla\mathbf u$ and $\mathbf n$; the
  velocity–level-set time-level mismatch (leia advects the interface **before**
  the momentum solve, so only $\mathbf u^{n},\mathbf u^{n-1}$ exist); the
  mismatch between the continuous strain formula and the strain the *discrete
  semi-Lagrangian operator actually applies*; band-gate clipping.
- **C — discrete zero-set preservation.** $\psi F$ vanishes on $\psi=0$ at the
  continuum. Discretely, multiplying cell values by a spatially varying
  positive $M_c$ and then re-fitting a quadratic moves the reconstructed zero
  contour by $O(|\nabla M|\,h^{2}/|\mathbf g|)$. Must be measured, not assumed —
  this repo's whole objection to redistancing is interface displacement.
- **D — coupled behaviour.** The actual goal: does a tighter band $|\nabla\psi|$
  reduce the *tangential* variation of the delivered $\kappa$, and does that
  change the $t_{\text{blow}}(N)$ **exponent**? A shifted prefactor is a tuning
  result; a changed exponent is the publishable one.

**Rule for the agent:** every proposed change states which of A/B/C/D it
targets. Changes targeting A alone are rejected after WP1 passes.

---

## 2. Work packages

### WP0 — the residual diagnostic $\delta_h$ and the offline noise-gain harness (do first; everything downstream scores on it)

Targets **B**. This WP can kill the method before a line of solver code is
written, which is its purpose.

**Definition (SL-specific, and computable here — unlike the note's general
$a^{\mathrm{eff}}$).** For the semi-Lagrangian update the effective strain of
the discrete transport operator is directly observable, because the
reconstruction supplies the gradient at both ends of the same characteristic:

$$a^{\mathrm{eff}}_c=-\frac{1}{\Delta t}\,
\log\!\frac{\,g^{A}_c\,}{\,g^{n}(\mathbf x_d)\,},\qquad
g^{n}(\mathbf x_d)=\bigl|\mathbf g_c+\mathbf H_c(\mathbf x_d-\mathbf x_c)\bigr|,$$

where $\mathbf g_c,\mathbf H_c$ are the fit coefficients of $\psi^n$ in the
arrival cell (`fitDerivatives`, exact for the quadratic fit — **zero extra
cost**) and $g^{A}_c$ is the fitted gradient magnitude after the advection
predictor. Then $\delta_c=\hat a_c-a^{\mathrm{eff}}_c$.

Tasks:

1. Function object or in-scheme diagnostic writing, per write time, band norms
   of $\hat a$, $a^{\mathrm{eff}}$ and $\delta$: `A_hat_max`, `A_eff_max`,
   `Delta_max`, `Delta_L2`, plus the clipped fraction. Band = the
   $\alpha\in(0,1)$ + $k$ face-ring band already used by the curvature path.
   Cheap (<1% step cost), halo-only, fit-free beyond `fitDerivatives`.
2. Run it **passively** (source inactive, `type none`) on: 2D reversed vortex
   $N=64,128,256$; 3D shear; the stationary droplet $N=64,128$; the translating
   droplet. Report $|\delta|_\infty\Delta t$ against the feasibility bound
   $1/(\pi W)$ for $W\in\{2,3\}$.
3. Offline noise-gain harness: extend the existing offline stencil machinery
   (`workflow/scripts/make_profile_mode_fig.py`, the one that reproduced solver
   $\kappa$ spikes to 0.00%) with the injector
   $\psi=d(\mathbf x)+\epsilon h\sin(\pi d/h+\phi)$ used by WP1 of the
   curvature plan, and measure the per-step gain of one source step on the 2h
   normal and 2h tangential modes over $\lambda\Delta t\in\{0.025,0.05,0.1,0.2,0.5,1\}$
   and $W\in\{1,2,3,4\}$.
4. **Acceptance:** the measured gain reproduces $1+\lambda\Delta t\,\pi W$ to
   within 20% (validating §0.3 and hence the design window); the tangential-mode
   gain is $1\pm10^{-3}$; parallel-consistent (serial ↔ np4 to ~1e-12, np4 ↔ np8).
5. **Gate:** if $|\delta|_\infty\Delta t>1/(\pi W)$ on the *kinematic* cases,
   stop and report — the method's feasibility condition fails on prescribed,
   noise-free velocity and no coupled run can rescue it. If it passes on
   kinematic cases but fails on the coupled ones, that is the *result*
   (parasitic $\nabla\mathbf u_h$ dominates $\delta$) and it retargets the
   program to the filtered-$\hat a$ arm of WP3.

### WP1 — `leiaTestCombinedSource`: exact unit tests before any solver run

Targets **A**. Standalone `leiaTest*` utility, repo pattern
(`applications/test/leiaTestMeanCurvature`), exits nonzero on mismatch.

1. **Source-only** ($\mathbf u=0$), planar $\psi_0=g_0x$: the exact solution
   stays linear in $x$ and its slope follows
   $g(t)=\bigl[1+(g_0^{-1}-1)e^{-\lambda t}\bigr]^{-1}$. Sweep
   $g_0\in\{0.2,0.6,2\}$, $\lambda\Delta t\in\{0.05,0.2,1\}$; check the
   exponential update reproduces the logistic map to the time-integration order
   and **never changes the sign of a stored value**.
2. **Affine velocity** $\mathbf v=A\mathbf x+\mathbf b$, planar interface,
   constant $\mathbf n$: $a=\mathbf n^{\mathsf T}A\mathbf n$ constant, and all
   four models have closed forms
   ($g=g_0e^{-at}$ / $g=g_0$ / $\dot g=g[(\beta-a)-g]$ / $\dot g=\lambda g(1-g)$).
   Cases: pure normal extension $A=\alpha\,\mathbf n\otimes\mathbf n$; rigid
   rotation $A^{\mathsf T}=-A$ (so $a=0$ — any drift is numerical); simple shear
   with a rotating normal; incompressible planar strain $\operatorname{tr}A=0$.
   **These are the tests that catch sign and transpose errors** in
   $(\nabla\mathbf u)_{ij}=\partial_i u_j$ (OpenFOAM convention) versus the
   note's $n^{\mathsf T}\nabla v\,n$.
3. **Zero-set preservation (C).** Frozen interface, source only, one step:
   measure the displacement of the reconstructed $\psi=0$ contour, per mesh, as
   a function of $\lambda\Delta t$ and $|\nabla M|$. Report against the GRL
   redistancing displacement that got redistancing retired.
4. **Guard branches:** $\chi$ flat core covers the full curvature stencil;
   $\epsilon_n$ regularisation; the fraction of band cells with
   $g_h\le c\,\epsilon_n$ (the note's masking diagnostic).
5. Acceptance: every closed form reproduced to the integrator's order; zero-set
   displacement quantified; nonzero exit on any mismatch.

### WP2 — the Eulerian `combined` model (cheap, one file, validates the continuum claims)

Targets **A**, and gives a controlled comparison against the two measured
halves in §0.1 in the framework where the note's derivation applies verbatim.

1. New RTS model `sdplsCombined` (`TypeName("combined")`) in
   `src/leiaLevelSet/sdplsSource/`, mirroring `sdplsBeta.C` exactly:
   `nonLinearPart(R, psi, U) = R + lambda*(1 - mag(grad(psi)))`, with `lambda`
   a required `get<scalar>` (explicit beats a hidden default — the measured
   `sdplsBeta` construction-failure lesson). ~80 lines. Add to `Make/files`.
   **Do not modify** `sdplsSource`, `sdplsR` or `sdplsBeta`.
2. Template token `SDPLS_LAMBDA` next to `SDPLS_BETA` in every `sdplsSource`
   block of `cases/*/system/fvSolution.template`, default in
   `cases/default.parameter`.
3. Study `config/combinedSourceStability.yaml`, copied from
   `config/sdplsStability.yaml`: `SDPLS_SOURCE: [noSource, R, beta, combined]`
   × `SOURCE_SCHEME: [explicit, simpleLinearImplicit, strictNegativeSpLinearImplicit]`
   × `SDPLS_LAMBDA` chosen from $\lambda\Delta t\in\{0.05,0.1,0.2\}$ at the
   study's CFL × CFL $\{0.1,0.5,1.0\}$, $N=128$, $T=2$.
   **Expected from §0.1:** `combined` + `explicit` inherits the `beta` +
   `explicit` $\Delta t$-independent divergence. That is a *prediction*; measure
   it. `strictNegativeSpLinearImplicit` should be stable by construction.
4. Extend to the benchVortex ladder ($T=2$ and $T=8$, $N$ up to 512) so the
   §0.1 table gains a `combined` row on identical footing.
5. **Acceptance to proceed to WP3:** at some $(\lambda,\text{discretization})$
   in the window of §0.3, `combined` band grad error is **below `noSource`**
   (i.e. $<1.46\times10^{-4}$ at $T=2$, $N=512$) with $E_{geom}$ and $E_{vol}$
   not worse than `noSource`. Anything else and the continuum argument does not
   survive the Eulerian discretisation, which must be understood before the SL
   path is attempted.

### WP3 — the semi-Lagrangian combined source (the production path)

Targets **B**, **C**. This is the work package that matters: the production
solver is `leiaSemiLagrangianLevelSetTwoPhaseFoam`, whose transport is an
**assignment** $\psi^{n+1}(\mathbf x_c)=\psi^{n}(\mathbf x_d)$, not a flux
update. The note's Eulerian FV/FE patterns do not apply; the characteristic form
does, and is simpler:

$$\psi^{n+1}(\mathbf x_c)=\psi^{n}(\mathbf x_d)\,
\exp\!\left(\int_{t^n}^{t^{n+1}}\! r\,\mathrm ds\right).$$

**Architecture (minimal diff, two-object discipline).** New RTS hierarchy
`slSource` in `src/leiaLevelSet/semiLagrangian/`:

```
class slSource
{
public:
    TypeName("none");
    declareRunTimeSelectionTable
    (
        autoPtr, slSource, Mesh, (const fvMesh& mesh), (mesh)
    );
    static autoPtr<slSource> New(const fvMesh& mesh);   // levelSet.semiLagrangian.source
    //- psi holds psi^A (post-advection) on entry, psi^{n+1} on exit.
    virtual void apply
    (
        volScalarField& psi,
        const volVectorField& Utraj,      // the SAME velocity the trajectory used
        const volTensorField& gradUtraj,  // the SAME gradU pointValueScheme built
        slReconstruction& recon,
        const scalar dt
    ) {}
};
```

composed and called by `slAdvection` — the composition root that already owns
`recon_`, `corrector_`, `scheme_` and is **not** itself an RTS model
(`slAdvection.C:60`), so composing there adds no model-class modification:

```
void Foam::slAdvection::advect(volScalarField& psi, const volVectorField& Unew,
                               const volVectorField& Uold)
{
    scheme_->advance(psi, Unew, Uold, recon_(), corrector_());
    source_->apply(psi, ..., recon_(), mesh_.time().deltaTValue());   // default none
}
```

Default `none` ⇒ **every existing case dictionary stays bit-identical** (verify
with a full-ladder regression, the `WP4`/K≡+0 pattern).

Models:

- `combined` — $r=\chi\bigl[\hat a+\lambda(1-g)\bigr]$, applied as
  $\psi\leftarrow e^{\Delta t\,r}\psi$ (sign-preserving for any
  $\Delta t\,r$; forward Euler can flip signs and create spurious zero
  crossings — do not offer it).
- `strainOnly` ($\lambda=0$) and `relaxationOnly` ($\hat a\equiv0$) — the two
  measured halves, on the *unchanged baseline value transport*, closing the gap
  the nSL article left open ("the next configuration in this line applies the
  strain factor on top of the unchanged baseline value transport — one change
  against the verified scheme").
- `multiplicativeCorrector` — WP6; the $\lambda\Delta t\to1$ arm.

**Quadrature arms** (a dictionary word, not separate classes):

| arm | $r$ evaluated at | order in $\Delta t$ | extra cost |
|---|---|---|---|
| `arrival` (**default**) | $\mathbf x_c$, from the fit of $\psi^{A}$ | 1 | one `recon.update()` |
| `departure` | $\mathbf x_d$, from the fit of $\psi^{n}$ via $\mathbf g_c+\mathbf H_c(\mathbf x_d-\mathbf x_c)$ | 1 | **zero** |
| `trapezoid` | mean of both | 2 | one `recon.update()` + feet exposure |

Start with `arrival`: the note's §"Operator splitting" is explicit that
"recomputation after the advection step is important — a single source
evaluation before advection can lag the cancellation by one time step and
increase $\delta_h$", and `arrival` is the post-advection evaluation. Add
`departure` second (free) and quantify the difference — it is a direct measure
of the source's temporal truncation contribution to $\delta$.

Mandatory implementation constraints (each has a measured reason in §0):

1. **Band gate is not optional.** $\chi$ from the existing `mollifier1` flat
   core, or the nSL band criterion $|\psi_c|/|\mathbf g|\le k_{\text{band}}\cdot
   \text{stencilRadius}$. Outside: $r\equiv0$ exactly, no taper leakage.
   Reason: unscoped strain gave $10^{6}$–$10^{34}$.
2. **$\hat a$ from the trajectory velocity's own $\nabla\mathbf u$**, i.e. the
   `gradU` already built at `pointValueScheme.C:353`, not a re-differentiation
   of $\mathbf U$. Reason: the note's checklist item 2 and the nSL precedent
   ($\nabla\mathbf u_n$ assembled analytically, never by differencing).
3. **$\mathbf n$ from the fit** ($\mathbf g/|\mathbf g|$ via `fitDerivatives`),
   never from `fvc::grad(psi)`. Reason: two-object discipline; the geometry path
   is the fit.
4. **Fit-order degradation ladder** for free: `fitDerivatives` returns 2/1/0;
   on 1 use $g$ with $\mathbf H=0$; on 0 **freeze the cell** ($r=0$).
5. **Clipping is a measurement, not a fix.** Any clip of $\hat a$, $r$ or the
   exponent must be counted and reported into $\delta$ (WP0). Silent clipping
   makes the method look stable while destroying the cancellation.
6. **Seam consistency.** Any per-cell quantity feeding a face value must be
   synchronised across processor seams; regression is np4 ↔ np8, not only
   serial ↔ np4.
7. **Guards centered and scale-free** (the absolute-coordinate determinant guard
   already zeroed the $N=512$ band once). Regression test per guard.

Acceptance: bit-identical baseline with `type none`; WP1 tests pass through the
solver path; step-cost delta measured and reported.

### WP4 — kinematic gates (BINDING METHODOLOGICAL RULE)

**Binding rule (user, 2026-08-07, carried from
`docs/plan-curvature-stabilization.md` §2/WP6):** $\psi$ advection is never
modified, suspended or specialised for a static test case. Any candidate must
be machinery that works identically for a deforming, translating interface, and
**promotion requires the moving-interface gates below**, not only a static
droplet score.

1. **2D reversed vortex**, $T=2$, $\sqrt2$ publication ladder
   $N\in\{32,45,64,90,128,181,256\}$, Co $\in\{\tfrac12,1\}$, source ACTIVE.
   Gate: shape and volume orders unchanged within fit tolerance against §0.5.
   **This is where the nSL `strain` mode died at $-1.24$** — it is the decisive
   gate, and it must be run before anything coupled.
2. **3D shear** and **3D deformation**, hex and cfMesh poly. Gate: §0.5 orders.
3. **Rigid translation corrugation gate** ($\sigma=0$, $N=64$, $R=6.4h$,
   $\mathbf U=(0.05,0,0)$): $m>4$ corrugation at $16h$ displacement must not
   exceed the `pointValue` baseline $0.209h$. §0.3(3) **predicts unchanged**;
   a change either way is information.
4. **Free-stream preservation**: $\max|\mathbf U-\mathbf U_{\text{trans}}|$
   stays $\sim3\times10^{-15}$ m/s.
5. **Profile gates** (the actual claim): band $\min/\max|\nabla\psi|$ must
   contract from the measured $[0.84,1.37]$ toward
   $[1-\eta/\lambda,\,1+\eta/\lambda]$ with $\eta=|\delta|_\infty$ from WP0 —
   and the measured tube must **agree with the predicted one**. That agreement,
   not the tube width itself, is the scientific result.
6. **Zero-set displacement (C)**: cumulative source-induced interface
   displacement over the run, separated from transport error by the
   frozen-interface test of WP1.3.

**Gate:** no coupled run until 1, 2 and 5 pass.

### WP5 — coupled capillary gates

Only after WP4.

1. **Stationary droplet**, uniform hex, water/air, $R=1$ mm, $N=64,128,256$,
   **horizon 0.3** (the 0.1 gate is measured too short to distinguish delay from
   closure). Arms: baseline / `strainOnly` / `relaxationOnly` / `combined` at
   two $\lambda$ / the WP3-filtered $\hat a$ variant. Scored on: onset
   ($\max|\mathbf U|>0.1$), $t_{\text{blow}}$, the WP0 A-spectrum,
   `minGradPsiBand`, $\kappa$ band error, $\Delta p$ (pure-phase probe,
   $\varepsilon=10^{-2}$ — **not** the $\alpha=\tfrac12$ partition), phase-volume
   drift, and $\delta_h$.
2. **Translating droplet** (`transISTKinematic`-type) and **oscillating
   droplet**: the moving-interface coupled gate required by the binding rule.
   No degradation relative to the static verdict.
3. **Exponent claim.** $t_{\text{blow}}(N)$, $N=64$–$512$, for the winning
   arm(s), refit against the measured onset law. **A changed exponent is the
   publishable result; a shifted prefactor is a tuning result.** Reference: the
   biharmonic band $\psi$-filter is measured as a ~5× fuse extension at $N=256$
   (blow-up $t=0.167$ vs unfiltered $0.035$) and is explicitly *a delay device,
   not a closure* — the combined source must be scored against that bar, not
   against the unfiltered baseline.
4. **Balanced-force non-regression**: exact constant $\kappa$ on a frozen circle
   still gives $\max|\mathbf U|\sim3\times10^{-11}$ m/s.
5. **Tangential-$\kappa$ diagnostic**: report $\|\nabla_{\!t}\kappa\|$ on active
   faces before/after. This is the mechanistic link the whole plan rests on; if
   the band tube tightens and $\|\nabla_{\!t}\kappa\|$ does not fall, the causal
   chain of §0.6 is wrong and that is the finding.

### WP6 — the multiplicative corrector as the $\lambda\Delta t\to1$ arm (falsification test)

Targets **B/C**. Zero extra implementation beyond a `target` word in the WP3
`combined` model:

- `logistic` (default): $g^{\text{tar}}=\mathcal G_{\Delta t}(g^{n})$, i.e. the
  finite-relaxation target — identical to `combined`.
- `immediate`: $g^{\text{tar}}=1$, $M=1/g^{A}$, i.e. $\lambda\Delta t=1$.

Run `immediate` on the WP4.1 vortex ladder and the WP4.3 translation rig.
**Prediction (§0.4): divergence at the nSL geometric-write-back rate
($\times1.7$ per 10 steps, $\Delta t$-independent).** Confirmation validates
§0.3 and retires the corrector as a separate route; refutation means §0.3 is
wrong and the whole design window must be re-derived. Either outcome is worth
the two runs.

### WP7 — evaluation matrix, article and deck data products

1. Arm table: baseline / `strainOnly` / `relaxationOnly` / `combined`($\lambda_1$,
   $\lambda_2$) / `combined`+filtered $\hat a$ / `immediate` / (reference)
   biharmonic filter — on the WP4 kinematic gates and the WP5 coupled gates.
2. Curated CSVs → `docs/combined-source-terms/cst-article/data/{tables,figures}`
   by script; register the theme in `workflow/scripts/paths.py::_THEMES` as
   `"combined-source-terms": ("combined-source-terms", "cst", "combined-source-terms.template.html")`
   when the first study runs.
3. Non-regression table of §0.5 reproduced for every promoted arm.
4. Article and deck updated from the curated data — never from a log.

---

## 3. Measured dead ends — DO NOT implement (with reasons)

Inherited from `docs/plan-curvature-stabilization.md` §3, plus the ones this
plan adds:

1. **Redistancing of any kind, any cadence.** Plane-based band rewrite is
   $O(h^{2}\kappa)$ **one-signed** and compounds; a redistancer rebuilds
   distance to whatever zero set it is handed, **including advection noise**,
   and the quadratic value fit has no maximum principle
   (`clipToStencilBounds` off). The whole point of $\psi F$ is that it never
   rewrites the band from geometry.
2. **Unscoped (whole-domain) source.** Measured: band gradient errors
   $10^{6}$–$10^{34}$. Flat-core band gate, always.
3. **Explicit lagging of the relaxation term.** `beta` + `explicit` diverges at
   every CFL, $\Delta t$-independently. Use the exponential update (SL) or an
   implicit/strict-negative-$S_p$ split (Eulerian).
4. **The immediate-target multiplicative corrector as a production route.**
   Same $\psi\,\delta g/g^{2}$ amplifier as the nSL geometric write-back, which
   diverged $\times1.7$ per 10 steps identically for three offset engines. Run
   it only as WP6's falsification arm.
5. **Normal extension of $\kappa$** (fetching $\kappa$ at the foot): no static
   gain ($O(h^{1.08})$/$O(h^{1.13})$), 40–85× $\kappa$ amplification in the
   coupled loop. Distance from the foot: yes. Curvature from the foot: no.
6. **Routing the source through the parallel-curve inverse.** In 3D the
   scalar-inverse error would be baked into the master field. The source
   operates on $\psi$ and $|\nabla\psi|$ only.
7. **Slope limiting of the geometry fit.** Barth–Jespersen collapses the order
   $3.0\to0.1$; Venkatakrishnan $3.0\to0.9$. Geometry fit stays unlimited with a
   conditioning fallback ladder.
8. **Velocity extension for pure advection**: worst accuracy at 12–27× cost
   (measured: `euler+VE:closestPoint` 6.63e+03 s vs `euler` 3.23e+02 s at $T=2$,
   $N=512$; 1.61e+04 vs 1.17e+03 at $T=8$).
9. **Basis rotation of the full quadratic fit**: mathematically a no-op.
10. **The K-less scalar parallel-surface inverse in 3D**: $O(h^{1.02})$ at raw
    magnitude on the sphere gate.
11. **$\epsilon_n$ as a substitute for controlling $g$.** Report the fraction of
    curvature-stencil cells with $g_h\le c\,\epsilon_n$; a non-negligible
    fraction means the regularisation is masking a loss of regularity.
12. **`nonLinearPart` returning a dimensionally naked expression.** Follow
    `sdplsBeta.C:72` — the `dimensioned<scalar>(dimless/dimTime, 1.0)` factor is
    load-bearing.

---

## 4. Repo conventions the agent must follow

- New models are runtime-selectable classes in the existing pluggable
  hierarchies; **never modify existing model classes to add behavior**. The
  Eulerian source is a new `sdplsSource` sibling; the SL source is a new
  hierarchy composed at `slAdvection`, whose default keeps every existing
  dictionary valid.
- Standalone physics tests as `leiaTest*` utilities **before** coupled runs.
- All runs config-driven (`config/*.yaml`); every figure/table regenerated from
  curated CSVs by script; write-time-only reporting cadence preserved so wall
  clocks stay comparable.
- Two-object discipline: the transport reconstruction and the geometry fit are
  separate objects with separate requirements. The source reads the **geometry**
  object.
- Degeneracy/conditioning guards centered and scale-free; regression test per
  guard.
- Seam consistency: any face value from rank-local fits must be synchronised
  (swap-and-average, `syncTools`); the regression is np4 ↔ np8 agreement of a
  **coupled** run.
- Every new field/diagnostic: serial ↔ np=4 agreement check.
- Dedicated `fvSchemes` keywords for source gradients (`gradPsiSdpls`,
  `gradUSdpls`); never silently couple to the advection scheme's limited
  $\nabla\psi$.
- Dictionary entries required, not defaulted, where a wrong silent default would
  produce a plausible-but-meaningless run (the `sdplsBeta` construction-failure
  lesson) — and verify an entry is actually **consumed** before building a
  hypothesis on it (the dead `limit` entry cost two studies).
- Build: `source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake`. All
  runs in WSL or on Lichtenberg, never Windows.

---

## 5. Open questions (tracked, not blocking)

1. **Is $\hat a$ worth filtering?** The note recommends evaluating $a$ on the
   reconstructed interface, tangentially filtering the **scalar**, and extending
   it along normals — rather than filtering $\mathbf D(\mathbf u)$ componentwise.
   In leia the interface reconstruction (`levelSetPlaneReconstruction`,
   `detrixheAslamPhaseIndicator`) already exists, so this is affordable; but
   filtering is a delay device in this repo's experience (WP3 of the curvature
   plan) and it moves the clipped part into $\delta$. Decide from WP0's
   $\delta$ decomposition, not a priori.
2. **Time level of $\hat a$.** The interface is advanced **before** the momentum
   solve, so $a$ is built from the AB2-extrapolated trajectory velocity, not
   from $\mathbf u^{n+1}$. The lag contributes to $\delta$; WP0's
   `arrival`-vs-`departure` comparison bounds it. Whether the $O(\Delta t)$ lag
   or the $O(h^{2})$ truncation dominates $\delta$ is unknown.
3. **$\operatorname{div}\mathbf u_h\ne0$ at cell centres.** $a=\mathbf n^{\mathsf T}
   \mathbf D\mathbf n$ needs no incompressibility, but the interpretation
   $a=-\operatorname{div}_\Sigma\mathbf u$ does. Report both and their
   difference as a solenoidality diagnostic.
4. **Should $\lambda$ be spatially varying?** $\lambda=\lambda(\mathbf x)$ only
   enters through $\Lambda(t)=\int\lambda\,\mathrm ds$; a resolution-indicator
   weighting would address the measured `R` over-flattening at $\kappa h\sim1$
   (§0.1) — the SDPLS deck's own open "mollify R by a resolution indicator".
   Do not add it before WP4 measures whether it is needed.
5. **3D anisotropy.** All noise-gain algebra in §0.3 is 1D in the normal
   coordinate. On a 3D body-diagonal stencil the $\zeta$-stations are 7 and the
   effective $k$ differs; whether $G\le1+\lambda\Delta t\,\pi W$ still bounds
   the gain is unverified.
6. **Interaction with the biharmonic band filter.** Both act on the band
   profile. Composed, do they subtract (filter removes the mode the source
   would amplify at the band edge) or compound? WP7 arm, not WP3.
7. **Contact lines.** The continuum source needs no new boundary condition
   merely because it is proportional to $\psi$, but the derivative and band
   stencils do. `cases/2Dcontactline-*` exist; out of scope until WP5 closes.
