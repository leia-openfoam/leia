# Attacking the parasitic-current instability, Shannon-style

Written from measured data only, no new runs. Every number below already exists in
this repository or in this session's measurements.

---

## 0. The evidence base, stated once

**Things that converge cleanly** (these are discretisation quantities):

| quantity | 2D | 3D |
|---|---|---|
| delivered non-gradient content at t=0 | +2.04, R² 0.9987 | +2.03, R² 0.9997 |
| band 2h mode amplitude at t=0 | +2.05, R² 0.9999 | +1.98, R² 0.9999 |
| max&#124;U&#124; at t=0 | +3.58, R² 1.0000 | +3.59 |
| div(U_cell) at t=0 | — | +3.46, +2.74 |

**Things that do not follow a power law at all** (these are growth-integrated):
2D end-state max&#124;U&#124; has local slopes 2.43 → 2.02 → 1.01 and a global fit of
+1.84 at R² 0.970; volume and shape global fits exceed 2, which is formally
impossible. A two-term model `A h² + B h` was fitted and **rejected** (requires
B < 0, misses N=512 by 112%).

**The per-step amplification is constant.** 3D ladder, growth rate of max&#124;U&#124;
fitted over t = 0.05..0.1:

| N_L | R/h | dt [s] | steps | r [1/s] | **r·dt** | e-folds r·T |
|---|---|---|---|---|---|---|
| 60 | 10.0 | 1.0861e-05 | 9207 | 40.47 | **4.395e-04** | 4.05 |
| 76 | 12.7 | 7.6186e-06 | 13126 | 52.32 | **3.986e-04** | 5.23 |
| 95 | 15.8 | 5.4514e-06 | 18344 | 90.55 | **4.936e-04** | 9.05 |

r rises 2.2× under refinement; r·dt is constant to ±24%. Since
`dt = CAPILLARY_DT_COEFF/nRef^1.5`, r ∝ 1/dt, and the e-folds accumulated by a fixed
physical time grow only because a finer mesh takes more steps to get there. The
earlier capillary-dt sweep independently measured ~90% of the growth rate as
dt-proportional at N = 128.

**Eliminated as the amplifier, each by measurement:**

- *The pressure–velocity coupling.* Imposed constant κ holds an m=2 perturbation at
  1.0000 over 3200 steps, Courant 1.2e-10, continuity 1e-12. Read independently
  against interFoam this session: same UEqn, same `phig`, same
  `p_rghEqn.flux() = rAUf·magSf·snGrad(p_rgh)`, same velocity update, `rAUf` cancels
  exactly, `grad(U)ᵗ` present in both, snGrad schemes matched on both fields.
- *Normal-direction ψ corrugation.* The WP0 band spectrum stays at 1.00× of its t=0
  value in every 3D arm **including the unstable one**, while max&#124;U&#124; grows 70×.
  The filter does what it was built to do and the instability proceeds anyway.
- *The Gaussian-curvature term K.* Removing it takes the delivered non-gradient
  content from order +2.03 to **+0.01** and 770× larger, and two arms diverge that
  previously reached t = 0.1. K is load-bearing, not the amplifier.
- *Solenoidality of the advecting field.* div(U_cell) converges at +2.7 to +3.5.
- *|∇ψ| drift.* Established downstream of the currents, not an independent driver.

**What the damage is.** Replacing the force by a pure discrete gradient drops
max&#124;U&#124; from 3.6e-5 to 1.3e-9. `σκ_f·snGrad(α)·|S_f|` is a discrete gradient of the
cell field `σκα` **if and only if κ is constant**, so the non-absorbable part is
`α_f·snGrad(κ_c)`. That is the one quantity the pressure solve cannot absorb, and the
only one with a clean second-order convergence.

**Bugs found this session.** The ψ filter was decomposition-dependent (53–205× worse
np-agreement than the unfiltered control) — fixed, and it **invalidates every
filtered result on record** pending the re-runs. The cell-centre inverse was applied
to 93.5% of cells that the curvature fill never populated — fixed, ~1e-7 effect on
the force. The t=0 Young–Laplace solve was unreferenced on a singular system and
inherited relTol 0.01, leaving the initial pressure ~1% off — fixed. Still open: that
solve uses `laplacian(p_rgh)` with unit coefficient while every step projects with
`laplacian(rAUf, p_rgh)`, and rAUf varies by roughly the density ratio, 839×.

---

## 0b. What the order parameter immediately showed (no new runs)

`workflow/scripts/per_step_gain.py` scores every existing arm on the per-step gain.
The estimator to quote is **gAvg = ln(u_end/u_0)/nSteps** -- no window, no fit, and it
predicts the end state by identity:

    u_end = u_0(h) * exp( gAvg(h) * nSteps(h) )

Verified: 3D N_L = 95 starts at 4.091e-05 with 7.45 e-folds -> 7.06e-02 against
7.04e-02 measured; 2D N = 512 starts at 4.672e-07 with 4.69 e-folds -> 5.08e-05
against 5.10e-05 measured.

Fitting the three factors against h over the ladders:

| | u_0 (initial error) | gAvg (per-step gain) | nSteps | e-folds |
|---|---|---|---|---|
| **2D** | h^+3.58 (R² 1.0000) | **h^+0.77** (R² 0.9975) | h^-1.50 | 1.02, 1.83, 2.87, 4.69 |
| **3D** | h^+3.59 (R² 0.9890) | **h^-1.90** (R² 0.8716) | h^-1.50 | 1.55, 2.60, 7.45 |

**The entire instability is the sign of one exponent.** The initial error converges at
h^3.6 in both dimensions and the step count grows as h^-1.5 in both. The only thing
that differs is the per-step gain: it CONVERGES in 2D (h^+0.8, so refinement wins) and
DIVERGES in 3D (h^-1.9, so refinement loses). That is the 2D/3D difference, stated as
one number.

Caveat on the 3D exponent: three points, R² 0.87, so -1.90 is not tightly determined.
The sign is, though: gAvg rises monotonically 1.688e-4 -> 1.984e-4 -> 4.061e-4.

This also CORRECTS the specification in section 5. The target is not "10x smaller g"
at one resolution -- it is **reverse the sign of the gain exponent**, from h^-1.9 to
positive. A constant-factor improvement in g buys a fixed number of levels and then
loses again.

**Falsifiable prediction for the N_L = 120 arm now running**, from extrapolating the
two fitted exponents: gAvg ~ 6.3e-4, e-folds ~ 16.5, u_0 ~ 1.8e-05, so
u_end ~ 2.7e+02 m/s, i.e. catastrophic divergence. If instead it lands near the
R/h = 15.8 arm's 7e-2, the h^-1.9 gain exponent is wrong, the 3D gain SATURATES, and
the target reverts to lowering a constant -- which would be much better news.

Also visible in the same table: theta = 0.05 has a HIGHER gAvg than theta = 0.2 at
every 2D resolution (5.93e-4 / 1.50e-4 / 1.96e-4 against 2.16e-4 / 1.37e-4 /
7.61e-5), which is the mechanism behind its known non-convergence, now expressed in
the order parameter rather than in survival times.

## 0c. First results under the plan (2026-08-19, post-fix)

### The seam fix reverses the 2D convergence claim

The 2D ladder re-run with the fixed filter, against the contaminated pre-fix numbers:

| N | pre-fix max&#124;U&#124; | POST-fix max&#124;U&#124; | change | pre-fix shape | POST-fix shape |
|---|---|---|---|---|---|
| 64 | 2.2400e-03 | **2.2433e-05** | −99.0% | 3.78e-05 | 1.6252e-06 |
| 128 | 4.1500e-04 | **2.8073e-05** | −93.2% | 5.24e-06 | 1.5318e-07 |
| 256 | 1.0200e-04 | **8.8277e-05** | −13.5% | 8.14e-08 | 5.8863e-08 |

The contamination **fell with refinement** (−99.0%, −93.2%, −13.5%), because at fixed
np = 8 a finer mesh puts proportionally less interface on processor seams. Removing it
therefore lowered the coarse values far more than the fine ones, which turns the
apparent convergence into growth: post-fix max&#124;U&#124; is 2.24e-5 → 2.81e-5 → 8.83e-5,
pairwise orders −0.32 and −1.65.

**The "first configuration to converge in every metric" result was an artefact of the
seam bug.** Shape and volume still converge cleanly (+3.41/+1.38 and +5.00/+4.09);
the parasitic current does not. That is qualitatively the same signature as interFoam
— interface converges, current does not — at four orders smaller magnitude.

In the order parameter the post-fix 2D gain **changes sign** with resolution:
gAvg = −7.606e-04 (N=64), −6.410e-05 (128), **+7.219e-05** (256). The loop damps at
coarse resolution and amplifies at fine, so there is a genuine stability boundary near
R/h ≈ 20 rather than a monotone trend. N = 512 still running.

### The fixed-dt test: prediction falsified, in the useful direction

Same three meshes, one dt = 5.4514e-6:

| N_L | R/h | native dt | gAvg | e-folds | fixed dt | gAvg | e-folds |
|---|---|---|---|---|---|---|---|
| 60 | 10.0 | 1.086e-05 | +1.688e-04 | +1.55 | 5.451e-06 | **+6.656e-05** | **+1.22** |
| 76 | 12.7 | 7.619e-06 | +1.984e-04 | +2.60 | 5.451e-06 | **−2.233e-05** | **−0.40** |
| 95 | 15.8 | 5.451e-06 | +4.061e-04 | +7.45 | 5.451e-06 | +2.984e-04 | +2.93 (53% done) |

The recorded prediction was that all three arms would accumulate the same ~9 e-folds,
so the coarse meshes would destabilise purely by taking more steps. **Falsified, and in
the opposite direction:** at the smaller step the coarse arms became *more* stable, and
N_L = 76 went from +2.60 e-folds to **−0.40**, i.e. from growth to decay, at unchanged h.

At fixed h (N_L = 60), halving the step gives gAvg ~ dt^+1.35 and e-folds ~ dt^+0.35.
So the per-step gain is **not** dt-independent, and the earlier "r·dt constant" across
the native ladder was h and dt covarying rather than a law. Reducing dt genuinely
helps — a temporal-splitting signature — which supports the mechanism but not the
earlier framing.

**The actionable part is the exponent.** e-folds ~ dt^0.35 means halving the
accumulated amplification requires dt smaller by 2^(1/0.35) = 7.2×, i.e. 7.2× the cost.
Brute-force step reduction is therefore **not** a practical fix; removing the splitting
is. That is what the `psiOuterCorrectorsGain3D` gate tests, and it is now the pivot of
the whole plan.

### Post-fix results, and what they retract (2026-08-19)

**Both headline claims from earlier in the session were seam-bug artefacts.** The
psi-filter fix (f83a1ab) changed every filtered number by 84-98%.

3D native-dt ladder, pre-fix against post-fix, same dt per arm and same horizon:

| N_L | R/h | | max&#124;U&#124; | volume | shape | min&#124;grad psi&#124; | gAvg |
|---|---|---|---|---|---|---|---|
| 60 | 10.0 | pre | 1.0035e-03 | 1.6314e-04 | 4.9531e-06 | 0.993619 | |
| | | **post** | **3.6708e-05** (−96.3%) | 1.3599e-04 | 5.8797e-07 (−88.1%) | 0.996648 | **−1.905e-04** |
| 76 | 12.7 | pre | 1.0574e-03 | 1.0514e-04 | 1.8848e-06 | 0.996480 | |
| | | **post** | **2.6153e-05** (−97.5%) | 9.8525e-05 | 3.5640e-07 (−81.1%) | 0.997886 | **−8.346e-05** |
| 95 | 15.8 | pre | 7.0362e-02 | 6.1514e-04 | 1.2074e-05 | 0.858594 | |
| | | **post** | **5.8169e-03** (−91.7%) | 4.2324e-04 | 1.8674e-06 (−84.5%) | 0.984926 | **+2.703e-04** |

Precisely what is retracted and what survives:

- **2D "converges in every metric" — WITHDRAWN.** Post-fix max&#124;U&#124; is
  2.24e-05 → 2.81e-05 → 8.83e-05, i.e. it *grows* (orders −0.32, −1.65). Shape and
  volume still converge (+3.41/+1.38, +5.00/+4.09).
- **3D "destabilises at R/h ~ 16" — the LOCATION survives, the SEVERITY does not.**
  The onset still sits between R/h = 12.7 and 15.8, but the finest arm ends at
  5.82e-03 rather than 7.04e-02 with min&#124;grad psi&#124; 0.985 rather than 0.859. It is a
  mild amplification, not a blow-up. Nothing in the post-fix 3D ladder diverges.
- **The gain changes sign with resolution in BOTH dimensions.** 2D gAvg −7.61e-04,
  −6.41e-05, **+7.22e-05** at N = 64/128/256; 3D gAvg −1.91e-04, −8.35e-05,
  **+2.70e-04** at R/h = 10.0/12.7/15.8. The coarse arms genuinely *damp*. So this is
  a stability boundary in resolution, not a monotone trend — near R/h ≈ 20 in 2D and
  ≈ 14 in 3D.

Cross-validation that the two 3D studies agree: N_L = 95's native dt IS 5.4514e-6, so
that arm appears in both the native and the fixed-dt study, and both give
max&#124;U&#124; = 5.8169e-03 — identical, as designed.

### The gate result: two candidates eliminated for ~128 core-h

`psiOuterCorrectorsGain3D`, quarter horizon, control inside the study:

| N_L | R/h | outerCorr | max&#124;U&#124; end | gAvg | e-folds |
|---|---|---|---|---|---|
| 60 | 10.0 | no | 2.7722e-04 | +1.164e-04 | +0.27 |
| 60 | 10.0 | **yes** | 2.7844e-04 (**+0.4%**) | +1.807e-04 | +0.42 |
| 95 | 15.8 | no | 2.3627e-04 | +3.828e-04 | +1.76 |
| 95 | 15.8 | **yes** | 2.3600e-04 (**−0.1%**) | +4.082e-04 | +1.87 |

Volume, shape and min&#124;grad psi&#124; agree to 3-4 digits, and the delivered non-gradient
content grows by the same factor either way (0.93x vs 0.92x at N_L=60, 2.33x vs 2.33x
at N_L=95). **Prediction (b) confirmed: the lag is exonerated.**

This eliminates TWO mechanisms in one test, because the switch changes both:

1. the **once-per-step force lag** — with it on, the interface pipeline re-runs on
   every outer iteration so sigma*kappa_f is no longer frozen against an implicit
   momentum LHS;
2. the **AB2 velocity extrapolation** — passes >= 2 discard `Uext + (dt/dt0)(Uext -
   UextOld)` and re-advect with the genuine current iterate.

Neither moves the answer by more than 0.4%.

**Consequence for Phase 3.** The semi-implicit / deferred-correction capillary force
(task #24) is the same class of fix — make the force implicit within the step — and
the gate says that class buys nothing here. It should be dropped from the plan unless
some other evidence revives it. That is ~128 core-h spent to avoid a much larger
implementation, which is exactly what the gate was built for.

**What remains.** dt reduction does help (e-folds ~ dt^0.35), so the continuous-time
limit appears stable, but neither of the two temporal candidates the switch tests is
responsible, and the coupling itself was already exonerated by the interFoam diff.
The sub-map gain decomposition (sec. 4) is now the priority: it is the only remaining
instrument that can localise a per-step gain to one arrow of the loop.

## 0d. The amplifier measured bare, and the two-factor law (2026-08-20)

`filterOffAmplifier3D` completed all four arms: the certified `L = 6R` box, `R/h`
= 10.0/12.7/15.8/20.0, `cellCentreInverse` delivery, **psi filter entirely off**,
quarter horizon `T = 0.025 s` so every arm is compared at the same *physical*
time. This is the run that separates the amplifier `A` from the filter's damping
`D`, because with the filter off the measured per-step rate *is* `A`.

| R/h | steps | u_0 (one-step kick) | max&#124;U&#124;(T) | A [1/step] | G = A*n [e-folds] | volume | shape |
|---|---|---|---|---|---|---|---|
| 10.0 | 2302 | 2.1203e-04 | 1.7249e-04 | **−8.965e-05** | **−0.21** | 1.402e-06 | 5.412e-07 |
| 12.7 | 3281 | 7.8220e-05 | 6.5617e-05 | **−5.355e-05** | **−0.18** | 1.232e-06 | 3.285e-07 |
| 15.8 | 4586 | 4.0848e-05 | 1.3838e-03 | **+7.681e-04** | **+3.52** | 1.795e-05 | 3.597e-07 |
| 20.0 | 6511 | 1.8911e-05 | 3.6447e-02 | **+1.162e-03** | **+7.56** | 5.450e-04 | 2.141e-06

**First conclusion: the filter is not the source.** With it removed entirely the
two fine arms still grow, at `A` = +7.68e-04 and +1.16e-03 per step — 3.5 and 7.6
e-folds over a quarter horizon. Whatever amplifies, it sits upstream of the
filter, and `A` changes sign between `R/h` = 12.7 and 15.8, the same place the
filtered ladder's root sits. The stability boundary is a property of `A`, not of
the filter.

### The two-factor law

Every arm's endpoint factorises exactly — this is an identity, not a fit:

```
max|U|(T)  =  u_0(h) * exp( G(h) ),        G(h) = A(h) * n(h) = A(h) * T/dt(h)
```

`u_0` is the first row of the metrics CSV: the velocity the interface acquires in
**one step** from the delivered curvature error. It should scale as the product of
the force error and the step over which it acts,

```
u_0  ~  (h^2 delivered-kappa error) * (h^1.5 capillary step)  =  h^3.5
```

and it does, in both dimensions and without tuning:

| pair | 3D, filter off | 2D, BDF2/upwind |
|---|---|---|
| coarse→mid | h^+4.22 | h^+3.62 |
| mid→fine | h^+2.91 | h^+3.51 |
| fine→finest | h^+3.30 | h^+3.69 |

The 2D triple (3.62, 3.51, 3.69) is tight enough to call `u_0 ~ h^3.5` measured.
So the *initial condition* of the instability converges at third-to-fourth order —
better than the h^2 delivery, because the capillary step shrinks too.

`G` is the exponent, and it goes the other way:

| | exponent of G in h | source |
|---|---|---|
| 3D, filter OFF | **h^−3.27** (G: +3.52 → +7.56) | `filterOffAmplifier3D` |
| 2D, BDF2/upwind | **h^−0.95, then h^−0.56** (G: +1.30 → +2.51 → +3.69) | `upwindConvection2D` |

**This is the sharpest statement of the problem so far.** Under refinement the
kick gets ~11x smaller per doubling while the number of e-folds acting on it gets
~1.5–2.1x bigger. Whether the parasitic current appears to *converge* or to
*diverge* is a race between a power law and an exponential, and the race is
decided by one number — the exponent `p` in `G ~ h^-p`:

```
d ln max|U|(T) / d ln h  =  3.5  −  p * G       currents still converge while  G < 3.5/p
```

- 2D, `p ~ 0.75`: the budget is `G < 4.7`. At `N = 512` and T = 0.025 s, `G` =
  3.69 — **just inside it**, which is exactly why the 2D ladder still shows
  max&#124;U&#124; falling (1.80e-05, the smallest on the ladder) while its growth rate is
  positive. Extrapolated to the full 0.1 s horizon at the same rate, `G ~ 14.8`
  and the sign of the trend reverses.
- 3D, `p = 3.27`: the budget is `G < 1.07`, and it was spent between `R/h` = 12.7
  and 15.8. Past it, refining by 1.26x makes the endpoint **26x worse**.

So 2D and 3D do not differ in *whether* the mechanism is unstable — both have `A
> 0` at fine h. They differ in `p`, by a factor of ~4.3. **That factor, not the
sign of `A`, is the quantity a fix has to change.** Any candidate can now be
scored on `p` at a quarter horizon instead of on `t_blow` at a full one, which is
4x cheaper per data point and immune to the e-fold-count spread that made
`t_blow` a poor score.

### The filter's damping has a sign flip

Matched against the θ = 0.2 arms of `psiOuterCorrectorsGain3D` (same box, same
horizon, same delivery) — and the initial kicks agree to 4 digits, 2.120e-04 vs
2.118e-04 at `R/h` = 10.0 and 4.085e-05 vs 4.061e-05 at 15.8, so the two arms
start from the same state and the endpoint ratio is the filter's whole effect:

| R/h | max&#124;U&#124;(T) filter OFF | max&#124;U&#124;(T) θ = 0.2 | G off → G on | filter verdict |
|---|---|---|---|---|
| 10.0 | 1.7249e-04 | 2.7722e-04 | −0.21 → **+0.27** | **1.61x WORSE** |
| 15.8 | 1.3838e-03 | 2.3627e-04 | +3.52 → **+1.76** | **5.86x BETTER** |

At `R/h` = 15.8 the filter removes half the growth rate. At `R/h` = 10.0, where
the amplifier is absent, it *creates* growth from a damped state. The biharmonic
band filter is therefore not a pure sink: it is a sink whose strength tracks the
corrugation content plus **its own weak source**, and the source is what is left
over where there is nothing to remove.

This answers the standing question of whether θ must be resolution-dependent:
**yes, and not as a tuning convenience.** At fixed θ the filter's own source term
sets a floor that dominates wherever `A < 0`, so θ must vanish with the
corrugation content it is meant to damp. It also explains the non-monotone θ
sweep — two terms of opposite sign scaling differently in θ cannot give a
monotone response.

### BDF2 and the convective scheme, on a matched window

`upwindConvection2D` ran the full 2D ladder twice, `upwind` against
`linearUpwind gradU`, both under BDF2. The two schemes agree to 4–5 significant
figures at every resolution — `u_0` identical to 5 digits, max&#124;U&#124; within 0.8% at
`N` = 64 and within 0.06% for `N >= 128`, volume and shape within 0.01%. **The
momentum convective scheme is inert on the stationary droplet**, now on a 4-point
2D ladder as well as the 3D pair. (This is the measurement that closed the
∇U-singularity hypothesis for the *stationary* case; it says nothing about
translating or oscillating interfaces, where the hypothesis stands.)

BDF2 vs Euler, re-scored **on a matched window** — the same step count, 1179 /
3333 / 9428, against `cellCentreInverseFilteredPostFix`:

| N | g Euler/linUpwind | g BDF2/upwind | Δg | Δmax&#124;U&#124; | Δvolume | Δshape |
|---|---|---|---|---|---|---|
| 64 | −2.600e-03 | −2.312e-03 | **+11.1%** | +40.5% | +1.1% | +1.2% |
| 128 | +3.790e-04 | +3.899e-04 | **+2.9%** | +3.7% | +0.4% | +0.4% |
| 256 | +2.744e-04 | +2.661e-04 | **−3.0%** | −7.5% | +0.4% | −0.2% |

Within ±3% and **sign-flipping** for `N >= 128`. `u_0` is identical to 4 digits
(OpenFOAM's `backward` starts up on Euler, so the first step is by construction
the same). The momentum time discretisation is not the amplifier, in either
dimension. Keeping BDF2 as the default is justified on formal grounds — the SL
foot-point trace is second-order in time and should not be fed a first-order
velocity — and it costs nothing measurable here.

> **A retraction of my own arithmetic.** An unmatched comparison of these same
> arms made BDF2 look 3.1x worse than Euler at `N` = 512. That was the horizon,
> not the scheme: `gAvg` is an endpoint estimator, the growth rate decays over the
> run, and a quarter-horizon window therefore reads higher than a full-horizon
> one for the *same* trace. Matched windows give ±3%. Any two arms compared on
> `gAvg` must have equal step counts.

### 2D volume error converges at fourth order

Worth recording because it is the one metric that is unambiguously good, and
because a single-metric view of this ladder would be misleading in the opposite
direction from usual — max&#124;U&#124; is non-monotone here while volume is clean:

| N | volume rel. error | order | shape (zero-set L2) | order |
|---|---|---|---|---|
| 64 | 3.059e-03 | — | 1.219e-06 | — |
| 128 | 1.442e-04 | **4.41** | 1.520e-07 | 3.00 |
| 256 | 8.933e-06 | **4.01** | 5.846e-08 | 1.38 |
| 512 | 4.978e-07 | **4.17** | 1.296e-08 | 2.17 |

## 0e. Popinet (2009) reframes the object we have been measuring (2026-08-23)

Popinet, *An accurate adaptive solver for surface-tension-driven interfacial
flows*, JCP 228:5838-5866 (the Gerris/Balsilisk surface-tension paper), section
6.2.1 and Figs 7-11. Read in full; the relevant content is not the height-function
curvature but his **interpretation of the stationary-droplet test**, which we had
not adopted and which our own data now confirms quantitatively.

### What he establishes

Balanced-force CSF plus height-function curvature recovers exact equilibrium --
velocity to round-off -- **independently of spatial resolution**, PROVIDED the
interface is given time to relax to its *numerical* equilibrium shape. The
initial transient is a train of oscillations of period ~0.4 T_sigma whose
wavelength, from T_lambda = sqrt(rho lambda^3/(pi sigma)) ~ 0.4 T_sigma, is
lambda ~ 0.8 D = pi D/4 -- consistent with the FOURFOLD symmetry of both his
velocity field (his Fig 9) and his curvature error along the circle (his Fig 6).

His reinterpretation, in his words: these oscillations are *not* parasitic
currents, they are "the result of a physically-consistent numerical solution of
the evolution of a perturbed initial interface shape", the perturbation being that
the exact circle imposed as initial condition is **not** the discrete equilibrium.
The resulting capillary waves are then "exponentially damped by viscosity -- again
in a physically-consistent manner -- on timescales of order T_v". True parasitic
currents, by contrast, "are continuously fed energy by an unphysical unbalance, do
not vanish with time". At La = infinity the oscillation neither grows nor decays:
constant amplitude, his Fig 8 top curve.

Two supporting measurements: at numerical equilibrium the **standard deviation of
the curvature along the interface** saturates at O(1e-11), i.e. round-off (his
Fig 10), which is condition 2 of his balanced-force argument; and the numerical
equilibrium shape converges to the exact circle at second order (his Fig 11).

And a methodological warning we have been violating: "test cases designed to
evaluate the accuracy of a given surface-tension implementation (for the
stationary droplet problem) need to make sure that simulations are run for
timescales comparable to T_v ... For shorter timescales the results obtained will
only reflect the accuracy of the initial curvature estimation."

His temporal scheme matters too. The volume fraction is **staggered in time**, c
at n +/- 1/2 against u at n, and the capillary force is evaluated at **n + 1/2** --
formally second order with no iteration. The explicit stability limit is
dt <= sqrt(rho h^3/(pi sigma)).

Finally, his *translating* droplet does NOT converge in time: the velocity norm
grows (his Fig 12) because advection error continually re-perturbs the shape, and
both velocity and shape errors are only ~first order (Figs 13, 14).

### Our oscillation is the m = 8 capillary wave

| | our stationary droplet |
|---|---|
| T_sigma = sqrt(rho D^3/sigma) | 10.48 ms |
| T_v(water) = D^2/nu | 4.0 s = 382 T_sigma |
| our horizon 0.025 s | 2.39 T_sigma = **0.6% of T_v** |
| our dt / Popinet's limit sqrt(rho h^3/(pi sigma)) | **0.164 at every R/h** |

The measured oscillation period of max&#124;U&#124; in the filter-off 2D twin, from an FFT
of the detrended log trace, is **resolution independent** across a factor three in
h -- 1.5626, 1.4704, 1.6667, 1.5626, 1.5625 ms at R/h = 10.0, 12.7, 15.8, 20.0,
31.7 -- against **1.455 ms predicted for m = 8** from Popinet's T_lambda with
lambda = 2 pi R/m. That is 7%, and inverting the measured period gives m = 7.6.
Independently, the azimuthal spectrum of &#124;U&#124; sampled on the ring r = 1.06 R is
dominated by **m = 8** at every resolution, with m = 4 and 16 harmonics at 16-22%
of the mean ring velocity, while the grid-scale mode number pi R/h runs 31 -> 79
over the same ladder. Two independent measurements of the same integer.

So the oscillation we have been calling an instability is the same object Popinet
identifies -- a capillary wave excited because the exact circle is not the discrete
equilibrium -- one harmonic higher than his (m = 8 against m = 4).

### The defect, restated

The physically correct envelope for that wave is viscous decay at
2 nu_water k^2 = **128 1/s**, i.e. 3.2 e-folds of DECAY over our 25 ms horizon.
The measured envelope rates of the same arms are

  R/h    10.0    12.7    15.8    20.0    25.0    31.7
  1/s    -5.3   -96.9    +5.7   +74.4  (+212)  +120.3

(the R/h = 25 entry is untrustworthy: its FFT returned the window length rather
than a period, so its detrending is wrong). Only R/h = 12.7 is near the physical
value. **The defect is therefore numerical ANTI-DAMPING of a physical capillary
oscillation, of order +100 to +250 1/s against a physical -128 1/s -- not a
spurious current fed by force imbalance.** That is a different quantity to fix,
with a different target: match 2 nu k^2, rather than drive max&#124;U&#124; to zero.

One contribution is quantifiable immediately. For an oscillator, a capillary force
evaluated explicitly at t^n gives amplification &#124;1 + i omega dt&#124;, i.e. anti-damping
+omega^2 dt/2 per unit time; with omega = 2 pi/1.5626 ms = 4021 rad/s that is
**+88 1/s at R/h = 10 falling to +22 1/s at R/h = 25**. Right sign, right order,
and exactly what Popinet's n + 1/2 staggering removes -- centring the force between
the two velocity levels makes the amplification factor exactly one, leaving
viscosity to do the damping.

It also explains the standing null result on `psiOuterCorrectors`. Iterating the
force toward n + 1 does not make the scheme neutral, it swaps anti-damping for
damping of the same O(omega^2 dt^2) magnitude -- which is 0.1-0.4% per step, the
size we measured. Neutrality needs the force at n + 1/2, not at n + 1.

Note also that because m = 8 is fixed and h-independent, omega is h-independent
while dt ~ h^1.5, so this particular anti-damping FALLS under refinement. That
predicts eventual restabilisation at fine enough h, and the 2D twin's per-step
gain peaking at R/h = 25 and falling at 31.7 (h^+3.69) is consistent with it.

### What we are missing, in order of cost

1. **The standard deviation of the delivered curvature along the interface.** We
   record `kErrL2Band` and `kErrLinfBand`, both errors against 1/R, which mix the
   absorbable uniform offset with the non-absorbable variation -- and the balanced
   force consumes only the variation. Popinet's Fig 10 plots exactly the variation
   and shows it must fall to round-off. This is the cheapest high-value metric we
   do not have.
2. **One arm run to a viscous timescale.** At 0.6% of T_v every number in this
   campaign describes a transient. A 2D R/h = 25 arm to t = 0.5 s is 48 T_sigma
   and 12% of T_v, 182k steps at 22500 cells -- affordable on a laptop. It answers
   the only question that matters: does the envelope turn over, i.e. does a
   numerical equilibrium exist for our delivery?
3. **Time-stagger the capillary force to n + 1/2.** Structural, and the one change
   with a derived reason to expect neutrality rather than a smaller prefactor.
4. **Test the fixed-point question directly** (already raised in sec. 5 of this
   plan, never measured): can any shape make our delivered kappa exactly constant?
   Popinet's height function has that property by construction; a quadratic WLS fit
   composed with the parallel-surface inverse may not.
5. **Set translating/oscillating expectations to first order.** Popinet's
   translating droplet grows in time and converges at ~first order; ours should not
   be held to second.
6. **Azimuthal spectrum of the delivered face-curvature error on the exact
   circle.** Popinet's height function gives a clean m = 4 pattern; we excite m = 8.
   If the static gate reproduces m = 8 in the delivered error, that is a specific
   signature of our delivery and it identifies which omega the anti-damping acts on.

### RETRACTION: our force is at n+1, not n, and the scheme is already symplectic

The mechanism attributed above -- "an explicitly evaluated capillary force
anti-damps at +omega^2 dt/2, which Popinet's n+1/2 staggering removes" -- is
WRONG, and so are the numbers derived from it (+22..88 1/s, the psiOuterCorrectors
explanation, and the prediction table originally written into
config/viscousHorizon2D.yaml).

Reading the solver settles it. `slAlphaEqn.H` advects psi IN PLACE, psi^n ->
psi^{n+1}, and the narrow band, phase indicator, alpha and the entire curvature
chain are rebuilt from the ADVECTED field before `UEqn.H:8` calls
`fSigma->faceSurfaceTensionForceFlux()`. The force is therefore built from
psi^{n+1}. For the linearised capillary oscillator that is

    eta^{n+1} = eta^n + dt u^n,      u^{n+1} = u^n - dt omega^2 eta^{n+1}

i.e. SYMPLECTIC EULER, whose amplification matrix has det(M) = 1 exactly, so
|lambda| = 1 while the eigenvalues remain complex -- which holds for
omega dt < 2. Computed per step:

| omega dt | (a) force at t^n | (b) force at t^{n+1} = OURS | (c) force at t^{n+1/2} = Popinet |
|---|---|---|---|
| 0.500 | 1.118034 | **1.000000** | 1.000000 |
| 1.047 | 1.447829 | **1.000000** | 1.000000 |
| 1.571 | 1.862268 | **1.000000** | 1.000000 |
| 2.000 | 2.236068 | **1.000000** | 1.000000 |
| 3.142 | 3.297296 | 7.743015 | 7.743015 |

Only (a), the fully explicit force, amplifies unconditionally. And (b) and (c)
are SPECTRALLY IDENTICAL: both have det = 1 and the same trace 2 - omega^2 dt^2,
hence the same eigenvalues, the same amplification AND the same phase error per
step. Centring the force at n+1/2 cannot by itself change linear stability, the
frequency, or the step limit. It changes the offset between the velocity and
interface samples -- the accuracy of the coupling -- and it makes our arrangement
match Popinet's.

WHAT SURVIVES. The measurements: the oscillation is the m = 8 capillary wave
(period and azimuthal spectrum agree independently), the physical envelope is
-128 1/s and we measure up to +212, the mode is mesh-locked, 83% of the curvature
error is non-absorbable variation, and our horizon is 0.6% of T_v.

THE MECHANISM THAT REPLACES IT. A symplectic integrator conserves only for a
force derived from a potential. Our delivered curvature carries an error that
depends on WHERE THE INTERFACE SITS RELATIVE TO THE MESH -- exactly the measured
m = 8 mesh-locked pattern -- so the force is not a gradient of any potential and
does net work around each oscillation cycle. Substituting
omega^2 -> omega^2[1 + eps cos(m theta(eta))] into the symplectic-Euler map
reproduces slow energy drift of either sign, where eps = 0 conserves. That is the
quantity a fix has to drive to zero: THE WORK DONE PER CYCLE BY THE
POSITION-DEPENDENT CURVATURE ERROR, not the time level of the force.

The step-limit measurement stands on its own and is unaffected: the practical
threshold is omega_grid dt ~ 1.0-1.3 against the omega dt < 2 bound of the linear
analysis, and Popinet's Eq. (18) is omega_grid dt = pi, already past it.

## 1. Cut it down

Eleven things are being tracked. **Two matter.**

Being tracked, and demonstrably not the order parameter:

| tracked | why it can be dropped |
|---|---|
| `t_blow` | already abandoned: e-fold count at blow-up varies 5.4–13.3 across the matrix |
| end-state max&#124;U&#124; as *the* metric | not a discretisation error; it is a growth integral, and its "order" is a category error |
| min&#124;∇ψ&#124; | measured downstream of the currents (it leaves 0.95 at t=0.087, the current passes 5e-3 at t=0.068) |
| A2h/A4h/A8h band spectrum | flat at 1.00× in the unstable 3D arm — currently uninformative there |
| div(U) vs div(φ) | converges at +2.7 to +3.5 |
| transport order (linear vs quadratic) | measured flat to 3–4 digits coupled; closed by the Δt^1.5 temporal cap |
| the K term | exonerated and load-bearing |
| 27 curvature deliveries × 11 extension branches | an option space, not a mechanism |
| filter strength θ | the damping knob, not the driver |
| volume and shape error | mandatory to report, but not the instability's order parameter |

The two that matter:

1. **g = r·dt, the per-step amplification** — dt-normalised, dimension-agnostic,
   measured constant to 24% across the 3D ladder, and dt-proportional per the
   earlier sweep. This is the order parameter.
2. **‖α_f·snGrad(κ_c)‖, the non-absorbable force content** — the only thing the
   projection cannot remove, and the only quantity with a clean +2.0 convergence.

Everything else is downstream, exonerated, or a knob. **Score on those two only.**

---

## 2. Look at problems already solved

Three neighbouring problems are solved, and this session measured two of them.

| solver | stable? | parasitic current converges? | why |
|---|---|---|---|
| interFoam | yes, everywhere | **no — it grows**: 2D 0.310 → 0.444 → 0.662 (orders −0.52, −0.57); Laplace jump 14% low and not improving | α bounded in [0,1] by MULES; `rhoPhi` built from the *same limited flux* that advanced α, so `ddt(ρ)+div(ρφ)=0` discretely |
| interFlow (geometric VOF, PLIC-RDF) | yes | volume exact to 12 digits; max&#124;U&#124; ≈ half interFoam's at N=64 (early, 20 steps) | bounded PLIC interface; curvature from a **reconstructed surface** — our premise |
| ours | no beyond R/h ≈ 16 in 3D | delivered force converges at +2.03 | unbounded ψ that drifts, no conservation law anchoring it |

**What the solved ones have in common: a bounded interface representation that cannot
drift, combined with whatever curvature you like.** interFoam pays for stability with
a non-convergent curvature. interFlow keeps a good curvature *and* a bounded
representation — and is the closest published analogue of what we want. The
balanced-force refined-level-set-grid line is the other precedent.

Shannon's two-small-jumps: stop trying to jump from "semi-Lagrangian level set with a
fitted curvature" straight to "stable and second order". Jump twice — (i) bound or
anchor the representation, (ii) keep the fitted curvature. That is literally the
interFlow architecture, and we now have it built and running for comparison.

---

## 3. Say the question a different way

The sunk cost is explicit: the campaign has produced 27 curvature deliveries. The
best of them, `cellCentreInverse`, achieved exactly what it was designed for — the
non-gradient content converges at +2.03 where every other cell field is +0.09 — **and
it did not buy stability.** That is the signal that the question is wrong, not the
answer.

Three reframings, in increasing sharpness:

- **A.** From "which delivery makes the non-gradient content converge?" (answered) to
  "what is the per-step gain of the closed loop, and what sets it?" This is an
  amplification question, and the loop-spectrum power-iteration harness already
  exists.
- **B.** "Why does a force converging at h² drive a response that does not converge?"
  Because the response accumulates over `N_steps ∝ h^-1.5` at a fixed per-step gain.
  So the real question is: **why is the per-step gain h-independent?** An
  h-independent per-step gain cannot come from a spatial truncation error. It comes
  from something that does not refine: an operator splitting, a once-per-step lag, a
  fixed relative error like the 1% initial pressure, or a solver tolerance.
- **C, the sharpest.** The exact solution is U = 0 for all time. So ask: **is the
  discrete equilibrium a fixed point of the one-step map at all?** If not, the
  per-step residual *is* g and everything else is bookkeeping. This is a static
  question — answerable in one step, not in 18,344.

Reframing B has an immediate consequence that **reverses a recommendation I made
earlier in this session**: I proposed building the tangential-structure diagnostic as
the next thing, on the grounds that K was exonerated and tangential structure was the
last dimension-specific candidate. But if the gain is h-independent and
dt-proportional, the driver is *temporal*, and a spatial mode diagnostic is the wrong
next instrument. Deprioritise it until the fixed-dt result is in.

---

## 4. Break it up and let yourself wander

The step is a composition of five sub-maps:

```
ψ ──fit──> κ_c ──deliver──> σκ_f ──force+project──> (φ, U) ──SL advect──> ψ
                                                        └──filter──┘
```

Measure the gain of **each sub-map in isolation** on the same perturbation, rather
than the gain of the whole loop. The harness for this exists (the quasi-static
power iteration). Sub-map gains are cheap — tens of steps, not tens of thousands —
and they localise the amplification to one arrow.

The wandering has already paid once and the result is under-used: the **gauge-freedom
measurement**. A δψ that vanishes on Γ is exactly neutral, power iteration converges
to λ = 1. That says the amplification lives entirely in the component of δψ that
*moves the zero set* — the normal displacement at the interface — and not in the
profile away from it. That is a large narrowing of the search space obtained from an
apparently unrelated poke, and it argues for decomposing δψ into (displacement of Γ)
⊕ (gauge) and measuring the gain of the displacement component alone.

---

## 5. Flip it

Assume the method **is** stable and second order at R/h = 40 in 3D. What must be true?

Working backwards from the numbers: at N_L = 120, t = 0.1 is 26,042 steps. For the
accumulated amplification to stay O(1) we need `g · N_steps ≲ 1`, i.e.

> **g ≲ 4e-5, against the measured 4.4e-4 — a 10× reduction in per-step gain.**

That is a *specification*, not an aspiration, and it is the most useful single number
this exercise produces. It tells us:

- a lever that improves the delivered force by 2 orders but leaves g alone is worthless
  for stability (which is exactly what `cellCentreInverse` did);
- a lever that reduces g by 10× wins even if it costs an order of accuracy;
- and the two candidate levers with a plausible 10× on g are both *temporal*:
  making the capillary force implicit or iterated (so the lag term vanishes at
  convergence), or adding a per-step damping that scales with g.

Second flip, on reframing C: assume the discrete equilibrium **is** a fixed point.
Then at equilibrium `σκ_f·snGrad(α)·|S_f|` must lie in the range of `snGrad` applied
to some cell field. For constant κ that field is `σκα`. So the admissible set of face
curvatures at equilibrium is

> `{ κ_f : σ κ_f snGrad(α) |S_f| ∈ range(snGrad) }`

and the **distance of our delivered κ_f from that set** is the cleanest possible order
parameter — static, per-step, no long run, and it subsumes the non-gradient content
measurement we already make. Note this is a *diagnostic*, not a scheme: the
potential form itself is closed as a scheme, and re-opening it is not proposed.

---

## 6. Make it bigger

Generalise from "make this configuration stable" to "score any configuration
cheaply".

1. **Make g a solver-computed column.** dt-normalised, dimension-agnostic, valid for
   2D and 3D, stationary and oscillating, any delivery. Every future run then scores
   itself on the order parameter instead of on t_blow.
2. **A gain gate before any long run.** Any new delivery, filter, or coupling must
   demonstrate a reduction in g over ~200 steps before it earns a ladder. This turns
   a 12-hour experiment into a two-minute one, and it is the difference between a job
   and an asset. It would have saved most of this campaign: every lever that bought
   prefactor and not exponent would have been rejected on day one.
3. **State the paper's claim at class level, not solver level.** "For balanced-force
   CSF with a reconstructed curvature, the parasitic current's growth is set by a
   per-step amplification of the non-absorbable force content `α_f·snGrad(κ_c)`; we
   measure that amplification and show it is h-independent under a once-per-step
   explicit coupling, so refinement at fixed capillary CFL accumulates more
   amplification, not less." That is a statement about the method class, and it is
   supported by interFoam and interFlow as the two comparison points.

---

## Execution order, by information per unit cost

**Phase 0 — restore a trustworthy baseline.** Already running: the four post-fix
studies. Nothing filtered can be interpreted until these land, because of the seam
bug.

**Phase 1 — cheap, decisive, no new code.**

1. **The fixed-dt ladder** (running). Three meshes at one dt. If the two coarse
   meshes destabilise purely by taking more, smaller steps, the mechanism is the
   once-per-step lag and reframing B is confirmed.
2. **`psiOuterCorrectors yes` at one resolution.** The switch exists and defaults off,
   so three outer correctors currently converge momentum and pressure against a force
   that cannot change. Predicted to change **g**, not merely t_blow — a sharper test
   than any survival comparison.
3. **The first-ten-step impulse**, free from the re-runs: the fixed Young–Laplace
   solve should reduce it. Quantifies how much of g is the ~1% initial pressure error.

**Phase 2 — small new measurement, high leverage.**

4. **g as a solver column** (item 6.1), plus the gain gate.
5. **Sub-map gain decomposition** (§4) on the existing power-iteration harness,
   restricted to the interface-displacement component of δψ that the gauge result
   already isolated.
6. **Distance from the admissible κ_f set** (§5), as a static per-step residual.

**Phase 3 — structural, only if Phase 1 confirms the lag.**

7. **Semi-implicit or iterated capillary force** (task #24, currently on hold at the
   user's direction). The specification from §5 is now quantitative: it must deliver
   g ≤ 4e-5. Its earlier licence rested on r₀ = 0, which was invalidated; the
   objection that it adds artificial σΔt interfacial dissipation stands and must be
   measured against the transport orders, not assumed away.
8. **The bounded-representation jump** (§2): anchor or bound the interface field the
   way interFoam and interFlow do, keeping the fitted curvature. The interFlow arm now
   built gives the reference numbers to judge whether it is worth the architectural
   cost.

**Deprioritised, with reason stated:** the tangential-structure diagnostic. It was my
previous recommendation; §3B argues it is the wrong instrument for an h-independent,
dt-proportional gain. Revisit if and only if the fixed-dt ladder shows the rates are
h-set after all.
