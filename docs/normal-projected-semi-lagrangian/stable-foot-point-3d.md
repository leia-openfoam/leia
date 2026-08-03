# Stable Foot-Point Algorithm for Implicit Surfaces in 3D

Self-contained description. Everything needed to implement the method is on this
page. (Companion to `normal-projected-semi-lagrangian.md`, Sec. 7.2, where this
algorithm is the offset engine of nSL Variant (a).)

## 1. Problem statement

Given:

* a scalar field $f:\mathbb R^3\to\mathbb R$, at least $C^2$ near its level sets,
  with $\nabla f\ne\mathbf 0$ in the region of interest,
* a target level $c\in\mathbb R$ defining the surface
  $\Sigma_c=\{\mathbf x\in\mathbb R^3: f(\mathbf x)=c\}$,
* a query point $\mathbf p\in\mathbb R^3$ in the vicinity of $\Sigma_c$,

find the foot point $\mathbf x^\Sigma$: the point of $\Sigma_c$ closest to
$\mathbf p$. Equivalently, $\mathbf x^\Sigma\in\Sigma_c$ such that
$\mathbf p-\mathbf x^\Sigma$ is parallel to $\nabla f(\mathbf x^\Sigma)$ (the
connecting segment is normal to the surface at the foot point, not at the query
point — this distinction is the entire difficulty).

Throughout, define the shifted residual $g(\mathbf x):=f(\mathbf x)-c$, so the
target surface is $g=0$. All formulas below use $g$; note $\nabla g=\nabla f$.

## 2. Building block: `surfacepoint` — pull a nearby point onto the surface

A point $\mathbf q$ close to $\Sigma_c$ but not exactly on it is brought onto the
surface by Newton iteration along the local gradient direction:

$$\mathbf q_{k+1}=\mathbf q_k-\frac{g(\mathbf q_k)}{|\nabla g(\mathbf q_k)|^2}\,\nabla g(\mathbf q_k)$$

repeated until $|\mathbf q_{k+1}-\mathbf q_k|$ is below tolerance.

This is a 1D Newton root-find for the scalar function
$t\mapsto g(\mathbf q+t\,\nabla g(\mathbf q))$, re-linearized at every step. It
converges fast when started near the surface, but it only zeroes the residual
$g$ — the resulting point is on the surface but is generally not the closest
point to any given $\mathbf p$.

## 3. The naive foot-point loop — and why it is unstable

The obvious approach alternates two projections. With $\mathbf p_i\in\Sigma_c$
the current surface iterate:

1. Foot point on the tangent plane. Project the query $\mathbf p$ onto the
   tangent plane of $\Sigma_c$ at $\mathbf p_i$:
   $$\mathbf q_i=\mathbf p-\frac{(\mathbf p-\mathbf p_i)\cdot\nabla g(\mathbf p_i)}{|\nabla g(\mathbf p_i)|^2}\,\nabla g(\mathbf p_i)$$
2. Return to the surface. $\mathbf p_{i+1}=\texttt{surfacepoint}(\mathbf q_i)$.

Repeat until $|\mathbf p_{i+1}-\mathbf p_i|<\varepsilon$.

This is a fixed-point iteration, not a Newton method for the constrained
closest-point problem. Each cycle uses only the tangent plane — a linear model
of the surface — and therefore ignores how the surface normal rotates as one
moves tangentially. The error contraction factor per cycle scales like
$|d|\,\kappa$, where $d$ is the signed distance of $\mathbf p$ from the surface
and $\kappa$ the local curvature. Consequences:

* convergence degrades as $|d|\kappa$ grows, and the iteration can stall or
  diverge once the query point is a significant fraction of the local radius of
  curvature away from the surface;
* on noisy or corrugated surfaces (locally large $\kappa$), the iterates can
  oscillate tangentially without settling, or settle to a point that is not the
  true foot point.

## 4. The stabilized algorithm

The fix uses no second derivatives. The key observation: after one naive cycle,
three points are available —

* $\mathbf p_i$: current surface point,
* $\mathbf q_i$: its tangent-plane foot point of $\mathbf p$,
* $\mathbf p_{i+1}$: the return of $\mathbf q_i$ to the surface.

Define the two step vectors

$$\mathbf f_1=\mathbf q_i-\mathbf p_i,\qquad \mathbf f_2=\mathbf p_{i+1}-\mathbf q_i .$$

$\mathbf f_1$ is the tangential (linear-model) step. $\mathbf f_2$ measures how
far the true surface pulled the linear prediction back — it is an empirical
curvature sample, obtained for free from the two projections already performed,
with first derivatives only.

These three points define a parabola — a local quadratic model of the surface
path from $\mathbf p_i$ toward the foot point:

$$\mathbf x(\alpha)=\mathbf p_i+\alpha\,\mathbf f_1+\alpha^2\,\mathbf f_2,
\qquad \mathbf x(0)=\mathbf p_i,\ \ \mathbf x(1)=\mathbf p_{i+1}.$$

Now find the $\alpha$ that brings $\mathbf x(\alpha)$ closest to the query point
$\mathbf p$, i.e. minimize $D(\alpha)=|\mathbf p-\mathbf x(\alpha)|^2$. Its
derivative is a cubic in $\alpha$; write it as
$\tfrac12 D'(\alpha) = -(a_0+a_1\alpha+a_2\alpha^2+a_3\alpha^3)$ with

$$a_0=(\mathbf p-\mathbf p_i)\cdot\mathbf f_1,\qquad
a_1=2\,\mathbf f_2\cdot(\mathbf p-\mathbf p_i)-|\mathbf f_1|^2,$$
$$a_2=-3\,\mathbf f_1\cdot\mathbf f_2,\qquad a_3=-2\,|\mathbf f_2|^2 .$$

One Newton step for the root of this cubic, started at $\alpha=1$ (the naive
result), gives

$$\alpha=1-\frac{a_0+a_1+a_2+a_3}{a_1+2a_2+3a_3}.$$

Accept the extrapolation only inside a guard interval $0<\alpha<\alpha_{\max}$
(rejecting it falls back to the plain cycle), then move along the parabola and
return to the surface:

$$\mathbf q_i\leftarrow\mathbf p_i+\alpha\,\mathbf f_1+\alpha^2\,\mathbf f_2,
\qquad \mathbf p_{i+1}\leftarrow\texttt{surfacepoint}(\mathbf q_i).$$

Loop until $|\mathbf p_{i+1}-\mathbf p_i|<\varepsilon$.

Suggested tolerances: $\varepsilon=10^{-6}$ (in units of the geometry),
$\alpha_{\max}=20$.

### 4.1 Full algorithm

```
input:  field g (= f - c), gradient grad_g, query point p
output: foot point on {g = 0}

p_i = surfacepoint(p)                     # initial on-surface point
repeat:
    # (1) foot point of p on the tangent plane at p_i
    lam = (p - p_i)·grad_g(p_i) / |grad_g(p_i)|^2
    q_i = p - lam * grad_g(p_i)

    # (2) back to the surface
    p_next = surfacepoint(q_i)

    f1 = q_i - p_i
    f2 = p_next - q_i

    # (3) quadratic (parabola) correction
    if |f1| > eps:
        a0 = (p - p_i)·f1
        a1 = 2*(f2·(p - p_i)) - |f1|^2
        a2 = -3*(f1·f2)
        a3 = -2*|f2|^2
        alpha = 1 - (a0 + a1 + a2 + a3)/(a1 + 2*a2 + 3*a3)

        if 0 < alpha < alpha_max:
            q_i    = p_i + alpha*f1 + alpha^2*f2
            p_next = surfacepoint(q_i)

    converged = |p_next - p_i| < eps
    p_i = p_next
until converged

return p_i
```

### 4.2 Why this works

The naive loop fails because a linear surface model cannot represent the
rotation of the normal along the surface. The parabola
$\mathbf p_i+\alpha\mathbf f_1+\alpha^2\mathbf f_2$ does carry that
information — its second-order coefficient $\mathbf f_2$ is the measured
deviation between the linear prediction and the actual surface over one step.
Minimizing the distance to $\mathbf p$ over this quadratic model, and only then
re-projecting, is effectively one Newton step on a curvature-aware model, at
the price of first derivatives only. The structure is an acceleration of the
fixed-point map: the sequence of iterates itself supplies the curvature
information that a full Newton method would obtain from the Hessian.

### 4.3 Practical notes for 3D

* Cost per cycle: two gradient-based `surfacepoint` solves (each a few scalar
  Newton steps), one $3$-vector tangent-plane projection, and a handful of dot
  products for the $\alpha$ correction. No linear systems, no Hessians.
* Degenerate guard: reject the correction if the denominator
  $a_1+2a_2+3a_3$ is near zero or if $\alpha$ leaves $(0,\alpha_{\max})$; the
  plain cycle is then used for that iteration.
* Starting point: any point in the vicinity of the surface works; starting from
  `surfacepoint(p)` (the gradient-projection of the query itself) is the
  natural choice.
* Sign/distance output: on convergence, the signed distance of $\mathbf p$ to
  the surface is $d=(\mathbf p-\mathbf x^\Sigma)\cdot\mathbf n$ with
  $\mathbf n=\nabla g(\mathbf x^\Sigma)/|\nabla g(\mathbf x^\Sigma)|$; the sign
  convention follows the orientation of $\nabla g$.
* Validity region: the closest point is unique only while $\mathbf p$ is closer
  to the surface than the local radius of curvature ($|d|\,\kappa_{\max}<1$);
  beyond that, any foot-point algorithm can jump between competing branches,
  and results should be treated accordingly.
