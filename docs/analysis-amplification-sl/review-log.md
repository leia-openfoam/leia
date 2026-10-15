# Review comments on the amplification report

Appended by `make review-amplification`. The annotated PDF is not
tracked, so this file is the durable record of each review round.

## 2026-10-15 17:16 -- 12 comment(s) from ms_reviewed.pdf

**[1] p.5 -- A compressing flow does not make phi grow. It steepens it.** (TomislavMaric)

> The point is what the exact solution does to the field as a whole. Along any trajectory ϕ is constant, by (4). So exact transport does not alter the collection of values that exist in the domain — it only moves them around, relabelling which cell holds which. For a closed or periodic domain this makes the exact update a rearrangement: every maximum, every minimum andil every norm of the field is p ...

emphasize

**[2] p.5 -- A compressing flow does not make phi grow. It steepens it.** (TomislavMaric)

> Equation (1) is not that equation. is a label carried by the fluid, not an amount of anything,l ϕ so a particle keeps its value whatever the flow around it does. What compression changes is howlf far apart the level sets lie: squeezing the fluid brings surfaces of different closer together and ϕ therefore raises ϕ|. The value does not grow; the gradient does. (That steepening is a real |∇ issue fo ...

The gradient norm transport equation for \phi states that the gradient norm rate of change is proportional to n^T (grad U) n 

where n = grad phi)/mag (grad phi)

or something in this direction, so even if the field is divergence free, only a velocity field that is constant in the normal direction of every hypersurface given by phi will not destroy mag(grad phi).

**[3] p.6 -- What was measured** (TomislavMaric)

> (14)

Add derivation steps and didactic explanation of how we get (14) from (13), express rho precisely.

**[4] p.6 -- What was measured** (TomislavMaric)

> over 1563 steps that factor is roughly 1000, so a truncation error far too small ρ = 1.00441 to see becomes a visible defect.

We do not know exact solution, or did we manufacture one here to get these amplification factors for a constant velocity and translating circle or something?

**[5] p.7 -- Step 1: the update is a weighted sum of neighbour values** (TomislavMaric)

> Start from a random field. Apply the rule, measure how much the size of the field changed, rescale it back to size one, and repeat two thousand times. The factor settles down to a single number: the largest factor by which the rule can stretch anything. That number is what Table 1 reports.

This needs to be written out also mathematically, in equations.

**[6] p.8 -- Step 3: from weights to a growth factor** (TomislavMaric)

> The weights on the right-hand side sum to one

How is this exactly ensured, describe.

**[7] p.9 -- Why each wave is treated separately** (TomislavMaric)

> Apply (16) to a wave and use ϕm+ = sj eiθ · (m+ sj ) eiθ · m eiθ · sj =

Write these steps explicitly, like in a book, making it easy to follow and to understand without looking up and doing substitutions in my mind.

**[8] p.9 -- Why each wave is treated separately** (TomislavMaric)

> that wave neither grows nor decays. Neutral.

So we want locally this, because this means locally we adhere to the material derivative equation? Connect these bullet points to the material derivative.

**[9] p.10 -- Why one dimension cannot answer the question** (TomislavMaric)

> That is exactly the Lax–Wendroff scheme, a classical method from 1960i whose properties are in every textbook. Our fit, in one dimension, reproduces it.

Very interesting.

**[10] p.10 -- Why one dimension cannot answer the question** (TomislavMaric)

> ≤ 1

This is stable, yes ,but smaller than one can introduce reduction of previous values?

**[11] p.10 -- Why one dimension cannot answer the question** (TomislavMaric)

> The script performs all five of these identities symbolically and refuses to report anything else if any of them fai

OK so we need to solve this equation for C, to get G(\theta) = 1 for all theta? This tells us C that we should then use locally to advect. How do we then do local time stepping so that the field is synchronously advected?

**[12] p.10 -- Why one dimension cannot answer the question** (TomislavMaric)

> (25)

visualize using tikz-pgf and two lines of 1D point values, one at t^n, another at t^{n+1} to see the connection.

## 2026-10-15 18:07 -- 10 comment(s) from ms_reviewed.pdf

**[1] p.11 -- Why each wave is treated separately** (TomislavMaric)

> ei θ · m ϕm = = cos(θ · m) + i sin(θ · m

Use a tikz-pgf schematic image in 1D and show the smallest wave, twice largest, then one more twice largest, and their superposition.

**[2] p.12 -- What the exact solution does to a wave, and what the three cases mean** (TomislavMaric)

> Co ˆu)

express this as well here, make descriptions didactic and self contained

**[3] p.13 -- The one-dimensional step, drawn** (TomislavMaric)

> is the target, not the ceilingi

Can we optimize weights for |G|=1 on a given mesh? Say hexahedral? Assume that there is some angle deviation and distance deviation w.r.t. perfect equidistant hexahedral mesh, which happens when unstructured hexahedral mesh is generated?

**[4] p.14 -- Can we choose C to make the scheme exact?** (TomislavMaric)

> Any velocity variation breaks it, and so does any mesh grading or refinement

and so does the capillary time step restriction

**[5] p.15 -- Why one dimension cannot answer the question** (TomislavMaric)

> There is a second quantity in circulation in this project, and (44) is the cleanest place to see whyi it misleads.

Should we then denote it here and discard it?

**[6] p.18 -- At the Courant numbers we advect at, the inflow cell stops transporting** (TomislavMaric)

> Extrapolating a quadratic is the textbook way to amplify: a parabola fitted on one side and continued past the data grows quickly, and the further you continue it the worse it gets. Equation (50) says exactly that, with the extrapolation distance measured in cell widths.

Emphasize.

**[7] p.18 -- At the Courant numbers we advect at, the inflow cell stops transporting** (TomislavMaric)

> What the boundary face value is really doing

Should we not apply the boundary condition there? But then again zero gradient is the worst possible boundary condition, because the gradient should not be zero at the inlet of \phi.. can we not extrapolate safely?

**[8] p.19 -- Is a smaller time step safer? No — and the curve has a closed form** (TomislavMaric)

> The freezing is strongest when the departure point coincides with the face, at — which is the Courant number of our 2D hexahedral advection case, Co = 0.5 including the convergence study that is the hex rung of the standing regression set.

OK so we use CFL=0.5 for now, it is still larger than the capillary time restriciton?

**[9] p.20 -- The face value's weight does not matter. Its presence does.** (TomislavMaric)

> Whether that materially affects those convergence orders is not established here.l The frozen cells are one layer at the inflow patch, and the interface in those studies spends mostl of its time away from that layer. But a measured order from a case whose inflow boundary is inert deserves to be checked rather than assumed, and the check is cheap: compare against and against a run at 0.45. stencilB ...

Can we not assume that we can extend curvature as constant in the inlet-normal direction, woudl that help ?

**[10] p.20 -- The face value's weight does not matter. Its presence does.** (TomislavMaric)

> which under zeroGradient

I think this might be a problem, because it is not zero-gradient, which means \phi is flat in the inlet-normal direction... this is inconsistent. What is the best practice for level set methods on the inlet?

