# VOID — produced on the closed-box `translatingDroplet2D` (before commit 440107f)

Every table in this directory was curated from runs of `cases/translatingDroplet2D`
made before 2026-09-02, when the case's `blockMeshDict` put **all four sides into one
`walls` patch**. The mesh therefore had no `inlet` and no `outlet`, while every field
in `0.org` declared both. OpenFOAM errors when a *mesh* patch is missing from a field
but silently ignores a *field* entry matching no mesh patch, so
`inlet { type fixedValue; value uniform (0.05 0 0); }` was parsed and discarded on
every run these tables are built from.

The case was a **closed slip-wall box**. Slip imposes `U.n = 0`, and on the x-normal
faces the normal *is* x, so the uniform stream the case initialises with is
incompatible with its own boundaries. The pressure projection annihilated it on step 1:

    first-step local continuity error   1.0e-05   (1.6e-12 by step 3)
    ambient mean Ux at t = 0.02        -2.0e-03   (should be +5.0e-02)
    droplet mean Ux                    +6.1e-02
    whole-domain mean Ux               -7.3e-07   (zero net flux, as a closed box must)

`maxMagUPrime = max|U - (U0,0,0)|` consequently measured the **annihilated free
stream** at roughly 2*U0 from the first step onward, not a spurious current.

Measured on the repaired mesh by `config/translatingFreeStreamGate2D`:

    max|U-U0| at step 1   1.031e-01 -> 2.159e-03   (48x)
    mean|U-U0|            5.011e-02 -> 1.160e-05   (4300x)
    maxMagU               0.05987   -> 0.05207 = U0
    step-1 continuity     1.00e-05  -> 2.32e-20

**Do not curate, cite or re-use anything here**, including metrics that look
unaffected by the boundary condition. Per CLAUDE.md ("A wrong setup voids its data"),
a setup that was wrong in one way is not assumed wrong in only that way. The
replacement studies are `config/translatingRepaired2D` and
`config/translatingRepairedEqualRho2D`.
