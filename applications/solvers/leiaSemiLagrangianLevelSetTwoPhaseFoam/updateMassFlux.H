// THE MASS FLUX rhoPhi and, for rhoLENT, the auxiliary density, on EVERY PIMPLE
// outer corrector. Shared by the semi-Lagrangian (slAlphaEqn.H) and the Eulerian
// (alphaEqn.H) two-phase solvers. Moved here VERBATIM from slAlphaEqn.H on
// 2026-09-27.
    // Mass flux and auxiliary density are updated on EVERY PIMPLE outer
    // iteration. This is load-bearing: pEqn changes phi^o, so rhoPhi and the
    // rhoLENT continuity update must be reassembled before the next UEqn.
    if (!freezeRhoPhi)
    {
        if (massFluxType == "interpolatedDensity")
        {
            rhoPhi = fvc::interpolate(rho)*phi;
        }
        else
        {
            rhoPhi = rhof*phi;
        }

        if (projectMassFlux)
        {
            // Correct the flux BEFORE anything reads it -- the auxiliary density
            // equation below and div(rhoPhi,U) must both see the corrected flux, or
            // the projection fixes a balance nobody uses.
            #include "massFluxProjection.H"
        }

        if (useRhoLENT)
        {
            #include "rhoLENTEqn.H"
        }
        else if (massResidualDiagnostic)
        {
            // Same residual, for a model that solves no density equation. Pure
            // diagnostic; inert unless massResidualDiagnostic is switched on.
            #include "massResidualDiag.H"
        }
    }
