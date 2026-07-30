/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::slAdvection

Description
    Cell-centred semi-Lagrangian advection of the level set:
    psi is constant along characteristics, so psi^{n+1}(x_c) = psi^n(x_d), with

      * the backward departure (foot) point from the explicit Taylor
        displacement integrator (Maric thesis Sec. 6.3.1, eqs. 165-167;
        local error O(dt^3) -- keep the 1/2 on the dt^2 term):

          x_d = x_c - u_c^{n+1} dt
                    + 1/2 [ (u_c^{n+1} - u_c^n)/dt
                            + (u_c^{n+1} . grad_c) u^{n+1} ] dt^2

      * the OLD field evaluated at x_d by a runtime-selectable slReconstruction
        (linearTaylor | linearWeightedLeastSquares | quadraticTaylor |
        quadraticWeightedLeastSquares), always the ARRIVAL cell's own polynomial
        => interpolation, not extrapolation, under CFL <= 1.

    No flux, no divergence, no linear solve. The characteristic feet are computed
    ONCE (the Taylor backtracking) and the OLD field is assembled at those feet by
    a runtime-selectable slCorrector (`correction`: direct | deferredCorrection),
    which composes with any runtime-selectable slReconstruction. Reads
    `levelSet { semiLagrangian { reconstruction; correction; CFLmax; } }`
    in fvSolution. grad u^{n+1} uses the LSQ scheme keyed
    `gradU` in fvSchemes (never Green-Gauss). The u^{n+1}, u^n fields are passed
    in, so the integrator works unchanged for a Navier-Stokes velocity; for the
    prescribed-velocity verification those fields are the analytic velocity
    sampled at cell centres (velocityModel::setVelocity), so no separate analytic
    path is required.

SourceFiles
    slAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef slAdvection_H
#define slAdvection_H

#include "fvMesh.H"
#include "fvSolution.H"
#include "volFields.H"
#include "Switch.H"
#include "autoPtr.H"
#include "slReconstruction.H"
#include "slCorrector.H"
#include "slScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class slAdvection Declaration
\*---------------------------------------------------------------------------*/

class slAdvection
{
    const fvMesh& mesh_;
    const fvSolution& fvSolution_;
    const dictionary& levelSetDict_;
    dictionary slDict_;

    //- The reconstruction of psi^n at the foot.
    autoPtr<slReconstruction> recon_;

    //- The strategy that assembles psi^{n+1} from the reconstruction at the feet
    //  (direct single pass | deferred correction). Selected by `correction`.
    autoPtr<slCorrector> corrector_;

    //- Convective CFL bound used for the foot-radius warning (the solver's
    //  adjustTimeStep/maxCo keeps the actual step under this).
    scalar CFLmax_;

    //- Reserved control: for prescribed-velocity verification the passed U is
    //  already the analytic cell-centre velocity, so the discrete-field path is
    //  used unconditionally; a future exact-gradient trajectory can honour it.
    Switch analyticVelocity_;

    //- Advection scheme (RTS): pointValue (default; the Lagrange interpolation
    //  psi^{n+1}(x_c) = psi^n(x_d), non-conservative) or fluxForm (conservative
    //  swept-face finite-volume semi-Lagrangian, CFL <= 1). Selected by
    //  levelSet.semiLagrangian.scheme.
    autoPtr<slScheme> scheme_;

public:

    //- Runtime type information
    ClassName("slAdvection");

    // Constructors

        explicit slAdvection(const fvMesh& mesh);

    // Selectors

        static autoPtr<slAdvection> New(const fvMesh& mesh);

    //- Destructor
    ~slAdvection() = default;

    // Member Functions

        const slReconstruction& reconstruction() const { return recon_(); }

        //- Advect one step: psi holds psi^n on entry, psi^{n+1} on exit.
        //  Unew = u^{n+1}, Uold = u^n at cell centres.
        //
        //  The foot is expanded about the ARRIVAL time, so Unew really must be the
        //  NEW-time velocity: passing u^n there instead leaves an O(dt^2 d_t u)
        //  foot error and drops the transport to first order in time (invisible in
        //  any steady or uniform test, where d_t u = 0 makes the two identical).
        //  A caller that cannot supply u^{n+1} -- e.g. a coupled solver advecting
        //  before its momentum solve -- must pass the Adams-Bashforth estimate
        //  u^n + (dt/dt0)(u^n - u^{n-1}) as Unew and u^n as Uold, which reduces the
        //  kernel to the equivalent departure-time-centred second-order trace
        //  x_d = x_c - dt/2 (3u^n - u^{n-1}) + 1/2 dt^2 (u.grad)u. See
        //  leiaSemiLagrangianLevelSetTwoPhaseFoam/slAlphaEqn.H.
        void advect
        (
            volScalarField& psi,
            const volVectorField& Unew,
            const volVectorField& Uold
        );

        //- Re-fit the CURRENT psi and fill kappa (1/length) with the interface mean
        //  curvature evaluated symbolically from the reconstruction (see
        //  slReconstruction::meanCurvature). Requires a quadratic reconstruction.
        void meanCurvature(const volScalarField& psi, volScalarField& kappa);

        //- As meanCurvature, but the SDF-assuming curvature kappa = Laplacian(psi)
        //  = tr(H) (see slReconstruction::meanCurvatureLaplacian). For comparison.
        void meanCurvatureLaplacian(const volScalarField& psi, volScalarField& kappa);

        //- As meanCurvature, but evaluated at the cell centre WITHOUT the foot-point
        //  normal extension (see slReconstruction::meanCurvatureNoExtension).
        void meanCurvatureNoExtension(const volScalarField& psi, volScalarField& kappa);

        //- As meanCurvature, but the foot found by a trust-region CLOSEST-POINT
        //  Newton on the cell's own quadratic (KKT system with the fit Hessian)
        //  instead of the gradient-projection flow; per-cell fallback to the
        //  cell-centre value (see slReconstruction::meanCurvatureClosestPoint).
        void meanCurvatureClosestPoint(const volScalarField& psi, volScalarField& kappa);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
