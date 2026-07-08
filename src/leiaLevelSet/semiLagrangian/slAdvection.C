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

\*---------------------------------------------------------------------------*/

#include "slAdvection.H"
#include "fvcGrad.H"
#include <cmath>   // std::isfinite

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slAdvection, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slAdvection::slAdvection(const fvMesh& mesh)
:
    mesh_(mesh),
    fvSolution_(mesh),
    levelSetDict_(fvSolution_.subDict("levelSet")),
    slDict_(levelSetDict_.subOrEmptyDict("semiLagrangian")),
    recon_(slReconstruction::New(mesh)),
    CFLmax_(slDict_.getOrDefault<scalar>("CFLmax", 1.0)),
    analyticVelocity_(slDict_.getOrDefault<Switch>("analyticVelocity", true))
{
    Info<< "slAdvection: CFLmax = " << CFLmax_
        << ", clipToStencilBounds = " << recon_->clipToStencilBounds()
        << endl;
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::slAdvection> Foam::slAdvection::New(const fvMesh& mesh)
{
    return autoPtr<slAdvection>(new slAdvection(mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::slAdvection::advect
(
    volScalarField& psi,
    const volVectorField& Unew,
    const volVectorField& Uold
)
{
    const scalar dt = mesh_.time().deltaTValue();
    const volVectorField& C = mesh_.C();

    // grad u^{n+1} at cell centres via least squares (fvSchemes key "gradU",
    // pointCellsLeastSquares); gradU[c] = d(U_j)/d(x_i) so (u.grad)u = (u & gradU).
    const volTensorField gradU(fvc::grad(Unew, "gradU"));

    // Prepare the reconstruction from psi^n (psi still holds the old field).
    recon_->update(psi);
    const slReconstruction& R = recon_();
    const bool clip = R.clipToStencilBounds();

    scalarField newPsi(mesh_.nCells());
    label nOutside = 0;
    label nNonFinite = 0;
    scalar maxRatio = 0;

    forAll(C, c)
    {
        const vector& uNew = Unew[c];
        const vector& uOld = Uold[c];

        // Taylor backward displacement (KEEP the 1/2 on the dt^2 term):
        //   x_d = x_c - u^{n+1} dt + 1/2 [ du/dt + (u^{n+1}.grad)u^{n+1} ] dt^2
        // du/dt ~ (u^{n+1} - u^n)/dt is first-order (it multiplies dt^2).
        const vector accel = (uNew - uOld)/dt + (uNew & gradU[c]);
        const point xd = C[c] - uNew*dt + 0.5*accel*dt*dt;

        // Foot-radius guard = the operational CFL<=1 check: the foot must stay
        // inside the point-neighbour hull, else evaluate() would extrapolate.
        const scalar disp = Foam::mag(xd - C[c]);
        const scalar radius = R.stencilRadius(c);
        if (radius > SMALL)
        {
            maxRatio = Foam::max(maxRatio, disp/radius);
            if (disp > radius)
            {
                ++nOutside;
            }
        }

        // evaluate() is the reconstruction (Barth-Jespersen slope-limited only
        // if limitSlope is on -- off by default because it costs convergence
        // order). Robustness: non-finite -> stencil mid; optional tight clip;
        // and an ALWAYS-ON generous value cap (mid +- 10*range) that only fires
        // on runaway growth, preventing overflow / SIGFPE without touching a
        // well-behaved (order-preserving) reconstruction.
        scalar v = R.evaluate(c, xd);
        scalar lo, hi;
        R.stencilRange(c, lo, hi);
        const scalar mid = 0.5*(lo + hi);
        if (!std::isfinite(v))
        {
            v = mid;
            ++nNonFinite;
        }
        else if (clip)
        {
            v = Foam::min(Foam::max(v, lo), hi);
        }
        else
        {
            const scalar cap = 10.0*Foam::max(hi - lo, SMALL);
            v = Foam::min(Foam::max(v, mid - cap), mid + cap);
        }
        newPsi[c] = v;
    }

    psi.primitiveFieldRef() = newPsi;
    psi.correctBoundaryConditions();

    // Optional post-advection fix-up (band model: re-extend psi outside the band
    // as a clean signed distance so freshly-entered band cells get a good value).
    recon_->postAdvect(psi);

    reduce(nOutside, sumOp<label>());
    reduce(nNonFinite, sumOp<label>());
    reduce(maxRatio, maxOp<scalar>());
    if (nNonFinite > 0)
    {
        WarningInFunction
            << "reconstruction produced a non-finite value in " << nNonFinite
            << " cells (reconstruction unstable at this CFL / resolution);"
            << " reset to the stencil mid-range." << endl;
    }
    if (nOutside > 0)
    {
        WarningInFunction
            << "semi-Lagrangian foot left the point-neighbour stencil in "
            << nOutside << " cells (max |x_d - x_c|/stencilRadius = "
            << maxRatio << "); CFL likely > 1 -- reduce maxCo." << endl;
    }
}

// ************************************************************************* //
