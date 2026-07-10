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
    corrector_(slCorrector::New(mesh, slDict_)),
    CFLmax_(slDict_.getOrDefault<scalar>("CFLmax", 1.0)),
    analyticVelocity_(slDict_.getOrDefault<Switch>("analyticVelocity", true))
{
    Info<< "slAdvection: CFLmax = " << CFLmax_
        << ", clipToStencilBounds = " << recon_->clipToStencilBounds()
        << ", correction = " << corrector_->type() << endl;
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

    // ------------------------------------------------------------------ //
    // Departure (foot) points: computed ONCE. The Taylor backward-foot
    // integrator is UNCHANGED -- it is only hoisted out of the reconstruction
    // loop below and cached, so the deferred-correction passes reuse the same
    // characteristic feet. (KEEP the 1/2 on the dt^2 term.)
    //   x_d = x_c - u^{n+1} dt + 1/2 [ du/dt + (u^{n+1}.grad)u^{n+1} ] dt^2
    // ------------------------------------------------------------------ //
    pointField feet(mesh_.nCells());
    forAll(C, c)
    {
        const vector& uNew = Unew[c];
        const vector& uOld = Uold[c];
        const vector accel = (uNew - uOld)/dt + (uNew & gradU[c]);
        feet[c] = C[c] - uNew*dt + 0.5*accel*dt*dt;
    }

    // ------------------------------------------------------------------ //
    // The selected correction strategy assembles psi^{n+1} in place from psi^n,
    // evaluating the reconstruction at the fixed feet and (for deferredCorrection)
    // rebuilding it from the current iterate. It owns the value cap / non-finite
    // reset / foot-radius guard. The backtracking above is not its concern.
    // ------------------------------------------------------------------ //
    corrector_->correct(psi, feet, recon_());

    // Optional post-advection fix-up (band model: re-extend psi outside the band
    // as a clean signed distance so freshly-entered band cells get a good value).
    recon_->postAdvect(psi);
}

// ************************************************************************* //
