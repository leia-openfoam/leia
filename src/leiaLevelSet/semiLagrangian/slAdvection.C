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
    analyticVelocity_(slDict_.getOrDefault<Switch>("analyticVelocity", true)),
    scheme_(slScheme::New(mesh))
{
    Info<< "slAdvection: scheme = " << scheme_->type()
        << ", CFLmax = " << CFLmax_
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
    // Delegate to the runtime-selected scheme (pointValue | fluxForm); both
    // reuse this object's reconstruction and correction strategy.
    scheme_->advance(psi, Unew, Uold, recon_(), corrector_());
}


void Foam::slAdvection::meanCurvature
(
    const volScalarField& psi,
    volScalarField& kappa
)
{
    // Re-fit the CURRENT psi so the stored per-cell quadratic reflects the new
    // interface, then let the (quadratic) reconstruction fill kappa symbolically.
    recon_->update(psi);
    recon_->meanCurvature(kappa);
}


void Foam::slAdvection::meanCurvatureLaplacian
(
    const volScalarField& psi,
    volScalarField& kappa
)
{
    recon_->update(psi);
    recon_->meanCurvatureLaplacian(kappa);
}


void Foam::slAdvection::meanCurvatureNoExtension
(
    const volScalarField& psi,
    volScalarField& kappa
)
{
    recon_->update(psi);
    recon_->meanCurvatureNoExtension(kappa);
}


void Foam::slAdvection::meanCurvatureClosestPoint
(
    const volScalarField& psi,
    volScalarField& kappa
)
{
    recon_->update(psi);
    recon_->meanCurvatureClosestPoint(kappa);
}

// ************************************************************************* //
