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

#include "directCorrector.H"
#include "addToRunTimeSelectionTable.H"
#include "slReconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directCorrector, 0);
    addToRunTimeSelectionTable(slCorrector, directCorrector, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directCorrector::directCorrector
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    slCorrector(mesh, dict)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::directCorrector::correct
(
    volScalarField& psi,
    const pointField& feet,
    slReconstruction& recon
) const
{
    recon.update(psi);                       // build the reconstruction from psi^n
    const slReconstruction& R = recon;
    const bool clip = R.clipToStencilBounds();

    footRadiusGuard(R, feet);

    scalarField newPsi(mesh_.nCells());
    label nNonFinite = 0;
    forAll(feet, c)
    {
        newPsi[c] = robustEvaluate(R, c, feet[c], clip, nNonFinite);
    }

    psi.primitiveFieldRef() = newPsi;
    psi.correctBoundaryConditions();

    warnNonFinite(nNonFinite);
}

// ************************************************************************* //
