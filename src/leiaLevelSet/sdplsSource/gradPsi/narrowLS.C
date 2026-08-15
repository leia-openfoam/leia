/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt
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

#include "narrowLS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- Look up the per-cell LLSQ plane normal, failing with the REQUIREMENT rather
//  than with "failed lookup of nc". Must be a free function called from the
//  member initializer list: a check in the constructor body runs too late,
//  because the reference member is bound first.
const Foam::volVectorField& checkedNc(const Foam::fvMesh& mesh)
{
    if (!mesh.foundObject<Foam::volVectorField>("nc"))
    {
        FatalErrorInFunction
            << "gradPsi type `narrowLS` needs the per-cell LLSQ plane normal "
            << "field `nc`, which is registered ONLY by phaseIndicator type "
            << "`geometric`. The selected phase indicator does not provide it."
            << nl
            << "Either set phaseIndicator to `geometric`, or use gradPsi type "
            << "`fvc`." << nl
            << "NOTE: switching the phase indicator also changes alpha, and so "
            << "the volume and shape errors -- a comparison that changes both "
            << "is not an ablation of the normal alone."
            << Foam::exit(Foam::FatalError);
    }
    return mesh.lookupObject<Foam::volVectorField>("nc");
}

} // End anonymous namespace

namespace Foam
{
    defineTypeNameAndDebug(narrowLS, false);
    addToRunTimeSelectionTable(gradPsi, narrowLS, Dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::narrowLS::narrowLS(const fvMesh& mesh)
    :
        gradPsi(mesh),
        narrowBand_(mesh.lookupObject<volScalarField>("NarrowBand")),
        nc_(checkedNc(mesh))
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

Foam::tmp<volVectorField> Foam::narrowLS::grad(const volScalarField& psi) const
{
    // Dedicated SDPLS gradient keyword (see gradPsi::grad).
    tmp<volVectorField> tgradpsi = fvc::grad(psi, "gradPsiSdpls");
    forAll(narrowBand_, cellID)
    {
        if (narrowBand_[cellID] == 1)
        {
            tgradpsi.ref()[cellID] = nc_[cellID];
        }
    }
    return tgradpsi;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
