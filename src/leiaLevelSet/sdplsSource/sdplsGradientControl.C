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

#include "sdplsGradientControl.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvSolution.H"

namespace Foam
{
    defineTypeNameAndDebug(sdplsGradientControl, false);
    addToRunTimeSelectionTable(sdplsSource, sdplsGradientControl, Dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sdplsGradientControl::sdplsGradientControl
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    sdplsSource(dict, mesh),
    law_(gradientControlLaw::New(dict.subDict("law")))
{
    if (!law_->active())
    {
        WarningInFunction
            << "sdplsSource gradientControl with law none and strain weight"
            << " none applies F = 0 everywhere; select noSource instead."
            << nl << endl;
    }

    const dictionary& fvs = static_cast<const fvSolution&>(mesh);
    if (law_->needsNormalStrain() && fvs.isDict("levelSet"))
    {
        const word extension =
            fvs.subDict("levelSet").subOrEmptyDict("velocityExtension")
               .getOrDefault<word>("type", "none");
        if (extension != "none")
        {
            WarningInFunction
                << "levelSet.velocityExtension selects '" << extension
                << "' and the strain weight is '" << law_->weight().type()
                << "': the extension already removes the normal strain from the"
                << " transport of |grad psi|, so this source cancels it twice."
                << nl << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::sdplsGradientControl::nonLinearPart
(
    const volScalarField& R,
    const volScalarField& psi,
    const volVectorField& U
) const
{
    const fvMesh& mesh = psi.mesh();
    const volVectorField gradPsi(grad(psi));

    // sigma = |symm(grad U)|_F, only for a law that reads it.
    tmp<volScalarField> tsigma;
    if (law_->needsStrainRate())
    {
        tsigma = mag(symm(fvc::grad(U, "gradUSdpls")));
    }

    tmp<volScalarField> tF
    (
        new volScalarField
        (
            IOobject
            (
                word(),
                fileName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 0.0)
        )
    );
    volScalarField& F = tF.ref();

    const gradientControlLaw& law = *law_;
    const bool withSigma = tsigma.valid();

    scalarField& Fi = F.primitiveFieldRef();
    forAll(Fi, celli)
    {
        Fi[celli] = law.F
        (
            magSqr(gradPsi[celli]),
            withSigma ? tsigma()[celli] : scalar(0),
            R[celli]
        );
    }

    // Boundary values for the written diagnostic field: the same law at the
    // patch values (the discretization reads the cell values only).
    volScalarField::Boundary& Fb = F.boundaryFieldRef();
    forAll(Fb, patchi)
    {
        const vectorField& gp = gradPsi.boundaryField()[patchi];
        const scalarField& Rp = R.boundaryField()[patchi];
        scalarField& Fp = Fb[patchi];
        forAll(Fp, facei)
        {
            Fp[facei] = law.F
            (
                magSqr(gp[facei]),
                withSigma ? tsigma().boundaryField()[patchi][facei] : scalar(0),
                Rp[facei]
            );
        }
    }

    return tF;
}


// ************************************************************************* //
