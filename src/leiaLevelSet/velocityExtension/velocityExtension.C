/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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

#include "velocityExtension.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityExtension, 0);
    defineRunTimeSelectionTable(velocityExtension, Mesh);

    // The base class is the no-op "none" model (identity correction).
    addToRunTimeSelectionTable(velocityExtension, velocityExtension, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityExtension::velocityExtension(const fvMesh& mesh)
:
    mesh_(mesh),
    fvSolution_(mesh),
    levelSetDict_(fvSolution_.subDict("levelSet")),
    velExtDict_(levelSetDict_.subOrEmptyDict("velocityExtension")),
    U_(mesh.lookupObject<volVectorField>
        (velExtDict_.getOrDefault<word>("U", "U"))),
    phi_(mesh.lookupObject<surfaceScalarField>
        (velExtDict_.getOrDefault<word>("flux", "phi"))),
    // zeroGradient BCs so the extension models can solve for Uext (U itself
    // carries non-solvable "calculated" BCs). Values are set from U in correct().
    Uext_
    (
        IOobject
        (
            "Uext",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Uext", U_.dimensions(), vector::zero),
        "zeroGradient"
    ),
    // Named "phi" (but NOT registered) so fvm::div(phiExt, psi) reuses the
    // existing div(phi,psi) scheme without clashing with the registered phi.
    phiExt_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false   // do not register
        ),
        phi_
    )
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::velocityExtension> Foam::velocityExtension::New
(
    const fvMesh& mesh
)
{
    const fvSolution& fvSolution(mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& velExtDict =
        levelSetDict.subOrEmptyDict("velocityExtension");
    const word modelType = velExtDict.getOrDefault<word>("type", "none");
    Info<< "Selecting velocityExtension " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "velocityExtension",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<velocityExtension>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::velocityExtension::correct()
{
    // "none": identity. Keep Uext/phiExt in sync with the base fields
    // (which may have been rescaled by oscillation).
    Uext_ == U_;
    phiExt_ == phi_;
}

// ************************************************************************* //
