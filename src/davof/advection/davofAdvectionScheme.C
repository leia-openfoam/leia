/*---------------------------------------------------------------------------*\
    leia :: davof :: davofAdvectionScheme (see davofAdvectionScheme.H)
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    Derived from TwoPhaseFlow's advectionSchemes base class:
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "davofAdvectionScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(davofAdvectionScheme, 0);
    defineRunTimeSelectionTable(davofAdvectionScheme, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofAdvectionScheme::davofAdvectionScheme
(
    const word& type,
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    IOdictionary
    (
        IOobject
        (
            "davofAdvectionScheme",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    advectionSchemesCoeffs_(alpha1.mesh().solverDict(alpha1.name())),
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    alphaPhi_
    (
        IOobject
        (
            "alphaPhi_",
            alpha1_.mesh().time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("zero", dimVol/dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::davofAdvectionScheme::modelDict() const
{
    return advectionSchemesCoeffs_;
}


Foam::dictionary& Foam::davofAdvectionScheme::modelDict()
{
    return advectionSchemesCoeffs_;
}


// ************************************************************************* //
