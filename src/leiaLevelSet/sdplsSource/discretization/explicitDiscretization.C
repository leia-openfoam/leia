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

#include "explicitDiscretization.H"
#include "addToRunTimeSelectionTable.H"
#include "sdplsSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitDiscretization, false);
    addToRunTimeSelectionTable(discretization, explicitDiscretization, Dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::explicitDiscretization::explicitDiscretization()
    :
        discretization()
{}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

Foam::tmp<scalarField> 
Foam::explicitDiscretization::
Sc(const volScalarField& nonLinearPart, const volScalarField& psi) const
{
    // OWNING copy, deliberately.
    //
    // The previous form was
    //     return (nonLinearPart * psi)().field();
    // which is a use-after-free: `nonLinearPart * psi` yields a
    // tmp<volScalarField> whose lifetime ends with the full-expression, while
    // `.field()` hands back a const reference into it. Binding that reference
    // to the returned tmp<scalarField> produces a NON-OWNING tmp pointing at
    // memory released before discretize() ever reads it, so the values then
    // multiplied by the cell volumes were whatever was left in the freed
    // block -- plausible-looking, allocator-dependent garbage.
    //
    // Caught by leiaTestSdplsSource: on an exactly affine field, `explicit`
    // disagreed with strictNegativeSpLinearImplicit in the regime where the
    // two have IDENTICAL Sc and Sp by construction (f_nl > 0), which no
    // legitimate discretization difference can produce.
    return tmp<scalarField>
    (
        new scalarField((nonLinearPart*psi)().primitiveField())
    );
}

Foam::tmp<scalarField> 
Foam::explicitDiscretization::Sp(const volScalarField& nonLinearPart) const
{
    return tmp<scalarField>(new scalarField(nonLinearPart.size(), 0.0));
}

// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
