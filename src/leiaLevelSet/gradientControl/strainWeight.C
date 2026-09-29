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

#include "strainWeight.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{
    defineTypeNameAndDebug(strainWeight, 0);
    defineRunTimeSelectionTable(strainWeight, dictionary);
    addToRunTimeSelectionTable(strainWeight, strainWeight, dictionary);

namespace strainWeights
{
    defineTypeNameAndDebug(full, 0);
    addToRunTimeSelectionTable(strainWeight, full, dictionary);

    defineTypeNameAndDebug(omega, 0);
    addToRunTimeSelectionTable(strainWeight, omega, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::strainWeight::strainWeight(const dictionary&)
{}


Foam::strainWeights::full::full(const dictionary& dict)
:
    strainWeight(dict)
{}


Foam::strainWeights::omega::omega(const dictionary& dict)
:
    strainWeight(dict),
    beta_(dict.get<scalar>("beta")),
    m_(dict.get<label>("m")),
    deltaS_(dict.get<scalar>("deltaS"))
{
    if (!(beta_ > 0) || m_ < 1 || !(deltaS_ > 0))
    {
        FatalIOErrorInFunction(dict)
            << "strainWeight omega needs beta > 0, m >= 1 and deltaS > 0; got"
            << " beta " << beta_ << ", m " << m_ << ", deltaS " << deltaS_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::strainWeight> Foam::strainWeight::New
(
    const dictionary& dict
)
{
    const word type = dict.getOrDefault<word>("type", strainWeight::typeName);

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "strainWeight", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<strainWeight>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::strainWeights::omega::w(const scalar q2) const
{
    const scalar s = (Foam::sqrt(Foam::max(q2, scalar(0))) - 1)/deltaS_;
    return 1 - Foam::exp(-beta_*Foam::pow(Foam::sqr(s), scalar(m_)));
}


// ************************************************************************* //
