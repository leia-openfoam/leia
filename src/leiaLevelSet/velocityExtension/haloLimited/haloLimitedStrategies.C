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

#include "haloLimitedStrategies.H"
#include "slReconstruction.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(extensionTravel, 0);
    defineRunTimeSelectionTable(extensionTravel, dictionary);
    defineTypeNameAndDebug(extensionWeight, 0);
    defineRunTimeSelectionTable(extensionWeight, dictionary);
    defineTypeNameAndDebug(extensionDirection, 0);
    defineRunTimeSelectionTable(extensionDirection, dictionary);
    defineTypeNameAndDebug(extensionSampler, 0);
    defineRunTimeSelectionTable(extensionSampler, dictionary);

namespace extensionTravels
{
    defineTypeNameAndDebug(capped, 0);
    addToRunTimeSelectionTable(extensionTravel, capped, dictionary);
}
namespace extensionWeights
{
    defineTypeNameAndDebug(fractionReached, 0);
    addToRunTimeSelectionTable(extensionWeight, fractionReached, dictionary);
}
namespace extensionDirections
{
    defineTypeNameAndDebug(levelSet, 0);
    addToRunTimeSelectionTable(extensionDirection, levelSet, dictionary);
}
namespace extensionSamplers
{
    defineTypeNameAndDebug(stencilFit, 0);
    addToRunTimeSelectionTable(extensionSampler, stencilFit, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extensionTravel>
Foam::extensionTravel::New(const dictionary& dict)
{
    const word type = dict.get<word>("type");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "extensionTravel", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<extensionTravel>(ctorPtr(dict));
}

Foam::autoPtr<Foam::extensionWeight>
Foam::extensionWeight::New(const dictionary& dict)
{
    const word type = dict.get<word>("type");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "extensionWeight", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<extensionWeight>(ctorPtr(dict));
}

Foam::autoPtr<Foam::extensionDirection>
Foam::extensionDirection::New(const dictionary& dict)
{
    const word type = dict.get<word>("type");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "extensionDirection", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<extensionDirection>(ctorPtr(dict));
}

Foam::autoPtr<Foam::extensionSampler>
Foam::extensionSampler::New(const fvMesh& mesh, const dictionary& dict)
{
    const word type = dict.get<word>("type");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "extensionSampler", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    return autoPtr<extensionSampler>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extensionTravels::capped::capped(const dictionary& dict)
:
    extensionTravel(dict),
    m_(dict.get<label>("m"))
{
    if (m_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "travel capped needs m >= 1; got " << m_ << exit(FatalIOError);
    }
}


Foam::extensionWeights::fractionReached::fractionReached(const dictionary& dict)
:
    extensionWeight(dict),
    beta_(dict.get<scalar>("beta"))
{
    if (!(beta_ > 0))
    {
        FatalIOErrorInFunction(dict)
            << "weight fractionReached needs beta > 0; got " << beta_
            << exit(FatalIOError);
    }
}


Foam::extensionDirections::levelSet::levelSet(const dictionary& dict)
:
    extensionDirection(dict)
{}


Foam::extensionSamplers::stencilFit::stencilFit
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    extensionSampler(mesh, dict),
    fit_(slReconstruction::New(mesh, "uncachedQuadraticWeightedLeastSquares"))
{}


Foam::extensionSamplers::stencilFit::~stencilFit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::extensionTravels::capped::S
(
    const scalar d,
    const scalar R
) const
{
    // S = d (1 + t^(2m))^(-1/(2m)), t = |d|/R. For t > 1 written as
    // sign(d) R (1 + t^(-2m))^(-1/(2m)): no power overflows, and |S| <= R holds
    // in floating point as well (R times a factor <= 1).
    const scalar t = Foam::mag(d)/R;
    const scalar e = scalar(-1)/(2*m_);
    if (t <= 1)
    {
        return d*Foam::pow(1 + Foam::pow(t, scalar(2*m_)), e);
    }
    return Foam::sign(d)*R*Foam::pow(1 + Foam::pow(1/t, scalar(2*m_)), e);
}


Foam::scalar Foam::extensionTravels::capped::fraction
(
    const scalar d,
    const scalar R
) const
{
    // c = S/d = (1 + t^(2m))^(-1/(2m)); c(0) = 1.
    const scalar t = Foam::mag(d)/R;
    const scalar e = scalar(-1)/(2*m_);
    if (t <= 1)
    {
        return Foam::pow(1 + Foam::pow(t, scalar(2*m_)), e);
    }
    return Foam::pow(1 + Foam::pow(1/t, scalar(2*m_)), e)/t;
}


Foam::scalar Foam::extensionDirections::levelSet::Q2(const scalar q2)
{
    if (q2 >= 0.25)
    {
        return q2;
    }
    const scalar r = 1 - 4*q2;
    return q2 + 0.25*r*r*r;
}


void Foam::extensionDirections::levelSet::evaluate
(
    const scalar psi,
    const vector& g,
    scalar& d,
    vector& e
) const
{
    const scalar s = Foam::sqrt(Q2(magSqr(g)));
    d = psi/s;
    e = g/s;
}


void Foam::extensionSamplers::stencilFit::update(const volScalarField& f)
{
    fit_->update(f);
}


Foam::scalar Foam::extensionSamplers::stencilFit::value
(
    const label celli,
    const point& x
) const
{
    return fit_->evaluateRaw(celli, x);
}


// ************************************************************************* //
