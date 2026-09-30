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

#include "gradientControlLaws.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

namespace Foam
{
namespace gradientControlLaws
{
    defineTypeNameAndDebug(linearQ, 0);
    addToRunTimeSelectionTable(gradientControlLaw, linearQ, dictionary);

    defineTypeNameAndDebug(linearZ, 0);
    addToRunTimeSelectionTable(gradientControlLaw, linearZ, dictionary);

    defineTypeNameAndDebug(cubicQ, 0);
    addToRunTimeSelectionTable(gradientControlLaw, cubicQ, dictionary);

    defineTypeNameAndDebug(cubicZ, 0);
    addToRunTimeSelectionTable(gradientControlLaw, cubicZ, dictionary);

    defineTypeNameAndDebug(twoThirdsZReg, 0);
    addToRunTimeSelectionTable(gradientControlLaw, twoThirdsZReg, dictionary);

    defineTypeNameAndDebug(saturatedLinearZ, 0);
    addToRunTimeSelectionTable(gradientControlLaw, saturatedLinearZ, dictionary);

    defineTypeNameAndDebug(softWall, 0);
    addToRunTimeSelectionTable(gradientControlLaw, softWall, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientControlLaws::muLaw::muLaw(const dictionary& lawDict)
:
    gradientControlLaw(lawDict),
    mu_(lawDict.get<scalar>("mu"))
{
    if (!(mu_ > 0))
    {
        FatalIOErrorInFunction(lawDict)
            << "gradient-control law " << lawDict.get<word>("type")
            << " needs mu > 0 [1/s]; got " << mu_ << exit(FatalIOError);
    }
}


Foam::gradientControlLaws::linearQ::linearQ(const dictionary& lawDict)
:
    muLaw(lawDict)
{}


Foam::gradientControlLaws::linearZ::linearZ(const dictionary& lawDict)
:
    muLaw(lawDict)
{}


Foam::gradientControlLaws::cubicQ::cubicQ(const dictionary& lawDict)
:
    muLaw(lawDict)
{}


Foam::gradientControlLaws::cubicZ::cubicZ(const dictionary& lawDict)
:
    muLaw(lawDict)
{}


Foam::gradientControlLaws::twoThirdsZReg::twoThirdsZReg
(
    const dictionary& lawDict
)
:
    muLaw(lawDict),
    eps_(lawDict.get<scalar>("eps"))
{
    if (!(eps_ > 0))
    {
        FatalIOErrorInFunction(lawDict)
            << "twoThirdsZReg needs eps > 0; got " << eps_ << exit(FatalIOError);
    }
}


Foam::gradientControlLaws::saturatedLinearZ::saturatedLinearZ
(
    const dictionary& lawDict
)
:
    muLaw(lawDict),
    c_(lawDict.get<scalar>("c"))
{
    if (!(c_ > 0))
    {
        FatalIOErrorInFunction(lawDict)
            << "saturatedLinearZ needs c > 0; got " << c_ << exit(FatalIOError);
    }
}


Foam::gradientControlLaws::softWall::softWall(const dictionary& lawDict)
:
    gradientControlLaw(lawDict),
    cKappa_(lawDict.get<scalar>("cKappa")),
    deltaS_(lawDict.get<scalar>("deltaS")),
    p_(lawDict.get<label>("p")),
    gamma_(lawDict.get<scalar>("gamma")),
    epsD_(lawDict.get<scalar>("epsD"))
{
    if (!(cKappa_ > 0) || !(deltaS_ > 0) || !(gamma_ > 0) || epsD_ < 0)
    {
        FatalIOErrorInFunction(lawDict)
            << "softWall needs cKappa > 0, deltaS > 0, gamma > 0 and"
            << " epsD >= 0; got cKappa " << cKappa_ << ", deltaS " << deltaS_
            << ", gamma " << gamma_ << ", epsD " << epsD_ << exit(FatalIOError);
    }
    if (p_ < 1 || p_ % 2 == 0)
    {
        FatalIOErrorInFunction(lawDict)
            << "softWall needs an odd exponent p >= 1, so that the wall"
            << " pushes q back from both sides; got p " << p_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gradientControlLaws::twoThirdsZReg::rate
(
    const scalar q2,
    const scalar
) const
{
    const scalar z = zeta(q2);
    return -mu_*z/Foam::pow(z*z + eps_*eps_, scalar(1)/6);
}


Foam::scalar Foam::gradientControlLaws::twoThirdsZReg::dRateDq
(
    const scalar q2,
    const scalar
) const
{
    // dG/dzeta = -mu (2/3 zeta^2 + eps^2)/(zeta^2 + eps^2)^(7/6); dzeta/dq = q.
    const scalar z = zeta(q2);
    const scalar r = z*z + eps_*eps_;
    return -mu_*q(q2)*(scalar(2)/3*z*z + eps_*eps_)/Foam::pow(r, scalar(7)/6);
}


Foam::scalar Foam::gradientControlLaws::softWall::ipow
(
    const scalar s,
    const label n
)
{
    scalar r = 1;
    for (label i = 0; i < n; ++i)
    {
        r *= s;
    }
    return r;
}


Foam::scalar Foam::gradientControlLaws::softWall::rate
(
    const scalar q2,
    const scalar sigma
) const
{
    const scalar s = (Foam::sqrt(Foam::max(q2, scalar(0))) - 1)/deltaS_;
    return -kappa(sigma)*Foam::tanh(gamma_*ipow(s, p_));
}


Foam::scalar Foam::gradientControlLaws::softWall::dRateDq
(
    const scalar q2,
    const scalar sigma
) const
{
    // d/ds tanh(gamma s^p) = gamma p s^(p-1) (1 - tanh^2); ds/dq = 1/deltaS.
    const scalar s = (Foam::sqrt(Foam::max(q2, scalar(0))) - 1)/deltaS_;
    const scalar t = Foam::tanh(gamma_*ipow(s, p_));
    return -kappa(sigma)*gamma_*p_*ipow(s, p_ - 1)*(1 - t*t)/deltaS_;
}


// ************************************************************************* //
