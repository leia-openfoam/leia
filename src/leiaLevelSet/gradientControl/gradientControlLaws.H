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

Namespace
    Foam::gradientControlLaws

Description
    The laws G of Foam::gradientControlLaw, with e = q - 1 and
    zeta = (q^2 - 1)/2 = e + e^2/2:

    \table
        type             | G(q)                                | coefficients
        linearQ          | -mu e                               | mu
        linearZ          | -mu zeta                            | mu
        cubicQ           | -mu e^3                             | mu
        cubicZ           | -mu zeta^3                          | mu
        twoThirdsZReg    | -mu zeta/(zeta^2 + eps^2)^(1/6)     | mu, eps
        saturatedLinearZ | -mu zeta/(1 + c|zeta|)              | mu, c
        softWall         | -K tanh[gamma ((q - 1)/deltaS)^p]   | cKappa, deltaS, p, gamma, epsD
    \endtable

    with K = cKappa sqrt(sigma^2 + epsD^2) and an odd integer p.

    The q- and z-laws are those of the companion dossier (Bothe, Maric, Soga):
    for the same mu the linear pair has the same slope -mu at q = 1, and the
    cubic pair the same zero slope. The linear z-law gives the exact logistic
    law z' = mu z (1 - z) for z = q^2 when the strain is cancelled. linearQ
    with mu is also the bounded-gradient equation of Bothe, Fricke and Soga
    (route S0 of the halo-limited dossier, target gradient 1).

    softWall is route S1a of the halo-limited dossier (D. Bothe, 2026-09-26):
    the rate is the local strain rate sigma = |D|_F times cKappa, so that
    the wall outruns the strain that drives q out of the band; epsD > 0 keeps
    a floor rate where D = 0. With deltaS = 0.08, p = 5 and
    gamma = artanh(0.9) the wall has 90 percent of its strength at
    |q - 1| = 0.08 and 99.97 percent at 0.10.

SourceFiles
    gradientControlLaws.C

\*---------------------------------------------------------------------------*/

#ifndef gradientControlLaws_H
#define gradientControlLaws_H

#include "gradientControlLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gradientControlLaws
{

/*---------------------------------------------------------------------------*\
                            Class muLaw Declaration
\*---------------------------------------------------------------------------*/

//- The common part of the laws G = mu g(q): the REQUIRED rate mu > 0.
//  Not selectable.
class muLaw
:
    public gradientControlLaw
{
protected:

        //- The rate [1/s]
        const scalar mu_;

        static scalar q(const scalar q2)
        {
            return Foam::sqrt(Foam::max(q2, scalar(0)));
        }

        static scalar zeta(const scalar q2)
        {
            return 0.5*(q2 - 1);
        }

public:

    explicit muLaw(const dictionary& lawDict);

    virtual ~muLaw() = default;

    virtual bool hasRate() const
    {
        return true;
    }
};


//- G = -mu (q - 1)
class linearQ
:
    public muLaw
{
public:

    TypeName("linearQ");

    explicit linearQ(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const
    {
        return -mu_*(q(q2) - 1);
    }

    virtual scalar dRateDq(const scalar, const scalar) const
    {
        return -mu_;
    }
};


//- G = -mu zeta, zeta = (q^2 - 1)/2
class linearZ
:
    public muLaw
{
public:

    TypeName("linearZ");

    explicit linearZ(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const
    {
        return -mu_*zeta(q2);
    }

    virtual scalar dRateDq(const scalar q2, const scalar) const
    {
        return -mu_*q(q2);
    }
};


//- G = -mu (q - 1)^3
class cubicQ
:
    public muLaw
{
public:

    TypeName("cubicQ");

    explicit cubicQ(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const
    {
        const scalar e = q(q2) - 1;
        return -mu_*e*e*e;
    }

    virtual scalar dRateDq(const scalar q2, const scalar) const
    {
        const scalar e = q(q2) - 1;
        return -3*mu_*e*e;
    }
};


//- G = -mu zeta^3
class cubicZ
:
    public muLaw
{
public:

    TypeName("cubicZ");

    explicit cubicZ(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const
    {
        const scalar z = zeta(q2);
        return -mu_*z*z*z;
    }

    virtual scalar dRateDq(const scalar q2, const scalar) const
    {
        const scalar z = zeta(q2);
        return -3*mu_*z*z*q(q2);
    }
};


//- G = -mu zeta/(zeta^2 + eps^2)^(1/6); the small-error slope is
//  mu eps^(-1/3), and |G| grows as mu |zeta|^(2/3) for |zeta| >> eps
class twoThirdsZReg
:
    public muLaw
{
    //- The regularization eps > 0
    const scalar eps_;

public:

    TypeName("twoThirdsZReg");

    explicit twoThirdsZReg(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const;

    virtual scalar dRateDq(const scalar q2, const scalar) const;
};


//- G = -mu zeta/(1 + c |zeta|); |G| < mu/c for every q
class saturatedLinearZ
:
    public muLaw
{
    //- The saturation c > 0
    const scalar c_;

public:

    TypeName("saturatedLinearZ");

    explicit saturatedLinearZ(const dictionary& lawDict);

    virtual scalar rate(const scalar q2, const scalar) const
    {
        const scalar z = zeta(q2);
        return -mu_*z/(1 + c_*Foam::mag(z));
    }

    virtual scalar dRateDq(const scalar q2, const scalar) const
    {
        const scalar d = 1 + c_*Foam::mag(zeta(q2));
        return -mu_*q(q2)/(d*d);
    }
};


/*---------------------------------------------------------------------------*\
                          Class softWall Declaration
\*---------------------------------------------------------------------------*/

//- G = -cKappa sqrt(sigma^2 + epsD^2) tanh[gamma ((q - 1)/deltaS)^p]
class softWall
:
    public gradientControlLaw
{
    // Private Data

        //- The ratio of the wall rate to the strain rate
        const scalar cKappa_;

        //- The shape width (not a threshold): deltaS = deltaA - eps_q
        const scalar deltaS_;

        //- The odd exponent p >= 1
        const label p_;

        //- The shape factor: tanh(gamma) is the strength at |q - 1| = deltaS
        const scalar gamma_;

        //- The floor rate [1/s] where the strain rate is zero
        const scalar epsD_;


    // Private Member Functions

        //- s^n for an integer n >= 0 (s^p keeps the sign of s for odd p)
        static scalar ipow(const scalar s, const label n);

        //- The rate cKappa sqrt(sigma^2 + epsD^2)
        scalar kappa(const scalar sigma) const
        {
            return cKappa_*Foam::sqrt(sigma*sigma + epsD_*epsD_);
        }

public:

    TypeName("softWall");

    explicit softWall(const dictionary& lawDict);

    virtual ~softWall() = default;

    virtual bool hasRate() const
    {
        return true;
    }

    virtual bool needsStrainRate() const
    {
        return true;
    }

    virtual scalar rate(const scalar q2, const scalar sigma) const;

    virtual scalar dRateDq(const scalar q2, const scalar sigma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace gradientControlLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
