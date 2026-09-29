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

Class
    Foam::strainWeight

Description
    The weight \f$ w(q) \f$ of the normal-strain term in the gradient-control
    coefficient
    \f{equation}{
        F = w(q)\, a + G(q),
        \qquad D_t \psi = \psi F,
        \qquad q = |\nabla\psi|,
        \qquad a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n},
    \f}
    (docs/plan-halo-limited-gradient-control.md, section 3.1). On the zero
    level set \f$ D_t q = q\,(F - a) \f$, so \f$ w = 1 \f$ cancels the strain
    exactly (the SDPLS source R) and \f$ w = 0 \f$ leaves it to the law
    \f$ G \f$ or to a velocity extension. A strategy of Foam::gradientControlLaw,
    selected from the \c strainWeight sub-dictionary of the law:
    \verbatim
        strainWeight
        {
            type    omega;      // none | full | omega
            beta    1;          // omega only
            m       2;          // omega only
            deltaS  0.08;       // omega only
        }
    \endverbatim
    The base class is the default \c none (\f$ w = 0 \f$). Every coefficient
    that a model reads is REQUIRED; there is no code default.

SourceFiles
    strainWeight.C

\*---------------------------------------------------------------------------*/

#ifndef strainWeight_H
#define strainWeight_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class strainWeight Declaration
\*---------------------------------------------------------------------------*/

class strainWeight
{
public:

    //- Runtime type information
    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        strainWeight,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from the strainWeight sub-dictionary
        explicit strainWeight(const dictionary& dict);

        //- No copy construct
        strainWeight(const strainWeight&) = delete;

        //- No copy assignment
        void operator=(const strainWeight&) = delete;


    // Selectors

        //- Select from the strainWeight sub-dictionary; a missing
        //  sub-dictionary or type selects none
        static autoPtr<strainWeight> New(const dictionary& dict);


    //- Destructor
    virtual ~strainWeight() = default;


    // Member Functions

        //- True if w is not identically zero, so that the consumer must
        //  compute the normal strain a
        virtual bool active() const
        {
            return false;
        }

        //- The weight w at q^2 = |grad psi|^2
        virtual scalar w(const scalar q2) const
        {
            return 0;
        }
};


namespace strainWeights
{

/*---------------------------------------------------------------------------*\
                            Class full Declaration
\*---------------------------------------------------------------------------*/

//- w = 1: the complete SDPLS strain cancellation (the source R).
class full
:
    public strainWeight
{
public:

    TypeName("full");

    explicit full(const dictionary& dict);

    virtual ~full() = default;

    virtual bool active() const
    {
        return true;
    }

    virtual scalar w(const scalar q2) const
    {
        return 1;
    }
};


/*---------------------------------------------------------------------------*\
                            Class omega Declaration
\*---------------------------------------------------------------------------*/

//- The directional weight of route S1c of the halo-limited dossier:
//  w = 1 - exp[-beta ((q - 1)/deltaS)^(2m)]. w(1) = 0, and w tends to 1 as
//  |q - 1| grows past deltaS, so only the dangerous stretching is cancelled.
class omega
:
    public strainWeight
{
    // Private Data

        //- The rate factor beta > 0
        const scalar beta_;

        //- The half exponent m >= 1 (the exponent is 2m, so w is even in q - 1)
        const label m_;

        //- The width deltaS > 0, shared with the soft wall
        const scalar deltaS_;

public:

    TypeName("omega");

    explicit omega(const dictionary& dict);

    virtual ~omega() = default;

    virtual bool active() const
    {
        return true;
    }

    virtual scalar w(const scalar q2) const;
};

} // End namespace strainWeights

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
