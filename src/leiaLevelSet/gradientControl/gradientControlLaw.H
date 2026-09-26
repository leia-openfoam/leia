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
    Foam::gradientControlLaw

Description
    The source law of the gradient-controlled level set
    (docs/plan-halo-limited-gradient-control.md, sections 3.1 and 4, D1):
    \f{equation}{
        D_t \psi = \psi F,
        \qquad
        F = w(q)\, a + G(q, \sigma),
        \qquad
        q = |\nabla\psi|,
    \f}
    with \f$ a = \mathbf{n}\cdot\mathbf{D}\cdot\mathbf{n} \f$ the normal strain,
    \f$ \sigma = |\mathbf{D}|_F \f$ the Frobenius norm of the strain-rate
    tensor, and \f$ w \f$ a Foam::strainWeight. On the zero level set
    \f$ D_t q = q\,(F - a) \f$, and \f$ G(1, \sigma) = 0 \f$ for every law,
    so a signed-distance level set is a fixed point.

    A law is pure scalar algebra: it evaluates \f$ G \f$ and
    \f$ \partial G/\partial q \f$ at one point. The consumers own the
    discretization -- sdplsGradientControl (Eulerian line, Phase D2) and
    slGradientControlSource (semi-Lagrangian line, Phase D3) -- and read the
    law from the \c law sub-dictionary of their source dictionary:
    \verbatim
        law
        {
            type        linearZ;    // none | linearQ | linearZ | cubicQ |
                                    // cubicZ | twoThirdsZReg |
                                    // saturatedLinearZ | softWall
            mu          10;         // [1/s]
            strainWeight { type none; }
        }
    \endverbatim
    The base class is the default \c none: \f$ G = 0 \f$. A law reads only its
    own coefficients, and every coefficient that it reads is REQUIRED (no code
    default), so a case cannot run a law with a value nobody chose.

    The q- and z-laws follow the companion dossier (Bothe, Maric, Soga), with
    \f$ e = q - 1 \f$ and \f$ \zeta = (q^2 - 1)/2 \f$. The normalization makes
    every pair agree to first order in \f$ e \f$ for the same \f$ \mu \f$.

SourceFiles
    gradientControlLaw.C

\*---------------------------------------------------------------------------*/

#ifndef gradientControlLaw_H
#define gradientControlLaw_H

#include "strainWeight.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class gradientControlLaw Declaration
\*---------------------------------------------------------------------------*/

class gradientControlLaw
{
protected:

    // Protected Data

        //- The weight w(q) of the normal-strain term
        autoPtr<strainWeight> strainWeight_;


public:

    //- Runtime type information
    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        gradientControlLaw,
        dictionary,
        (const dictionary& lawDict),
        (lawDict)
    );


    // Constructors

        //- Construct from the law sub-dictionary; builds the strain weight
        //  from its strainWeight sub-dictionary
        explicit gradientControlLaw(const dictionary& lawDict);

        //- No copy construct
        gradientControlLaw(const gradientControlLaw&) = delete;

        //- No copy assignment
        void operator=(const gradientControlLaw&) = delete;


    // Selectors

        //- Select from the law sub-dictionary (the entry type is REQUIRED)
        static autoPtr<gradientControlLaw> New(const dictionary& lawDict);


    //- Destructor
    virtual ~gradientControlLaw() = default;


    // Member Functions

        //- True if G is not identically zero
        virtual bool hasRate() const
        {
            return false;
        }

        //- True if rate() reads sigma = |D|_F; the consumer computes the
        //  strain-rate tensor only then (or for the strain weight)
        virtual bool needsStrainRate() const
        {
            return false;
        }

        //- G(q^2, sigma) [1/s]; zero at q^2 = 1 for every law
        virtual scalar rate(const scalar q2, const scalar sigma) const
        {
            return 0;
        }

        //- dG/dq at (q^2, sigma) [1/s], for the stiffness number
        //  Co_B = dt |q dG/dq| of the source update
        virtual scalar dRateDq(const scalar q2, const scalar sigma) const
        {
            return 0;
        }

        //- The strain weight
        const strainWeight& weight() const
        {
            return *strainWeight_;
        }

        //- True if F is not identically zero
        bool active() const
        {
            return hasRate() || strainWeight_->active();
        }

        //- True if F reads the normal strain a
        bool needsNormalStrain() const
        {
            return strainWeight_->active();
        }

        //- F = w(q) a + G(q, sigma) [1/s]
        scalar F(const scalar q2, const scalar sigma, const scalar a) const
        {
            return strainWeight_->w(q2)*a + rate(q2, sigma);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
