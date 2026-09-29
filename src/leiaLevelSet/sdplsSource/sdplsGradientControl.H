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
    Foam::sdplsGradientControl

Description
    The gradient-control source of the Eulerian line
    (docs/plan-halo-limited-gradient-control.md, D2):
    \f{equation}{
        D_t \psi = \psi F,
        \qquad
        F = w(q)\, a + G(q, \sigma),
    \f}
    with \f$ q = |\nabla\psi| \f$ from the \c gradPsi strategy of
    Foam::sdplsSource, \f$ a = \mathbf{n}\cdot\nabla\mathbf{u}\cdot\mathbf{n} \f$
    the field R of the base class, and \f$ \sigma = |\mathrm{symm}(\nabla
    \mathbf{u})|_F \f$ (scheme \c gradUSdpls), computed only when the law
    reads it. The law and its strain weight are a Foam::gradientControlLaw
    (libleiaGradientControl) from the \c law sub-dictionary:
    \verbatim
        sdplsSource
        {
            type            gradientControl;
            discretization  simpleLinearImplicit;
            gradPsi         fvc;
            mollifier       { type none; }
            law
            {
                type    linearZ;
                mu      10;
                strainWeight { type none; }
            }
        }
    \endverbatim
    F is the nonLinearPart() of the base class, so the discretization and the
    mollifier strategies apply unchanged. Special cases: law none with strain
    weight full is the source R; law linearQ with mu 1 and strain weight none
    is the source beta with beta = 1 (the bounded-gradient equation);
    law linearQ with strain weight full is the combined source of
    docs/plan-combined-source-terms.md.

    A velocity extension already removes the normal strain from the transport
    of q. With an extension and a strain weight other than none, the strain is
    cancelled twice; the model warns at construction.

SourceFiles
    sdplsGradientControl.C

\*---------------------------------------------------------------------------*/

#ifndef sdplsGradientControl_H
#define sdplsGradientControl_H

#include "sdplsSource.H"
#include "gradientControlLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sdplsGradientControl Declaration
\*---------------------------------------------------------------------------*/

class sdplsGradientControl
:
    public sdplsSource
{
    // Private Data

        //- The law G and its strain weight w
        autoPtr<gradientControlLaw> law_;


public:

    //- Runtime type information
    TypeName("gradientControl");


    // Constructors

        //- Construct from the sdplsSource dictionary
        sdplsGradientControl(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~sdplsGradientControl() = default;


    // Member Functions

        //- The law
        const gradientControlLaw& law() const
        {
            return *law_;
        }


protected:

    // Protected Member Functions

        //- F = w(q) R + G(q, sigma)
        virtual tmp<volScalarField> nonLinearPart
        (
            const volScalarField& R,
            const volScalarField& psi,
            const volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
