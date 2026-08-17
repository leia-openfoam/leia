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
    Foam::sdplsRdiv

Description
    The strain-cancelling SDPLS source in CONSERVATIVE (divergence) form.

    Two exact identities turn the algebraic source into a surface divergence:
    \f$ \nabla_{\!s}\psi = 0 \f$ (the surface gradient of a level set vanishes,
    since \f$ \nabla\psi \parallel \nhat \f$) and, for incompressible flow,
    \f$ a = \nhat\cdot\nabla\mathbf{v}\cdot\nhat = -\nabla_{\!s}\cdot\mathbf{v} \f$
    (the rate of interfacial area generation). Together:
    \f$ S = a\,\psi = -\nabla_{\!s}\cdot(\psi\mathbf{v}) \f$, no approximation.

    Writing \f$ \nabla_{\!s}\cdot \f$ as the volumetric divergence minus its
    normal--normal projection and expanding the remainder with
    \f$ \mathbf{w} \equiv (\nhat\cdot\mathbf{v})\,\nhat \f$ (the normal
    velocity), the equation this source assembles is, relative to the solver's
    base transport \f$ \partial_t\psi + \nabla\cdot(\phi\psi)
    - \psi\nabla\cdot\phi = S_{\mathrm{matrix}} \f$:

    \f{equation}{
        S_{\mathrm{matrix}} \;=\;
            \nabla\cdot(\phi_w\,\psi) \;-\; \nabla\cdot(\phi\,\psi)
            \;-\; \psi\,\bigl(\nabla\cdot\mathbf{w} - a\bigr),
    \f}

    all terms implicit in \f$ \psi \f$ (two \c fvm::div and one \c fvm::Sp),
    with only \f$ \nhat \f$ lagged -- the same Picard structure as the
    algebraic source, already measured converged at three passes. At the
    continuum, \f$ \nabla\cdot(\phi\psi)=\mathbf{w}\cdot\nabla\psi+\psi\nabla
    \cdot\mathbf{v} \f$ on a level set, and the expression collapses to
    \f$ \psi\,(a - \nabla\cdot\mathbf{v}) \f$ -- the algebraic source plus the
    same compressibility correction the base equation already carries.

    WHY. The coupled stationary-droplet studies measured the algebraic source
    steadily DRAINING the droplet (relative volume error \f$ -1 \f$ by
    \f$ t=0.1 \f$ at every resolution) while sourceless transport conserves to
    \f$ 10^{-5} \f$; the drain survives the topological cut-off and the
    production (fitted) curvature, the CSF sign pairing is verified, the
    assembly matches the verified kinematic form, and the strain is recomputed
    every assembly. The failure mode is a conservation failure, and this form
    removes it structurally: what enters through \c fvm::div fluxes cannot
    drain.

    THE TRADE, stated up front: the algebraic source vanishes identically on
    \f$ \psi=0 \f$ -- pointwise-zero interface displacement. Here that property
    holds analytically but not cell-by-cell: discrete fluxes into an interface
    cell need not cancel. Conservation is bought at the price of that pointwise
    guarantee, which is exactly what the coupled measurement says is the right
    trade.

    The mollifier is deliberately NOT composed with this form: a cut-off
    multiplying a flux breaks the conservation this class exists to provide.
    A dictionary selecting one gets a warning, not silent acceptance.

    Selected as:
    \verbatim
        sdplsSource
        {
            type            Rdiv;
            gradPsi         fvc;      // the normal's source, as for R
            // discretization is IGNORED: the assembly is fixed by the form.
        }
    \endverbatim

SourceFiles
    sdplsRdiv.C

\*---------------------------------------------------------------------------*/

#ifndef sdplsRdiv_H
#define sdplsRdiv_H

#include "sdplsR.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sdplsRdiv Declaration
\*---------------------------------------------------------------------------*/

class sdplsRdiv
:
    public sdplsR
{

public:

    // Static Data Members

        TypeName("Rdiv");

    // Constructors

        //- Construct from components
        sdplsRdiv(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~sdplsRdiv() = default;

    // Member Functions

        //- The conservative assembly (see Description). Overrides the
        //  algebraic-coefficient path entirely; discretization() is unused.
        virtual tmp<fvScalarMatrix>
        fvmsdplsSource(const volScalarField& psi, const volVectorField& U);
};

} // End namespace Foam

#endif

// ************************************************************************* //
