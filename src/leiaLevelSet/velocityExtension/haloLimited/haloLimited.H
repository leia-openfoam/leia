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
    Foam::velocityExtensions::haloLimited

Description
    The halo-limited directional velocity extension (halo-limited dossier,
    D. Bothe, 2026-09-26; docs/plan-halo-limited-gradient-control.md, D4):
    \f{equation}{
        \mathbf{u}_H(\mathbf{x}) = [1 - w]\,\mathbf{u}(\mathbf{x})
                                 + w\,\mathbf{u}(\mathbf{Y}),
        \qquad
        \mathbf{Y} = \mathbf{x} - S_R(d)\,\mathbf{e},
        \qquad
        w = w(S_R(d)/d),
    \f}
    with (d, e) from the direction strategy, S_R from the travel strategy
    (|S_R| < R = radiusCells h), w from the weight strategy and u(Y) from the
    sampler strategy. On the interface Y = x, w = 1 and the normal derivative
    of Y vanishes, so n.grad(u_H).n = 0 there: the extension removes the normal
    strain from the transport of |grad psi| without a closest-point search, and
    no sample point leaves the radius R.

    Flux, CORRECTION form (the pressure-corrected flux stays where w = 0):
    \f{equation}{
        \Phi^H_f = \Phi_f
                 + w_f\,[\mathbf{u}_h(\mathbf{Y}_f) - \mathbf{u}_h(\mathbf{x}_f)]
                   \cdot\mathbf{S}_f,
    \f}
    where u_h is the mean of the owner's and the neighbour's local models, both
    evaluated at the same two points; the constant part of each model cancels.
    Physical patches keep phi. On a coupled face each side evaluates its own
    cell's model and the two corrections are averaged in one orientation, so
    both sides carry exactly opposite fluxes. Uext_c = (1 - w_c) U_c + w_c
    u_h(Y_c) with the cell's own model (for the cellCentred trace and for
    diagnostics).

    The cell size h_c is the smallest centre distance to a face neighbour,
    coupled neighbours included (the interfaceExtension definition); a face
    uses min(h_owner, h_neighbour). Neither depends on the decomposition.

    Fields (written): hlWeight, hlDistance (d), hlReach (|Y - x|/R, always
    < 1), hlCorrection (the face flux correction). Per call the log reports
    max hlReach and the volume-weighted L2 of div(phiExt) - div(phi)
    (collective reductions on every rank).

    \verbatim
        velocityExtension
        {
            type          haloLimited;
            radiusCells   1;                                // 0 < R/h <= 1
            travel        { type capped;          m    2; }
            weight        { type fractionReached; beta 1; }
            direction     { type levelSet; }
            sampler       { type stencilFit; }
        }
    \endverbatim
    The gradient of psi uses the scheme gradPsiExtension.

SourceFiles
    haloLimited.C

\*---------------------------------------------------------------------------*/

#ifndef haloLimited_H
#define haloLimited_H

#include "velocityExtension.H"
#include "haloLimitedStrategies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityExtensions
{

/*---------------------------------------------------------------------------*\
                         Class haloLimited Declaration
\*---------------------------------------------------------------------------*/

class haloLimited
:
    public velocityExtension
{
    // Private Data

        //- The level-set field
        const volScalarField& psi_;

        //- The sampling radius in cell sizes, 0 < radiusCells <= 1
        const scalar radiusCells_;

        autoPtr<extensionTravel> travel_;
        autoPtr<extensionWeight> weight_;
        autoPtr<extensionDirection> direction_;
        autoPtr<extensionSampler> sampler_;

        //- The cell size h_c (a field, so that coupled patches carry the
        //  neighbour's value)
        volScalarField cellSize_;

        volScalarField hlWeight_;
        volScalarField hlDistance_;
        volScalarField hlReach_;
        surfaceScalarField hlCorrection_;


    // Private Member Functions

        void computeCellSize();


public:

    //- Runtime type information
    TypeName("haloLimited");


    // Constructors

        explicit haloLimited(const fvMesh& mesh);


    //- Destructor
    virtual ~haloLimited() = default;


    // Member Functions

        //- Update Uext and phiExt from U, phi and psi
        virtual void correct();

        //- The sampling radius in cell sizes
        scalar radiusCells() const
        {
            return radiusCells_;
        }

        const extensionTravel& travelStrategy() const
        {
            return *travel_;
        }

        const extensionWeight& weightStrategy() const
        {
            return *weight_;
        }

        const extensionDirection& directionStrategy() const
        {
            return *direction_;
        }

        //- The cell size h_c
        const volScalarField& cellSize() const
        {
            return cellSize_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace velocityExtensions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
