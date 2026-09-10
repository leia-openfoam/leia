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
    Foam::stencilBoundsValueBound

Description
    The quasi-monotone clip of the reconstructed value to the arrival cell's
    stencil [min, max]. Reproduces `clipToStencilBounds true` bit for bit,
    including `clipRegion` and `clipKeepExtrema`, both of which stay on
    slReconstruction and are read through its accessors.

    FALSIFIED AS A FIX, and kept so its measurements stay reproducible. Two
    numbers say why. On Popinet's 2D hexahedral translating droplet at N = 64 the
    global clip moves the volume error +30 % and the shape error +5 %, because it
    flattens the level set's own extrema: the apex of the distance cone at the
    droplet centre and the box corners. The exemption written to protect those
    extrema then removes the bound from the cells that need it -- gate G4
    measured 59.2 % of the acting cells to be stencil extrema themselves, and the
    polyhedral failure returned at step 506 against the control's 527.

    The defect is structural, not a tuning problem. A monotone bound cannot
    represent an extremum, and a signed distance field has genuine interior
    extrema on its medial axis. See lipschitzConeValueBound for the bound that
    admits them.

SourceFiles
    stencilBoundsValueBound.C

\*---------------------------------------------------------------------------*/

#ifndef stencilBoundsValueBound_H
#define stencilBoundsValueBound_H

#include "slValueBound.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class stencilBoundsValueBound Declaration
\*---------------------------------------------------------------------------*/

class stencilBoundsValueBound
:
    public slValueBound
{
public:

    //- Runtime type information
    TypeName("stencilBounds");

    //- Construct from mesh and the semiLagrangian sub-dict.
    stencilBoundsValueBound(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~stencilBoundsValueBound() = default;


    // Member Functions

        virtual void interval
        (
            const slReconstruction& recon,
            const label c,
            const point& foot,
            const bool eligible,
            const scalar stencilLo,
            const scalar stencilHi,
            slBoundResult& r
        ) const override;
};


} // End namespace Foam

#endif

// ************************************************************************* //
