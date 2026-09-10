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

#include "stencilBoundsValueBound.H"
#include "slReconstruction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stencilBoundsValueBound, 0);
    addToRunTimeSelectionTable(slValueBound, stencilBoundsValueBound, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stencilBoundsValueBound::stencilBoundsValueBound
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    slValueBound(mesh, dict)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::stencilBoundsValueBound::interval
(
    const slReconstruction& recon,
    const label c,
    const point& foot,
    const bool eligible,
    const scalar stencilLo,
    const scalar stencilHi,
    slBoundResult& r
) const
{
    bool bound = eligible;

    // A quasi-monotone bound cannot represent an extremum. Where psi_c is already
    // the extremum of its own stencil, stencilLo (or stencilHi) IS psi_c, so every
    // reconstructed value beyond it is pulled back to psi_c and the extremum is
    // flattened -- at every step, for as long as the run lasts. The level set has
    // such extrema by construction: the apex of the distance cone at the droplet
    // centre, and the box corner farthest from the interface. The fit undershoots
    // at the apex because a smooth quadratic cannot follow a non-differentiable
    // minimum, so the clip fires there on EVERY mesh, hexahedral meshes included.
    //
    // The test is exact in floating point and carries no coefficient: the range is
    // the min and max over a stencil that CONTAINS psi_c, so psi_c == stencilLo
    // holds bitwise when the cell is the minimum. A FLAT stencil is not an
    // extremum -- the clip must keep enforcing exactness on a constant field.
    //
    // MEASURED, gate G4: the exemption and the bound are in direct conflict.
    // 59.2 % of the cells the bound must act on are themselves stencil extrema,
    // so the exemption removes the bound where it is needed and the failure
    // returns. Kept because the arms that measured it must stay reproducible.
    if (bound && recon.clipKeepExtrema() && stencilLo != stencilHi)
    {
        const scalar psiC = recon.stencilCellValue(c);
        if (psiC == stencilLo || psiC == stencilHi)
        {
            bound = false;
        }
    }

    if (bound)
    {
        // The clip is a BOUND, not a strength: it carries no coefficient and it
        // changes the value ONLY when the fit puts it outside the stencil range,
        // which is exactly when the update creates a new extremum. A cell whose
        // reconstruction stays inside its stencil bounds is bit-unchanged.
        r.lo = stencilLo;
        r.hi = stencilHi;
        r.enforced = true;
    }
    else
    {
        runawayCap(stencilLo, stencilHi, r);
    }
}

// ************************************************************************* //
