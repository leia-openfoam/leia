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
    Foam::noValueBound

Description
    No bound on the reconstructed value, beyond the runaway cap that predates
    the slValueBound family (stencil mid +- 10*stencil range). Reproduces
    `clipToStencilBounds false` bit for bit, and it is what the sentinel
    `valueBound fromClipSwitch` resolves to in that case.

    This is the default configuration of the method, and it is NOT stable: the
    frozen-velocity power iteration measures rho(B) > 1 on every mesh. It is the
    default because every bound tried so far costs more than it buys, so the
    unbounded operator remains the reference the candidates are scored against.

    Reported as inert(), so the corrector skips the bound diagnostics and the
    transport-spectrum probe may treat the operator as linear.

SourceFiles
    noValueBound.C

\*---------------------------------------------------------------------------*/

#ifndef noValueBound_H
#define noValueBound_H

#include "slValueBound.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noValueBound Declaration
\*---------------------------------------------------------------------------*/

class noValueBound
:
    public slValueBound
{
public:

    //- Runtime type information
    TypeName("none");

    //- Construct from mesh and the semiLagrangian sub-dict.
    noValueBound(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~noValueBound() = default;


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

        virtual bool inert() const override
        {
            return true;
        }
};


} // End namespace Foam

#endif

// ************************************************************************* //
