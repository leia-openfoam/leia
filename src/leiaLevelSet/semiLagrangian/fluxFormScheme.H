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
    Foam::fluxFormScheme

Description
    Conservative (flux-form) semi-Lagrangian, CFL <= 1. The face psi-flux is the
    reconstructed psi^n integrated over the swept prism; at CFL <= 1 the prism
    lies in the face-upwind cell, so F_f = phi_f * psiTilde with psiTilde the
    UPWIND cell's reconstruction at the swept-prism centroid x_f - 1/2 u_f dt.
    psi^{n+1} = psi^n - dt div(F) telescopes over internal faces -> discrete
    conservation of the integral of psi (exact for a divergence-free flux). The
    corrector is unused. CFL > 1 would need the multi-cell preimage integral
    (documented scope limit; not implemented).

SourceFiles
    fluxFormScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fluxFormScheme_H
#define fluxFormScheme_H

#include "slScheme.H"

namespace Foam
{

class fluxFormScheme
:
    public slScheme
{
public:

    TypeName("fluxForm");

    explicit fluxFormScheme(const fvMesh& mesh);

    virtual ~fluxFormScheme() = default;

    virtual void advance
    (
        volScalarField& psi,
        const volVectorField& Unew,
        const volVectorField& Uold,
        slReconstruction& recon,
        slCorrector& corrector
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
