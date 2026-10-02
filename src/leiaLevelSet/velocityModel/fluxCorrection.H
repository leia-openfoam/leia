/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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

Global
    Foam::correctFlux

Description
    Helmholtz/projection enforcement of a divergence-free volumetric flux:
    solves laplacian(p) = div(phi) and sets phi -= grad p (= pEqn.flux()).
    Moved out of the solver header so both the solver (-fluxCorrection) and the
    velocityExtension models (projectFlux) can use it. Since 2026-10-01 the
    prescribed velocity models use it too (velocityModel.fluxProjection
    helmholtz, the polyhedral meshes): there the face-centre flux
    U(x_f) & S_f of an analytically solenoidal field is not discretely
    divergence-free.

    ONE reference cell for the whole decomposed mesh. The pure-Neumann
    Laplacian is pinned at a single cell, on the master rank only. Until
    2026-10-01 every rank pinned its own local cell 0, so a run on N ranks
    solved a problem with N pinned cells: the projected flux kept a
    divergence at N - 1 cells and depended on the decomposition.

    fluxDivergenceReport prints three norms of div(phi): the maximum and the
    volume-weighted mean of |div(phi)| [1/s], and the largest non-cancelling
    fraction |sum_f phi_f| / sum_f |phi_f| of a cell [-]. Collective: call it
    on every rank.

SourceFiles
    fluxCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrection_H
#define fluxCorrection_H

#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Project phi onto a divergence-free space (in place).
void correctFlux(surfaceScalarField& phi);

//- Print max and mean |div(phi)| and the largest non-cancelling flux fraction.
void fluxDivergenceReport(const surfaceScalarField& phi, const char* label);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
