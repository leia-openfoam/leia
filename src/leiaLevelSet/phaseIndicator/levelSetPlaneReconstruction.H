/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::leastSquaresPlaneCoeffs

Description
    Linear least-squares reconstruction of the level-set field in a single
    cell:

    \f[
        \psi^l_c(\mathbf{x}) = \mathbf{n}_c \cdot \mathbf{x} + d_c
    \f]

    The fit uses the cell-centre values of \c psi at \c cellI, its
    cell-neighbours, and (in parallel) the neighbours across coupled patches.
    Returned as the coefficient list \c {n_x, n_y, n_z, d}.

    Shared by geometricPhaseIndicator and detrixheAslamPhaseIndicator so both
    use the identical reconstruction; only the volume-fraction step differs.

SourceFiles
    levelSetPlaneReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetPlaneReconstruction_H
#define levelSetPlaneReconstruction_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "scalarList.H"
#include "scalarField.H"
#include "vectorField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Neighbour \c psi and cell-centre values across coupled (processor) patches.
//
//  \c patchNeighbourField() is a *paired* MPI exchange: every rank sharing a
//  processor patch must call it the same number of times, in the same order.
//  Calling it per narrow-band cell (inside the reconstruction loop) diverges
//  across ranks -- each rank's narrow-band cells on a boundary differ -- which
//  deadlocks or aborts with MPI_ERR_TRUNCATE. This evaluates it ONCE for every
//  coupled patch up front (a matched call on all ranks); the per-cell
//  reconstruction then only indexes the cached fields by [patchI][patchFaceI].
//  Non-coupled patches are left unset.
class coupledFaceNeighbours
{
    PtrList<scalarField> psiNei_;
    PtrList<vectorField> CNei_;

public:

    coupledFaceNeighbours(const fvMesh& mesh, const volScalarField& psi);

    //- True if patchI is a coupled patch with cached neighbour data.
    bool valid(const label patchI) const { return psiNei_.set(patchI); }

    //- Neighbour-side psi on coupled patchI (indexed by patch-local face).
    const scalarField& psi(const label patchI) const { return psiNei_[patchI]; }

    //- Neighbour-side cell centres on coupled patchI (patch-local face).
    const vectorField& C(const label patchI) const { return CNei_[patchI]; }
};


//- Least-squares plane coefficients {n_x, n_y, n_z, d} for psi in cellI.
//  Coupled-patch neighbours are supplied precomputed (see coupledFaceNeighbours)
//  so this routine performs no MPI and is safe to call per cell.
scalarList leastSquaresPlaneCoeffs
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const label cellI,
    const coupledFaceNeighbours& coupledNei
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
