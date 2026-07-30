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

Namespace function
    Foam::computePlaneAnchors

Description
    Signed-distance Dirichlet anchors for the plane-anchored redistancer
    models, built from the SAME linear least-squares level-set planes
    (leastSquaresPlaneCoeffs) that the Detrixhe-Aslam / geometric phase
    indicators and the geometric face-area mass flux read. Anchoring the
    redistanced field to these planes realigns psi to the discrete interface
    the rest of the solver already sees (the level-set analog of
    geometricVoF's reconstructed distance function, plicRDF), instead of
    displacing it -- the classical danger of reinitialization.

    Anchor set:
      - NARROW-BAND cells (NarrowBand == 1, the psi sign-change cells): the
        cell's own plane, NORMALIZED (the raw LSQ gradient is not unit --
        it deviates by exactly the |grad psi| drift being corrected):
        dist_c = (n & x_c + d)/|n_raw|. The plane foot
        x_c - dist_c*nHat is also returned (seed for the foot-wave model).
      - FIRST-RING cells (band face-neighbours outside the band): the
        normal-direction-weighted average of the adjacent band cells'
        plane distances evaluated at the ring cell centre, with plicRDF
        weights w = sqr(dHat & nHat), dHat the unit vector toward the
        donor's plane FOOT (Scheufler & Roenby 2019). One ring only: the
        plane extrapolation error grows as d^2*kappa/2, and anchors are
        hard constraints the fill step cannot correct.

    Parallel: the coupledFaceNeighbours cache (matched MPI call, constructed
    unconditionally on all ranks) makes the band-cell planes correct across
    processor boundaries; ring donors across processor boundaries are
    SKIPPED in this version (the affected ring cells are simply not
    anchored; the Eikonal/wave fill compensates). Serial runs are the
    validation target.

SourceFiles
    planeAnchors.C

\*---------------------------------------------------------------------------*/

#ifndef planeAnchors_H
#define planeAnchors_H

#include "fvMesh.H"
#include "volFields.H"

namespace Foam
{

//- Anchor data for plane-anchored redistancing
struct planeAnchorData
{
    //- All anchor cells (band + first ring), duplicate-free
    labelList cells;

    //- SIGNED plane distance per anchor cell (same order as cells)
    scalarField dist;

    //- Subset: band anchor cells (own plane available)
    labelList bandCells;

    //- SIGNED plane distance per band anchor cell
    scalarField bandDist;

    //- Plane foot point per band anchor cell (x_c - dist*nHat)
    vectorField bandFoot;

    //- Unit plane normal per band anchor cell (donor-plane wave payload:
    //  a receiving cell evaluates |(x - foot) & nHat|, the CONTINUOUS
    //  distance to the donor plane, instead of the scalloped distance to
    //  the discrete foot-point cloud)
    vectorField bandNHat;
};


//- Build the anchors from the current psi and narrow band (0/1) field.
//  curvatureCorrection: remove the chord-plane bias of the least-squares
//  fit (~ <t^2>*kappa/2, one-signed interface displacement on curved
//  interfaces) by residual regression over the fit stencil; exactly inert
//  on planar data. First-ring cells anchor the SINGLE best-aligned donor
//  plane's evaluation (selection, not averaging).
planeAnchorData computePlaneAnchors
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand,
    const bool curvatureCorrection = true
);

} // End namespace Foam

#endif

// ************************************************************************* //
