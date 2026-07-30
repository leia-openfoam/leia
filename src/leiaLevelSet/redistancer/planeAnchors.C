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

#include "planeAnchors.H"
#include "levelSetPlaneReconstruction.H"
#include "DynamicList.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

//- Chord-plane (curvature) bias of the least-squares plane fit: for a
//  smooth interface, psi over the stencil is locally
//  psi ~ |n|*(s + kappa*t^2/2) with t the tangential offset -- the fitted
//  plane CONSTANT absorbs the stencil mean of the quadratic term, which
//  displaces the reconstructed interface one-sidedly by ~<t^2>*kappa/2
//  (measured: E_geom == E_vol exactly on the circle test). Regressing the
//  fit RESIDUALS against the centered tangential square over the SAME
//  stencil (cell + face neighbours + coupled) recovers the quadratic
//  coefficient c = |n|*kappa/2 without extra communication; subtracting
//  c*<t^2> from the plane constant removes the bias. Exactly zero on
//  planar data (residuals vanish) -- machine exactness is preserved.
scalar chordBiasCorrection
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const label cellI,
    const coupledFaceNeighbours& coupledNei,
    const vector& nHat,
    const vector& nRaw,
    const scalar dRaw
)
{
    const auto& C = mesh.C();
    const auto& Nc = mesh.cellCells()[cellI];

    // Gather the fit stencil (mirrors leastSquaresPlaneCoeffs exactly).
    DynamicList<vector> xs(Nc.size() + 7);
    DynamicList<scalar> ps(Nc.size() + 7);
    xs.append(C[cellI]); ps.append(psi[cellI]);
    for (const label k : Nc)
    {
        xs.append(C[k]); ps.append(psi[k]);
    }
    const auto& pbm = mesh.boundaryMesh();
    for (const label faceJ : mesh.cells()[cellI])
    {
        if (mesh.isInternalFace(faceJ)) continue;
        const label patchL = pbm.whichPatch(faceJ);
        if (patchL < 0 || !coupledNei.valid(patchL)) continue;
        const label faceP = faceJ - pbm[patchL].start();
        xs.append(coupledNei.C(patchL)[faceP]);
        ps.append(coupledNei.psi(patchL)[faceP]);
    }

    // Stencil centroid; centered tangential offset squares.
    vector xbar(Zero);
    forAll(xs, i) { xbar += xs[i]; }
    xbar /= xs.size();

    scalarField t2(xs.size());
    scalar t2mean = 0;
    forAll(xs, i)
    {
        const vector r = xs[i] - xbar;
        t2[i] = magSqr(r) - sqr(r & nHat);
        t2mean += t2[i];
    }
    t2mean /= xs.size();

    // Residual regression r_k = c*(t2_k - <t2>).
    scalar num = 0, den = 0;
    forAll(xs, i)
    {
        const scalar res = ps[i] - ((nRaw & xs[i]) + dRaw);
        const scalar dt2 = t2[i] - t2mean;
        num += res*dt2;
        den += sqr(dt2);
    }
    if (den < VSMALL)
    {
        return 0;   // no tangential spread (e.g. 1D): nothing to correct
    }
    return (num/den)*t2mean;   // c*<t2>: the constant absorbed by the fit
}

} // End anonymous namespace


planeAnchorData computePlaneAnchors
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand,
    const bool curvatureCorrection
)
{
    const vectorField& C = mesh.cellCentres();
    const labelListList& cellCells = mesh.cellCells();

    // Matched MPI call: construct unconditionally on every rank.
    const coupledFaceNeighbours coupledNei(mesh, psi);

    // --- Pass 1: band-cell planes (normalized) ------------------------------
    boolList hasPlane(mesh.nCells(), false);
    scalarField cellDist(mesh.nCells(), 0.0);   // signed plane distance
    vectorField cellNHat(mesh.nCells(), Zero);
    vectorField cellFoot(mesh.nCells(), Zero);

    DynamicList<label> bandCells(mesh.nCells()/10);

    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] != 1) continue;

        const scalarList pc =
            leastSquaresPlaneCoeffs(mesh, psi, cellI, coupledNei);

        const vector nc(pc[0], pc[1], pc[2]);
        const scalar nmag = mag(nc);
        if (nmag < SMALL)
        {
            // Degenerate reconstruction: this band cell is not anchored.
            continue;
        }

        // Signed-distance plane: dist(x) = (n & x) + d with |n| = 1. The
        // normalization is essential -- the raw LSQ gradient magnitude
        // deviates from 1 by exactly the |grad psi| drift being corrected.
        const vector n = nc/nmag;

        // Curvature correction of the chord-plane bias (see helper above):
        // subtract the quadratic term the plane constant absorbed.
        scalar dRawCorr = pc[3];
        if (curvatureCorrection)
        {
            dRawCorr -=
                chordBiasCorrection(mesh, psi, cellI, coupledNei, n, nc, pc[3]);
        }
        const scalar d = dRawCorr/nmag;

        const scalar distC = (n & C[cellI]) + d;

        hasPlane[cellI] = true;
        cellDist[cellI] = distC;
        cellNHat[cellI] = n;
        cellFoot[cellI] = C[cellI] - distC*n;
        bandCells.append(cellI);
    }

    // --- Pass 2: first-ring cells (weighted donor-plane average) ------------
    // Loop ring cells directly (one visit each): outside the band, with at
    // least one anchored band face-neighbour.
    DynamicList<label> ringCells(bandCells.size());
    DynamicList<scalar> ringDist(bandCells.size());

    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] == 1) continue;

        const labelList& nbrs = cellCells[cellI];

        // SINGLE best-aligned donor (alignment of the direction toward the
        // donor's plane foot with the donor normal, Scheufler & Roenby
        // 2019) -- used to SELECT, not to average: a weighted average over
        // planes that disagree by O(h^2*kappa) injects O(h*kappa) gradient
        // noise between neighbouring anchors, while a single donor plane
        // keeps |grad psi| = 1 exactly within each donor region (same
        // donor-plane evaluation the far-field wave performs).
        scalar wBest = -1;
        scalar dBest = 0;
        for (const label nbrI : nbrs)
        {
            if (!hasPlane[nbrI]) continue;

            const vector& n = cellNHat[nbrI];
            const scalar dist =
                (n & C[cellI]) + (cellDist[nbrI] - (n & C[nbrI]));

            const vector toFoot = cellFoot[nbrI] - C[cellI];
            const scalar toFootMag = mag(toFoot);
            const scalar w =
                (toFootMag > SMALL)
              ? sqr((toFoot/toFootMag) & n) + SMALL
              : 1.0;

            if (w > wBest)
            {
                wBest = w;
                dBest = dist;
            }
        }

        if (wBest > 0)
        {
            ringCells.append(cellI);
            ringDist.append(dBest);
        }
    }

    // --- Assemble ------------------------------------------------------------
    planeAnchorData anchors;

    anchors.bandCells = labelList(bandCells);
    anchors.bandDist.setSize(anchors.bandCells.size());
    anchors.bandFoot.setSize(anchors.bandCells.size());
    anchors.bandNHat.setSize(anchors.bandCells.size());
    forAll(anchors.bandCells, i)
    {
        const label cellI = anchors.bandCells[i];
        anchors.bandDist[i] = cellDist[cellI];
        anchors.bandFoot[i] = cellFoot[cellI];
        anchors.bandNHat[i] = cellNHat[cellI];
    }

    anchors.cells.setSize(anchors.bandCells.size() + ringCells.size());
    anchors.dist.setSize(anchors.cells.size());
    forAll(anchors.bandCells, i)
    {
        anchors.cells[i] = anchors.bandCells[i];
        anchors.dist[i] = anchors.bandDist[i];
    }
    forAll(ringCells, i)
    {
        anchors.cells[anchors.bandCells.size() + i] = ringCells[i];
        anchors.dist[anchors.bandCells.size() + i] = ringDist[i];
    }

    return anchors;
}

} // End namespace Foam

// ************************************************************************* //
