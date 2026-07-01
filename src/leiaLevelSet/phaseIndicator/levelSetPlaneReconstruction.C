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

\*---------------------------------------------------------------------------*/

#include "levelSetPlaneReconstruction.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "coupledFvPatch.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::coupledFaceNeighbours::coupledFaceNeighbours
(
    const fvMesh& mesh,
    const volScalarField& psi
)
:
    psiNei_(mesh.boundary().size()),
    CNei_(mesh.boundary().size())
{
    const auto& patches = mesh.boundary();
    const auto& psiBf = psi.boundaryField();
    const auto& CBf = mesh.C().boundaryField();

    // One matched patchNeighbourField() call per coupled patch on every rank.
    forAll(patches, patchI)
    {
        if (isA<coupledFvPatch>(patches[patchI]))
        {
            psiNei_.set(patchI, psiBf[patchI].patchNeighbourField().ptr());
            CNei_.set(patchI, CBf[patchI].patchNeighbourField().ptr());
        }
    }
}


Foam::scalarList Foam::leastSquaresPlaneCoeffs
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const label cellI,
    const coupledFaceNeighbours& coupledNei
)
{
    // Assemble the 4x4 LLSQ linear system for psi^l(x) = nc & x + dc.
    simpleMatrix LLSQ
    (
        4  /* size 4x4 */,
        0. /* init coeff value*/,
        0. /* init source value*/
    );
    scalarList planeCoeffs(4, 0);
    scalarList LLSQsource(4, 0);

    const auto& cellCells = mesh.cellCells();
    const auto& cellCenters = mesh.C();

    // NOTE: Nc doesn't contain cellI in OpenFOAM: cellI contributions
    // are handled additionally.
    const auto& Nc = cellCells[cellI];

    // equations 0,1,2
    // \frac{\partial e^{lsq}_c} {\partial n_{c,cmpt}} = 0
    for (char row = 0; row < 3; ++row)
    {
        // Contributions from cellI : not in Nc
        for (char col = 0; col < 3; ++col)
            LLSQ(row, col) +=
                cellCenters[cellI][col]*cellCenters[cellI][row];
        LLSQ(row, 3) += cellCenters[cellI][row];
        LLSQsource[row] += psi[cellI]*cellCenters[cellI][row];

        // Contributions from Nc
        forAll(Nc, cellK)
        {
            for (char col = 0; col < 3; ++col)
                LLSQ(row, col) +=
                    cellCenters[Nc[cellK]][col]*cellCenters[Nc[cellK]][row];
            LLSQ(row, 3) += cellCenters[Nc[cellK]][row];
            LLSQsource[row] += psi[Nc[cellK]]*cellCenters[Nc[cellK]][row];
        }
    }
    // equation 3
    // \frac{\partial e^{lsq}_c}{\partial d_c} = 0
    for (char col = 0; col < 3; ++col)
        LLSQ(3, col) += cellCenters[cellI][col];
    LLSQsource[3] += psi[cellI];

    forAll(Nc, cellK)
    {
        for (char col = 0; col < 3; ++col)
            LLSQ(3, col) += cellCenters[Nc[cellK]][col];
        LLSQsource[3] += psi[Nc[cellK]];
    }
    // - dc coefficient is |Nc|: in OpenFOAM it is |Nc| + 1 (for cellI)
    LLSQ(3, 3) = Nc.size() + 1;

    // Contributions from face-adjacent cells across coupled (processor)
    // boundaries. The neighbour psi / cell-centre values are precomputed once
    // per reconstruction pass (coupledFaceNeighbours) -- no MPI happens here, so
    // this is safe to call for whichever narrow-band cells this rank owns.
    const auto& pbm = mesh.boundaryMesh();
    const auto& nBandCell = mesh.cells()[cellI];
    forAll(nBandCell, faceI)
    {
        const label faceJ = nBandCell[faceI];
        if (mesh.isInternalFace(faceJ))
        {
            continue;
        }

        // Owning boundary patch of faceJ, and its patch-local face index.
        const label patchL = pbm.whichPatch(faceJ);
        if ((patchL < 0) || !coupledNei.valid(patchL))
        {
            continue;   // physical boundary (not coupled): no neighbour cell
        }
        const label faceP = faceJ - pbm[patchL].start();

        const vector& cellCenter = coupledNei.C(patchL)[faceP];
        const scalar  psiValue   = coupledNei.psi(patchL)[faceP];

        for (char row = 0; row < 3; ++row)
        {
            for (char col = 0; col < 3; ++col)
                LLSQ(row, col) += cellCenter[col]*cellCenter[row];
            LLSQ(row, 3) += cellCenter[row];
            LLSQsource[row] += psiValue*cellCenter[row];
        }
        for (char col = 0; col < 3; ++col)
            LLSQ(3, col) += cellCenter[col];
        LLSQ(3, 3) += 1;
        LLSQsource[3] += psiValue;
    }

    // Handle empty (2D/1D) directions: a one-cell-thick mesh has no variation
    // along an empty direction, so the corresponding normal component is
    // undetermined and the 4x4 system is singular. Pin nc_cmpt = 0 there and
    // decouple it, reducing the fit to the active directions. mesh.geometricD()
    // is +1 for active directions and -1 for empty ones.
    const Vector<label> geomD = mesh.geometricD();
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (geomD[cmpt] < 0)
        {
            for (label k = 0; k < 4; ++k)
            {
                LLSQ(cmpt, k) = 0;
                LLSQ(k, cmpt) = 0;
            }
            LLSQ(cmpt, cmpt) = 1;
            LLSQsource[cmpt] = 0;
        }
    }

    scalarField& source = LLSQ.source();
    source = LLSQsource;
    planeCoeffs = LLSQ.solve(); // TODO: Improve this. Gauss substitution. TM.

    return planeCoeffs;
}

// ************************************************************************* //
