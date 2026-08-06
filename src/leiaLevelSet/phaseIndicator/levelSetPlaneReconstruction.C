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

#include <cmath>

// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- 3x3 determinant by cofactor expansion (no divisions).
inline scalar det3x3
(
    scalar a, scalar b, scalar c,
    scalar d, scalar e, scalar f,
    scalar g, scalar h, scalar i
)
{
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}

} // End anonymous namespace
} // End namespace Foam

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

    // Guard against a singular normal-equations matrix. On general polyhedra a
    // narrow-band cell can have a face-neighbour cloud whose centres are (near)
    // coplanar, making the 4x4 LLSQ system rank-deficient; simpleMatrix::solve()
    // then divides by a ~0 pivot and returns non-finite coefficients, which blow
    // up downstream (nc/|nc| = inf/inf = NaN -> SIGFPE). Detect this with the
    // Hadamard ratio |det|/prod(diag) (~1 for well-spread rows, -> 0 for a
    // coplanar/degenerate cloud) and return a flat (zero) plane; the callers
    // treat |n| ~ 0 as "degenerate" and fall back to the sign of psi.
    //
    // The test matrix is assembled in CENTERED coordinates (displacements from
    // the cell centre), NOT from LLSQ itself: with the absolute coordinates of
    // the solve matrix the ratio shrinks like (h/|x|)^4 as the mesh refines --
    // the columns (x, y, 1) become collinear -- and a fixed threshold misfires
    // on perfectly regular stencils (measured: every N=512 band cell of the
    // 0.01 m stationary-droplet case returned the flat-plane fallback, killing
    // the plane-based phase indicators and curvature models at fine h).
    // Centered displacements make the ratio a pure shape measure of the
    // neighbour cloud: scale-free, ~O(1e-1) on hexahedral stencils at any N,
    // -> 0 exactly for the coplanar clouds the guard exists to catch. The
    // SOLVE still uses the absolute system, so every cell that passed before
    // keeps byte-identical coefficients.
    scalar Mc[4][4] = {{0}};
    {
        const vector& xc = cellCenters[cellI];

        auto addSample = [&](const vector& x)
        {
            const vector dx = x - xc;
            for (char row = 0; row < 3; ++row)
            {
                for (char col = 0; col < 3; ++col)
                {
                    Mc[row][col] += dx[col]*dx[row];
                }
                Mc[row][3] += dx[row];
                Mc[3][row] += dx[row];
            }
            Mc[3][3] += 1;
        };

        addSample(cellCenters[cellI]);
        forAll(Nc, cellK)
        {
            addSample(cellCenters[Nc[cellK]]);
        }
        forAll(nBandCell, faceI)
        {
            const label faceJ = nBandCell[faceI];
            if (mesh.isInternalFace(faceJ))
            {
                continue;
            }
            const label patchL = pbm.whichPatch(faceJ);
            if ((patchL < 0) || !coupledNei.valid(patchL))
            {
                continue;
            }
            addSample(coupledNei.C(patchL)[faceJ - pbm[patchL].start()]);
        }

        // Same empty-direction pinning as the solve matrix.
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            if (geomD[cmpt] < 0)
            {
                for (label k = 0; k < 4; ++k)
                {
                    Mc[cmpt][k] = 0;
                    Mc[k][cmpt] = 0;
                }
                Mc[cmpt][cmpt] = 1;
            }
        }

        // Scale away the mixed length dimensions so the ratio compares like
        // with like: row/column 0..2 carry one length factor each relative to
        // row/column 3. Normalising by the mean squared displacement makes
        // every entry O(1) on a healthy stencil.
        scalar h2 = 0;
        label nS = 0;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            if (geomD[cmpt] > 0) { h2 += Mc[cmpt][cmpt]; ++nS; }
        }
        h2 = (nS > 0 && Mc[3][3] > 0) ? h2/(nS*Mc[3][3]) : 0;
        if (h2 <= VSMALL)
        {
            return scalarList(4, scalar(0));   // all samples coincide
        }
        const scalar hs = Foam::sqrt(h2);
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            if (geomD[cmpt] < 0) { continue; }
            for (label k = 0; k < 4; ++k)
            {
                Mc[cmpt][k] /= hs;
                Mc[k][cmpt] /= hs;
            }
        }
    }

    const scalar diagProd =
        Mc[0][0]*Mc[1][1]*Mc[2][2]*Mc[3][3];
    const scalar detLLSQ =
        Mc[0][0]*det3x3(Mc[1][1],Mc[1][2],Mc[1][3],
                        Mc[2][1],Mc[2][2],Mc[2][3],
                        Mc[3][1],Mc[3][2],Mc[3][3])
      - Mc[0][1]*det3x3(Mc[1][0],Mc[1][2],Mc[1][3],
                        Mc[2][0],Mc[2][2],Mc[2][3],
                        Mc[3][0],Mc[3][2],Mc[3][3])
      + Mc[0][2]*det3x3(Mc[1][0],Mc[1][1],Mc[1][3],
                        Mc[2][0],Mc[2][1],Mc[2][3],
                        Mc[3][0],Mc[3][1],Mc[3][3])
      - Mc[0][3]*det3x3(Mc[1][0],Mc[1][1],Mc[1][2],
                        Mc[2][0],Mc[2][1],Mc[2][2],
                        Mc[3][0],Mc[3][1],Mc[3][2]);
    if (diagProd <= VSMALL || mag(detLLSQ) < 1e-10*mag(diagProd))
    {
        return scalarList(4, scalar(0));   // degenerate -> flat plane
    }

    // KNOWN LIMIT: the shape test certifies the stencil cloud, not the
    // conditioning of the ABSOLUTE-coordinate solve below, which degrades
    // like (|x|/h)^2. Fine for domains near the origin (relative error
    // ~1e-6 at |x|/h = 1e5); a mesh sitting ~1e7 cell sizes away from the
    // coordinate origin would get finite-but-inaccurate normals that the
    // isfinite backstop cannot catch. The clean fix -- assemble and solve in
    // centered coordinates and shift d back -- changes every result at
    // rounding level, so it must land together with a revalidation of the
    // phase-indicator gates, not as a drive-by.
    planeCoeffs = LLSQ.solve(); // TODO: Improve this. Gauss substitution. TM.

    // Belt-and-suspenders: a near-singular (but above-threshold) solve can still
    // return non-finite coefficients; treat those as degenerate too.
    for (const scalar c : planeCoeffs)
    {
        if (!std::isfinite(c))
        {
            return scalarList(4, scalar(0));
        }
    }

    return planeCoeffs;
}

// ************************************************************************* //
