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

#include "bandQuadraticWLSQReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "QRMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bandQuadraticWLSQReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        bandQuadraticWLSQReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bandQuadraticWLSQReconstruction::bandQuadraticWLSQReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    activeDirs_(),
    ncoeffFull_(0),
    cellSize_(),
    band_(),
    coeffs_(),
    ncoeff_(),
    nLayersBand_(slDict_.getOrDefault<label>("nLayersBand", 3)),
    bandGuard_(slDict_.getOrDefault<label>("bandGuard", 2)),
    initialized_(false)
{}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::bandQuadraticWLSQReconstruction::computeCellSize()
{
    cellSize_.setSize(mesh_.nCells(), GREAT);
    const surfaceScalarField& dc = mesh_.deltaCoeffs();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, f)
    {
        const scalar d = 1.0/dc[f];
        cellSize_[own[f]] = Foam::min(cellSize_[own[f]], d);
        cellSize_[nei[f]] = Foam::min(cellSize_[nei[f]], d);
    }
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const labelUList& fc = mesh_.boundary()[patchI].faceCells();
            const fvsPatchScalarField& pdc = dc.boundaryField()[patchI];
            forAll(fc, i)
            {
                cellSize_[fc[i]] = Foam::min(cellSize_[fc[i]], 1.0/pdc[i]);
            }
        }
    }
}


void Foam::bandQuadraticWLSQReconstruction::setupBasis()
{
    const Vector<label> gd = mesh_.geometricD();
    DynamicList<label> dirs(3);
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (gd[cmpt] == 1)
        {
            dirs.append(cmpt);
        }
    }
    activeDirs_ = dirs;
    const label nd = activeDirs_.size();
    // nd linear + nd pure-quadratic + nd*(nd-1)/2 cross.
    ncoeffFull_ = 2*nd + (nd*(nd - 1))/2;
}


void Foam::bandQuadraticWLSQReconstruction::basis
(
    const vector& d,
    const label ncoeff,
    scalarList& b
) const
{
    b.setSize(ncoeff);
    const label nd = activeDirs_.size();
    forAll(activeDirs_, a)
    {
        b[a] = d[activeDirs_[a]];
    }
    if (ncoeff == nd)
    {
        return;   // linear fallback
    }
    label k = nd;
    forAll(activeDirs_, a)
    {
        const scalar da = d[activeDirs_[a]];
        b[k++] = 0.5*da*da;
    }
    for (label a = 0; a < nd; ++a)
    {
        for (label bcol = a + 1; bcol < nd; ++bcol)
        {
            b[k++] = d[activeDirs_[a]]*d[activeDirs_[bcol]];
        }
    }
}


Foam::label Foam::bandQuadraticWLSQReconstruction::fitCell
(
    const List<vector>& X,
    const List<scalar>& s,
    scalarList& coeffs
) const
{
    const label nd = activeDirs_.size();
    const label nNbr = X.size() - 1;
    const point xc = X[0];
    const scalar psiC = s[0];

    label ncoeff = ncoeffFull_;
    if (nNbr < ncoeffFull_)
    {
        ncoeff = (nNbr >= nd) ? nd : 0;
    }
    coeffs.setSize(ncoeff);
    coeffs = 0.0;
    if (ncoeff == 0)
    {
        return 0;
    }

    scalarRectangularMatrix A(nNbr, ncoeff, Zero);
    scalarList bw(nNbr);
    scalarList b;
    for (label i = 0; i < nNbr; ++i)
    {
        const vector d = X[i + 1] - xc;
        const scalar wi = 1.0/Foam::max(Foam::mag(d), SMALL);
        basis(d, ncoeff, b);
        for (label k = 0; k < ncoeff; ++k)
        {
            A(i, k) = wi*b[k];
        }
        bw[i] = wi*(s[i + 1] - psiC);
    }
    const scalarRectangularMatrix Minv(MatrixTools::pinv(A));
    for (label k = 0; k < ncoeff; ++k)
    {
        scalar acc = 0;
        for (label i = 0; i < nNbr; ++i)
        {
            acc += Minv(k, i)*bw[i];
        }
        coeffs[k] = acc;
    }
    return ncoeff;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::bandQuadraticWLSQReconstruction::update(const volScalarField& psiOld)
{
    collectStencil(psiOld);

    if (!initialized_)
    {
        setupBasis();
        computeCellSize();
        initialized_ = true;
    }

    const label nCells = mesh_.nCells();
    band_.setSize(nCells);
    ncoeff_.setSize(nCells);
    coeffs_.setSize(nCells);

    const scalar bw = scalar(nLayersBand_ + bandGuard_);
    label nBand = 0;
    forAll(band_, c)
    {
        band_[c] = (Foam::mag(psiOld[c]) <= bw*cellSize_[c]);
        if (band_[c])
        {
            ncoeff_[c] = fitCell(stencilC_[c], stencilPsi_[c], coeffs_[c]);
            ++nBand;
        }
        else
        {
            ncoeff_[c] = 0;
            coeffs_[c].clear();
        }
    }

    reduce(nBand, sumOp<label>());
    Info<< "bandQuadraticWLSQ: quadratic fit in " << nBand
        << " band cells (far field re-extended)" << endl;

    computeLimiters();
}


Foam::scalar Foam::bandQuadraticWLSQReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];

    if (band_[c] && ncoeff_[c] > 0)
    {
        const vector d = x - mesh_.C()[c];
        scalarList b;
        basis(d, ncoeff_[c], b);
        scalar val = psiC;
        for (label k = 0; k < ncoeff_[c]; ++k)
        {
            val += coeffs_[c][k]*b[k];
        }
        return val;
    }
    return psiC;   // far-field: placeholder; overwritten by postAdvect().
}


void Foam::bandQuadraticWLSQReconstruction::postAdvect(volScalarField& psi)
{
    const scalar hMax = gMax(cellSize_);
    const label nSweeps = nLayersBand_ + bandGuard_ + 2;
    const scalar farInit = 1e3*hMax;

    scalarField& p = psi.primitiveFieldRef();
    forAll(p, c)
    {
        if (!band_[c])
        {
            const scalar sgn = (p[c] >= 0) ? 1.0 : -1.0;
            p[c] = sgn*farInit;
        }
    }
    psi.correctBoundaryConditions();

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        List<List<scalar>> sPsi;
        stencil_.collectData(psi, sPsi);

        const scalarField pIn(psi.primitiveField());
        scalarField pOut(pIn);
        forAll(band_, c)
        {
            if (band_[c])
            {
                continue;
            }
            const List<scalar>& s = sPsi[c];
            const List<vector>& X = stencilC_[c];
            const label n = Foam::min(s.size(), X.size());
            if (n < 2)
            {
                continue;
            }
            const scalar sgn = (pIn[c] >= 0) ? 1.0 : -1.0;
            scalar best = Foam::mag(pIn[c]);
            for (label i = 1; i < n; ++i)
            {
                const scalar cand =
                    Foam::mag(s[i]) + Foam::mag(X[i] - X[0]);
                best = Foam::min(best, cand);
            }
            pOut[c] = sgn*best;
        }
        psi.primitiveFieldRef() = pOut;
        psi.correctBoundaryConditions();
    }
}

// ************************************************************************* //
