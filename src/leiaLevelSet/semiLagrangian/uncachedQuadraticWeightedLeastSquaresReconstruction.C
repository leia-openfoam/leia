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

#include "uncachedQuadraticWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uncachedQuadraticWeightedLeastSquaresReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        uncachedQuadraticWeightedLeastSquaresReconstruction,
        Mesh
    );
}

namespace
{
// In-place Cholesky solve of the SPD system A x = b, A stored row-major n x n
// (n <= 9), b overwritten with x. Returns false if A is not positive-definite
// (near-singular / rank-deficient stencil) so the caller can fall back.
static bool choleskySolve(Foam::scalar* A, Foam::scalar* b, const Foam::label n)
{
    using Foam::scalar;
    using Foam::label;
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= Foam::SMALL) { return false; }          // not SPD
        d = Foam::sqrt(d);
        A[j*n + j] = d;
        for (label i = j + 1; i < n; ++i)
        {
            scalar s = A[i*n + j];
            for (label k = 0; k < j; ++k) { s -= A[i*n + k]*A[j*n + k]; }
            A[i*n + j] = s/d;
        }
    }
    // forward solve L y = b
    for (label i = 0; i < n; ++i)
    {
        scalar s = b[i];
        for (label k = 0; k < i; ++k) { s -= A[i*n + k]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    // backward solve L^T x = y
    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = b[i];
        for (label k = i + 1; k < n; ++k) { s -= A[k*n + i]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    return true;
}
} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::
uncachedQuadraticWeightedLeastSquaresReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    activeDirs_(),
    ncoeffFull_(0),
    nNbr_(),
    maxNbr_(0),
    ncoeff_(),
    coeffsFlat_(),
    built_(false),
    ridgeEps_(slDict_.getOrDefault<scalar>("ridgeEps", 0))
{}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

// Identical basis to quadraticWeightedLeastSquares (constant-free).
void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::basis
(
    const vector& d,
    const label ncoeff,
    scalar* b
) const
{
    const label nd = activeDirs_.size();
    forAll(activeDirs_, a) { b[a] = d[activeDirs_[a]]; }
    if (ncoeff == nd) { return; }        // linear fallback basis

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


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::setupBasis()
{
    const Vector<label> gd = mesh_.geometricD();
    DynamicList<label> dirs(3);
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (gd[cmpt] == 1) { dirs.append(cmpt); }
    }
    activeDirs_ = dirs;
    const label nd = activeDirs_.size();
    ncoeffFull_ = 2*nd + (nd*(nd - 1))/2;
}


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::build()
{
    setupBasis();
    const label nd = activeDirs_.size();

    const label nCells = mesh_.nCells();
    nNbr_.setSize(nCells);
    ncoeff_.setSize(nCells);
    maxNbr_ = 0;

    label nFallback = 0;
    forAll(stencilC_, c)
    {
        const label nNbr = stencilC_[c].size() - 1;   // skip self
        nNbr_[c] = nNbr;
        maxNbr_ = Foam::max(maxNbr_, nNbr);

        label ncoeff = ncoeffFull_;
        if (nNbr < ncoeffFull_)
        {
            ncoeff = (nNbr >= nd) ? nd : 0;   // linear fallback, else constant
            ++nFallback;
        }
        ncoeff_[c] = ncoeff;
    }

    reduce(nFallback, sumOp<label>());
    if (nFallback > 0)
    {
        Info<< "uncachedQuadraticWeightedLeastSquares: " << nFallback
            << " cells with under-determined stencils fell back to a linear fit"
            << endl;
    }
    built_ = true;
    // NB: do NOT releaseStencilCentres() -- the fit is reassembled every step.
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::update
(
    const volScalarField& psiOld
)
{
    collectStencil(psiOld);       // fills stencilPsi_ (+ stencilC_ once); sets ptr

    if (!built_) { build(); }

    coeffsFlat_.setSize(mesh_.nCells()*ncoeffFull_);   // one contiguous buffer

    scalar A[81];      // ncoeff x ncoeff normal-equations matrix (<= 9x9)
    scalar g[9];       // rhs / solution
    scalar brow[9];    // one basis row

    label nSingular = 0;
    forAll(stencilPsi_, c)
    {
        const label nc = ncoeff_[c];
        if (nc == 0) { continue; }    // constant reconstruction (= psi_c)

        const List<vector>& C = stencilC_[c];    // [0] = arrival cell c
        const List<scalar>& s = stencilPsi_[c];
        const label nNbr = nNbr_[c];
        const point xc = C[0];
        const scalar psiC = s[0];

        // Assemble the SYMMETRIC weighted normal equations M = A^T W^2 A,
        // g = A^T W^2 ds, in one sweep over the stencil. Same 1/|d| IDW weight
        // as quadraticWeightedLeastSquares -> the same weighted LS minimiser.
        for (label k = 0; k < nc; ++k)
        {
            g[k] = 0;
            for (label l = 0; l < nc; ++l) { A[k*nc + l] = 0; }
        }
        for (label i = 0; i < nNbr; ++i)
        {
            const vector d = C[i + 1] - xc;
            const scalar w = 1.0/Foam::max(Foam::mag(d), SMALL);
            const scalar w2 = w*w;
            basis(d, nc, brow);
            const scalar ds = s[i + 1] - psiC;
            for (label k = 0; k < nc; ++k)
            {
                const scalar wb = w2*brow[k];
                g[k] += wb*ds;
                for (label l = 0; l < nc; ++l) { A[k*nc + l] += wb*brow[l]; }
            }
        }

        if (ridgeEps_ > 0)
        {
            scalar tr = 0;
            for (label k = 0; k < nc; ++k) { tr += A[k*nc + k]; }
            const scalar r = ridgeEps_*tr/nc;
            for (label k = 0; k < nc; ++k) { A[k*nc + k] += r; }
        }

        const bool ok = choleskySolve(A, g, nc);
        if (!ok) { ++nSingular; }
        scalar* cf = &coeffsFlat_[c*ncoeffFull_];
        for (label k = 0; k < nc; ++k) { cf[k] = ok ? g[k] : 0; }
    }

    reduce(nSingular, sumOp<label>());
    if (nSingular > 0)
    {
        WarningInFunction
            << "uncachedQuadraticWeightedLeastSquares: normal equations not SPD in "
            << nSingular << " cells (near-degenerate stencil); fell back to the"
            << " constant reconstruction there. Consider a small ridgeEps." << endl;
    }

    computeLimiters();
}


Foam::scalar Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];
    const label nc = ncoeff_[c];
    if (nc == 0) { return psiC; }

    const vector d = x - mesh_.C()[c];
    scalar b[9];
    basis(d, nc, b);
    const scalar* cf = &coeffsFlat_[c*ncoeffFull_];
    scalar val = psiC;
    for (label k = 0; k < nc; ++k) { val += cf[k]*b[k]; }
    return val;
}

// ************************************************************************* //
