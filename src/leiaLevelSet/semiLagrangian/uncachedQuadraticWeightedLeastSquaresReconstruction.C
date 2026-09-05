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

// Cholesky factorisation of the SPD matrix A (row-major n x n, n <= 9, destroyed)
// returning its SMALLEST pivot sqrt(d_j); 0 if the factorisation breaks down. On a
// Jacobi-scaled (unit-diagonal) matrix every pivot is <= 1 and the smallest one
// measures how nearly the last basis directions are linear combinations of the
// earlier ones: min pivot ~ 1/sqrt(condition number).
static Foam::scalar minCholeskyPivot(Foam::scalar* A, const Foam::label n)
{
    using Foam::scalar;
    using Foam::label;
    scalar minPiv = Foam::GREAT;
    for (label j = 0; j < n; ++j)
    {
        scalar d = A[j*n + j];
        for (label k = 0; k < j; ++k) { d -= A[j*n + k]*A[j*n + k]; }
        if (d <= 0) { return 0; }
        d = Foam::sqrt(d);
        minPiv = Foam::min(minPiv, d);
        A[j*n + j] = d;
        for (label i = j + 1; i < n; ++i)
        {
            scalar s = A[i*n + j];
            for (label k = 0; k < j; ++k) { s -= A[i*n + k]*A[j*n + k]; }
            A[i*n + j] = s/d;
        }
    }
    return minPiv;
}


// Householder-QR least squares of the m x n (m >= n, n <= 9) system Q R x = y,
// with the DESIGN matrix Ad stored row-major (m rows, n cols) and rhs y (length
// m); solution written to x (length n). Solves min ||Ad x - y||_2 WITHOUT
// forming the normal equations (no condition-number squaring). Returns false
// only if a pivot is numerically zero (rank-deficient below n) so the caller
// can fall back. m <= maxNbr (~26), n <= 9 -> tiny dense factorisation.
static bool householderQRSolve
(
    Foam::scalar* Ad, Foam::scalar* y, Foam::scalar* x,
    const Foam::label m, const Foam::label n
)
{
    using Foam::scalar;
    using Foam::label;
    for (label k = 0; k < n; ++k)
    {
        // Householder vector for column k (rows k..m-1).
        scalar normx = 0;
        for (label i = k; i < m; ++i) { normx += Ad[i*n + k]*Ad[i*n + k]; }
        normx = Foam::sqrt(normx);
        if (normx <= Foam::SMALL) { return false; }        // rank-deficient
        const scalar akk = Ad[k*n + k];
        const scalar alpha = (akk >= 0) ? -normx : normx;
        // v = column - alpha e_k ; store v in-place in Ad[k..m-1][k]
        Ad[k*n + k] = akk - alpha;
        scalar vtv = Ad[k*n + k]*Ad[k*n + k];
        for (label i = k + 1; i < m; ++i) { vtv += Ad[i*n + k]*Ad[i*n + k]; }
        if (vtv <= Foam::SMALL) { Ad[k*n + k] = akk; continue; }
        // apply H = I - 2 v v^T / (v^T v) to the trailing columns of Ad and to y
        for (label j = k + 1; j < n; ++j)
        {
            scalar s = 0;
            for (label i = k; i < m; ++i) { s += Ad[i*n + k]*Ad[i*n + j]; }
            s = 2.0*s/vtv;
            for (label i = k; i < m; ++i) { Ad[i*n + j] -= s*Ad[i*n + k]; }
        }
        scalar sy = 0;
        for (label i = k; i < m; ++i) { sy += Ad[i*n + k]*y[i]; }
        sy = 2.0*sy/vtv;
        for (label i = k; i < m; ++i) { y[i] -= sy*Ad[i*n + k]; }
        Ad[k*n + k] = alpha;                                 // R diagonal
    }
    // back-substitute R x = (Q^T y)[0..n-1]
    for (label i = n - 1; i >= 0; --i)
    {
        scalar s = y[i];
        for (label j = i + 1; j < n; ++j) { s -= Ad[i*n + j]*x[j]; }
        if (Foam::mag(Ad[i*n + i]) <= Foam::SMALL) { return false; }
        x[i] = s/Ad[i*n + i];
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
    ridgeEps_(slDict_.getOrDefault<scalar>("ridgeEps", 0)),
    quadPivotTol_(slDict_.getOrDefault<scalar>("quadraticPivotTol", 0.3)),
    writeFitOrder_(slDict_.getOrDefault<bool>("writeFitOrder", false)),
    fit_(slDict_.getOrDefault<word>("fit", "normalEquations")),
    curvatureNewtonIters_(slDict_.getOrDefault<label>("curvatureNewtonIters", 3)),
    closestPointIters_(slDict_.getOrDefault<label>("closestPointNewtonIters", 10)),
    offsetCorrection_(slDict_.getOrDefault<word>("offsetCorrection", "none")),
    fpTolRel_(slDict_.getOrDefault<scalar>("footPointTolRel", 1e-6)),
    fpAlphaMax_(slDict_.getOrDefault<scalar>("footPointAlphaMax", 20)),
    fpCycles_(slDict_.getOrDefault<label>("footPointCycles", 10)),
    fpSurfIters_(slDict_.getOrDefault<label>("footPointSurfIters", 20)),
    offsetBeta_(slDict_.getOrDefault<scalar>("offsetBeta", 0.25))
{
    if
    (
        fpTolRel_ <= 0 || fpAlphaMax_ <= 0 || fpCycles_ < 1 || fpSurfIters_ < 1
    )
    {
        FatalIOErrorInFunction(slDict_)
            << "footPoint controls must be positive "
            << "(footPointTolRel " << fpTolRel_
            << ", footPointAlphaMax " << fpAlphaMax_
            << ", footPointCycles " << fpCycles_
            << ", footPointSurfIters " << fpSurfIters_ << ")"
            << exit(FatalIOError);
    }
    if (fit_ != "normalEquations" && fit_ != "householderQR")
    {
        FatalIOErrorInFunction(slDict_)
            << "fit must be normalEquations or householderQR, got '"
            << fit_ << "'" << exit(FatalIOError);
    }
    if
    (
        offsetCorrection_ != "none"
     && offsetCorrection_ != "psi"
     && offsetCorrection_ != "psiOverGradPsi"
     && offsetCorrection_ != "quadraticRoot"
    )
    {
        FatalIOErrorInFunction(slDict_)
            << "offsetCorrection must be none, psi, psiOverGradPsi or "
            << "quadraticRoot, got '" << offsetCorrection_ << "'"
            << exit(FatalIOError);
    }
    if (offsetBeta_ <= 0 || offsetBeta_ >= 1)
    {
        FatalIOErrorInFunction(slDict_)
            << "offsetBeta must lie in (0,1), got " << offsetBeta_
            << exit(FatalIOError);
    }
    Info<< "uncachedQuadraticWeightedLeastSquares: fit = " << fit_
        << ", offsetCorrection = " << offsetCorrection_ << endl;
}

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


// grad P(d) = g_lin + H d, read directly from the fitted coefficients (same
// coeff<->basis map as basis()): linear cf[0..nd-1], diagonal Hessian
// cf[nd..2nd-1], off-diagonal (cross) Hessian cf[2nd..].
void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::gradFromCoeffs
(
    const scalar* cf,
    const label nc,
    const vector& d,
    vector& g
) const
{
    g = Zero;
    const label nd = activeDirs_.size();

    for (label a = 0; a < nd; ++a) { g[activeDirs_[a]] += cf[a]; }
    if (nc == nd) { return; }             // linear fit: constant gradient

    for (label a = 0; a < nd; ++a)
    {
        g[activeDirs_[a]] += cf[nd + a]*d[activeDirs_[a]];
    }
    label k = 2*nd;
    for (label a = 0; a < nd; ++a)
    {
        for (label bcol = a + 1; bcol < nd; ++bcol)
        {
            g[activeDirs_[a]]    += cf[k]*d[activeDirs_[bcol]];
            g[activeDirs_[bcol]] += cf[k]*d[activeDirs_[a]];
            ++k;
        }
    }
}


// Symmetric Hessian H of the quadratic fit (constant over the cell). The 0.5 on
// the diagonal basis term means d^2 P/dx_a^2 = cf[nd+a]; cross terms are cf[2nd..].
void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::hessianFromCoeffs
(
    const scalar* cf,
    const label nc,
    symmTensor& H
) const
{
    H = Zero;
    const label nd = activeDirs_.size();
    if (nc < ncoeffFull_) { return; }     // linear/constant fit: no curvature

    scalar Hm[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (label a = 0; a < nd; ++a)
    {
        Hm[activeDirs_[a]][activeDirs_[a]] = cf[nd + a];
    }
    label k = 2*nd;
    for (label a = 0; a < nd; ++a)
    {
        for (label bcol = a + 1; bcol < nd; ++bcol)
        {
            const label i = activeDirs_[a];
            const label j = activeDirs_[bcol];
            Hm[i][j] = cf[k];
            Hm[j][i] = cf[k];
            ++k;
        }
    }
    H = symmTensor
    (
        Hm[0][0], Hm[0][1], Hm[0][2],
                  Hm[1][1], Hm[1][2],
                            Hm[2][2]
    );
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
    for (label c = 0; c < nCells; ++c)
    {
        const label nNbr = stencilSize(c) - 1;   // skip self
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

    // GEOMETRIC conditioning of the quadratic fit (2026-09-05). A stencil can be
    // large enough for the quadratic and still not span it: cfMesh's boundary-layer
    // and size-transition polyhedra offer 10 point-neighbours for 9 coefficients,
    // and their weighted normal matrices had condition numbers of 1e7-1e12 with
    // curvature coefficients of 1e3-1e8 that the absolute pivot test in
    // choleskySolve (d <= SMALL on entries of order h^4) never sees. Such a fit is
    // exact at d = 0 and garbage at any displacement -- invisible on the
    // stationary droplet (|u| ~ 1e-5 m/s, foot displacement 1e-10) and on the
    // kinematic gates (u = 0 on the walls, where these cells live), catastrophic on
    // a translating polyhedral case (the interface moved 12 U dt in one step, blow-up
    // at step 3, measured on popinetTranslating3D). The decision is made ONCE per
    // mesh from geometry alone: assemble the weighted normal matrix of the full
    // quadratic, Jacobi-scale it to unit diagonal, take its smallest Cholesky pivot,
    // and below quadraticPivotTol use the LINEAR fit, whose coefficients are the
    // first nd of the same buffer (every consumer honours ncoeff_[c]).
    // quadraticPivotTol 0 reproduces the former behaviour exactly.
    label nIllCond = 0;
    label nIllCondLin = 0;
    label nMarginal = 0;
    scalar minPivotAll = GREAT;
    scalarField pivotField(writeFitOrder_ ? nCells : 0, -1.0);   // quadratic pivot per cell (diagnostic)
    if (quadPivotTol_ > 0)
    {
        scalar A[81];
        scalar brow[9];
        scalar s[9];
        for (label c = 0; c < nCells; ++c)
        {
            if (ncoeff_[c] != ncoeffFull_) { continue; }
            const label nc = ncoeffFull_;
            const label nNbr = nNbr_[c];
            const point xc = stencilC(c, 0);
            for (label k = 0; k < nc*nc; ++k) { A[k] = 0; }
            for (label i = 0; i < nNbr; ++i)
            {
                const vector d = stencilC(c, i + 1) - xc;
                const scalar w2 = 1.0/Foam::max(magSqr(d), SMALL);   // (1/|d|)^2
                basis(d, nc, brow);
                for (label k = 0; k < nc; ++k)
                {
                    const scalar wb = w2*brow[k];
                    for (label l = 0; l < nc; ++l) { A[k*nc + l] += wb*brow[l]; }
                }
            }
            scalar minPiv = 0;
            bool scaled = true;
            for (label k = 0; k < nc; ++k)
            {
                const scalar akk = A[k*nc + k];
                if (!(akk > 0)) { scaled = false; break; }
                s[k] = 1.0/Foam::sqrt(akk);
            }
            if (scaled)
            {
                for (label k = 0; k < nc; ++k)
                {
                    for (label l = 0; l < nc; ++l) { A[k*nc + l] *= s[k]*s[l]; }
                }
                minPiv = minCholeskyPivot(A, nc);
            }
            minPivotAll = Foam::min(minPivotAll, minPiv);
            if (writeFitOrder_) { pivotField[c] = minPiv; }
            if (minPiv < quadPivotTol_)
            {
                ncoeff_[c] = nd;
                ++nIllCond;
            }
            else if (minPiv < 1.5*quadPivotTol_)
            {
                ++nMarginal;   // kept quadratic, but close to the line: reported
            }
        }
        // The same test on the LINEAR fit of every cell that is not quadratic (the
        // demoted ones and the under-determined ones): a stencil whose neighbours are
        // nearly coplanar leaves the gradient component normal to that plane free, and
        // a garbage gradient costs coefficient x |u dt| at FIRST power at the foot.
        // Measured 2026-09-05 on the translating polyhedral droplet: two boundary-layer
        // cells survived the quadratic demotion and were wrong by 2.4e-2 after three
        // steps. Such a cell keeps its own value (constant reconstruction).
        for (label c = 0; c < nCells; ++c)
        {
            if (ncoeff_[c] != nd) { continue; }
            const label nc = nd;
            const label nNbr = nNbr_[c];
            const point xc = stencilC(c, 0);
            for (label k = 0; k < nc*nc; ++k) { A[k] = 0; }
            for (label i = 0; i < nNbr; ++i)
            {
                const vector d = stencilC(c, i + 1) - xc;
                const scalar w2 = 1.0/Foam::max(magSqr(d), SMALL);
                basis(d, nc, brow);
                for (label k = 0; k < nc; ++k)
                {
                    const scalar wb = w2*brow[k];
                    for (label l = 0; l < nc; ++l) { A[k*nc + l] += wb*brow[l]; }
                }
            }
            scalar minPiv = 0;
            bool scaled = true;
            for (label k = 0; k < nc; ++k)
            {
                const scalar akk = A[k*nc + k];
                if (!(akk > 0)) { scaled = false; break; }
                s[k] = 1.0/Foam::sqrt(akk);
            }
            if (scaled)
            {
                for (label k = 0; k < nc; ++k)
                {
                    for (label l = 0; l < nc; ++l) { A[k*nc + l] *= s[k]*s[l]; }
                }
                minPiv = minCholeskyPivot(A, nc);
            }
            if (minPiv < quadPivotTol_)
            {
                ncoeff_[c] = 0;
                ++nIllCondLin;
            }
        }
    }
    reduce(nIllCond, sumOp<label>());
    reduce(nIllCondLin, sumOp<label>());
    reduce(nMarginal, sumOp<label>());
    reduce(minPivotAll, minOp<scalar>());
    if (quadPivotTol_ > 0)
    {
        Info<< "uncachedQuadraticWeightedLeastSquares: " << nIllCond
            << " cells with ill-conditioned quadratic stencils (scaled Cholesky pivot < "
            << quadPivotTol_ << ", smallest " << minPivotAll
            << ") use a linear fit; " << nIllCondLin
            << " of those with an ill-conditioned linear stencil too keep their cell value; "
            << nMarginal << " quadratic cells within a factor 1.5 above the tolerance" << endl;
    }
    if (writeFitOrder_)
    {
        volScalarField fitOrder
        (
            IOobject("slFitOrder", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar(dimless, 0)
        );
        forAll(ncoeff_, c)
        {
            fitOrder[c] = (ncoeff_[c] == ncoeffFull_) ? 2 : (ncoeff_[c] == nd ? 1 : 0);
        }
        fitOrder.write();
        volScalarField fitPivot
        (
            IOobject("slFitPivot", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar(dimless, -1)
        );
        forAll(pivotField, c) { fitPivot[c] = pivotField[c]; }
        fitPivot.write();
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

    const bool useQR = (fit_ == "householderQR");

    scalar A[81];      // ncoeff x ncoeff normal-equations matrix (<= 9x9)
    scalar g[9];       // rhs / solution
    scalar brow[9];    // one basis row
    // QR scratch: weighted design rows (m x nc, row-major) and weighted rhs.
    List<scalar> Ad(useQR ? maxNbr_*ncoeffFull_ : 0);
    List<scalar> yq(useQR ? maxNbr_ : 0);

    label nSingular = 0;
    forAll(stencilPsi_, c)
    {
        const label nc = ncoeff_[c];
        if (nc == 0) { continue; }    // constant reconstruction (= psi_c)

        const List<scalar>& s = stencilPsi_[c];
        const label nNbr = nNbr_[c];
        const point xc = stencilC(c, 0);         // arrival cell centre (== mesh_.C()[c])
        const scalar psiC = s[0];
        scalar* cf = &coeffsFlat_[c*ncoeffFull_];

        if (useQR)
        {
            // Assemble the WEIGHTED DESIGN system W A x = W ds directly and solve
            // by Householder QR -- no normal-equations squaring. Same 1/|d| IDW
            // weight, so the minimiser is identical to the Cholesky path on a
            // well-conditioned cell but far better conditioned near-degenerate.
            for (label i = 0; i < nNbr; ++i)
            {
                const vector d = stencilC(c, i + 1) - xc;
                const scalar w = 1.0/Foam::max(Foam::mag(d), SMALL);
                basis(d, nc, brow);
                for (label k = 0; k < nc; ++k) { Ad[i*nc + k] = w*brow[k]; }
                yq[i] = w*(s[i + 1] - psiC);
            }
            const bool ok =
                (nNbr >= nc)
             && householderQRSolve(Ad.data(), yq.data(), g, nNbr, nc);
            if (!ok) { ++nSingular; }
            for (label k = 0; k < nc; ++k) { cf[k] = ok ? g[k] : 0; }
            continue;
        }

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
            const vector d = stencilC(c, i + 1) - xc;
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
        for (label k = 0; k < nc; ++k) { cf[k] = ok ? g[k] : 0; }
    }

    reduce(nSingular, sumOp<label>());
    if (nSingular > 0)
    {
        WarningInFunction
            << "uncachedQuadraticWeightedLeastSquares: least-squares fit failed in "
            << nSingular << " cells (near-degenerate stencil); fell back to the"
            << " constant reconstruction there."
            << (useQR ? "" : " Consider a small ridgeEps or fit householderQR.")
            << endl;
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


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::meanCurvature
(
    volScalarField& kappa
) const
{
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    scalar bb[9];
    forAll(k, c)
    {
        const label nc = ncoeff_[c];
        if (nc < ncoeffFull_) { continue; }   // need the full quadratic (Hessian)

        const scalar psiC = (*psiOldPtr_)[c];
        const scalar* cf = &coeffsFlat_[c*ncoeffFull_];

        // Restrict to the near-interface band: kappa in far cells multiplies a
        // zero snGrad(alpha) in the CSF force, so it is never used.
        // ISO-AGNOSTIC gate: the distance to the interface is estimated as
        // |psi|/|grad psi| (first order, epsilon-guarded) -- the raw |psi| value is
        // a distance ONLY for a signed-distance field; on an algebraic level set
        // (|grad psi| ~ 1e3) the raw test fails in EVERY cell and silently zeroes
        // kappa (measured on the equal-axes-ellipsoid droplet).
        {
            vector dG(Zero), gG(Zero);
            gradFromCoeffs(cf, nc, dG, gG);
            if (Foam::mag(psiC)/max(Foam::mag(gG), SMALL) > 3.0*stencilRadius(c))
            {
                continue;
            }
        }

        symmTensor H;
        hessianFromCoeffs(cf, nc, H);

        // Foot-point Newton from the cell centre onto the zero set {P = 0}. Cells
        // sharing an interface-normal ray converge to the same foot, so evaluating
        // kappa THERE makes it constant along the normal (normal-constant
        // extension) -- purely symbolic and local, no field sweep.
        vector delta(Zero);
        for (label it = 0; it < curvatureNewtonIters_; ++it)
        {
            basis(delta, nc, bb);
            scalar p = psiC;
            for (label m = 0; m < nc; ++m) { p += cf[m]*bb[m]; }
            vector gr(Zero);
            gradFromCoeffs(cf, nc, delta, gr);
            const scalar g2 = gr & gr;
            if (g2 < SMALL) { break; }
            delta -= (p/g2)*gr;
        }

        // kappa = div(grad psi/|grad psi|) = (tr(H)|g|^2 - g.(H.g))/|g|^3, evaluated
        // at the foot point (same sign/definition as traceGradGradPsiSnGradAlpha).
        vector gf(Zero);
        gradFromCoeffs(cf, nc, delta, gf);
        const scalar gm = Foam::mag(gf);
        if (gm < SMALL) { continue; }

        const scalar kd = (tr(H)*gm*gm - (gf & (H & gf)))/(gm*gm*gm);
        k[c] =
            (offsetCorrection_ == "none")
          ? kd
          : offsetCorrected(kd, curvatureOffset(psiC, gf, gm, H));
    }

    kappa.correctBoundaryConditions();
}


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::meanCurvatureLaplacian
(
    volScalarField& kappa
) const
{
    // Curvature under the signed-distance assumption: kappa = Laplacian(psi) =
    // tr(H). The Hessian is constant over a cell for a quadratic fit, so no
    // foot-point projection is needed. Equal to the full div(grad psi/|grad psi|)
    // of meanCurvature() only when |grad psi| = 1 (and the second derivative along
    // the normal vanishes) -- i.e. for a true signed distance.
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        const label nc = ncoeff_[c];
        if (nc < ncoeffFull_) { continue; }   // need the full quadratic (Hessian)

        const scalar psiC = (*psiOldPtr_)[c];
        const scalar* cf = &coeffsFlat_[c*ncoeffFull_];
        // Iso-agnostic band gate (distance ~ |psi|/|grad psi|; see meanCurvature).
        {
            vector dG(Zero), gG(Zero);
            gradFromCoeffs(cf, nc, dG, gG);
            if (Foam::mag(psiC)/max(Foam::mag(gG), SMALL) > 3.0*stencilRadius(c))
            {
                continue;
            }
        }
        symmTensor H;
        hessianFromCoeffs(cf, nc, H);
        vector d0(Zero), g0(Zero);
        gradFromCoeffs(cf, nc, d0, g0);
        const scalar gm0 = max(Foam::mag(g0), SMALL);
        k[c] =
            (offsetCorrection_ == "none")
          ? tr(H)
          : offsetCorrected(tr(H), curvatureOffset(psiC, g0, gm0, H));
    }

    kappa.correctBoundaryConditions();
}


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::meanCurvatureNoExtension
(
    volScalarField& kappa
) const
{
    // Identical closed form to meanCurvature -- kappa = (tr(H)|g|^2 - g.(H.g))/|g|^3
    // -- but the gradient g is taken AT THE CELL CENTRE (delta = 0), skipping the
    // foot-point Newton projection. The curvature is therefore NOT extended constant
    // along the interface normal; this isolates the effect of that extension.
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    forAll(k, c)
    {
        const label nc = ncoeff_[c];
        if (nc < ncoeffFull_) { continue; }

        const scalar psiC = (*psiOldPtr_)[c];
        const scalar* cf = &coeffsFlat_[c*ncoeffFull_];

        symmTensor H;
        hessianFromCoeffs(cf, nc, H);

        // Gradient at the cell centre: g = grad P(0) = the linear coefficients.
        vector d0(Zero);
        vector gc(Zero);
        gradFromCoeffs(cf, nc, d0, gc);
        const scalar gm = Foam::mag(gc);
        if (gm < SMALL) { continue; }

        // Iso-agnostic band gate (distance ~ |psi|/|grad psi|; see meanCurvature).
        if (Foam::mag(psiC)/gm > 3.0*stencilRadius(c)) { continue; }

        const scalar kd = (tr(H)*gm*gm - (gc & (H & gc)))/(gm*gm*gm);
        k[c] =
            (offsetCorrection_ == "none")
          ? kd
          : offsetCorrected(kd, curvatureOffset(psiC, gc, gm, H));
    }

    kappa.correctBoundaryConditions();
}


void Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::
meanCurvatureClosestPoint
(
    volScalarField& kappa
) const
{
    // Symbolic kappa at the interface foot of the cell's OWN quadratic, the
    // foot found as the TRUE closest point: min |delta|^2 s.t. P(delta) = 0.
    // KKT: delta = lambda*gradP(delta), P(delta) = 0, solved by Newton on
    // F(delta, lambda) with the analytic Jacobian
    //     [ I - lambda*H   -gradP(delta) ]
    //     [ gradP(delta)^T       0       ]
    // (H = the constant fit Hessian). This differs from meanCurvature()'s
    // gradient-projection flow (delta -= P/|gradP|^2 gradP), which follows
    // curving gradient lines and lands NEAR but not AT the closest point.
    // Guards: trust region |delta| <= stencilRadius(c) (a capped foot means
    // the zero set lies outside the trusted fit region -- outer band cells),
    // pivot/degeneracy checks on the 4x4 solve, non-finite checks; every
    // guard failure falls back to the cell-centre (NoExtension) value, so
    // only the force-support layer is foot-evaluated.
    scalarField& k = kappa.primitiveFieldRef();
    k = 0;

    scalar bb[9];
    forAll(k, c)
    {
        const label nc = ncoeff_[c];
        if (nc < ncoeffFull_) { continue; }

        const scalar psiC = (*psiOldPtr_)[c];
        const scalar* cf = &coeffsFlat_[c*ncoeffFull_];

        symmTensor H;
        hessianFromCoeffs(cf, nc, H);
        const tensor Ht(H);

        // Cell-centre gradient; the fallback (NoExtension) value.
        vector d0(Zero);
        vector g0(Zero);
        gradFromCoeffs(cf, nc, d0, g0);
        const scalar gm0 = Foam::mag(g0);
        if (gm0 < SMALL) { continue; }

        // Iso-agnostic band gate (distance ~ |psi|/|grad psi|; see meanCurvature).
        if (Foam::mag(psiC)/gm0 > 3.0*stencilRadius(c)) { continue; }

        const scalar kCentre =
            (tr(H)*gm0*gm0 - (g0 & (H & g0)))/(gm0*gm0*gm0);
        k[c] = kCentre;

        const scalar rTrust = stencilRadius(c);
        const scalar tolP = 1e-8*rTrust;

        // Initial guess: one linearized projection step.
        scalar lambda = -psiC/(gm0*gm0);
        vector delta = lambda*g0;

        bool ok = false;
        bool failed = false;
        for (label it = 0; it < closestPointIters_; ++it)
        {
            basis(delta, nc, bb);
            scalar p = psiC;
            for (label m = 0; m < nc; ++m) { p += cf[m]*bb[m]; }

            vector gr(Zero);
            gradFromCoeffs(cf, nc, delta, gr);
            if (magSqr(gr) < SMALL) { failed = true; break; }

            const vector rV = delta - lambda*gr;
            if (Foam::mag(p) < tolP && Foam::mag(rV) < 1e-8*rTrust)
            {
                ok = true;
                break;
            }

            // Newton system A*[ddelta; dlambda] = -[rV; p], 4x4 Gaussian
            // elimination with partial pivoting and a pivot guard (empty
            // mesh directions give clean identity rows).
            scalar A[4][5];
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    A[i][j] = (i == j ? 1.0 : 0.0) - lambda*Ht(i, j);
                }
                A[i][3] = -gr[i];
                A[i][4] = -rV[i];
                A[3][i] = gr[i];
            }
            A[3][3] = 0;
            A[3][4] = -p;

            for (label col = 0; col < 4 && !failed; ++col)
            {
                label piv = col;
                for (label r = col + 1; r < 4; ++r)
                {
                    if (Foam::mag(A[r][col]) > Foam::mag(A[piv][col]))
                    {
                        piv = r;
                    }
                }
                if (Foam::mag(A[piv][col]) < VSMALL) { failed = true; break; }
                if (piv != col)
                {
                    for (label j = col; j < 5; ++j)
                    {
                        const scalar t = A[col][j];
                        A[col][j] = A[piv][j];
                        A[piv][j] = t;
                    }
                }
                for (label r = col + 1; r < 4; ++r)
                {
                    const scalar f = A[r][col]/A[col][col];
                    for (label j = col; j < 5; ++j) { A[r][j] -= f*A[col][j]; }
                }
            }
            if (failed) { break; }

            scalar sol[4];
            for (label i = 3; i >= 0; --i)
            {
                scalar s = A[i][4];
                for (label j = i + 1; j < 4; ++j) { s -= A[i][j]*sol[j]; }
                sol[i] = s/A[i][i];
            }
            if
            (
                !std::isfinite(sol[0]) || !std::isfinite(sol[1])
             || !std::isfinite(sol[2]) || !std::isfinite(sol[3])
            )
            {
                failed = true;
                break;
            }

            delta += vector(sol[0], sol[1], sol[2]);
            lambda += sol[3];

            // Trust region: a capped step means the foot is outside the
            // trusted fit region -- reject and keep the centre value.
            if (Foam::mag(delta) > rTrust) { failed = true; break; }
        }

        if (failed || !ok) { continue; }   // k[c] stays at kCentre

        vector gf(Zero);
        gradFromCoeffs(cf, nc, delta, gf);
        const scalar gm = Foam::mag(gf);
        if (gm < SMALL) { continue; }

        const scalar kFoot = (tr(H)*gm*gm - (gf & (H & gf)))/(gm*gm*gm);
        if (std::isfinite(kFoot)) { k[c] = kFoot; }
    }

    kappa.correctBoundaryConditions();
}


Foam::label
Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::fitDerivatives
(
    const label c,
    vector& g,
    symmTensor& H
) const
{
    g = Zero;
    H = Zero;

    const label nc = ncoeff_[c];
    if (nc == 0) { return 0; }                        // constant fallback

    const scalar* cf = &coeffsFlat_[c*ncoeffFull_];

    // Gradient at the cell centre = the linear coefficients (d = 0).
    const vector d0(Zero);
    gradFromCoeffs(cf, nc, d0, g);

    if (nc < ncoeffFull_) { return 1; }               // linear fallback: H = 0

    hessianFromCoeffs(cf, nc, H);
    return 2;
}


Foam::scalar
Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::signedOffset
(
    const label c,
    bool& fallback
) const
{
    vector g(Zero);
    symmTensor H(Zero);
    const label order = fitDerivatives(c, g, H);

    const scalar psiC = (*psiOldPtr_)[c];
    const scalar gm = Foam::mag(g);

    if (order == 0 || gm < SMALL)
    {
        // No usable normal ray: report the fallback branch; the caller (the
        // nSL scheme freezes such cells before asking) decides what to do.
        fallback = true;
        return 0;
    }

    const vector n = g/gm;
    const scalar hnn = n & (H & n);                   // 0 for a linear fit

    return offsetDistance(psiC, gm, hnn, fallback);
}


Foam::scalar
Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::footPointDistance
(
    const label c,
    const point& p,
    const scalar level,
    bool& ok
) const
{
    // Stabilized foot-point on the cell's own quadratic model
    // (stable-foot-point-3d.md). The shifted residual is
    // G(d) = R_c(x_c + d) - level with grad G(d) = g + H d, both read from the
    // fit coefficients -- pure polynomial arithmetic, no field access. All
    // positions are displacements d from the cell centre.
    ok = false;

    const label nc = ncoeff_[c];
    if (nc == 0) { return 0; }

    const scalar* cf = &coeffsFlat_[c*ncoeffFull_];
    const scalar psiC = (*psiOldPtr_)[c];
    const point& xc = mesh_.C()[c];
    const scalar radius = stencilRadius(c);
    if (radius < SMALL) { return 0; }

    const scalar eps = fpTolRel_*radius;
    const vector dp = p - xc;                          // query displacement

    // G(d) = psiC + sum_k cf_k b_k(d) - level.
    auto G = [&](const vector& d) -> scalar
    {
        scalar b[9];
        basis(d, nc, b);
        scalar v = psiC - level;
        for (label k = 0; k < nc; ++k) { v += cf[k]*b[k]; }
        return v;
    };

    // surfacepoint: Newton along the local gradient onto {G = 0} (note Sec. 2).
    // Returns false on a degenerate gradient or non-convergence.
    auto surfacepoint = [&](vector& d) -> bool
    {
        for (label k = 0; k < fpSurfIters_; ++k)
        {
            vector gr(Zero);
            gradFromCoeffs(cf, nc, d, gr);
            const scalar g2 = gr & gr;
            if (g2 < SMALL) { return false; }

            const vector step = -(G(d)/g2)*gr;
            d += step;
            if (Foam::mag(step) < eps) { return true; }
        }
        return false;
    };

    // Initial on-surface point: surfacepoint of the query itself (Sec. 4.3).
    vector dpi = dp;
    if (!surfacepoint(dpi)) { return 0; }

    bool converged = false;
    for (label cycle = 0; cycle < fpCycles_ && !converged; ++cycle)
    {
        // (1) foot point of the query on the tangent plane at p_i.
        vector gi(Zero);
        gradFromCoeffs(cf, nc, dpi, gi);
        const scalar gi2 = gi & gi;
        if (gi2 < SMALL) { return 0; }

        const scalar lam = ((dp - dpi) & gi)/gi2;
        vector dqi = dp - lam*gi;

        // (2) back to the surface.
        vector dpn = dqi;
        if (!surfacepoint(dpn)) { return 0; }

        // (3) parabola correction from the two step vectors (Sec. 4).
        const vector f1 = dqi - dpi;
        const vector f2 = dpn - dqi;

        if (Foam::mag(f1) > eps)
        {
            const vector rel = dp - dpi;               // p - p_i
            const scalar a0 = rel & f1;
            const scalar a1 = 2.0*(f2 & rel) - magSqr(f1);
            const scalar a2 = -3.0*(f1 & f2);
            const scalar a3 = -2.0*magSqr(f2);
            const scalar den = a1 + 2.0*a2 + 3.0*a3;

            if (Foam::mag(den) > SMALL)
            {
                const scalar alpha = 1.0 - (a0 + a1 + a2 + a3)/den;

                if (alpha > 0 && alpha < fpAlphaMax_)
                {
                    dqi = dpi + alpha*f1 + sqr(alpha)*f2;
                    dpn = dqi;
                    if (!surfacepoint(dpn)) { return 0; }
                }
                // else: keep the plain-cycle dpn (degenerate guard, Sec. 4.3)
            }
        }

        converged = (Foam::mag(dpn - dpi) < eps);
        dpi = dpn;
    }

    // Trust region: the model is only meaningful where the stencil sampled it.
    if (!converged || Foam::mag(dpi) > radius) { return 0; }

    vector gS(Zero);
    gradFromCoeffs(cf, nc, dpi, gS);
    const scalar gSm = Foam::mag(gS);
    if (gSm < SMALL) { return 0; }

    // Signed distance, sign along the model gradient (positive on the +psi
    // side), consistent with signedOffset's convention.
    ok = true;
    return ((dp - dpi) & gS)/gSm;
}


Foam::scalar
Foam::uncachedQuadraticWeightedLeastSquaresReconstruction::curvatureAtFootPoint
(
    const label c,
    const point& p,
    bool& ok
) const
{
    // Fit curvature AT the closest point of p on the cell's own quadratic
    // zero set. The foot search below DUPLICATES footPointDistance on purpose:
    // that function sits inside every measured production path (the parallel-
    // surface delivery, nSL Variant (a)) and must stay byte-untouched; a
    // shared-helper refactor would re-generate its code and forfeit the
    // bit-identity gates. Keep the two in sync when the algorithm changes.
    ok = false;

    const label nc = ncoeff_[c];
    if (nc == 0) { return 0; }

    const scalar* cf = &coeffsFlat_[c*ncoeffFull_];
    const scalar psiC = (*psiOldPtr_)[c];
    const point& xc = mesh_.C()[c];
    const scalar radius = stencilRadius(c);
    if (radius < SMALL) { return 0; }

    const scalar eps = fpTolRel_*radius;
    const vector dp = p - xc;

    auto G = [&](const vector& d) -> scalar
    {
        scalar b[9];
        basis(d, nc, b);
        scalar v = psiC;                                // level = 0: zero set
        for (label k = 0; k < nc; ++k) { v += cf[k]*b[k]; }
        return v;
    };

    auto surfacepoint = [&](vector& d) -> bool
    {
        for (label k = 0; k < fpSurfIters_; ++k)
        {
            vector gr(Zero);
            gradFromCoeffs(cf, nc, d, gr);
            const scalar g2 = gr & gr;
            if (g2 < SMALL) { return false; }

            const vector step = -(G(d)/g2)*gr;
            d += step;
            if (Foam::mag(step) < eps) { return true; }
        }
        return false;
    };

    vector dpi = dp;
    if (!surfacepoint(dpi)) { return 0; }

    bool converged = false;
    for (label cycle = 0; cycle < fpCycles_ && !converged; ++cycle)
    {
        vector gi(Zero);
        gradFromCoeffs(cf, nc, dpi, gi);
        const scalar gi2 = gi & gi;
        if (gi2 < SMALL) { return 0; }

        const scalar lam = ((dp - dpi) & gi)/gi2;
        vector dqi = dp - lam*gi;

        vector dpn = dqi;
        if (!surfacepoint(dpn)) { return 0; }

        const vector f1 = dqi - dpi;
        const vector f2 = dpn - dqi;

        if (Foam::mag(f1) > eps)
        {
            const vector rel = dp - dpi;
            const scalar a0 = rel & f1;
            const scalar a1 = 2.0*(f2 & rel) - magSqr(f1);
            const scalar a2 = -3.0*(f1 & f2);
            const scalar a3 = -2.0*magSqr(f2);
            const scalar den = a1 + 2.0*a2 + 3.0*a3;

            if (Foam::mag(den) > SMALL)
            {
                const scalar alpha = 1.0 - (a0 + a1 + a2 + a3)/den;

                if (alpha > 0 && alpha < fpAlphaMax_)
                {
                    dqi = dpi + alpha*f1 + sqr(alpha)*f2;
                    dpn = dqi;
                    if (!surfacepoint(dpn)) { return 0; }
                }
            }
        }

        converged = (Foam::mag(dpn - dpi) < eps);
        dpi = dpn;
    }

    if (!converged || Foam::mag(dpi) > radius) { return 0; }

    // kappa = div(grad R/|grad R|) of the fit AT the foot: gradient there from
    // the coefficients, the (constant) quadratic Hessian likewise.
    vector g(Zero);
    gradFromCoeffs(cf, nc, dpi, g);
    const scalar gm = Foam::mag(g);
    if (gm < SMALL) { return 0; }

    symmTensor H(Zero);
    hessianFromCoeffs(cf, nc, H);

    ok = true;
    return (tr(H)*gm*gm - (g & (H & g)))/(gm*gm*gm);
}

// ************************************************************************* //
