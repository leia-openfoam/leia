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

#include "signedDistanceQuadraticWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(signedDistanceQuadraticWeightedLeastSquaresReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        signedDistanceQuadraticWeightedLeastSquaresReconstruction,
        Mesh
    );
    // Short selection alias.
    addNamedToRunTimeSelectionTable
    (
        slReconstruction,
        signedDistanceQuadraticWeightedLeastSquaresReconstruction,
        Mesh,
        sdQuadraticWLSQ
    );
}

namespace
{
// In-place Cholesky solve of the SPD system A x = b, A stored row-major n x n
// (n <= 9), b overwritten with x. Returns false if A is not positive-definite.
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
    for (label i = 0; i < n; ++i)                         // forward solve L y = b
    {
        scalar s = b[i];
        for (label k = 0; k < i; ++k) { s -= A[i*n + k]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    for (label i = n - 1; i >= 0; --i)                    // backward solve L^T x = y
    {
        scalar s = b[i];
        for (label k = i + 1; k < n; ++k) { s -= A[k*n + i]*b[k]; }
        b[i] = s/A[i*n + i];
    }
    return true;
}
} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::
signedDistanceQuadraticWeightedLeastSquaresReconstruction
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
    normalStrainRescale_(slDict_.getOrDefault<Switch>("normalStrainRescale", true)),
    footpointNewtonIters_(slDict_.getOrDefault<label>("footpointNewtonIters", 3)),
    reprojectBandOnly_(slDict_.getOrDefault<Switch>("reprojectBandOnly", true)),
    reprojBandRadii_(slDict_.getOrDefault<scalar>("reprojBandRadii", 3.0)),
    gradLo_(slDict_.getOrDefault<scalar>("gradLo", 0.4)),
    gradHi_(slDict_.getOrDefault<scalar>("gradHi", 2.5)),
    gradU_(),
    dt_(0),
    haveStrain_(false)
{
    // A reprojected distance is not a stencil value: the slope limiter and the
    // stencil-range clamp (which assume evaluateRaw reproduces psi_c in range)
    // are incompatible with this model. Force them off so evaluate() passes
    // evaluateRaw() through unchanged.
    limitSlope_ = false;
    clipToStencilBounds_ = false;
}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::basis
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


void Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::gradFromCoeffs
(
    const scalar* cf,
    const label nc,
    const vector& d,
    vector& g
) const
{
    // P(d) = psi_c + sum_k cf[k] b_k(d); grad P is the analytic derivative of the
    // constant-free basis w.r.t. the active-direction components.
    g = Zero;
    const label nd = activeDirs_.size();

    // linear part: cf[0..nd-1] multiply d_a
    for (label a = 0; a < nd; ++a) { g[activeDirs_[a]] += cf[a]; }
    if (nc == nd) { return; }             // linear fit: constant gradient

    // pure quadratic: cf[nd + a] multiplies 0.5 d_a^2 -> derivative cf[nd+a] d_a
    for (label a = 0; a < nd; ++a)
    {
        g[activeDirs_[a]] += cf[nd + a]*d[activeDirs_[a]];
    }
    // cross terms: cf[k] multiplies d_a d_b -> d/d a += cf d_b, d/d b += cf d_a
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


void Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::setupBasis()
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


void Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::build()
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
            ncoeff = (nNbr >= nd) ? nd : 0;       // linear fallback, else constant
            ++nFallback;
        }
        ncoeff_[c] = ncoeff;
    }

    reduce(nFallback, sumOp<label>());
    if (nFallback > 0)
    {
        Info<< "signedDistanceQuadraticWeightedLeastSquares: " << nFallback
            << " cells with under-determined stencils fell back to a linear fit"
            << endl;
    }
    built_ = true;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::update
(
    const volScalarField& psiOld
)
{
    collectStencil(psiOld);       // fills stencilPsi_ (+ stencilC_ once); sets ptr

    if (!built_) { build(); }

    coeffsFlat_.setSize(mesh_.nCells()*ncoeffFull_);

    scalar A[81];      // ncoeff x ncoeff normal-equations matrix (<= 9x9)
    scalar g[9];       // rhs / solution
    scalar brow[9];    // one basis row

    label nSingular = 0;
    forAll(stencilPsi_, c)
    {
        const label nc = ncoeff_[c];
        if (nc == 0) { continue; }    // constant reconstruction (= psi_c)

        const List<scalar>& s = stencilPsi_[c];
        const label nNbr = nNbr_[c];
        const point xc = stencilC(c, 0);
        const scalar psiC = s[0];

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
        scalar* cf = &coeffsFlat_[c*ncoeffFull_];
        for (label k = 0; k < nc; ++k) { cf[k] = ok ? g[k] : 0; }
    }

    reduce(nSingular, sumOp<label>());
    if (nSingular > 0)
    {
        WarningInFunction
            << "signedDistanceQuadraticWeightedLeastSquares: normal equations not"
            << " SPD in " << nSingular << " cells (near-degenerate stencil); fell"
            << " back to the constant reconstruction there." << endl;
    }

    // --- strain data for the (1 + eps_nn*dt) rescaling ---------------------- //
    // U in the registry is the new-time velocity (u^{n+1}) at this point -- the
    // same field slAdvection uses for the foot -- so grad U here == its gradU.
    haveStrain_ = false;
    if (normalStrainRescale_ && mesh_.foundObject<volVectorField>("U"))
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        gradU_ = fvc::grad(U, "gradU")().primitiveField();
        dt_ = mesh_.time().deltaTValue();
        haveStrain_ = true;
    }

    computeLimiters();
}


Foam::scalar Foam::signedDistanceQuadraticWeightedLeastSquaresReconstruction::evaluateRaw
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

    // value P0 = P(x_d)
    scalar P0 = psiC;
    for (label k = 0; k < nc; ++k) { P0 += cf[k]*b[k]; }

    // Band-only (the spec is band-only): the reprojection is a near-interface
    // distance; far from the interface (large |P0|) keep the plain value-copy, so a
    // noisy far-field gradient cannot feed garbage back into the band.
    if (reprojectBandOnly_ && mag(P0) > reprojBandRadii_*stencilRadius(c))
    {
        return P0;
    }

    // gradient g = grad P(x_d)
    vector gvec(Zero);
    gradFromCoeffs(cf, nc, d, gvec);
    const scalar gmag = mag(gvec);

    // Trust region: only reproject where the fitted |grad P| is near the SDF value 1;
    // outside [gradLo, gradHi] the P0/|grad P| division is untrustworthy (would blow
    // up / flatten), so fall back to the raw value P0.
    if (gmag < gradLo_ || gmag > gradHi_) { return P0; }

    const vector n = gvec/(gmag + SMALL);   // defensive epsilon guard (gmag is
                                            // trust-region-bounded above, but +SMALL
                                            // is negligible and avoids any 0/0)

    // Signed distance from the foot to the ZERO set of the local quadratic, by a Newton
    // flow to the surface using the reconstruction's value and symbolic gradient
    // grad P(delta) = g + H*delta at each step -- NOT a single gradient-descent step.
    // delta starts at the foot displacement (x_d - x_c) and converges onto {P = 0};
    // footpointNewtonIters_ = 1 reproduces the linear P0/|grad P|.
    vector delta = d;
    for (label it = 0; it < footpointNewtonIters_; ++it)
    {
        scalar bb[9];
        basis(delta, nc, bb);
        scalar p = psiC;
        for (label k = 0; k < nc; ++k) { p += cf[k]*bb[k]; }
        vector gr(Zero);
        gradFromCoeffs(cf, nc, delta, gr);
        const scalar g2 = gr & gr;
        if (g2 < SMALL) { break; }
        delta -= (p/g2)*gr;
    }
    scalar d1 = Foam::sign(P0)*Foam::mag(delta - d);   // distance foot -> quadric surface

    // normal-strain rescaling: old-interface distance at the foot -> new-interface
    // distance at the arrival cell. D at the arrival cell c; sign is +.
    if (haveStrain_)
    {
        const tensor gu = gradU_[c];
        const tensor D = 0.5*(gu + gu.T());
        const scalar epsNN = n & (D & n);
        d1 *= (1.0 + epsNN*dt_);
    }

    return d1;
}

// ************************************************************************* //
