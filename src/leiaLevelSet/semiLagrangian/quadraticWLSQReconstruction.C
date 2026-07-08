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

#include "quadraticWLSQReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "QRMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadraticWLSQReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        quadraticWLSQReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadraticWLSQReconstruction::quadraticWLSQReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    activeDirs_(),
    ncoeffFull_(0),
    MW_(),
    nNbr_(),
    maxNbr_(0),
    ncoeff_(),
    coeffs_(),
    built_(false)
{}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::quadraticWLSQReconstruction::basis
(
    const vector& d,
    const label ncoeff,
    scalar* b
) const
{
    const label nd = activeDirs_.size();

    // Linear part: d along each active direction.
    forAll(activeDirs_, a)
    {
        b[a] = d[activeDirs_[a]];
    }
    if (ncoeff == nd)
    {
        return;   // linear fallback basis
    }

    // Quadratic part: 1/2 d_a^2 (pure) then d_a d_b (cross, a < b).
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


void Foam::quadraticWLSQReconstruction::setupBasis()
{
    // Active (non-empty) directions from the mesh: geometricD() is +1 for a
    // used direction, -1 for an empty (2D) one.
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


void Foam::quadraticWLSQReconstruction::build()
{
    setupBasis();
    const label nd = activeDirs_.size();

    const label nCells = mesh_.nCells();
    MW_.setSize(nCells);
    nNbr_.setSize(nCells);
    ncoeff_.setSize(nCells);
    maxNbr_ = 0;

    label nFallback = 0;
    scalarList brow(ncoeffFull_);   // reused basis buffer (one row of A)
    forAll(stencilC_, c)
    {
        const List<vector>& C = stencilC_[c];     // [0] = arrival cell c
        const label nNbr = C.size() - 1;          // skip self
        const point xc = C[0];
        nNbr_[c] = nNbr;
        maxNbr_ = Foam::max(maxNbr_, nNbr);

        // Full quadratic if over-determined, else linear fallback, else const.
        label ncoeff = ncoeffFull_;
        if (nNbr < ncoeffFull_)
        {
            ncoeff = (nNbr >= nd) ? nd : 0;
            ++nFallback;
        }
        ncoeff_[c] = ncoeff;

        if (ncoeff == 0)
        {
            MW_[c].clear();   // constant reconstruction (= psi_c)
            continue;
        }

        // Weighted design matrix A_{ik} = w_i b_k(d_i); Minv = pinv(A).
        scalarRectangularMatrix A(nNbr, ncoeff, Zero);
        scalarList wloc(nNbr);
        for (label i = 0; i < nNbr; ++i)
        {
            const vector d = C[i + 1] - xc;
            const scalar wi = 1.0/Foam::max(Foam::mag(d), SMALL);
            wloc[i] = wi;
            basis(d, ncoeff, brow.data());
            for (label k = 0; k < ncoeff; ++k)
            {
                A(i, k) = wi*brow[k];
            }
        }
        const scalarRectangularMatrix Minv(MatrixTools::pinv(A));

        // Fold the IDW weight into the map and store flat row-major so the
        // per-step update is a contiguous mat-vec: coeffs_k = sum_i MW[k,i] ds_i,
        // ds_i = psi_i - psi_c (weight already applied).
        MW_[c].setSize(ncoeff*nNbr);
        for (label k = 0; k < ncoeff; ++k)
        {
            scalar* row = MW_[c].data() + k*nNbr;
            for (label i = 0; i < nNbr; ++i)
            {
                row[i] = Minv(k, i)*wloc[i];
            }
        }
    }

    reduce(nFallback, sumOp<label>());
    if (nFallback > 0)
    {
        Info<< "quadraticWLSQ: " << nFallback
            << " cells with under-determined stencils fell back to a linear fit"
            << endl;
    }
    built_ = true;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::quadraticWLSQReconstruction::update(const volScalarField& psiOld)
{
    collectStencil(psiOld);   // fills stencilPsi_ (+ stencilC_ once); sets ptr

    if (!built_)
    {
        build();
    }

    coeffs_.setSize(mesh_.nCells());
    scalarList ds(maxNbr_);   // reused neighbour-difference scratch (no per-cell alloc)
    forAll(stencilPsi_, c)
    {
        const label ncoeff = ncoeff_[c];
        coeffs_[c].setSize(ncoeff);
        if (ncoeff == 0)
        {
            continue;   // constant reconstruction
        }

        const List<scalar>& s = stencilPsi_[c];   // [0] = psi_c
        const scalar psiC = s[0];
        const label nNbr = nNbr_[c];
        const scalar* sp = s.cdata();

        // Neighbour differences (weight is folded into MW_).
        scalar* dsp = ds.data();
        for (label i = 0; i < nNbr; ++i)
        {
            dsp[i] = sp[i + 1] - psiC;
        }
        // coeffs_k = sum_i MW[k,i] ds_i. Contiguous rows, folded weight, no
        // per-cell allocation. (A manual 4-accumulator unroll was measured to
        // REGRESS here: the stencil is short (~8-12 neighbours), so the unroll
        // setup/combine outweighs any pipelining gain and interferes with the
        // -O3 heuristics. The plain loop is both faster and clearer.)
        const scalar* MW = MW_[c].cdata();
        scalar* cf = coeffs_[c].data();
        for (label k = 0; k < ncoeff; ++k)
        {
            const scalar* row = MW + k*nNbr;
            scalar acc = 0;
            for (label i = 0; i < nNbr; ++i)
            {
                acc += row[i]*dsp[i];
            }
            cf[k] = acc;
        }
    }
    computeLimiters();
}


Foam::scalar Foam::quadraticWLSQReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];
    const label ncoeff = ncoeff_[c];
    if (ncoeff == 0)
    {
        return psiC;
    }

    const vector d = x - mesh_.C()[c];
    scalar b[9];   // max basis size (3D full quadratic); no allocation
    basis(d, ncoeff, b);
    const scalar* cf = coeffs_[c].cdata();
    scalar val = psiC;
    for (label k = 0; k < ncoeff; ++k)
    {
        val += cf[k]*b[k];
    }
    return val;
}

// ************************************************************************* //
