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

#include "defectCorrectedIDWReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(defectCorrectedIDWReconstruction, 0);
    addToRunTimeSelectionTable
    (
        slReconstruction,
        defectCorrectedIDWReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::defectCorrectedIDWReconstruction::defectCorrectedIDWReconstruction
(
    const fvMesh& mesh
)
:
    slReconstruction(mesh),
    gradPsi_
    (
        IOobject
        (
            "slIdwGradPsi",   // unique name: several models may coexist on one mesh
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("slIdwGradPsi", dimless, vector::zero),
        "zeroGradient"
    ),
    stencilGrad_(),
    power_(slDict_.getOrDefault<scalar>("idwPower", 2)),
    delta_(),
    m_(),
    maxM_(0),
    built_(false),
    nIters_(slDict_.getOrDefault<label>("dcIters", 3)),
    maxIters_(slDict_.getOrDefault<label>("dcMaxIters", 10)),
    dcTol_(slDict_.getOrDefault<scalar>("dcTol", 1e-3)),
    omega_(Foam::min(scalar(1), slDict_.getOrDefault<scalar>("dcRelax", 0.8))),
    regEps_(slDict_.getOrDefault<scalar>("dcRegEps", 0.1))
{
    maxIters_ = Foam::max(maxIters_, nIters_);
    Info<< "defectCorrectedIDW: iterative nodal-defect correction, p=" << power_
        << " dcIters=" << nIters_ << " dcMaxIters=" << maxIters_
        << " dcRelax=" << omega_ << " dcRegEps=" << regEps_
        << " dcTol=" << dcTol_ << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::defectCorrectedIDWReconstruction::update(const volScalarField& psiOld)
{
    // Cell-point-cell LSQ gradient (fvSchemes key "gradPsi" = pointCellsLeastSquares,
    // linear-exact). Values-only assignment: psi is a signed distance, so a
    // dimensioned operator= on the dimless helper field would abort.
    const volVectorField g(fvc::grad(psiOld, "gradPsi"));
    gradPsi_.primitiveFieldRef() = g.primitiveField();
    gradPsi_.correctBoundaryConditions();

    // psiOld + stencil centres (+ cached radius); also sets psiOldPtr_.
    collectStencil(psiOld);

    // grad(psi) over the SAME haloed CPC stencil (incl. remote neighbours),
    // aligned with stencilC_/stencilPsi_ for the per-node Taylor lift.
    stencil_.collectData(gradPsi_, stencilGrad_);

    // One-off sizing (mesh static within a run).
    if (!built_)
    {
        const label nc = mesh_.nCells();
        m_.setSize(nc);
        delta_.setSize(nc);
        maxM_ = 0;
        for (label c = 0; c < nc; ++c)
        {
            m_[c] = stencilSize(c);
            delta_[c].setSize(m_[c], 0.0);
            maxM_ = Foam::max(maxM_, m_[c]);
        }
        built_ = true;
    }

    // Reused per-cell work buffers (avoid per-cell allocation at 128^3).
    scalarField Aflat(maxM_*maxM_);   // row-major regularised row-stochastic blend
    scalarField Adiag(maxM_);
    scalarField b(maxM_);

    label nFallback = 0;
    List<vector> X(maxM_);   // reused per-cell stencil centres (from mesh_.C()/tail)

    forAll(stencilPsi_, c)
    {
        const label m = m_[c];
        scalarList& d = delta_[c];
        d = 0.0;                               // reset the defects each step
        if (m < 2) { continue; }               // degenerate stencil -> single pass

        for (label i = 0; i < m; ++i) { X[i] = stencilC(c, i); }  // [0] == cell c
        const List<scalar>& P = stencilPsi_[c];
        const List<vector>& G = stencilGrad_[c];

        const scalar hc2  = Foam::max(radius_[c]*radius_[c], VSMALL);
        const scalar eps2 = regEps_*regEps_*hc2;
        const scalar pe   = 0.5*power_;

        // Regularised, row-stochastic blend matrix A_{DC} = wreg_DC / sum_E wreg_DE.
        // wreg(r) = (h_c^2/(|r|^2 + eps^2 h_c^2))^{p/2}; diagonal = (1/regEps^2)^{p/2}.
        bool dominant = true;
        for (label D = 0; D < m; ++D)
        {
            scalar rowSum = 0;
            scalar* row = &Aflat[D*m];
            for (label C = 0; C < m; ++C)
            {
                const scalar r2 = Foam::magSqr(X[D] - X[C]);
                const scalar w  = Foam::pow(hc2/(r2 + eps2), pe);
                row[C] = w;
                rowSum += w;
            }
            const scalar inv = 1.0/Foam::max(rowSum, VSMALL);
            for (label C = 0; C < m; ++C) row[C] *= inv;
            Adiag[D] = row[D];
            if (Adiag[D] <= 0.5) { dominant = false; }
        }

        // Weak diagonal dominance (near-degenerate/irregular stencil): the
        // Gauss-Seidel contraction is not guaranteed -> keep delta = 0 (the
        // single-pass blend) for this cell rather than iterate a non-contraction.
        if (!dominant) { ++nFallback; continue; }

        // Fixed gradient-defect load b_D = sum_C A_{DC} grad_C.(x_D - x_C).
        for (label D = 0; D < m; ++D)
        {
            const scalar* row = &Aflat[D*m];
            scalar acc = 0;
            for (label C = 0; C < m; ++C)
            {
                acc += row[C]*(G[C] & (X[D] - X[C]));
            }
            b[D] = acc;
        }

        scalar lo, hi;
        stencilRange(c, lo, hi);
        const scalar scale = (hi - lo) + VSMALL;

        // Damped Gauss-Seidel sweeps for the nodal defects; arrival node pinned.
        for (label it = 0; it < maxIters_; ++it)
        {
            scalar rmax = 0;
            d[0] = 0.0;                          // pin self-defect -> preserves contract
            for (label D = 1; D < m; ++D)
            {
                const scalar* row = &Aflat[D*m];
                scalar Su = b[D];
                for (label C = 0; C < m; ++C) { Su += row[C]*(P[C] + d[C]); }
                const scalar rD = P[D] - Su;     // nodal residual
                d[D] += omega_*rD/Adiag[D];
                rmax = Foam::max(rmax, Foam::mag(rD));
            }
            if (it + 1 >= nIters_ && rmax <= dcTol_*scale) { break; }
        }
        d[0] = 0.0;                              // re-pin (safety)
    }

    reduce(nFallback, sumOp<label>());
    if (nFallback > 0)
    {
        Info<< "defectCorrectedIDW: single-pass fallback (weak diagonal dominance) in "
            << nFallback << " cells" << endl;
    }

    computeLimiters();
}


Foam::scalar Foam::defectCorrectedIDWReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const List<scalar>& P = stencilPsi_[c];
    const List<vector>& G = stencilGrad_[c];
    const scalarList&   d = delta_[c];          // cached nodal defects (d[0] == 0)

    scalar sumW = 0;
    scalar phid = 0;
    const label m = stencilSize(c);
    for (label i = 0; i < m; ++i)
    {
        const vector r = x - stencilC(c, i);    // [0] == arrival cell c
        // Corrected (effective-nodal-value) Taylor value from stencil cell i.
        const scalar taylor = P[i] + d[i] + (G[i] & r);
        const scalar dist = Foam::mag(r);
        if (dist < SMALL)
        {
            // x coincides with this cell centre (in particular x == x_c, i == 0,
            // d[0] == 0): return the exact value -> evaluateRaw(c, x_c) == psi_c.
            return taylor;
        }
        const scalar w = 1.0/Foam::pow(dist, power_);
        phid += w*taylor;
        sumW += w;
    }
    return phid/Foam::max(sumW, SMALL);
}

// ************************************************************************* //
