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

#include "normalProjectedScheme.H"
#include "slReconstruction.H"
#include "slCorrector.H"
#include "fvSolution.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(normalProjectedScheme, 0);
    addToRunTimeSelectionTable(slScheme, normalProjectedScheme, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalProjectedScheme::normalProjectedScheme(const fvMesh& mesh)
:
    slScheme(mesh),
    minGradPsi_(0.05),
    haveHistory_
    (
        // Restart detection: a previous run wrote uProjOld.nSL (AUTO_WRITE
        // below), so its presence in the start-time directory means a valid
        // previous-step projection exists and the AB2 startup is skipped.
        IOobject
        (
            "uProjOld.nSL",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ).typeHeaderOk<volVectorField>(true)
    ),
    uProjOld_
    (
        IOobject
        (
            "uProjOld.nSL",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        "zeroGradient"
    )
{
    const fvSolution& fvSol(mesh);
    const dictionary& slDict =
        fvSol.subDict("levelSet").subOrEmptyDict("semiLagrangian");

    minGradPsi_ = slDict.getOrDefault<scalar>("minGradPsi", 0.05);

    if (minGradPsi_ <= 0)
    {
        FatalIOErrorInFunction(slDict)
            << "minGradPsi must be positive, got " << minGradPsi_
            << exit(FatalIOError);
    }

    Info<< "normalProjected SL scheme: minGradPsi = " << minGradPsi_
        << ", previous-step projection "
        << (haveHistory_ ? "read from disk" : "not found (AB2 startup)")
        << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::normalProjectedScheme::advance
(
    volScalarField& psi,
    const volVectorField& Unew,
    const volVectorField& Uold,
    slReconstruction& recon,
    slCorrector& corrector
)
{
    // Unew is deliberately unused: the departure-centred trace needs only the
    // time-n velocity (Uold at every call site) and the stored projection of
    // the previous step -- the only centring whose normals exist at their own
    // time level (see the class Description). The corrector is unused: there
    // is no reconstruction evaluation at the foot to correct, and its
    // stencil-bounds clip would be wrong for a renormalized value.

    // Fits of psi^n. The coefficients are increments about the LIVE psi[c],
    // so psi must not be modified until the per-cell loop has read everything
    // it needs -- assemble into scratch fields and assign at the end.
    recon.update(psi);

    const scalar dt = mesh_.time().deltaTValue();
    const scalar dt0 = mesh_.time().deltaT0Value();
    // The temporal difference u_n^n - u_n^{n-1} spans the PREVIOUS step while
    // the trace spans the CURRENT one; the ratio keeps the Adams-Bashforth
    // extrapolation time-centred under a varying step (identical role to the
    // dtRatio of the coupled solver's call site).
    const scalar r = (dt0 > SMALL) ? dt/dt0 : 1.0;

    const volVectorField& C = mesh_.C();

    // grad(u^n) for the analytic assembly of grad(u_n). OpenFOAM convention
    // throughout: (grad u)_ij = d u_j / d x_i, so (a & gradU) is (a . grad) u.
    const volTensorField gradU(fvc::grad(Uold, "gradU"));

    scalarField newPsi(mesh_.nCells());
    vectorField uProj(mesh_.nCells(), Zero);

    label nFrozen = 0;
    label nFallback = 0;
    label nOutside = 0;
    label nNonFinite = 0;
    scalar maxRatio = 0;

    forAll(newPsi, c)
    {
        const scalar psiC = psi[c];

        vector g(Zero);
        symmTensor H(Zero);
        const label order = recon.fitDerivatives(c, g, H);
        const scalar gm = Foam::mag(g);

        // Freeze: no usable normal ray (constant-fit fallback, the clamped
        // far-field plateau, or a skeleton/ridge cell). uProj stays Zero,
        // which the AB2 restart below treats as "no history" on reactivation.
        if (order == 0 || gm < minGradPsi_)
        {
            newPsi[c] = psiC;
            ++nFrozen;
            continue;
        }

        const vector n = g/gm;
        const scalar vn = Uold[c] & n;
        const vector un = vn*n;
        uProj[c] = un;

        // Previous-step projection (previous step's normals, stored). A cell
        // frozen last step stored exactly Zero: restart its AB2 from the
        // current projection (one first-order step for that cell, standard
        // Adams-Bashforth startup).
        vector unOld = haveHistory_ ? uProjOld_[c] : un;
        if (magSqr(unOld) == 0) { unOld = un; }

        // Analytic grad(u_n) = grad(v_n) (x) n + v_n grad(n), assembled from
        // the fit -- never by differencing the projected field. In OpenFOAM
        // index convention ((grad u)_ij = d_i u_j):
        //   d_i n_j   = (H_ij - (H.n)_i n_j)/|g|
        //   d_i v_n   = (gradU . n)_i + (gradN . u)_i
        //   d_i(u_n)_j = (grad v_n)_i n_j + v_n (grad n)_ij
        const tensor gradN = (tensor(H) - (H & n)*n)/gm;
        const vector gradVn = (gradU[c] & n) + (gradN & Uold[c]);
        const tensor gradUn = gradVn*n + vn*gradN;

        // (u_n . grad) u_n -- same contraction as the baseline's convective
        // term (pointValueScheme: uNew & gradU[c]).
        const vector conv = un & gradUn;

        // Departure-centred second-order foot along u_n; at r = 1 the linear
        // displacement is exactly -(dt/2)(3 u_n^n - u_n^{n-1}) = the
        // time-centred midpoint speed.
        const point xd =
            C[c] - dt*(un + 0.5*r*(un - unOld)) + 0.5*dt*dt*conv;

        // Foot-radius guard (inline replica of slCorrector::footRadiusGuard,
        // which is protected there): warn-only, no clamp.
        const scalar disp = Foam::mag(xd - C[c]);
        const scalar radius = recon.stencilRadius(c);
        if (radius > SMALL)
        {
            maxRatio = Foam::max(maxRatio, disp/radius);
            if (disp > radius) { ++nOutside; }
        }

        // Geometric offset of the cell from the fit's zero set (stable
        // quadratic root; first-order fallback) plus the normal displacement
        // of the trace. Zero-preservation of the conversion guarantees the
        // renormalization never moves the front.
        bool fb = false;
        const scalar dC = recon.signedOffset(c, fb);
        if (fb) { ++nFallback; }

        const scalar delta = (xd - C[c]) & n;

        scalar v = dC + delta;
        if (!std::isfinite(v))
        {
            v = psiC;
            ++nNonFinite;
        }
        newPsi[c] = v;
    }

    psi.primitiveFieldRef() = newPsi;
    psi.correctBoundaryConditions();

    // Optional band-model far-field fix-up (no-op for the standard models).
    recon.postAdvect(psi);

    // Store THIS step's projections (this step's normals) as the next step's
    // u_n^{n-1} -- the design note's storage requirement (Sec. 4.2).
    uProjOld_.primitiveFieldRef() = uProj;
    uProjOld_.correctBoundaryConditions();
    haveHistory_ = true;

    reduce(nFrozen, sumOp<label>());
    reduce(nFallback, sumOp<label>());
    reduce(nOutside, sumOp<label>());
    reduce(nNonFinite, sumOp<label>());
    reduce(maxRatio, maxOp<scalar>());

    Info<< "normalProjected SL: frozen " << nFrozen
        << ", offset fallback " << nFallback
        << ", feet outside stencil " << nOutside
        << " (max |x_d - x_c|/radius = " << maxRatio << ")"
        << ", non-finite reset " << nNonFinite << endl;

    if (nOutside > 0)
    {
        WarningInFunction
            << "normal-projected foot left the stencil in " << nOutside
            << " cells (max ratio " << maxRatio
            << "); CFL likely > 1 -- reduce maxCo." << endl;
    }
}

// ************************************************************************* //
