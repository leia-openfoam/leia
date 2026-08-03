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
    bandRadii_(1.0),
    renormalization_("strain"),
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
    bandRadii_ = slDict.getOrDefault<scalar>("bandRadii", 1.0);
    renormalization_ =
        slDict.getOrDefault<word>("renormalization", "strain");

    if
    (
        renormalization_ != "geometric"
     && renormalization_ != "strain"
     && renormalization_ != "none"
    )
    {
        FatalIOErrorInFunction(slDict)
            << "renormalization must be geometric, strain or none, got '"
            << renormalization_ << "'" << exit(FatalIOError);
    }

    if (minGradPsi_ <= 0 || bandRadii_ <= 0)
    {
        FatalIOErrorInFunction(slDict)
            << "minGradPsi and bandRadii must be positive, got "
            << minGradPsi_ << ", " << bandRadii_ << exit(FatalIOError);
    }
    Info<< "normalProjected SL scheme: minGradPsi = " << minGradPsi_
        << ", bandRadii = " << bandRadii_
        << ", renormalization = " << renormalization_
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
    label nFar = 0;
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

        const scalar radius = recon.stencilRadius(c);

        // Foot-radius guard (inline replica of slCorrector::footRadiusGuard,
        // which is protected there): warn-only, no clamp.
        const scalar disp = Foam::mag(xd - C[c]);
        if (radius > SMALL)
        {
            maxRatio = Foam::max(maxRatio, disp/radius);
            if (disp > radius) { ++nOutside; }
        }

        // Geometric offset of the cell from the fit's zero set, plus the
        // normal displacement of the trace. Zero-preservation of the
        // conversion guarantees the renormalization never moves the front.
        //
        // Conversion ladder (see the header): the quadratic root is meaningful
        // only for BAND-SCALE offsets -- its psi_c*h_nn term amplifies fitted-
        // Hessian noise linearly in the offset -- so it is gated to
        // |psi_c|/|g| <= bandRadii*stencilRadius (the curvature path's
        // iso-agnostic band gate). Outside, the first-order psi_c/|g| is exact
        // on a clean signed distance and Hessian-free; it needs |g| inside the
        // trust region, else the cell is geometrically degenerate (cone tip,
        // skeleton) and freezes for the step.
        const scalar ratio = Foam::mag(psiC)/gm;
        const scalar delta = (xd - C[c]) & n;

        scalar v;
        if
        (
            renormalization_ == "geometric"
         && ratio <= bandRadii_*radius
        )
        {
            // Band: renormalizing geometric update d_c + delta. MEASURED
            // UNSTABLE as a per-step write-back (see the header); retained as
            // a selectable mode for the record and for future damped variants.
            bool fb = false;
            const scalar dC = recon.signedOffset(c, fb);
            if (fb) { ++nFallback; }
            v = dC + delta;
        }
        else if
        (
            renormalization_ == "strain"
         && ratio <= bandRadii_*radius
        )
        {
            // BAND ONLY -- the gate is essential: the ODE below describes the
            // stretching of a band-scale DISTANCE whose ray foot experiences
            // eps_nn. Applied to the large far-field psi with the LOCAL
            // eps_nn, the factor compounds exp(int eps dt) over the run --
            // measured on the reversed 2D vortex as an exponential far-field
            // blow-up (band gradient error 1e6..1e34 across the ladder)
            // before this gate existed.

            // Raw ray transport + the exact gradient-magnitude ODE
            // D|grad psi|/Dt = -|grad psi| eps_nn, integrated one step:
            // distances stretch by (1 + eps_nn dt) under normal strain
            // (design note Sec. 7.1). eps_nn comes from grad(u) -- smooth,
            // solver-supplied -- and the fit enters only through n, with the
            // whole correction multiplied by dt: noise gain O(dt), unlike the
            // geometric write-back.
            const tensor& gU = gradU[c];
            const scalar epsNN = n & (symm(gU) & n);
            v = (psiC + delta*gm)*(1.0 + epsNN*dt);
            ++nFar;
        }
        else
        {
            // Far field: RAW first-order transport along the normal ray,
            // psi^{n+1} = psi(x_c + delta n) = psi_c + delta |g| + O(delta^2).
            // Deliberately NOT renormalized: any conversion divides by the
            // fitted |g| and so amplifies fit error linearly in the offset
            // (measured ~30 h one-step errors at boundary-ring cells with
            // one-sided stencils). Here fit error enters only through
            // delta*|g| -- a per-step transport increment of ~1e-2 h -- and
            // a plateau cell (|g| ~ 0) freezes automatically. The far field
            // keeps its own |grad psi| like the baseline scheme; only the
            // band is renormalized, and only the band feeds the physics.
            v = psiC + delta*gm;
            ++nFar;
        }
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
    reduce(nFar, sumOp<label>());
    reduce(nFallback, sumOp<label>());
    reduce(nOutside, sumOp<label>());
    reduce(nNonFinite, sumOp<label>());
    reduce(maxRatio, maxOp<scalar>());

    Info<< "normalProjected SL: frozen " << nFrozen
        << ", far-field first-order " << nFar
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
