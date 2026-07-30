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

#include "integralSurfaceTension.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(integralSurfaceTension, false);
addToRunTimeSelectionTable(surfaceTensionForce, integralSurfaceTension, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

integralSurfaceTension::integralSurfaceTension(const fvMesh& mesh)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    psi_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("levelSet", "psi")
        )
    ),
    kappa_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("kappa", "kappa")
        )
    ),
    curvatureSource_
    (
        surfTensionDict_.getOrDefault<word>("curvatureSource", "field")
    ),
    constantCurvature_
    (
        surfTensionDict_.getOrDefault<scalar>("curvature", 0.0)
    ),
    // Default ON: the psi-derivative interface tangents are divided by the local
    // |grad psi|, SMOOTHED (3 Jacobi passes -- kills the few-percent grid-scale noise
    // of central differences that otherwise re-enters the stress as a vortical driver)
    // and interpolated CONSISTENTLY to each crossing (every stress entry touching a
    // crossing uses the same G). This makes the scheme ISO-VALUE-AGNOSTIC: MEASURED on
    // the static N=32 droplet, the SDF circle completes with decaying currents and
    // Laplace 73.6/72.74, and the ALGEBRAIC equal-axes-ellipsoid circle (|grad psi|
    // ~1000, no redistancing) ALSO completes with decaying currents (0.26 -> 5.3e-3)
    // and Laplace 71-75.5 -- where the unnormalised scheme gives pLaplace ~1.5e5 and
    // blows up at t~1.9e-4. Set false only to reproduce the historical SDF-only form.
    normalizeGrad_(surfTensionDict_.getOrDefault<bool>("normalizeGradient", true)),
    addressingBuilt_(false),
    nx_(0), ny_(0), dx_(0),
    dirX_(-1), dirY_(-1)
{
    if
    (
        curvatureSource_ != "field"
     && curvatureSource_ != "constant"
     && curvatureSource_ != "none"
    )
    {
        FatalIOErrorInFunction(surfTensionDict_)
            << "Unknown curvatureSource '" << curvatureSource_ << "'. "
            << "Valid choices are field, constant and none."
            << exit(FatalIOError);
    }

    Info<< "integralSurfaceTension (CST, Abu-Al-Saud/Popinet/Tchelepi 2018): "
        << "2D uniform-hex prototype; curvatureSource=" << curvatureSource_;
    if (curvatureSource_ == "constant")
    {
        Info<< ", curvature=" << constantCurvature_ << " 1/m";
    }
    Info<< endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void integralSurfaceTension::buildAddressing() const
{
    // Two non-empty global directions
    const Vector<label> gd = mesh_.geometricD();
    DynamicList<label> dirs(2);
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        if (gd[cmpt] == 1) { dirs.append(cmpt); }
    }
    if (dirs.size() != 2)
    {
        FatalErrorInFunction
            << "integralSurfaceTension: prototype requires a 2D mesh "
            << "(exactly one empty direction)." << abort(FatalError);
    }
    dirX_ = dirs[0]; dirY_ = dirs[1];

    const vectorField& C = mesh_.cellCentres();
    scalar xmin = GREAT, xmax = -GREAT, ymin = GREAT, ymax = -GREAT;
    forAll(C, c)
    {
        xmin = min(xmin, C[c][dirX_]); xmax = max(xmax, C[c][dirX_]);
        ymin = min(ymin, C[c][dirY_]); ymax = max(ymax, C[c][dirY_]);
    }

    // uniform spacing from the smallest centre-to-centre distance
    scalar d = GREAT;
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, f) { d = min(d, mag(C[own[f]] - C[nei[f]])); }
    dx_ = d;

    nx_ = label(std::lround((xmax - xmin)/dx_)) + 1;
    ny_ = label(std::lround((ymax - ymin)/dx_)) + 1;

    if (nx_*ny_ != mesh_.nCells())
    {
        FatalErrorInFunction
            << "integralSurfaceTension: mesh is not a uniform single-block 2D "
            << "hex grid (nx*ny = " << nx_ << "*" << ny_ << " != nCells = "
            << mesh_.nCells() << ")." << abort(FatalError);
    }

    cellI_.setSize(mesh_.nCells());
    cellJ_.setSize(mesh_.nCells());
    ijCell_.setSize(mesh_.nCells(), -1);
    forAll(C, c)
    {
        const label i = label(std::lround((C[c][dirX_] - xmin)/dx_));
        const label j = label(std::lround((C[c][dirY_] - ymin)/dx_));
        cellI_[c] = i; cellJ_[c] = j;
        ijCell_[i + j*nx_] = c;
    }
    forAll(ijCell_, k)
    {
        if (ijCell_[k] < 0)
        {
            FatalErrorInFunction
                << "integralSurfaceTension: structured (i,j) addressing failed "
                << "(non-uniform mesh?)." << abort(FatalError);
        }
    }
    addressingBuilt_ = true;

    Info<< "integralSurfaceTension: structured block " << nx_ << " x " << ny_
        << ", dx = " << dx_ << endl;
}


tmp<surfaceScalarField>
integralSurfaceTension::calcFaceSurfaceTensionForceFlux() const
{
    if (!addressingBuilt_) { buildAddressing(); }

    const scalar gamma = sigma_.value();
    const scalar dx = dx_;
    const scalarField& psi = psi_.primitiveField();
    const scalarField& kap = kappa_.primitiveField();

    auto P = [&](const label i, const label j) -> scalar
    {
        const label c = at(i, j);
        return (c >= 0) ? psi[c] : GREAT;    // GREAT: never a sign change
    };
    auto K = [&](const label i, const label j) -> scalar
    {
        if (curvatureSource_ == "none")
        {
            return 0.0;
        }
        if (curvatureSource_ == "constant")
        {
            return constantCurvature_;
        }
        const label c = at(i, j);
        return (c >= 0) ? kap[c] : 0.0;
    };

    // Local |grad psi| for the ISO-AGNOSTIC tangent normalisation. The stress below
    // builds interface tangents from psi-DERIVATIVES, which are UNIT tangents only for
    // a signed distance (|grad psi| = 1). Dividing by |grad psi| recovers the true unit
    // tangent for any smooth level set / algebraic implicit surface.
    //   CRITICAL: use the SAME structured central differences as the stress stencil,
    //   NOT fvc::grad (Gauss). On a true SDF the Gauss |grad psi| ~ 0.94 (6% low near a
    //   curved interface), and dividing by an INCONSISTENT gradient perturbs the stress
    //   out of balance -> the well-balanced property is lost and the drop blows up. The
    //   structured central difference is 1 + O(h^2) on an SDF (consistent), so the
    //   normalisation is a genuine near no-op there and the correction on a non-SDF.
    scalarField gpsi(mesh_.nCells(), 1.0);
    if (normalizeGrad_)
    {
        forAll(psi, c)
        {
            const label i = cellI_[c], j = cellJ_[c];
            const scalar pxm = P(i - 1, j), pxp = P(i + 1, j);
            const scalar pym = P(i, j - 1), pyp = P(i, j + 1);
            const scalar gx = (pxp < GREAT/2 && pxm < GREAT/2) ? (pxp - pxm)/(2*dx) : 0;
            const scalar gy = (pyp < GREAT/2 && pym < GREAT/2) ? (pyp - pym)/(2*dx) : 0;
            gpsi[c] = max(Foam::sqrt(gx*gx + gy*gy), SMALL);
        }
        // SMOOTH the normalisation field. |grad psi| is a smooth O(1) SCALE factor
        // (not interface data), but central differences of the discrete psi carry
        // few-percent GRID-SCALE noise; dividing the stress by a noisy G injects a
        // grid-scale tangential stress oscillation -- the measured vortical driver
        // that destabilised the first (unsmoothed) retrofit on the SDF regression.
        // Three Jacobi passes kill the grid-scale noise without biasing the scale.
        for (label s = 0; s < 3; ++s)
        {
            const scalarField g0(gpsi);
            forAll(gpsi, c)
            {
                const label i = cellI_[c], j = cellJ_[c];
                scalar sum = 0; label nn = 0;
                const label cxm = at(i - 1, j), cxp = at(i + 1, j);
                const label cym = at(i, j - 1), cyp = at(i, j + 1);
                if (cxm >= 0) { sum += g0[cxm]; ++nn; }
                if (cxp >= 0) { sum += g0[cxp]; ++nn; }
                if (cym >= 0) { sum += g0[cym]; ++nn; }
                if (cyp >= 0) { sum += g0[cyp]; ++nn; }
                if (nn > 0) { gpsi[c] = 0.5*g0[c] + 0.5*sum/nn; }
            }
        }
    }
    auto G = [&](const label i, const label j) -> scalar
    {
        if (!normalizeGrad_) { return 1.0; }
        const label c = at(i, j);
        return (c >= 0) ? gpsi[c] : 1.0;
    };

    // ---- diagonal stress components at cell centres (Algorithm 1) ----------
    // sigmaXX(i,j): crossings of the VERTICAL centre-line segment of cell (i,j)
    // (x-momentum control-volume face); tangent x-component from vertical psi
    // differences; jump correction -sign(psi_c)*kappa~*(1/2 - xi).
    scalarField sigmaXX(mesh_.nCells(), 0.0);
    scalarField sigmaYY(mesh_.nCells(), 0.0);

    forAll(psi, c)
    {
        const label i = cellI_[c], j = cellJ_[c];
        const scalar pc = psi[c];

        // sigma^xx: k = -1, +1 over vertical neighbours (i, j+k)
        for (label k = -1; k <= 1; k += 2)
        {
            const scalar pk = P(i, j + k);
            if (pk >= GREAT/2) { continue; }              // block edge
            if (pc*(pc + pk) > 0) { continue; }           // no crossing in half
            if (mag(pc - pk) < VSMALL) { continue; }      // degenerate segment
            const scalar pm = P(i, j - 1), pp = P(i, j + 1);
            if (pm >= GREAT/2 || pp >= GREAT/2) { continue; }
            const scalar xi = pc/(pc - pk);
            const scalar gN = G(i, j) + xi*(G(i, j + k) - G(i, j));   // |grad psi| @ crossing
            const scalar tx = (0.5*(pp - pm) + k*xi*(pm - 2.0*pc + pp))/dx/max(gN, SMALL);
            const scalar kt = K(i, j) + xi*(K(i, j + k) - K(i, j));
            sigmaXX[c] += gamma*(mag(tx)/dx - sign(pc)*kt*(0.5 - xi));
        }

        // sigma^yy: k over horizontal neighbours (i+k, j)
        for (label k = -1; k <= 1; k += 2)
        {
            const scalar pk = P(i + k, j);
            if (pk >= GREAT/2) { continue; }
            if (pc*(pc + pk) > 0) { continue; }
            if (mag(pc - pk) < VSMALL) { continue; }
            const scalar pm = P(i - 1, j), pp = P(i + 1, j);
            if (pm >= GREAT/2 || pp >= GREAT/2) { continue; }
            const scalar xi = pc/(pc - pk);
            const scalar gN = G(i, j) + xi*(G(i + k, j) - G(i, j));   // |grad psi| @ crossing
            const scalar ty = (0.5*(pp - pm) + k*xi*(pm - 2.0*pc + pp))/dx/max(gN, SMALL);
            const scalar kt = K(i, j) + xi*(K(i + k, j) - K(i, j));
            sigmaYY[c] += gamma*(mag(ty)/dx - sign(pc)*kt*(0.5 - xi));
        }
    }

    // ---- off-diagonal components at corners (Algorithm 2) ------------------
    // Corner (i-1/2, j-1/2) is indexed by (i,j), i in [0..nx], j in [0..ny]
    // (only interior corners with all four cells present are non-zero).
    // sigmaXY: x-momentum flux through the horizontal CV face; the interface
    // crossing is sought on the x-segment between the vertical-face midpoints.
    const label ncx = nx_ + 1, ncy = ny_ + 1;
    scalarField sigmaXY(ncx*ncy, 0.0);
    scalarField sigmaYX(ncx*ncy, 0.0);

    for (label j = 1; j < ny_; ++j)
    {
        for (label i = 1; i < nx_; ++i)
        {
            const scalar pmm = P(i - 1, j - 1);   // (i-1, j-1)
            const scalar pmp = P(i - 1, j);       // (i-1, j)
            const scalar ppm = P(i, j - 1);       // (i,   j-1)
            const scalar ppp = P(i, j);           // (i,   j)
            if (max(max(mag(pmm), mag(pmp)), max(mag(ppm), mag(ppp)))
                >= GREAT/2) { continue; }

            const label cidx = i + j*ncx;

            // |grad psi| at the SEGMENT ENDPOINTS (face midpoints) for the iso-agnostic
            // tangent normalisation; interpolated TO THE CROSSING with the same xi as
            // the tangent, so every stress entry touching a crossing is scaled by the
            // SAME consistent G (a flat corner average breaks the pairwise contact-
            // force telescoping of the scheme).
            const scalar gLm = 0.5*(G(i - 1, j - 1) + G(i - 1, j));   // left  x-face
            const scalar gRp = 0.5*(G(i, j - 1)     + G(i, j));       // right x-face
            const scalar gBm = 0.5*(G(i - 1, j - 1) + G(i, j - 1));   // bottom y-face
            const scalar gTp = 0.5*(G(i - 1, j)     + G(i, j));       // top    y-face

            // sigma^xy at corner (i-1/2, j-1/2)
            {
                const scalar phiL = pmp + pmm;      // left end (x2 the average)
                const scalar phiR = ppp + ppm;      // right end
                if (phiL*phiR <= 0 && mag(phiL - phiR) > VSMALL)
                {
                    const scalar xi = phiL/(phiL - phiR);
                    const scalar gX = max(gLm + xi*(gRp - gLm), SMALL);
                    const scalar tx =
                        (pmp - pmm + xi*(ppp - pmp + pmm - ppm))/dx/gX;
                    sigmaXY[cidx] = -gamma*sign(ppp + ppm)*tx/dx;
                }
            }

            // sigma^yx at the same corner (rotated indices)
            {
                const scalar phiB = ppm + pmm;      // bottom end
                const scalar phiT = ppp + pmp;      // top end
                if (phiB*phiT <= 0 && mag(phiB - phiT) > VSMALL)
                {
                    const scalar xi = phiB/(phiB - phiT);
                    const scalar gY = max(gBm + xi*(gTp - gBm), SMALL);
                    const scalar ty =
                        (ppm - pmm + xi*(ppp - ppm + pmm - pmp))/dx/gY;
                    sigmaYX[cidx] = -gamma*sign(ppp + pmp)*ty/dx;
                }
            }
        }
    }

    // ---- face-normal assembly: fSigma_f = [-div(sigma)] . n_f --------------
    tmp<surfaceScalarField> tfSigma
    (
        new surfaceScalarField
        (
            IOobject
            (
                "fSigmaCST",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "fSigmaCST",
                sigma_.dimensions()/dimLength/dimLength,   // N/m^3
                Zero
            )
        )
    );
    surfaceScalarField& fSigma = tfSigma.ref();
    scalarField& fIn = fSigma.primitiveFieldRef();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& C = mesh_.cellCentres();

    forAll(own, f)
    {
        const label cO = own[f], cN = nei[f];
        const vector d = C[cN] - C[cO];

        if (mag(d[dirX_]) > mag(d[dirY_]))
        {
            // x-face between (i-1, j) and (i, j); i = the +x cell's index
            const label sgn = (d[dirX_] > 0) ? 1 : -1;
            const label cP = (sgn > 0) ? cO : cN;    // -x cell (i-1, j)
            const label cU = (sgn > 0) ? cN : cO;    // +x cell (i,   j)
            const label i = cellI_[cU], j = cellJ_[cU];
            const scalar dSxx = (sigmaXX[cU] - sigmaXX[cP])/dx;
            const scalar dSxy =
                (sigmaXY[i + (j + 1)*ncx] - sigmaXY[i + j*ncx])/dx;
            // Overall sign: the paper's printed equilibrium relation
            // (App. B: p_i - p_{i-1} + dSigma = 0) delivers the INVERTED
            // Laplace jump on a signed-distance circle (hand-evaluated at the
            // poles, and measured: pLaplace = -72.7 instead of +72.7, with an
            // anti-restoring capillary instability growing at the grid-
            // capillary frequency). Their Eq. (8) uses +div(sigma) while their
            // stress identity states -div(pI + sigma) -- the physical overall
            // sign is +dSigma here (jump p_water = p_air + gamma*kappa, and
            // capillary waves restored, not amplified).
            fIn[f] = sgn*(dSxx + dSxy);
        }
        else
        {
            // y-face between (i, j-1) and (i, j)
            const label sgn = (d[dirY_] > 0) ? 1 : -1;
            const label cP = (sgn > 0) ? cO : cN;    // -y cell (i, j-1)
            const label cU = (sgn > 0) ? cN : cO;    // +y cell (i, j)
            const label i = cellI_[cU], j = cellJ_[cU];
            const scalar dSyy = (sigmaYY[cU] - sigmaYY[cP])/dx;
            const scalar dSyx =
                (sigmaYX[(i + 1) + j*ncx] - sigmaYX[i + j*ncx])/dx;
            fIn[f] = sgn*(dSyy + dSyx);   // sign: see x-face comment
        }
    }

    // Boundary faces: zero (the interface must stay away from boundaries in
    // this prototype; the droplet case satisfies this by construction).

    // Sanity report (cheap; also traps NaN/Inf before they reach the pressure
    // solver as an unexplained GAMG FPE).
    scalar fmin = GREAT, fmax = -GREAT;
    label nBad = 0;
    forAll(fIn, f)
    {
        if (!std::isfinite(fIn[f])) { ++nBad; }
        else { fmin = min(fmin, fIn[f]); fmax = max(fmax, fIn[f]); }
    }
    Info<< "integralSurfaceTension: face force density in [" << fmin << ", "
        << fmax << "] N/m^3";
    if (nBad > 0)
    {
        Info<< "  (" << nBad << " NON-FINITE FACE VALUES!)";
    }
    Info<< endl;

    // fIn is signed in the owner-to-neighbour face direction and is therefore
    // an oriented face quantity. Preserve that metadata through the area
    // integration so the API validator and pressure-flux algebra agree.
    tfSigma.ref().setOriented();
    return tfSigma*mesh_.magSf();
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
