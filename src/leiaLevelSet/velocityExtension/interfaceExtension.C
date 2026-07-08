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

#include "interfaceExtension.H"
#include "fluxCorrection.H"
#include "interpolationCellPoint.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceExtension::interfaceExtension(const fvMesh& mesh)
:
    velocityExtension(mesh),
    psi_(mesh.lookupObject<volScalarField>
        (velExtDict_.getOrDefault<word>("levelSet", "psi"))),
    alpha_(mesh.lookupObject<volScalarField>
        (velExtDict_.getOrDefault<word>("alpha", "alpha"))),
    nLayers_(velExtDict_.getOrDefault<label>("nLayers", 3)),
    nDescent_(velExtDict_.getOrDefault<label>("nDescent", 3)),
    nAnchorLayers_(velExtDict_.getOrDefault<label>("nAnchorLayers", 0)),
    projectFlux_(velExtDict_.getOrDefault<bool>("projectFlux", false)),
    maxScale_(velExtDict_.getOrDefault<scalar>("maxScale", 5.0)),
    band_
    (
        IOobject
        (
            "velExtBand",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    nHat_
    (
        IOobject
        (
            "velExtNHat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless, vector::zero)
    ),
    seedCells_()
{
    // Anchors must stay strictly inside the application band.
    if (nAnchorLayers_ > nLayers_ - 1)
    {
        WarningInFunction
            << "nAnchorLayers = " << nAnchorLayers_
            << " clamped to nLayers - 1 = " << nLayers_ - 1 << endl;
        nAnchorLayers_ = nLayers_ - 1;
    }
    if (nAnchorLayers_ < 0)
    {
        nAnchorLayers_ = 0;
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceExtension::computeNormals()
{
    tmp<volVectorField> tg = fvc::grad(psi_);
    const volVectorField& g = tg();
    nHat_ = g/(mag(g) + dimensionedScalar("e", g.dimensions(), SMALL));
}


void Foam::interfaceExtension::computeCellSize()
{
    if (cellSize_.size() == mesh_.nCells())
    {
        return;   // mesh is static; compute once
    }

    cellSize_.setSize(mesh_.nCells(), GREAT);

    // Minimum face-neighbour centre distance = 1/deltaCoeffs. Only internal
    // and coupled faces contribute -- empty (2D) directions drop out, so this
    // is the in-plane h for 2D meshes of arbitrary thickness (cbrt(V) is not).
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


void Foam::interfaceExtension::buildBand()
{
    band_ == dimensionedScalar("0", dimless, 0);

    // Seed = the single interface cell layer: cells cut by psi = 0, i.e. with a
    // fractional volume fraction 0 < alpha < 1. These carry the Dirichlet
    // interface velocity (the only cells pinned in the extension solve).
    labelHashSet seed;
    forAll(alpha_, c)
    {
        if (alpha_[c] > SMALL && alpha_[c] < 1 - SMALL)
        {
            band_[c] = 1;
            seed.insert(c);
        }
    }

    // Anchor layers: dilate the SEED set by nAnchorLayers cell-neighbour
    // layers. Every anchor cell is foot-pointed (buildSeedVelocity) and pinned
    // as Dirichlet data exactly like the interface layer -- otherwise the PDE
    // stencil of an interface-adjacent cell mixes solved and pinned values and
    // Uext develops kinks at the interface layer's edge.
    const labelListList& cc = mesh_.cellCells();
    {
        labelList afront(seed.toc());
        for (label layer = 0; layer < nAnchorLayers_; ++layer)
        {
            labelHashSet next;
            forAll(afront, i)
            {
                const labelList& nbrs = cc[afront[i]];
                forAll(nbrs, k)
                {
                    const label nb = nbrs[k];
                    if (!seed.found(nb))
                    {
                        band_[nb] = 1;
                        seed.insert(nb);
                        next.insert(nb);
                    }
                }
            }
            afront = next.toc();
        }
    }
    seedCells_ = seed.toc();

    // Application band: |psi|-BASED, one guard layer wider than the nominal
    // nLayers band -- band_ = {|psi_c| <= (nLayers+1)*cellSize_c} (union the
    // seed staircase). A cell-neighbour-hop dilation from the staircase seed
    // layer under-covers the |psi| shell at corners (cells with |psi| < 3h can
    // be > 3 hops away), leaving raw-velocity cells INSIDE the shell -- the
    // measured O(1) ring at the band edge. The +1 guard layer keeps the
    // gradient stencils of every |psi| <= nLayers*h cell inside the assigned
    // region (the Uext -> U fade starts beyond it, see updateFlux).
    forAll(psi_, c)
    {
        if (Foam::mag(psi_[c]) <= (nLayers_ + 1)*cellSize_[c])
        {
            band_[c] = 1;
        }
    }
}


void Foam::interfaceExtension::buildSeedVelocity()
{
    interpolationCellPoint<scalar> psiInterp(psi_);
    interpolationCellPoint<vector> Uinterp(U_);
    const auto& C = mesh_.C();

    forAll(seedCells_, i)
    {
        const label c = seedCells_[i];

        // Gradient descent from the cell centre to the interface (psi = 0)
        // along the normal, re-evaluating psi at the moving point each step
        // (Newton iteration; reduces to the foot point for a signed distance).
        //
        // psi is not guaranteed to be a signed distance (no redistancing under
        // strong stretching), so a step can overshoot and leave the cell -- or
        // the mesh. Re-localise the point each step with findCell and only ever
        // interpolate inside the cell that contains it; if a step would leave
        // the mesh, keep the last valid point. This avoids the degenerate-tet
        // SIGFPE in interpolationCellPoint::findTetrahedron.
        point x = C[c];
        label cx = c;
        for (label it = 0; it < nDescent_; ++it)
        {
            const point xNew = x - psiInterp.interpolate(x, cx)*nHat_[c];
            const label cn = mesh_.findCell(xNew);
            if (cn < 0)
            {
                break;
            }
            x = xNew;
            cx = cn;
        }

        // IDW (cell-point) interpolation of U at the interface foot point.
        Uext_[c] = Uinterp.interpolate(x, cx);
    }
}


void Foam::interfaceExtension::clipExtended()
{
    // Runaway safety net: bound the extended velocity (whole domain -- the
    // seamless flux advects psi with Uext everywhere) by a multiple of the
    // global peak base-velocity magnitude. Prescribed-flux verification is
    // bounded (|U| ~ O(1)); an extension that overshoots this is diverging.
    if (maxScale_ <= 0)
    {
        return;
    }
    const scalar Umax = maxScale_*max(mag(U_)).value();
    if (Umax <= SMALL)
    {
        return;
    }
    forAll(Uext_, c)
    {
        const scalar m = mag(Uext_[c]);
        if (m > Umax)
        {
            Uext_[c] *= Umax/m;
        }
    }
}


void Foam::interfaceExtension::updateFlux()
{
    // WHOLE-DOMAIN extension flux -- deliberately no band seam. A hard switch
    // from the extension flux to the base flux across one face is a velocity
    // discontinuity that writes irreversible |grad psi| kinks just outside the
    // band; under flow reversal they return onto Sigma (measured at T=2,
    // N=128: band L2 of ||grad psi|-1| grew ~10x through the seam).
    //
    // Smoothly FADE the extension into the base velocity outside the band:
    // w = 1 for |psi| <= L (band half-width, L = nLayers*cellSize), cosine
    // ramp to w = 0 by 2L. pseudoTime fades naturally (its march starts from
    // Uext == U and has finite reach); anisotropicDiffusion is whole-domain
    // and NEEDS the fade -- advecting far-field psi with the approximate
    // extension distorts it globally and feeds back through
    // n = grad(psi)/|grad(psi)| (measured: global gradient error 0.18 -> 1.0
    // without the fade). The blended Uext stays C^1-smooth, so the flux has
    // no discontinuity anywhere.
    forAll(Uext_, c)
    {
        // Fade starts beyond the assigned band (nLayers + 1 guard layer), so
        // the gradient stencils of all |psi| <= nLayers*h cells see only
        // extended values.
        const scalar L = (nLayers_ + 1)*cellSize_[c];
        const scalar a = (Foam::mag(psi_[c]) - L)/Foam::max(L, SMALL);
        if (a >= 1.0)
        {
            Uext_[c] = U_[c];
        }
        else if (a > 0.0)
        {
            const scalar w = 0.5*(1.0 + Foam::cos(M_PI*a));
            Uext_[c] = w*Uext_[c] + (1.0 - w)*U_[c];
        }
    }
    Uext_.correctBoundaryConditions();
    const surfaceScalarField phiU(fvc::interpolate(Uext_) & mesh_.Sf());
    phiExt_.primitiveFieldRef() = phiU.primitiveField();

    // Keep the base boundary fluxes: the domain boundary keeps its base
    // (e.g. impermeable) behaviour; the extension only matters near Sigma.
    phiExt_.boundaryFieldRef() = phi_.boundaryField();

    if (projectFlux_)
    {
        correctFlux(phiExt_);
    }
}


void Foam::interfaceExtension::correct()
{
    Uext_ == U_;            // start from the current base velocity
    computeCellSize();
    computeNormals();
    buildBand();
    buildSeedVelocity();    // Dirichlet data at the interface layer
    extend();               // model-specific propagation through the band
    clipExtended();         // safety cap against runaway feedback
    updateFlux();
}

// ************************************************************************* //
