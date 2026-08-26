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
    fadeMode_(velExtDict_.getOrDefault<word>("fadeMode", "levelSet")),
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

    // ---- DILATE, ACROSS PROCESSOR BOUNDARIES, RECORDING THE LAYER INDEX.
    //
    // One dilation now serves both consumers. The anchor set (Dirichlet data)
    // is {layer <= nAnchorLayers_} and the application band is
    // {layer <= nBand}, so the two can no longer disagree about where the
    // interface region is.
    //
    // THE PARALLEL DEFECT THIS FIXES: the previous loop used mesh_.cellCells()
    // -- LOCAL topology -- with no coupled-patch exchange, so a cell one layer
    // across a processor boundary was never reached, the band was torn at every
    // boundary, and the result depended on the DECOMPOSITION. The marker is now
    // carried in a volScalarField whose coupled patches are evaluated once per
    // pass, exactly as interfaceBandMollifier and neighbourNarrowBand do. The
    // exchange has to happen once PER LAYER: layer k+1 on this rank may only be
    // reachable through a cell that layer k placed on the far side.
    //
    // EXTENT. nBand = nLayers_ + 1 reproduces the historical full-strength
    // region L = (nLayers_+1)*cellSize_, and nFade = 2*nBand reproduces the
    // outer edge of the historical taper at 2L -- the layer index simply
    // replaces |psi|/cellSize as the coordinate.
    const label nBand = nLayers_ + 1;
    const label nFade = 2*nBand;

    layer_.setSize(mesh_.nCells());
    layer_ = -1;
    forAll(alpha_, c)
    {
        if (alpha_[c] > SMALL && alpha_[c] < 1 - SMALL)
        {
            layer_[c] = 0;
        }
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    for (label pass = 1; pass <= nFade; ++pass)
    {
        volScalarField inBand
        (
            IOobject
            (
                word(), fileName(), mesh_,
                IOobject::NO_READ, IOobject::NO_WRITE, false
            ),
            mesh_,
            dimensionedScalar(dimless, 0.0),
            "zeroGradient"
        );
        forAll(layer_, cellI)
        {
            inBand[cellI] = (layer_[cellI] >= 0) ? 1.0 : 0.0;
        }
        inBand.correctBoundaryConditions();

        labelField next(layer_);

        forAll(own, faceI)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            if (layer_[o] >= 0 && next[n] < 0) next[n] = pass;
            if (layer_[n] >= 0 && next[o] < 0) next[o] = pass;
        }

        // patchNeighbourField() is the membership on the OTHER side of a
        // coupled patch -- the neighbour CELL value, which is what must be read
        // here. Without this the dilation stops at the boundary.
        forAll(inBand.boundaryField(), patchI)
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            if (!isA<coupledFvPatch>(patch)) continue;

            const scalarField neiField
            (
                inBand.boundaryField()[patchI].patchNeighbourField()
            );
            const labelUList& fc = patch.faceCells();
            forAll(neiField, i)
            {
                if (neiField[i] > 0.5 && next[fc[i]] < 0)
                {
                    next[fc[i]] = pass;
                }
            }
        }

        layer_ = next;
    }

    // Anchor set: the seed plus nAnchorLayers_ dilations, i.e. exactly the
    // cells the previous separate anchor loop produced -- now derived from the
    // same, parallel-correct, dilation.
    labelHashSet anchors;
    forAll(layer_, c)
    {
        if (layer_[c] >= 0 && layer_[c] <= nAnchorLayers_)
        {
            anchors.insert(c);
        }
    }
    seedCells_ = anchors.toc();

    // Application band.
    //
    // In levelSet mode the band must reproduce the historical one EXACTLY in
    // serial, or every existing study silently changes: that band was the
    // ANCHOR staircase (seed dilated nAnchorLayers_ times) unioned with the
    // |psi| shell below, NOT the full nBand dilation. Measured: taking
    // {layer <= nBand} instead moved band mean |d(psi)/dx| on the 1D gate from
    // 0.942131 to 0.949671 at nLayers = 3 -- small, but it is a change to
    // published behaviour with no measurement asking for it.
    //
    // In layer mode the full nBand dilation IS the band, because there is no
    // |psi| shell to union with.
    const label bandExtent = (fadeMode_ == "levelSet") ? nAnchorLayers_ : nBand;
    forAll(layer_, c)
    {
        if (layer_[c] >= 0 && layer_[c] <= bandExtent)
        {
            band_[c] = 1;
        }
    }

    if (fadeMode_ == "levelSet")
    {
        // THE |psi| BAND, still the default. Its own comment recorded that a
        // hop dilation
        // "under-covers the |psi| shell at corners", leaving raw-velocity cells
        // inside the shell. That was a genuine defect of a HOP BAND FEEDING A
        // |psi| FADE: two different metrics disagreeing about the same region.
        // With the fade on the layer index the mismatch cannot arise, because
        // band and fade are then measured in the same coordinate.
        forAll(psi_, c)
        {
            if (Foam::mag(psi_[c]) <= (nLayers_ + 1)*cellSize_[c])
            {
                band_[c] = 1;
            }
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
    // THE FADE COORDINATE. Historically this was |psi|/L with
    // L = (nLayers_+1)*cellSize_, which is a DISTANCE only while
    // |grad psi| = 1 -- the property this class exists to maintain. As psi
    // flattens the |psi| = L contour drifts outward through the mesh, so the
    // fade boundary moves and the blend reintroduces the raw velocity inside
    // whatever region is being measured.
    //
    // MEASURED, cases/1Dstretch, band mean |d(psi)/dx| at t = 1 against the
    // exact 1, MEASUREMENT window fixed at |psi| <= 3h, only the extension
    // width varied:
    //     nLayers      3          8         32        200
    //     |dpsi/dx|  0.942131  0.999111  1.000000  1.000000
    // The entire gradient error was this boundary. Widen it until it leaves the
    // measurement window and the model reproduces the rigid translation
    // exactly -- which is what it should do for a normal-constant Uext.
    //
    // The layer index is purely topological and cannot drift, so the fade sits
    // where it is meant to at the DEFAULT nLayers.
    //
    // FLAT ON THE SEED LAYER, deliberately, following interfaceBandMollifier:
    // a taper with a kink at the interface would add variation of Uext across
    // the two cells straddling psi = 0, which is exactly where it does most
    // damage to the zero contour.
    const label nBand = nLayers_ + 1;
    const label nFade = 2*nBand;

    forAll(Uext_, c)
    {
        scalar a;   // 0 at the inner edge of the fade, 1 at the outer edge

        if (fadeMode_ == "levelSet")
        {
            const scalar L = (nLayers_ + 1)*cellSize_[c];
            a = (Foam::mag(psi_[c]) - L)/Foam::max(L, SMALL);
        }
        else
        {
            const label l = layer_[c];
            // Outside the dilated region entirely: raw velocity.
            a = (l < 0) ? 1.0
              : scalar(l - nBand)/scalar(Foam::max(nFade - nBand, 1));
        }

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
