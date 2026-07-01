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
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceExtension::computeNormals()
{
    tmp<volVectorField> tg = fvc::grad(psi_);
    const volVectorField& g = tg();
    nHat_ = g/(mag(g) + dimensionedScalar("e", g.dimensions(), SMALL));
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
    seedCells_ = seed.toc();

    // Application band = interface layer dilated by nLayers cell-neighbour
    // layers. Uext is used for the advecting flux here; the extension itself is
    // solved over the whole domain.
    const labelListList& cc = mesh_.cellCells();
    labelList front(seedCells_);
    for (label layer = 0; layer < nLayers_; ++layer)
    {
        labelHashSet next;
        forAll(front, i)
        {
            const labelList& nbrs = cc[front[i]];
            forAll(nbrs, k)
            {
                const label nb = nbrs[k];
                if (band_[nb] < 0.5)
                {
                    band_[nb] = 1;
                    next.insert(nb);
                }
            }
        }
        front = next.toc();
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
    // Runaway safety net: bound the extended velocity in the band by a multiple
    // of the global peak base-velocity magnitude. Prescribed-flux verification
    // is bounded (|U| ~ O(1)); an extension that overshoots this is diverging.
    if (maxScale_ <= 0)
    {
        return;
    }
    const scalar Umax = maxScale_*max(mag(U_)).value();
    if (Umax <= SMALL)
    {
        return;
    }
    forAll(band_, c)
    {
        if (band_[c] > 0.5)
        {
            const scalar m = mag(Uext_[c]);
            if (m > Umax)
            {
                Uext_[c] *= Umax/m;
            }
        }
    }
}


void Foam::interfaceExtension::updateFlux()
{
    // Outside the band, keep the base velocity.
    forAll(Uext_, c)
    {
        if (band_[c] < 0.5)
        {
            Uext_[c] = U_[c];
        }
    }
    Uext_.correctBoundaryConditions();

    // Advection flux: base flux everywhere, Uext-derived flux on band faces.
    surfaceScalarField phiU(fvc::interpolate(Uext_) & mesh_.Sf());
    phiExt_ == phi_;

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, f)
    {
        if (band_[own[f]] > 0.5 || band_[nei[f]] > 0.5)
        {
            phiExt_[f] = phiU[f];
        }
    }

    auto& phiExtBf = phiExt_.boundaryFieldRef();
    const auto& phiUBf = phiU.boundaryField();
    forAll(mesh_.boundary(), patchI)
    {
        const labelUList& fc = mesh_.boundary()[patchI].faceCells();
        forAll(fc, i)
        {
            if (band_[fc[i]] > 0.5)
            {
                phiExtBf[patchI][i] = phiUBf[patchI][i];
            }
        }
    }

    if (projectFlux_)
    {
        correctFlux(phiExt_);
    }
}


void Foam::interfaceExtension::correct()
{
    Uext_ == U_;            // start from the current base velocity
    computeNormals();
    buildBand();
    buildSeedVelocity();    // Dirichlet data at the interface layer
    extend();               // model-specific propagation through the band
    clipExtended();         // safety cap against runaway feedback
    updateFlux();
}

// ************************************************************************* //
