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

#include "projectionBlend.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "rescaledDistance.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineTypeNameAndDebug(projectionBlend, 0);
    addToRunTimeSelectionTable(velocityExtension, projectionBlend, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::projectionBlend::projectionBlend(const fvMesh& mesh)
:
    velocityExtension(mesh),
    psi_
    (
        mesh.lookupObject<volScalarField>
        (
            velExtDict_.getOrDefault<word>("levelSet", "psi")
        )
    ),
    order_(velExtDict_.getOrDefault<word>("distance", "second") == "first" ? 1 : 2),
    w1_(velExtDict_.getOrDefault<scalar>("w1", 1.5)),
    w2_(velExtDict_.getOrDefault<scalar>("w2", 4.0)),
    gradScheme_(velExtDict_.getOrDefault<word>("gradScheme", "rescaleGradPsi")),
    hessianScheme_(velExtDict_.getOrDefault<word>("hessianScheme", "rescaleHessian")),
    gradUScheme_(velExtDict_.getOrDefault<word>("gradUScheme", "gradUSdpls")),
    gradFloor_(velExtDict_.getOrDefault<scalar>("gradFloor", 1e-3)),
    h_(mesh.nCells(), 0)
{
    const word dist(velExtDict_.getOrDefault<word>("distance", "second"));
    if (dist != "first" && dist != "second")
    {
        FatalIOErrorInFunction(velExtDict_)
            << "projectionBlend: distance must be 'first' or 'second', got '"
            << dist << "'" << exit(FatalIOError);
    }
    if (!(w2_ > w1_ && w1_ >= 0))
    {
        FatalIOErrorInFunction(velExtDict_)
            << "projectionBlend: need 0 <= w1 < w2, got w1 = " << w1_
            << ", w2 = " << w2_ << exit(FatalIOError);
    }

    // Cell size from the volume over the geometric (non-empty) directions.
    const Vector<label>& geomD = mesh.geometricD();
    const vector span = mesh.bounds().span();
    scalar emptyExtent = 1;
    label nD = 0;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (geomD[d] == 1)
        {
            ++nD;
        }
        else
        {
            emptyExtent *= span[d];
        }
    }
    forAll(h_, celli)
    {
        h_[celli] = Foam::pow(mesh.V()[celli]/emptyExtent, 1.0/nD);
    }

    Info<< "projectionBlend velocity extension: distance " << dist
        << ", w1 " << w1_ << " h, w2 " << w2_ << " h, grad " << gradScheme_
        << ", Hessian " << hessianScheme_ << ", grad U " << gradUScheme_
        << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::projectionBlend::correct()
{
    const fvMesh& mesh = mesh_;

    volVectorField gradPsi(fvc::grad(psi_, gradScheme_));
    gradPsi.correctBoundaryConditions();

    tmp<volTensorField> tH;
    if (order_ == 2)
    {
        tH = fvc::grad(gradPsi, hessianScheme_);
    }

    const volTensorField gradU(fvc::grad(U_, gradUScheme_));

    vectorField& ue = Uext_.primitiveFieldRef();

    label nProjected = 0;
    forAll(ue, celli)
    {
        const vector& u = U_[celli];
        ue[celli] = u;

        const vector& g = gradPsi[celli];
        const scalar q = mag(g);
        const scalar dnn =
            (order_ == 2 && q > 0) ? ((g/q) & (tH()[celli] & (g/q))) : scalar(0);
        scalar dist = 0;
        bool negDisc = false;
        if
        (
            !rescaledDistance::estimate
            (
                order_, psi_[celli], q, dnn, gradFloor_, false, 0, dist, negDisc
            )
        )
        {
            continue;
        }

        const scalar a = mag(dist)/h_[celli];
        scalar delta = 1;
        if (a <= w1_)
        {
            delta = 0;
        }
        else if (a < w2_)
        {
            const scalar s = (a - w1_)/(w2_ - w1_);
            delta = 0.5*(1 - Foam::cos(constant::mathematical::pi*s));
        }
        if (delta >= 1)
        {
            continue;
        }

        // U at the projected point: a first-order Taylor step from x_c.
        const vector dx = -dist*(g/q);
        const vector uP = u + (dx & gradU[celli]);
        ue[celli] = delta*u + (1 - delta)*uP;
        ++nProjected;

    }
    Uext_.correctBoundaryConditions();

    // Flux: interpolated on internal and coupled faces, the prescribed flux
    // on physical boundary faces.
    phiExt_ = fvc::interpolate(Uext_) & mesh.Sf();
    surfaceScalarField::Boundary& pb = phiExt_.boundaryFieldRef();
    forAll(pb, patchi)
    {
        if (!pb[patchi].coupled())
        {
            pb[patchi] = phi_.boundaryField()[patchi];
        }
    }

    if (debug)
    {
        Info<< "projectionBlend: " << returnReduce(nProjected, sumOp<label>())
            << " cells projected" << endl;
    }
}


// ************************************************************************* //
