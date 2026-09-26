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

#include "slGradientControlSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

namespace Foam
{
    defineTypeNameAndDebug(slGradientControlSource, 0);
    addToRunTimeSelectionTable(slSource, slGradientControlSource, Mesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slGradientControlSource::slGradientControlSource(const fvMesh& mesh)
:
    slSource(mesh),
    bandCells_(dict_.get<scalar>("bandCells")),
    law_(gradientControlLaw::New(dict_.subDict("law"))),
    h_(mesh.nCells(), GREAT),
    F_
    (
        IOobject
        (
            "slSourceF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, 0)
    ),
    clamp_
    (
        IOobject
        (
            "slSourceClamp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    nClampTotal_(0)
{
    if (!(bandCells_ > 0))
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source gradientControl needs bandCells > 0; got "
            << bandCells_ << exit(FatalIOError);
    }
    if (!law_->active())
    {
        WarningInFunction
            << "semiLagrangian.source gradientControl with law none and strain"
            << " weight none applies F = 0 everywhere; select none instead."
            << nl << endl;
    }
    computeCellSize();
    Info<< "slGradientControlSource: band |psi|/|grad psi| <= " << bandCells_
        << " h, law " << law_->type() << ", strain weight "
        << law_->weight().type() << ", |dt F| clamp " << maxExponent << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slGradientControlSource::computeCellSize()
{
    const volVectorField& C = mesh_.C();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    h_ = GREAT;
    forAll(nei, facei)
    {
        const scalar d = mag(C[nei[facei]] - C[own[facei]]);
        h_[own[facei]] = min(h_[own[facei]], d);
        h_[nei[facei]] = min(h_[nei[facei]], d);
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const vectorField Cn(C.boundaryField()[patchi].patchNeighbourField());
        forAll(fc, i)
        {
            h_[fc[i]] = min(h_[fc[i]], mag(Cn[i] - C[fc[i]]));
        }
    }
    // A cell without a face neighbour (a one-cell mesh) keeps no band.
    forAll(h_, celli)
    {
        if (h_[celli] >= GREAT)
        {
            h_[celli] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slGradientControlSource::apply
(
    volScalarField& psi,
    const volVectorField& Utrajectory,
    slReconstruction&,
    const scalar dt
)
{
    const gradientControlLaw& law = *law_;
    const volVectorField g(fvc::grad(psi, "gradPsiSource"));

    tmp<volSymmTensorField> tD;
    if (law.needsStrainRate() || law.needsNormalStrain())
    {
        tD = symm(fvc::grad(Utrajectory, "gradUSource"));
    }

    scalarField& psiI = psi.primitiveFieldRef();
    scalarField& FI = F_.primitiveFieldRef();
    scalarField& clampI = clamp_.primitiveFieldRef();
    FI = 0;
    clampI = 0;

    label nBand = 0, nClamp = 0;
    scalar maxAbsExp = 0;
    forAll(psiI, celli)
    {
        const scalar q2 = magSqr(g[celli]);
        if (!(q2 > VSMALL))
        {
            continue;
        }
        const scalar q = Foam::sqrt(q2);
        if (mag(psiI[celli]) > bandCells_*h_[celli]*q)
        {
            continue;
        }
        ++nBand;

        scalar sigma = 0, a = 0;
        if (tD.valid())
        {
            const symmTensor& D = tD()[celli];
            sigma = mag(D);
            const vector n = g[celli]/q;
            a = n & D & n;
        }
        const scalar F = law.F(q2, sigma, a);
        scalar x = dt*F;
        if (mag(x) > maxExponent)
        {
            x = sign(x)*maxExponent;
            clampI[celli] = 1;
            ++nClamp;
        }
        maxAbsExp = max(maxAbsExp, mag(x));
        FI[celli] = F;
        psiI[celli] *= Foam::exp(x);
    }

    // Processor and physical patch values of the new psi (the halo of the
    // next step's stencils).
    psi.correctBoundaryConditions();
    F_.correctBoundaryConditions();
    clamp_.correctBoundaryConditions();

    // Collective reductions on every rank, outside any master guard.
    const label nBandAll = returnReduce(nBand, sumOp<label>());
    const label nClampAll = returnReduce(nClamp, sumOp<label>());
    const scalar maxAll = returnReduce(maxAbsExp, maxOp<scalar>());
    nClampTotal_ += nClampAll;

    Info<< "slGradientControlSource: band cells " << nBandAll
        << ", clamped " << nClampAll << " (run total " << nClampTotal_
        << "), max |dt F| " << maxAll << endl;
}


// ************************************************************************* //
