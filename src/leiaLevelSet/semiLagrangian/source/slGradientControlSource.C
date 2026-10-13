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
#include "slReconstruction.H"

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
    qScheme_(dict_.getOrDefault<word>("qScheme", "centred")),
    anchorInterface_(dict_.getOrDefault<bool>("anchorInterface", false)),
    strainNormal_(dict_.getOrDefault<word>("strainNormal", "psi")),
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
    if (qScheme_ != "centred" && qScheme_ != "upwindFit")
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source gradientControl qScheme must be centred"
            << " or upwindFit; got " << qScheme_ << exit(FatalIOError);
    }
    if (strainNormal_ != "psi" && strainNormal_ != "fixedX")
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source gradientControl strainNormal must be psi"
            << " or fixedX; got " << strainNormal_ << exit(FatalIOError);
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
        << law_->weight().type() << ", q " << qScheme_
        << (anchorInterface_ ? ", interface cells anchored (F = 0)" : "")
        << (strainNormal_ == "fixedX" ? ", strain direction fixed to e_x" : "")
        << ", |dt F| clamp " << maxExponent << endl;
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


void Foam::slGradientControlSource::upwindFitGradient
(
    volScalarField& psi,
    slReconstruction& fit,
    vectorField& gUp,
    boolList& valid
) const
{
    // Pass 1: the fit of every cell on the current psi; its gradient and
    // Hessian at the cell centre, exchanged across coupled patches.
    psi.correctBoundaryConditions();
    fit.update(psi);

    volVectorField gFit
    (
        IOobject
        (
            "slSourceFitGrad", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, Zero)
    );
    volSymmTensorField HFit
    (
        IOobject
        (
            "slSourceFitHess", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedSymmTensor(dimless/sqr(dimLength), Zero)
    );
    vectorField& gFitI = gFit.primitiveFieldRef();
    symmTensorField& HFitI = HFit.primitiveFieldRef();
    forAll(gFitI, celli)
    {
        fit.fitDerivatives(celli, gFitI[celli], HFitI[celli]);
    }
    gFit.correctBoundaryConditions();
    HFit.correctBoundaryConditions();

    // Pass 2: the upwind-weighted mean of the neighbours' fit gradients,
    // each evaluated at the centre of c.
    const scalarField& psiI = psi.primitiveField();
    const volVectorField& C = mesh_.C();
    vectorField acc(mesh_.nCells(), Zero);
    scalarField wsum(mesh_.nCells(), Zero);

    auto add = [&]
    (
        const label c,
        const scalar psiU,
        const vector& gU,
        const symmTensor& HU,
        const vector& xcMinusXu
    )
    {
        const scalar s = (psiI[c] > 0 ? 1.0 : (psiI[c] < 0 ? -1.0 : 0.0));
        const scalar w = s*(psiI[c] - psiU);
        if (w > 0)
        {
            acc[c] += w*(gU + (HU & xcMinusXu));
            wsum[c] += w;
        }
    };

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(nei, facei)
    {
        const label o = own[facei];
        const label u = nei[facei];
        const vector dou(C[o] - C[u]);
        add(o, psiI[u], gFitI[u], HFitI[u], dou);
        add(u, psiI[o], gFitI[o], HFitI[o], -dou);
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField psiN(psi.boundaryField()[patchi].patchNeighbourField());
        const vectorField gN(gFit.boundaryField()[patchi].patchNeighbourField());
        const symmTensorField HN(HFit.boundaryField()[patchi].patchNeighbourField());
        // delta(): from the centre of the patch cell to the neighbour centre
        const vectorField delta(p.delta());
        forAll(fc, i)
        {
            add(fc[i], psiN[i], gN[i], HN[i], -delta[i]);
        }
    }

    gUp.setSize(mesh_.nCells());
    valid.setSize(mesh_.nCells());
    forAll(gUp, celli)
    {
        valid[celli] = (wsum[celli] > 0);
        gUp[celli] = valid[celli] ? acc[celli]/wsum[celli] : vector(Zero);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slGradientControlSource::apply
(
    volScalarField& psi,
    const volVectorField& Utrajectory,
    slReconstruction& geometryFit,
    const scalar dt
)
{
    const gradientControlLaw& law = *law_;

    // qScheme centred: g = fvc::grad (the inert default, unchanged).
    // qScheme upwindFit: g from the fits of the upwind face neighbours.
    const bool upwind = (qScheme_ == "upwindFit");
    vectorField gUp;
    boolList upValid;
    tmp<volVectorField> tg;
    if (upwind)
    {
        upwindFitGradient(psi, geometryFit, gUp, upValid);
    }
    else
    {
        tg = fvc::grad(psi, "gradPsiSource");
    }

    // anchorInterface: the cells with a face neighbour of the opposite sign,
    // from the arrival field before any value changes.
    boolList anchored;
    if (anchorInterface_)
    {
        psi.correctBoundaryConditions();
        const scalarField& p0 = psi.primitiveField();
        anchored.setSize(mesh_.nCells(), false);
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        forAll(nei, facei)
        {
            if (p0[own[facei]]*p0[nei[facei]] <= 0)
            {
                anchored[own[facei]] = true;
                anchored[nei[facei]] = true;
            }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& pp = mesh_.boundary()[patchi];
            if (!pp.coupled())
            {
                continue;
            }
            const labelUList& fc = pp.faceCells();
            const scalarField pN(psi.boundaryField()[patchi].patchNeighbourField());
            forAll(fc, i)
            {
                if (p0[fc[i]]*pN[i] <= 0)
                {
                    anchored[fc[i]] = true;
                }
            }
        }
    }

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
        if (upwind && !upValid[celli])
        {
            continue;
        }
        if (anchorInterface_ && anchored[celli])
        {
            continue;
        }
        const vector& gc = (upwind ? gUp[celli] : tg()[celli]);
        const scalar q2 = magSqr(gc);
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
            const vector n =
                (strainNormal_ == "fixedX" ? vector(1, 0, 0) : gc/q);
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
