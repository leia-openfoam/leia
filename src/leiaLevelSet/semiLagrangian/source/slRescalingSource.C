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

#include "slRescalingSource.H"
#include "addToRunTimeSelectionTable.H"
#include "slReconstruction.H"
#include "rescaledDistance.H"

namespace Foam
{
    defineTypeNameAndDebug(slRescalingSource, 0);
    addToRunTimeSelectionTable(slSource, slRescalingSource, Mesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slRescalingSource::slRescalingSource(const fvMesh& mesh)
:
    slSource(mesh),
    order_(2),
    trigger_(dict_.getOrDefault<word>("trigger", "gradientDeviation")),
    theta_(dict_.getOrDefault<scalar>("theta", 0.1)),
    maxEvents_(dict_.getOrDefault<label>("maxEvents", -1)),
    bandLayers_(dict_.getOrDefault<label>("bandLayers", 3)),
    gradFloor_(dict_.getOrDefault<scalar>("gradFloor", 1e-3)),
    nEvents_(0),
    nCalls_(0),
    nRescaledTotal_(0),
    nSkippedTotal_(0)
{
    const word orderName(dict_.getOrDefault<word>("order", "second"));
    if (orderName == "first")
    {
        order_ = 1;
    }
    else if (orderName == "second")
    {
        order_ = 2;
    }
    else
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source rescaling order must be first or"
            << " second; got " << orderName << exit(FatalIOError);
    }
    if (trigger_ != "gradientDeviation" && trigger_ != "everyStep")
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source rescaling trigger must be"
            << " gradientDeviation or everyStep; got " << trigger_
            << exit(FatalIOError);
    }
    if (!(theta_ > 0))
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source rescaling needs theta > 0; got "
            << theta_ << exit(FatalIOError);
    }
    if (bandLayers_ < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source rescaling needs bandLayers >= 0; got "
            << bandLayers_ << exit(FatalIOError);
    }
    if (!(gradFloor_ > 0))
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source rescaling needs gradFloor > 0; got "
            << gradFloor_ << exit(FatalIOError);
    }

    Info<< "slRescalingSource: order " << orderName << ", trigger "
        << trigger_ << ", theta " << theta_ << ", maxEvents " << maxEvents_
        << ", region NarrowBand + " << bandLayers_ << " layers, gradFloor "
        << gradFloor_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::slRescalingSource::~slRescalingSource()
{
    // The totals are already reduced over the ranks, and Info prints on the
    // master only. The destructor makes no collective call.
    Info<< "slRescalingSource: total " << nEvents_ << " events in "
        << nCalls_ << " calls, " << nRescaledTotal_
        << " cell values rescaled, " << nSkippedTotal_ << " skipped" << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slRescalingSource::regionMask(boolList& mask) const
{
    if (!mesh_.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "semiLagrangian.source rescaling needs the registered field"
            << " NarrowBand: set levelSet.narrowBand.type to signChange in"
            << " fvSolution." << exit(FatalError);
    }
    const volScalarField& nb =
        mesh_.lookupObject<volScalarField>("NarrowBand");

    // A marker field, so that each layer reads the neighbour values across
    // the coupled patches. The physical patches are zeroGradient. The
    // pattern is redistancer::topologicalMask (libleiaRedistancer, which
    // this library does not link).
    volScalarField m
    (
        IOobject
        (
            "slRescalingRegion",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        word("zeroGradient")
    );
    scalarField& mI = m.primitiveFieldRef();
    const scalarField& nbI = nb.primitiveField();
    forAll(mI, celli)
    {
        mI[celli] = (nbI[celli] > 0.5) ? 1.0 : 0.0;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    for (label layer = 0; layer < bandLayers_; ++layer)
    {
        // The exchange of the marker of the previous layer.
        m.correctBoundaryConditions();

        // The new cells are collected first and assigned after the loops:
        // one pass adds exactly one layer on every rank.
        boolList add(mesh_.nCells(), false);
        forAll(nei, facei)
        {
            if (mI[own[facei]] > 0.5 || mI[nei[facei]] > 0.5)
            {
                add[own[facei]] = true;
                add[nei[facei]] = true;
            }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& p = mesh_.boundary()[patchi];
            if (!p.coupled())
            {
                continue;
            }
            const scalarField nbr
            (
                m.boundaryField()[patchi].patchNeighbourField()
            );
            const labelUList& fc = p.faceCells();
            forAll(fc, i)
            {
                if (nbr[i] > 0.5)
                {
                    add[fc[i]] = true;
                }
            }
        }
        forAll(add, celli)
        {
            if (add[celli])
            {
                mI[celli] = 1.0;
            }
        }
    }

    mask.setSize(mesh_.nCells());
    forAll(mask, celli)
    {
        mask[celli] = (mI[celli] > 0.5);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slRescalingSource::apply
(
    volScalarField& psi,
    const volVectorField&,
    slReconstruction& geometryFit,
    const scalar
)
{
    ++nCalls_;

    // The event limit. nEvents_ is equal on every rank, so every rank
    // returns here in the same call.
    if (maxEvents_ >= 0 && nEvents_ >= maxEvents_)
    {
        return;
    }

    // The fit of the arrival field. The coupled-patch values must be
    // current before the fit reads its stencils.
    psi.correctBoundaryConditions();
    geometryFit.update(psi);

    boolList inRegion;
    regionMask(inRegion);

    // Pass 1: the new values from the arrival field. No value is assigned.
    scalarField& psiI = psi.primitiveFieldRef();
    scalarField psiNew(psiI);
    boolList changed(mesh_.nCells(), false);

    label nRegion = 0, nRescaled = 0, nSkipped = 0, nNegDisc = 0;
    scalar maxDev = 0;
    vector g(Zero);
    symmTensor H(Zero);
    forAll(psiI, celli)
    {
        if (!inRegion[celli])
        {
            continue;
        }
        ++nRegion;

        const label fitOrder = geometryFit.fitDerivatives(celli, g, H);
        if (fitOrder < 1)
        {
            // The constant fallback: the fit has no gradient.
            ++nSkipped;
            continue;
        }
        const scalar q = mag(g);
        maxDev = max(maxDev, mag(q - 1));

        // The linear fallback has H = 0: the first-order map.
        const label cellOrder = (fitOrder >= 2 ? order_ : 1);
        scalar dnn = 0;
        if (cellOrder == 2 && q >= gradFloor_)
        {
            const vector n(g/q);
            dnn = n & H & n;
        }

        scalar dist = psiI[celli];
        bool negDisc = false;
        if
        (
            !rescaledDistance::estimate
            (
                cellOrder, psiI[celli], q, dnn, gradFloor_, false, 0,
                dist, negDisc
            )
        )
        {
            // |g| < gradFloor: the cell keeps its value.
            ++nSkipped;
            continue;
        }
        if (negDisc)
        {
            ++nNegDisc;
        }
        psiNew[celli] = dist;
        changed[celli] = true;
        ++nRescaled;
    }

    // The trigger. Collective reductions on every rank, outside any master
    // guard. The decision is then equal on every rank.
    reduce(maxDev, maxOp<scalar>());
    reduce(nRegion, sumOp<label>());
    reduce(nRescaled, sumOp<label>());
    reduce(nSkipped, sumOp<label>());
    reduce(nNegDisc, sumOp<label>());

    const bool fire = (trigger_ == "everyStep") || (maxDev > theta_);
    if (!fire)
    {
        return;
    }

    // Pass 2: the assignment, then the patch values of the new psi.
    forAll(psiI, celli)
    {
        if (changed[celli])
        {
            psiI[celli] = psiNew[celli];
        }
    }
    psi.correctBoundaryConditions();

    ++nEvents_;
    nRescaledTotal_ += nRescaled;
    nSkippedTotal_ += nSkipped;

    Info<< "slRescalingSource: event " << nEvents_ << " at time index "
        << mesh_.time().timeIndex() << ", max||g| - 1| " << maxDev
        << ", region cells " << nRegion << ", rescaled " << nRescaled
        << ", skipped " << nSkipped << ", negative discriminant "
        << nNegDisc << endl;
}


// ************************************************************************* //
