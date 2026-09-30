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

#include "haloLimited.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "coupledPolyPatch.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "syncTools.H"

namespace Foam
{
namespace velocityExtensions
{
    defineTypeNameAndDebug(haloLimited, 0);
    addToRunTimeSelectionTable(velocityExtension, haloLimited, Mesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityExtensions::haloLimited::haloLimited(const fvMesh& mesh)
:
    velocityExtension(mesh),
    psi_
    (
        mesh.lookupObject<volScalarField>
        (
            velExtDict_.getOrDefault<word>("levelSet", "psi")
        )
    ),
    radiusCells_(velExtDict_.get<scalar>("radiusCells")),
    travel_(extensionTravel::New(velExtDict_.subDict("travel"))),
    weight_(extensionWeight::New(velExtDict_.subDict("weight"))),
    direction_(extensionDirection::New(velExtDict_.subDict("direction"))),
    sampler_(extensionSampler::New(mesh, velExtDict_.subDict("sampler"))),
    cellSize_
    (
        IOobject
        (
            "hlCellSize",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimLength, 0),
        calculatedFvPatchScalarField::typeName
    ),
    hlWeight_
    (
        IOobject("hlWeight", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    hlDistance_
    (
        IOobject("hlDistance", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimLength, 0)
    ),
    hlReach_
    (
        IOobject("hlReach", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    hlCorrection_
    (
        IOobject("hlCorrection", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(phi_.dimensions(), 0)
    )
{
    if (!(radiusCells_ > 0) || radiusCells_ > 1)
    {
        FatalIOErrorInFunction(velExtDict_)
            << "haloLimited needs 0 < radiusCells <= 1 (the sampler evaluates"
            << " the local model of the face's own cells); got " << radiusCells_
            << exit(FatalIOError);
    }
    computeCellSize();
    Info<< "haloLimited: radius " << radiusCells_ << " h, travel "
        << travel_->type() << ", weight " << weight_->type() << ", direction "
        << direction_->type() << ", sampler " << sampler_->type() << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::velocityExtensions::haloLimited::computeCellSize()
{
    // The smallest face-neighbour centre distance (1/deltaCoeffs); internal and
    // coupled faces only, so empty (2D) directions drop out (interfaceExtension).
    scalarField& h = cellSize_.primitiveFieldRef();
    h = GREAT;
    const surfaceScalarField& dc = mesh_.deltaCoeffs();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, f)
    {
        const scalar d = 1.0/dc[f];
        h[own[f]] = min(h[own[f]], d);
        h[nei[f]] = min(h[nei[f]], d);
    }
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled())
        {
            const labelUList& fc = mesh_.boundary()[patchi].faceCells();
            const fvsPatchScalarField& pdc = dc.boundaryField()[patchi];
            forAll(fc, i)
            {
                h[fc[i]] = min(h[fc[i]], 1.0/pdc[i]);
            }
        }
    }
    forAll(h, celli)
    {
        if (h[celli] >= GREAT)
        {
            h[celli] = 0;   // no face neighbour: no sampling
        }
    }
    // Coupled patches: the neighbour cell's h, for the face size.
    cellSize_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityExtensions::haloLimited::correct()
{
    const fvMesh& mesh = mesh_;
    const label nInt = mesh.nInternalFaces();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& C = mesh.C().primitiveField();
    const vectorField& Cf = mesh.Cf().primitiveField();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& h = cellSize_.primitiveField();

    const volVectorField g(fvc::grad(psi_, "gradPsiExtension"));
    const surfaceScalarField psif(linearInterpolate(psi_));
    const surfaceVectorField gf(linearInterpolate(g));

    // The sample point Y, the weight w, the distance d and the reach |Y - x|/R.
    auto sample = [&]
    (
        const scalar psiv, const vector& gv, const scalar R, const point& x,
        point& Y, scalar& w, scalar& d, scalar& reach
    )
    {
        vector e;
        direction_->evaluate(psiv, gv, d, e);
        if (!(R > 0))
        {
            Y = x; w = 0; reach = 0;
            return;
        }
        const scalar S = travel_->S(d, R);
        w = weight_->w(travel_->fraction(d, R));
        Y = x - S*e;
        reach = mag(S*e)/R;
    };

    // ---- geometry: cells, internal faces, coupled boundary faces ----------
    const label nCells = mesh.nCells();
    pointField Yc(nCells);
    scalarField wc(nCells), dc(nCells), rc(nCells);
    forAll(Yc, c)
    {
        sample(psi_[c], g[c], radiusCells_*h[c], C[c], Yc[c], wc[c], dc[c], rc[c]);
    }

    pointField Yf(nInt);
    scalarField wf(nInt);
    scalar maxReach = gMax(rc);
    for (label f = 0; f < nInt; ++f)
    {
        scalar d, r;
        const scalar hf = min(h[own[f]], h[nei[f]]);
        sample(psif[f], gf[f], radiusCells_*hf, Cf[f], Yf[f], wf[f], d, r);
        maxReach = max(maxReach, r);
    }

    const label nPatches = mesh.boundary().size();
    List<pointField> Yb(nPatches);
    List<scalarField> wb(nPatches);
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField& hN = cellSize_.boundaryField()[patchi];
        const scalarField& pp = psif.boundaryField()[patchi];
        const vectorField& gp = gf.boundaryField()[patchi];
        const vectorField& Cp = mesh.Cf().boundaryField()[patchi];
        Yb[patchi].setSize(p.size());
        wb[patchi].setSize(p.size());
        forAll(fc, i)
        {
            scalar d, r;
            const scalar hf = min(h[fc[i]], hN[i]);
            sample(pp[i], gp[i], radiusCells_*hf, Cp[i], Yb[patchi][i], wb[patchi][i], d, r);
            maxReach = max(maxReach, r);
        }
    }
    reduce(maxReach, maxOp<scalar>());

    // ---- sampling, one velocity component at a time ------------------------
    // dUf: u_h(Y_f) - u_h(x_f) (mean of the owner's and the neighbour's model);
    // dUb: the same with the local cell's model only (coupled faces);
    // UYc: u_h(Y_c) with the cell's own model.
    vectorField dUf(nInt, Zero);
    List<vectorField> dUb(nPatches);
    forAll(Yb, patchi)
    {
        dUb[patchi].setSize(Yb[patchi].size(), Zero);
    }
    vectorField UYc(U_.primitiveField());

    const Vector<label>& solD = mesh.solutionD();
    for (direction k = 0; k < vector::nComponents; ++k)
    {
        if (solD[k] < 0)
        {
            continue;   // empty direction: the component stays that of U
        }
        const volScalarField Uk(U_.component(k));
        sampler_->update(Uk);

        for (label f = 0; f < nInt; ++f)
        {
            const label P = own[f], N = nei[f];
            dUf[f][k] =
                0.5
               *(
                    (sampler_->value(P, Yf[f]) - sampler_->value(P, Cf[f]))
                  + (sampler_->value(N, Yf[f]) - sampler_->value(N, Cf[f]))
                );
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!p.coupled())
            {
                continue;
            }
            const labelUList& fc = p.faceCells();
            const vectorField& Cp = mesh.Cf().boundaryField()[patchi];
            forAll(fc, i)
            {
                dUb[patchi][i][k] =
                    sampler_->value(fc[i], Yb[patchi][i])
                  - sampler_->value(fc[i], Cp[i]);
            }
        }
        forAll(UYc, c)
        {
            UYc[c][k] = sampler_->value(c, Yc[c]);
        }
    }

    // ---- the flux correction -------------------------------------------------
    surfaceScalarField& corr = hlCorrection_;
    scalarField& ci = corr.primitiveFieldRef();
    for (label f = 0; f < nInt; ++f)
    {
        ci[f] = wf[f]*(dUf[f] & Sf[f]);
    }

    // Coupled faces: each side's correction in the owner side's orientation,
    // summed across the coupling and halved, so that the two sides carry
    // exactly opposite values.
    scalarField bcorr(mesh.nBoundaryFaces(), Zero);
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (!p.coupled())
        {
            continue;
        }
        const scalar s =
            refCast<const coupledPolyPatch>(p.patch()).owner() ? 1 : -1;
        const vectorField& Sp = mesh.Sf().boundaryField()[patchi];
        const label b0 = p.start() - nInt;
        forAll(Sp, i)
        {
            bcorr[b0 + i] = s*wb[patchi][i]*(dUb[patchi][i] & Sp[i]);
        }
    }
    syncTools::syncBoundaryFaceList(mesh, bcorr, plusEqOp<scalar>());

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        fvsPatchScalarField& cp = corr.boundaryFieldRef()[patchi];
        cp = 0;   // physical patches keep phi
        if (!p.coupled())
        {
            continue;
        }
        const scalar s =
            refCast<const coupledPolyPatch>(p.patch()).owner() ? 1 : -1;
        const label b0 = p.start() - nInt;
        forAll(cp, i)
        {
            cp[i] = s*0.5*bcorr[b0 + i];
        }
    }

    phiExt_ == phi_ + corr;

    // ---- the cell velocity -----------------------------------------------------
    Uext_ == U_;
    vectorField& Ue = Uext_.primitiveFieldRef();
    forAll(Ue, c)
    {
        Ue[c] = U_[c] + wc[c]*(UYc[c] - U_[c]);   // = (1 - w) U + w u_h(Y)
    }
    Uext_.correctBoundaryConditions();
    // Physical patches carry U (zeroGradient would copy the extended cell value).
    volVectorField::Boundary& Ub = Uext_.boundaryFieldRef();
    forAll(Ub, patchi)
    {
        if (!Ub[patchi].coupled())
        {
            Ub[patchi] == U_.boundaryField()[patchi];
        }
    }

    // ---- diagnostics (collective, on every rank) --------------------------------
    hlWeight_.primitiveFieldRef() = wc;
    hlDistance_.primitiveFieldRef() = dc;
    hlReach_.primitiveFieldRef() = rc;
    hlWeight_.correctBoundaryConditions();
    hlDistance_.correctBoundaryConditions();
    hlReach_.correctBoundaryConditions();

    const scalarField divCorr(fvc::div(corr)().primitiveField());
    const scalarField& V = mesh.V().field();
    const scalar l2 =
        Foam::sqrt(gSum(V*sqr(divCorr))/max(gSum(V), VSMALL));

    Info<< "haloLimited: max |Y - x|/R " << maxReach
        << ", L2 of div(phiExt - phi) " << l2 << endl;
}


// ************************************************************************* //
