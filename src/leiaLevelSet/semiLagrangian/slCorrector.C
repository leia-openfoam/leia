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

#include "slCorrector.H"
#include "slReconstruction.H"
#include "coupledFvPatch.H"
#include "Switch.H"
#include <cmath>   // std::isfinite

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slCorrector, 0);
    defineRunTimeSelectionTable(slCorrector, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slCorrector::slCorrector(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    dict_(dict),
    bound_
    (
        slValueBound::New
        (
            mesh,
            dict,
            // The legacy Switch, read from the SAME sub-dict and with the SAME
            // default that slReconstruction uses, so the sentinel resolves
            // identically on both sides and no case changes behaviour.
            dict.getOrDefault<Switch>("clipToStencilBounds", false)
        )
    )
{
    bound_->printBanner();
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::slCorrector> Foam::slCorrector::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType = dict.getOrDefault<word>("correction", "direct");
    Info<< "Selecting slCorrector " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "slCorrector",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<slCorrector>(ctorPtr(mesh, dict));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::slCorrector::robustEvaluate
(
    const slReconstruction& recon,
    const label c,
    const point& foot,
    const bool eligible,
    slBoundTally& tally
) const
{
    scalar v = recon.evaluate(c, foot);
    scalar lo, hi;
    recon.stencilRange(c, lo, hi);
    const scalar mid = 0.5*(lo + hi);

    if (!std::isfinite(v))
    {
        v = mid;
        ++tally.nNonFinite;
        return v;
    }

    // The bound decides ONE thing: the interval. It has no side effects, it may
    // not reduce, and it gets the stencil range passed in so no bound recomputes
    // it. `none` returns the runaway cap, which reproduces the pre-family path
    // exactly: min(max(v, mid - cap), mid + cap) with cap = 10*max(hi - lo, SMALL).
    slBoundResult r;
    bound_->interval(recon, c, foot, eligible, lo, hi, r);

    const scalar vb = Foam::min(Foam::max(v, r.lo), r.hi);

    if (r.enforced)
    {
        // A bound is a BOUND, not a strength: it changes the value ONLY when the
        // fit puts it outside the admissible interval. A cell whose reconstruction
        // stays inside is bit-unchanged. The runaway cap is NOT counted as a
        // firing (r.enforced is false there), so the counters keep the meaning
        // they had before the bound became selectable.
        if (vb != v)
        {
            tally.fired[c] = true;
            tally.delta[c] = vb - v;
        }
        const scalar d = Foam::max(Foam::mag(foot - mesh_.C()[c]), SMALL);
        tally.slack[c] = Foam::min(r.hi - v, v - r.lo)/d;
    }

    if (r.inadmissible)
    {
        tally.inadmissible[c] = true;
    }

    return vb;
}


void Foam::slCorrector::buildClipMask
(
    const volScalarField& psi,
    const slReconstruction& recon,
    boolList& clipCell
) const
{
    clipCell.setSize(mesh_.nCells(), false);

    if (bound_->inert())
    {
        return;                                  // no bound acts anywhere
    }

    clipCell = true;

    if (recon.clipRegion() == "all")
    {
        return;                                  // the measured GLOBAL clip
    }

    // ------------------------------------------------------------------
    // clipRegion outsideBand: withhold the clip from the SIGN-CHANGE NARROW
    // BAND, and apply it in the far field.
    //
    // WHY the band must be excluded. A quasi-monotone clip cannot preserve an
    // extremum, and the level set has legitimate extrema. The one that matters
    // for the reported metrics sits AT the interface: the fit there must be free
    // to overshoot its stencil bounds, or the interface metrics move. MEASURED
    // on Popinet's 2D hexahedral translating droplet at N = 64, the global clip
    // costs +30 % volume error and +5 % shape error. Excluding the band removes
    // that cost by construction, because the clip then cannot touch a cell that
    // the interface passes through.
    //
    // WHY the far field needs it. The polyhedral far-field defect grows a false
    // zero set from the one-sided small cells of a cfMesh boundary slab, where
    // the fit's amplification bound Lambda reaches 1.26 against 1.05 on hex
    // (slFitAmplification). The clip is the only bound that removes the growth
    // without a coefficient: Lambda = 1 means the update is a convex combination
    // of the stencil values, and the clip enforces exactly that property.
    //
    // WHY this form is admissible. The criterion is the repository's own
    // signChangeNarrowBand: a cell is in the band when psi changes sign (or
    // vanishes) across one of its faces. It is stated per cell over face
    // neighbours, so it assumes no mesh structure, no interface shape and no
    // mode; it uses a compact stencil; and it names no mesh type. The coupled
    // patches are handled below, so a rank boundary classifies its cells exactly
    // as the serial run does.
    //
    // psi with GUARANTEED-CURRENT coupled-patch values. On a processor boundary
    // the far side comes from patchNeighbourField(), which returns whatever was
    // last exchanged; deriving the criterion from a synced COPY removes the
    // dependence on whatever the caller did beforehand. Without it the band --
    // and therefore the answer -- depends on the DECOMPOSITION, silently. This
    // is the same guarantee narrowBand::psiSynced() makes, and the same class of
    // defect this repository has hit in the psi filter and the band dilation.
    volScalarField psiSync(psi);
    psiSync.correctBoundaryConditions();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, facei)
    {
        if (psiSync[own[facei]]*psiSync[nei[facei]] <= 0)
        {
            clipCell[own[facei]] = false;
            clipCell[nei[facei]] = false;
        }
    }

    const volScalarField::Boundary& psiBf = psiSync.boundaryField();
    const labelUList& faceOwner = mesh_.faceOwner();
    forAll(psiBf, patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];

        if (!isA<coupledFvPatch>(p))
        {
            continue;      // a physical patch has no cell on the far side
        }

        const tmp<scalarField> tnbr = psiBf[patchi].patchNeighbourField();
        const scalarField& nbr = tnbr();

        forAll(nbr, facei)
        {
            const label c = faceOwner[facei + p.start()];

            if (psiSync[c]*nbr[facei] <= 0)
            {
                clipCell[c] = false;
            }
        }
    }
}


void Foam::slCorrector::reportClipActivity
(
    const slReconstruction& recon,
    const boolList& clipCell,
    const slBoundTally& tally
) const
{
    if (bound_->inert())
    {
        return;
    }

    const boolList& firedCell = tally.fired;

    if (firedEver_.size() != firedCell.size())
    {
        firedEver_.setSize(firedCell.size(), false);
    }

    label nEligible = 0;
    label nFired = 0;
    label nInadmissible = 0;
    scalar maxDelta = 0;
    scalar minSlack = 0;
    forAll(clipCell, c)
    {
        if (clipCell[c]) { ++nEligible; }
        if (firedCell[c]) { ++nFired; firedEver_[c] = true; }
        if (tally.inadmissible[c]) { ++nInadmissible; }
        maxDelta = Foam::max(maxDelta, Foam::mag(tally.delta[c]));
        minSlack = Foam::min(minSlack, tally.slack[c]);
    }
    // Every reduction happens AFTER the cell loop and on every rank. A collective
    // inside a loop over cells deadlocks as soon as two ranks hold different cell
    // counts, and that has cost this campaign a full run.
    reduce(nEligible, sumOp<label>());
    reduce(nFired, sumOp<label>());
    reduce(nInadmissible, sumOp<label>());
    reduce(maxDelta, maxOp<scalar>());
    reduce(minSlack, minOp<scalar>());
    nInadmissibleTotal_ += nInadmissible;

    // Accumulate EVERY step. The write-time sample alone is not evidence: the
    // inflowOnly gate reported 0 bounded at every write while its metric CSV had
    // already left the baseline at step 20.
    nClipFiredTotal_ += nFired;
    if (nFired > 0)
    {
        ++nClipStepsFired_;
        if (firstFireIndex_ < 0)
        {
            firstFireIndex_ = mesh_.time().timeIndex();
        }
    }

    if (!mesh_.time().writeTime())
    {
        return;
    }

    // KEEP THIS LINE'S SHAPE. workflow/scripts/make_clip_region_gate_table.py
    // parses it, so the bound-independent numbers stay where they are and the
    // new ones go on their own line below.
    Info<< "slCorrector: quasi-monotone clip (clipRegion = "
        << recon.clipRegion() << "): " << nEligible
        << " eligible cells, " << nFired << " bounded this step, "
        << nClipFiredTotal_ << " cell-steps bounded in total over "
        << nClipStepsFired_ << " steps, first at step " << firstFireIndex_
        << "; writing slClipEligible, slClipFired and slClipFiredEver" << endl;

    Info<< "slCorrector: valueBound " << bound_->type()
        << ": max |delta| " << maxDelta
        << ", min slack/|d| " << minSlack
        << ", " << nInadmissible << " inadmissible this step, "
        << nInadmissibleTotal_ << " in total"
        << "; writing slBoundDelta, slBoundSlack and slBoundInadmissible" << endl;

    // The count says the clip acted. Only the FIELD says where, and where is what
    // decides whether the region rule is right: a firing inside the band is a
    // different defect from a firing in the far field, and the two need opposite
    // remedies.
    volScalarField eligible
    (
        IOobject("slClipEligible", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    volScalarField fired
    (
        IOobject("slClipFired", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    volScalarField firedEver
    (
        IOobject("slClipFiredEver", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    // The bound's own read-outs. slBoundDelta is the interface damage in the
    // units of psi, which is how the monotone clip's +30 % volume cost was
    // traced to six cells. slBoundInadmissible is the eikonal-drift measurement:
    // an empty interval means the OLD field already violates the Lipschitz
    // condition over that stencil.
    volScalarField boundDelta
    (
        IOobject("slBoundDelta", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    volScalarField boundSlack
    (
        IOobject("slBoundSlack", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    volScalarField boundInadmissible
    (
        IOobject("slBoundInadmissible", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar(dimless, 0)
    );
    forAll(clipCell, c)
    {
        eligible[c] = clipCell[c] ? 1.0 : 0.0;
        fired[c] = firedCell[c] ? 1.0 : 0.0;
        firedEver[c] = firedEver_[c] ? 1.0 : 0.0;
        boundDelta[c] = tally.delta[c];
        boundSlack[c] = tally.slack[c];
        boundInadmissible[c] = tally.inadmissible[c] ? 1.0 : 0.0;
    }
    eligible.write();
    fired.write();
    firedEver.write();
    boundDelta.write();
    boundSlack.write();
    boundInadmissible.write();
}


void Foam::slCorrector::footRadiusGuard
(
    const slReconstruction& recon,
    const pointField& feet
) const
{
    const volVectorField& C = mesh_.C();
    label nOutside = 0;
    scalar maxRatio = 0;
    forAll(feet, c)
    {
        const scalar disp = Foam::mag(feet[c] - C[c]);
        const scalar radius = recon.stencilRadius(c);
        if (radius > SMALL)
        {
            maxRatio = Foam::max(maxRatio, disp/radius);
            if (disp > radius) { ++nOutside; }
        }
    }
    reduce(nOutside, sumOp<label>());
    reduce(maxRatio, maxOp<scalar>());
    if (nOutside > 0)
    {
        WarningInFunction
            << "semi-Lagrangian foot left the point-neighbour stencil in "
            << nOutside << " cells (max |x_d - x_c|/stencilRadius = "
            << maxRatio << "); CFL likely > 1 -- reduce maxCo." << endl;
    }
}


void Foam::slCorrector::warnNonFinite(label nNonFinite) const
{
    reduce(nNonFinite, sumOp<label>());
    if (nNonFinite > 0)
    {
        WarningInFunction
            << "reconstruction produced a non-finite value in " << nNonFinite
            << " cells (reconstruction unstable at this CFL / resolution);"
            << " reset to the stencil mid-range." << endl;
    }
}

// ************************************************************************* //
