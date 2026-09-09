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
    dict_(dict)
{}

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
    const bool clip,
    label& nNonFinite,
    boolList& firedCell
) const
{
    scalar v = recon.evaluate(c, foot);
    scalar lo, hi;
    recon.stencilRange(c, lo, hi);
    const scalar mid = 0.5*(lo + hi);

    if (!std::isfinite(v))
    {
        v = mid;
        ++nNonFinite;
        return v;
    }

    bool bound = clip;

    // A quasi-monotone bound cannot represent an extremum. Where psi_c is already
    // the extremum of its own stencil, lo (or hi) IS psi_c, so every reconstructed
    // value beyond it is pulled back to psi_c and the extremum is flattened -- at
    // every step, for as long as the run lasts. The level set has such extrema by
    // construction: the apex of the distance cone at the droplet centre, and the
    // box corner farthest from the interface. The fit undershoots at the apex
    // because a smooth quadratic cannot follow a non-differentiable minimum, so
    // the clip fires there on EVERY mesh, hexahedral meshes included.
    //
    // The test is exact in floating point and carries no coefficient: lo and hi
    // are the min and max over a stencil that CONTAINS psi_c, so psi_c == lo holds
    // bitwise when the cell is the minimum. A FLAT stencil (lo == hi) is not an
    // extremum -- the clip must keep enforcing exactness on a constant field.
    //
    // What stays bounded is the SPURIOUS extremum: the fit putting the value
    // outside the range of a cell that was NOT an extremum. That is the polyhedral
    // far-field defect the clip exists for.
    if (bound && recon.clipKeepExtrema() && lo != hi)
    {
        const scalar psiC = recon.stencilCellValue(c);
        if (psiC == lo || psiC == hi)
        {
            bound = false;
        }
    }

    if (bound)
    {
        // The clip is a BOUND, not a strength: it carries no coefficient and it
        // changes the value ONLY when the fit puts it outside the stencil range,
        // which is exactly when the update creates a new extremum. A cell whose
        // reconstruction stays inside its stencil bounds is bit-unchanged.
        const scalar vClipped = Foam::min(Foam::max(v, lo), hi);
        if (vClipped != v)
        {
            firedCell[c] = true;
            v = vClipped;
        }
    }
    else
    {
        const scalar cap = 10.0*Foam::max(hi - lo, SMALL);
        v = Foam::min(Foam::max(v, mid - cap), mid + cap);
    }
    return v;
}


void Foam::slCorrector::buildClipMask
(
    const volScalarField& psi,
    const slReconstruction& recon,
    boolList& clipCell
) const
{
    clipCell.setSize(mesh_.nCells(), false);

    if (!recon.clipToStencilBounds())
    {
        return;                                  // the clip is off everywhere
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
    const boolList& firedCell
) const
{
    if (!recon.clipToStencilBounds())
    {
        return;
    }

    if (firedEver_.size() != firedCell.size())
    {
        firedEver_.setSize(firedCell.size(), false);
    }

    label nEligible = 0;
    label nFired = 0;
    forAll(clipCell, c)
    {
        if (clipCell[c]) { ++nEligible; }
        if (firedCell[c]) { ++nFired; firedEver_[c] = true; }
    }
    reduce(nEligible, sumOp<label>());
    reduce(nFired, sumOp<label>());

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

    Info<< "slCorrector: quasi-monotone clip (clipRegion = "
        << recon.clipRegion() << "): " << nEligible
        << " eligible cells, " << nFired << " bounded this step, "
        << nClipFiredTotal_ << " cell-steps bounded in total over "
        << nClipStepsFired_ << " steps, first at step " << firstFireIndex_
        << "; writing slClipEligible, slClipFired and slClipFiredEver" << endl;

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
    forAll(clipCell, c)
    {
        eligible[c] = clipCell[c] ? 1.0 : 0.0;
        fired[c] = firedCell[c] ? 1.0 : 0.0;
        firedEver[c] = firedEver_[c] ? 1.0 : 0.0;
    }
    eligible.write();
    fired.write();
    firedEver.write();
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
