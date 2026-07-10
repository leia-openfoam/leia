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

#include "slReconstruction.H"
#include "addToRunTimeSelectionTable.H"
#include "centredCPCCellToCellStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slReconstruction, 0);
    defineRunTimeSelectionTable(slReconstruction, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slReconstruction::slReconstruction(const fvMesh& mesh)
:
    mesh_(mesh),
    fvSolution_(mesh),
    levelSetDict_(fvSolution_.subDict("levelSet")),
    slDict_(levelSetDict_.subOrEmptyDict("semiLagrangian")),
    clipToStencilBounds_
    (
        slDict_.getOrDefault<Switch>("clipToStencilBounds", false)
    ),
    // MEASURED: Barth-Jespersen over the wide point-neighbour stencil with the
    // IDW-weighted quadratic fit over-restricts (spurious overshoot at the
    // far/diagonal neighbours) and collapses the convergence order (quadratic
    // 3.0 -> 0.1). Off by default -- the value cap in slAdvection prevents
    // overflow WITHOUT hurting order. Opt in only if oscillations demand it
    // (a nearest-neighbour or Venkatakrishnan variant would be needed to keep
    // second order).
    limitSlope_(slDict_.getOrDefault<Switch>("limitSlope", false)),
    phi_(mesh.nCells(), 1.0),
    stencil_(centredCPCCellToCellStencilObject::New(mesh)),
    stencilC_(),
    haveCentres_(false),
    radius_(),
    stencilPsi_(),
    psiOldPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::slReconstruction> Foam::slReconstruction::New
(
    const fvMesh& mesh
)
{
    const fvSolution& fvSolution(mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");
    const dictionary& slDict = levelSetDict.subOrEmptyDict("semiLagrangian");
    const word modelType =
        slDict.getOrDefault<word>("reconstruction", "quadraticWeightedLeastSquares");
    Info<< "Selecting slReconstruction " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "slReconstruction",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<slReconstruction>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::slReconstruction::collectStencil(const volScalarField& psiOld)
{
    // One paired halo exchange (matched on every rank). Gathers psiOld over
    // each cell's compact CPC stencil, including remote (processor/cyclic)
    // point neighbours; stencil_[c][0] is the arrival cell c itself.
    stencil_.collectData(psiOld, stencilPsi_);

    if (!haveCentres_)
    {
        stencil_.collectData(mesh_.C(), stencilC_);
        // Cache the per-cell interpolation radius NOW so stencilRadius() no
        // longer reads the centres; caching models can then free stencilC_
        // (releaseStencilCentres) after build() -- ~26 vectors/cell is GBs at 128^3.
        radius_.setSize(mesh_.nCells());
        forAll(stencilC_, c)
        {
            const List<vector>& C = stencilC_[c];
            scalar r = 0;
            for (label i = 1; i < C.size(); ++i)
            {
                r = Foam::max(r, Foam::mag(C[i] - C[0]));
            }
            radius_[c] = r;
        }
        haveCentres_ = true;
    }
    psiOldPtr_ = &psiOld;
}


void Foam::slReconstruction::stencilRange
(
    const label c,
    scalar& lo,
    scalar& hi
) const
{
    const List<scalar>& s = stencilPsi_[c];
    lo = hi = s[0];
    forAll(s, i)
    {
        lo = Foam::min(lo, s[i]);
        hi = Foam::max(hi, s[i]);
    }
}


Foam::scalar Foam::slReconstruction::stencilRadius(const label c) const
{
    // Cached in collectStencil() so this stays valid after the centres are freed.
    return radius_[c];
}


void Foam::slReconstruction::releaseStencilCentres()
{
    // The slope limiter (computeLimiters) still needs the centres; keep them if on.
    if (limitSlope_)
    {
        return;
    }
    stencilC_.clear();   // frees the per-cell centre lists (radius_ is retained)
}


void Foam::slReconstruction::computeLimiters()
{
    phi_.setSize(mesh_.nCells());
    if (!limitSlope_)
    {
        phi_ = 1.0;
        return;
    }

    // Barth-Jespersen: scale each cell's reconstruction increment so the
    // UNLIMITED value at every stencil-neighbour centre stays within the stencil
    // min/max of psiOld. phi in [0,1]; phi = 1 on smooth, consistent data.
    forAll(phi_, c)
    {
        const List<scalar>& s = stencilPsi_[c];
        const List<vector>& X = stencilC_[c];
        const scalar psiC = s[0];
        scalar lo = s[0], hi = s[0];
        forAll(s, i)
        {
            lo = Foam::min(lo, s[i]);
            hi = Foam::max(hi, s[i]);
        }

        scalar phi = 1.0;
        for (label i = 1; i < X.size(); ++i)
        {
            const scalar d = evaluateRaw(c, X[i]) - psiC;   // unlimited increment
            if (d > VSMALL)
            {
                phi = Foam::min(phi, Foam::min(scalar(1), (hi - psiC)/d));
            }
            else if (d < -VSMALL)
            {
                phi = Foam::min(phi, Foam::min(scalar(1), (lo - psiC)/d));
            }
        }
        phi_[c] = Foam::max(scalar(0), phi);
    }
}


Foam::scalar Foam::slReconstruction::evaluate
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];
    return psiC + phi_[c]*(evaluateRaw(c, x) - psiC);
}

// ************************************************************************* //
