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
#include "centredCFCCellToCellStencilObject.H"

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
    // Backward compatible: an explicit slopeLimiter wins; otherwise the legacy
    // limitSlope Switch maps true->barthJespersen, false->none.
    limiterType_
    (
        slDict_.getOrDefault<word>
        (
            "slopeLimiter",
            limitSlope_ ? word("barthJespersen") : word("none")
        )
    ),
    venkK_(slDict_.getOrDefault<scalar>("venkK", 5.0)),
    phi_(mesh.nCells(), 1.0),
    // Departure-foot reconstruction stencil. Default cell-POINT-cell (CPC): a hex
    // cell has only 6 face-neighbours < the 9 needed for a 3D quadratic value fit,
    // so structured meshes require the wider point-neighbour stencil. General
    // polyhedra (~12-16 faces) already over-determine the quadratic with the
    // compact cell-FACE-cell (CFC) stencil, which is ~3-4x smaller -> much cheaper
    // per-cell fit and better conditioned (no sprawl into the coplanar far field).
    // Selectable per run: levelSet { semiLagrangian { stencil point | face; } }.
    stencil_
    (
        slDict_.getOrDefault<word>("stencil", "point") == "face"
      ? static_cast<const extendedCentredCellToCellStencil&>
            (centredCFCCellToCellStencilObject::New(mesh))
      : static_cast<const extendedCentredCellToCellStencil&>
            (centredCPCCellToCellStencilObject::New(mesh))
    ),
    nLocal_(mesh.nCells()),
    centreTail_(),
    haveCentres_(false),
    radius_(),
    stencilPsi_(),
    psiOldPtr_(nullptr)
{
    if (limiterType_ != "none"
     && limiterType_ != "barthJespersen"
     && limiterType_ != "venkatakrishnan")
    {
        FatalIOErrorInFunction(slDict_)
            << "slopeLimiter must be none, barthJespersen or venkatakrishnan, "
            << "got '" << limiterType_ << "'" << exit(FatalIOError);
    }

    Info<< "slReconstruction: departure-foot stencil = "
        << slDict_.getOrDefault<word>("stencil", "point")
        << " (point = cell-point-cell, face = cell-face-cell)"
        << ", slopeLimiter = " << limiterType_;
    if (limiterType_ == "venkatakrishnan")
    {
        Info<< " (venkK = " << venkK_ << ")";
    }
    Info<< endl;
}

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

    return New(mesh, modelType);
}


Foam::autoPtr<Foam::slReconstruction> Foam::slReconstruction::New
(
    const fvMesh& mesh,
    const word& modelType
)
{
    Info<< "Selecting slReconstruction " << modelType << endl;

    auto* ctorPtr = MeshConstructorTable(modelType);

    if (!ctorPtr)
    {
        const fvSolution& fvSolution(mesh);

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
        // Build the boundary/halo centre tail ONCE; local centres are read from
        // mesh_.C() on demand via stencilC(), so they are not duplicated per
        // stencil entry. Cache the per-cell interpolation radius via the accessor.
        buildCentreTail();
        radius_.setSize(mesh_.nCells());
        forAll(radius_, c)
        {
            const label n = stencilSize(c);
            const point x0 = stencilC(c, 0);
            scalar r = 0;
            for (label i = 1; i < n; ++i)
            {
                r = Foam::max(r, Foam::mag(stencilC(c, i) - x0));
            }
            radius_[c] = r;
        }
        haveCentres_ = true;
    }
    psiOldPtr_ = &psiOld;
}


void Foam::slReconstruction::buildCentreTail()
{
    // Replicate extendedCellToFaceStencil::collectData's compact layout for the
    // cell centres, then keep only the [nLocal_ ..) tail (boundary faces +
    // processor halo). Internal-cell entries are served from mesh_.C() directly,
    // so stencilC(c,i) == the former collectData(mesh_.C())[c][i] for every entry.
    const mapDistribute& map = stencil_.map();
    const volVectorField& C = mesh_.C();

    List<vector> flat(map.constructSize(), Zero);
    forAll(C, celli)
    {
        flat[celli] = C[celli];                    // internal cells [0 .. nCells)
    }
    forAll(C.boundaryField(), patchi)              // boundary-face centres
    {
        const fvPatchField<vector>& pf = C.boundaryField()[patchi];
        label nc = pf.patch().start() - mesh_.nInternalFaces() + mesh_.nCells();
        forAll(pf, i)
        {
            flat[nc++] = pf[i];
        }
    }
    map.distribute(flat);                          // fill processor-halo entries

    nLocal_ = mesh_.nCells();
    centreTail_.setSize(flat.size() - nLocal_);
    forAll(centreTail_, i)
    {
        centreTail_[i] = flat[nLocal_ + i];
    }
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


void Foam::slReconstruction::computeLimiters()
{
    phi_.setSize(mesh_.nCells());
    if (limiterType_ == "none")
    {
        phi_ = 1.0;
        return;
    }

    const bool venk = (limiterType_ == "venkatakrishnan");

    // Both variants scale each cell's reconstruction increment so the UNLIMITED
    // value at every stencil-neighbour centre stays within the stencil min/max
    // of psiOld; phi in [0,1], phi = 1 => unlimited.
    //  - barthJespersen: hard cap phi_i = min(1, Delta/d).
    //  - venkatakrishnan: SMOOTH cap that -> 1 in smooth regions (epsilon^2 =
    //    (venkK*h)^3 with h = the cell's interpolation radius), so it does not
    //    collapse the convergence order the way the hard cap does.
    forAll(phi_, c)
    {
        const List<scalar>& s = stencilPsi_[c];
        const scalar psiC = s[0];
        scalar lo = s[0], hi = s[0];
        forAll(s, i)
        {
            lo = Foam::min(lo, s[i]);
            hi = Foam::max(hi, s[i]);
        }

        const scalar eps2 =
            venk ? Foam::pow3(venkK_*radius_[c]) : 0.0;

        scalar phi = 1.0;
        const label n = stencilSize(c);
        for (label i = 1; i < n; ++i)
        {
            const scalar d = evaluateRaw(c, stencilC(c, i)) - psiC;
            scalar phiI = 1.0;
            if (d > VSMALL)
            {
                const scalar D = hi - psiC;                 // allowed + increment
                phiI = venk
                  ? ((D*D + eps2)*d + 2.0*d*d*D)/(D*D + 2.0*d*d + D*d + eps2)/d
                  : Foam::min(scalar(1), D/d);
            }
            else if (d < -VSMALL)
            {
                const scalar D = lo - psiC;                 // allowed - increment
                phiI = venk
                  ? ((D*D + eps2)*d + 2.0*d*d*D)/(D*D + 2.0*d*d + D*d + eps2)/d
                  : Foam::min(scalar(1), D/d);
            }
            phi = Foam::min(phi, phiI);
        }
        phi_[c] = Foam::min(scalar(1), Foam::max(scalar(0), phi));
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
