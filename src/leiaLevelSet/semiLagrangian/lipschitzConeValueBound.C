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

#include "lipschitzConeValueBound.H"
#include "slReconstruction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lipschitzConeValueBound, 0);
    addToRunTimeSelectionTable(slValueBound, lipschitzConeValueBound, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lipschitzConeValueBound::lipschitzConeValueBound
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    slValueBound(mesh, dict),
    lipschitzMode_(dict.getOrDefault<word>("lipschitzMode", "unity")),
    L0_(dict.getOrDefault<scalar>("lipschitzConstant", 1.0)),
    onInadmissible_(dict.getOrDefault<word>("onInadmissible", "cellOnly")),
    useBoundaryFaces_(dict.getOrDefault<Switch>("coneBoundaryFaces", false)),
    measureL_(lipschitzMode_ == "stencil"),
    recoverCellOnly_(onInadmissible_ == "cellOnly")
{
    if (lipschitzMode_ != "unity" && lipschitzMode_ != "stencil")
    {
        FatalIOErrorInFunction(dict)
            << "lipschitzMode must be unity or stencil, got '"
            << lipschitzMode_ << "'" << exit(FatalIOError);
    }

    if (onInadmissible_ != "cellOnly" && onInadmissible_ != "none")
    {
        FatalIOErrorInFunction(dict)
            << "onInadmissible must be cellOnly or none, got '"
            << onInadmissible_ << "'" << exit(FatalIOError);
    }

    // A NEGATIVE or zero L makes the interval collapse onto the stencil values
    // and every cell inadmissible. That is not a configuration, it is a typo.
    if (L0_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "lipschitzConstant must be positive, got " << L0_
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::lipschitzConeValueBound::interval
(
    const slReconstruction& recon,
    const label c,
    const point& foot,
    const bool eligible,
    const scalar stencilLo,
    const scalar stencilHi,
    slBoundResult& r
) const
{
    if (!eligible)
    {
        runawayCap(stencilLo, stencilHi, r);
        return;
    }

    // The eikonal value 1 is the floor, never the ceiling. For an exact distance
    // field the measured owner-to-neighbour slope is at most 1, and it is much
    // less than 1 in a direction across the gradient, so taking the measured
    // value ALONE would give a bound tighter than the truth and would clip
    // correct values. max() trusts the eikonal value and widens only where the
    // data says the field is steeper -- which is exactly the drift this mode is
    // meant to expose.
    const scalar L =
        measureL_
      ? Foam::max(L0_, recon.stencilLipschitz(c))
      : L0_;

    scalar lo, hi;
    const bool ok =
        recon.stencilConeRange(c, foot, L, useBoundaryFaces_, lo, hi);

    if (ok)
    {
        r.lo = lo;
        r.hi = hi;
        r.enforced = true;
        return;
    }

    // EMPTY interval. The old data is not L-Lipschitz over this stencil, so the
    // bound has no mathematical meaning here. Widening L until it fits, or
    // swapping lo and hi, would destroy that meaning; the recovery path is
    // documented instead, and the cell is counted so the crossing fraction stays
    // a measurement of the eikonal drift.
    r.inadmissible = true;

    if (recoverCellOnly_)
    {
        // The owner-only interval. It is never empty, and it keeps the one bound
        // that is proved: |psi_c^{n+1} - psi_c^n| <= L|delta_c|.
        const scalar psiC = recon.stencilCellValue(c);
        const scalar span = L*Foam::mag(foot - mesh_.C()[c]);
        r.lo = psiC - span;
        r.hi = psiC + span;
        r.enforced = true;
    }
    else
    {
        runawayCap(stencilLo, stencilHi, r);
    }
}


void Foam::lipschitzConeValueBound::printBanner() const
{
    Info<< "    lipschitzCone: mode " << lipschitzMode_
        << ", L0 " << L0_
        << ", onInadmissible " << onInadmissible_
        << ", boundaryFaces " << useBoundaryFaces_ << endl;

    if (L0_ != 1.0)
    {
        Info<< "    WARNING: lipschitzConstant " << L0_ << " != 1 is a TUNED"
            << " coefficient; the coefficient-free claim does not hold." << endl;
    }
}

// ************************************************************************* //
