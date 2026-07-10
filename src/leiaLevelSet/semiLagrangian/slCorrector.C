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
    label& nNonFinite
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
    }
    else if (clip)
    {
        v = Foam::min(Foam::max(v, lo), hi);
    }
    else
    {
        const scalar cap = 10.0*Foam::max(hi - lo, SMALL);
        v = Foam::min(Foam::max(v, mid - cap), mid + cap);
    }
    return v;
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
