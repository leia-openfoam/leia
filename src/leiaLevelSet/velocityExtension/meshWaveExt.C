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

#include "meshWaveExt.H"
#include "addToRunTimeSelectionTable.H"
#include "MeshWave.H"
#include "wallPointData.H"
#include "interpolationCellPoint.H"
#include "surfaceInterpolate.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshWaveExt, 0);
    addToRunTimeSelectionTable(velocityExtension, meshWaveExt, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWaveExt::meshWaveExt(const fvMesh& mesh)
:
    interfaceExtension(mesh)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::meshWaveExt::extend()
{
    // Seed every psi-sign-change internal face with the interface point
    // nearest to it (face foot point) and the interface velocity there; the
    // wave then delivers the nearest seed's payload to every cell -- a
    // discrete, distance-ordered closest-point extension on arbitrary
    // polyhedral meshes, parallel-native via FaceCellWave.
    interpolationCellPoint<vector> Uinterp(U_);

    const surfaceScalarField psiF(fvc::interpolate(psi_));
    const surfaceVectorField nF(fvc::interpolate(nHat_));
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    DynamicList<label> seedFaces(seedCells_.size()*2);
    DynamicList<wallPointData<vector>> seedInfo(seedCells_.size()*2);

    forAll(own, f)
    {
        if (psi_[own[f]]*psi_[nei[f]] < 0)
        {
            const point& xf = mesh_.Cf()[f];
            const vector n = nF[f]/(Foam::mag(nF[f]) + SMALL);
            const point fp = xf - psiF[f]*n;   // face foot point on Sigma

            // Interface velocity at the foot point (owner-cell localized;
            // falls back to the owner-cell value if the foot point left the
            // processor-local mesh).
            const label cf = mesh_.findCell(fp);
            const vector Ufp =
                (cf >= 0) ? Uinterp.interpolate(fp, cf) : U_[own[f]];

            seedFaces.append(f);
            seedInfo.append
            (
                wallPointData<vector>(fp, Ufp, magSqr(xf - fp))
            );
        }
    }

    // Whole-mesh wave (cheap); only band cells consume the payload -- the
    // seamless fade handles the transition to the base velocity beyond it.
    MeshWave<wallPointData<vector>> wave
    (
        mesh_,
        seedFaces,
        seedInfo,
        mesh_.globalData().nTotalCells() + 1
    );

    const List<wallPointData<vector>>& cellInfo = wave.allCellInfo();
    int dummyTd = 0;
    forAll(band_, c)
    {
        if (band_[c] > 0.5 && cellInfo[c].valid(dummyTd))
        {
            Uext_[c] = cellInfo[c].data();
        }
    }
}

// ************************************************************************* //
