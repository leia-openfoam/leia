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

#include "planeFootWaveRedistancer.H"
#include "planeAnchors.H"
#include "addToRunTimeSelectionTable.H"
#include "MeshWave.H"
#include "wallPointData.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(planeFootWaveRedistancer, false);
addToRunTimeSelectionTable(redistancer, planeFootWaveRedistancer, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

planeFootWaveRedistancer::planeFootWaveRedistancer(const fvMesh& mesh)
:
    redistancer(mesh),
    mesh_(mesh)
{
    Info<< "planeFootWaveRedistancer (geometric closest-foot MeshWave)"
        << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool planeFootWaveRedistancer::doRedistance(volScalarField& psi)
{
    if (!mesh_.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "planeFootWave redistancing requires the registered "
            << "'NarrowBand' field: set levelSet.narrowBand.type to "
            << "signChange (or a type derived from it) in fvSolution."
            << exit(FatalError);
    }

    const volScalarField& narrowBand =
        mesh_.lookupObject<volScalarField>("NarrowBand");

    const planeAnchorData anchors =
        computePlaneAnchors
        (
            mesh_, psi, narrowBand,
            redistDict_.getOrDefault<bool>("curvatureCorrection", true)
        );

    if (returnReduce(anchors.bandCells.size(), sumOp<label>()) == 0)
    {
        Info<< "planeFootWaveRedistancer: no interface anchors, skipping."
            << endl;
        return false;
    }

    const scalarField psiOld(psi.primitiveField());
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();

    // Seed every face of every band anchor cell with that cell's DONOR
    // PLANE: origin = plane foot point (the wave's nearest-donor ordering
    // metric), payload = unit plane normal. A face shared by two anchor
    // cells would be seeded twice and FaceCellWave::setFaceInfo blindly
    // overwrites -- deduplicate, keeping the seed with the smaller distSqr.
    Map<label> faceToSeed(8*anchors.bandCells.size());
    DynamicList<label> seedFaces(8*anchors.bandCells.size());
    DynamicList<wallPointData<vector>> seedInfo(8*anchors.bandCells.size());

    forAll(anchors.bandCells, i)
    {
        const label cellI = anchors.bandCells[i];
        const point& foot = anchors.bandFoot[i];
        const vector& nHat = anchors.bandNHat[i];

        for (const label faceI : mesh_.cells()[cellI])
        {
            const scalar dist2 = magSqr(Cf[faceI] - foot);

            const auto iter = faceToSeed.cfind(faceI);
            if (iter.good())
            {
                if (dist2 < seedInfo[iter.val()].distSqr())
                {
                    seedInfo[iter.val()] =
                        wallPointData<vector>(foot, nHat, dist2);
                }
            }
            else
            {
                faceToSeed.insert(faceI, seedFaces.size());
                seedFaces.append(faceI);
                seedInfo.append(wallPointData<vector>(foot, nHat, dist2));
            }
        }
    }

    // Whole-mesh nearest-donor wave (parallel-native FaceCellWave). The
    // wave ORDERS donors by distance to the foot point; the received VALUE
    // is the CONTINUOUS distance to the donor plane (geometricVoF/plicRDF's
    // trick): evaluating |(x - foot) & nHat| instead of |x - foot| removes
    // the scalloping floor of the discrete foot-point cloud -- measured on
    // the static gate, point-cloud distances RAISED the band deviation an
    // event is supposed to lower.
    MeshWave<wallPointData<vector>> wave
    (
        mesh_,
        seedFaces,
        seedInfo,
        mesh_.globalData().nTotalCells() + 1
    );

    const List<wallPointData<vector>>& cellInfo = wave.allCellInfo();

    // Unsigned donor-plane distance, then re-sign from the pre-redistance
    // field; anchor cells (band + first ring) keep their (signed) plane
    // distances exactly.
    scalarField& psiIn = psi.primitiveFieldRef();
    int dummyTd = 0;
    forAll(psiIn, celli)
    {
        const scalar d =
            cellInfo[celli].valid(dummyTd)
          ? mag((C[celli] - cellInfo[celli].origin()) & cellInfo[celli].data())
          : mag(psiOld[celli]);   // disconnected region: keep the old value

        psiIn[celli] = (psiOld[celli] >= 0) ? d : -d;
    }
    forAll(anchors.cells, i)
    {
        psiIn[anchors.cells[i]] = anchors.dist[i];
    }
    psi.correctBoundaryConditions();

    Info<< "planeFootWaveRedistancer: " << anchors.bandCells.size()
        << " foot seeds over " << seedFaces.size() << " faces" << endl;

    return true;
}

} // End namespace Foam

// ************************************************************************* //
