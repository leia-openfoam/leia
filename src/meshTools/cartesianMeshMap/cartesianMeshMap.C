/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Julian Reitzel
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

#include "cartesianMeshMap.H"


// * * * * * * * * * * * * * * * Free Functions * * * * * * * * * * * * * * //

bool Foam::isMeshBlockMesh(const Foam::polyMesh& mesh, const Foam::blockMesh& blockMesh)
{
    // The mesh is considered a blockMesh if its points are equal to the one of
    // the blockMesh created from the blockMeshDict

    const pointField& meshPoints = mesh.points();
    const pointField& blockMeshPoints = blockMesh.points();
    std::vector<point> mp(meshPoints.cbegin(), meshPoints.cend());
    std::vector<point> bmp(blockMeshPoints.cbegin(), blockMeshPoints.cend());
    std::sort(mp.begin(), mp.end());
    std::sort(bmp.begin(), bmp.end());


    std::vector<point> difference{};
    std::set_difference
    (
        mp.cbegin(), 
        mp.cend(), 
        bmp.cbegin(),
        bmp.cend(),
        std::inserter(difference, difference.begin())
    );
    return difference.empty();
}


bool Foam::isMeshCartesian(const Foam::fvMesh& mesh)
{
    // equidistant
    const surfaceScalarField magDeltas = mag(mesh.delta());
    scalar reference_delta = magDeltas[0];
    bool equidistant = std::all_of(magDeltas.cbegin(), magDeltas.cend(), 
        [reference_delta](scalar delta){return delta == reference_delta;});

    // orthogonal
    bool orthogonal = !mesh.checkFaceOrthogonality();

    // alignment with coordinate axes
    bool alignment = false;
    if (orthogonal)
    {
        const cellShape& cellShape0 = mesh.cellShapes()[0];
        const faceList facesCell0 = cellShape0.faces();
        alignment = true;
        vector nx(1,0,0);
        vector ny(0,1,0);
        vector nz(0,0,1);
        for (const face& f : facesCell0)
        {
            vector nf = f.unitNormal(mesh.points());
            if (!(nf == nx || nf == ny || nf == nz ))
            {
                alignment = false;
            }
        }
    }

    // Check
    return orthogonal && equidistant && alignment;
}


bool Foam::isBlockMeshOneBlock(const Foam::blockMesh& blockMesh)
{
    if (blockMesh.size() == 1) // Single block
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
const Foam::blockMesh& Foam::cartesianMeshMap::CreateBlockMeshFromDict()
{
    if (!blockMeshPtr_)
    {
        word regionName(polyMesh::defaultRegion);

        const word dictName("blockMeshDict");

        autoPtr<IOdictionary> meshDictPtr;

        {
            fileName dictPath;
            const word& regionDir = polyMesh::regionName(regionName);

            // Assume dictionary is to be found in the system directory
            dictPath = mesh().time().system()/regionDir/dictName;

            IOobject meshDictIO
            (
                dictPath,
                mesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
            {
                FatalErrorInFunction
                    << meshDictIO.objectPath() << nl
                    << exit(FatalError);
            }

            Info<< "Creating block mesh from "
                << meshDictIO.objectRelPath() << endl;

            meshDictPtr = autoPtr<IOdictionary>::New(meshDictIO);
        }

        const IOdictionary& meshDict = *meshDictPtr;

        blockMeshPtr_ = autoPtr<blockMesh>(new blockMesh(meshDict, regionName));
    }

    return *blockMeshPtr_;
}

const Foam::blockMesh& Foam::cartesianMeshMap::getBlockMesh() const
{
    if (!blockMeshPtr_)
    {
        FatalErrorInFunction
        << "blockMeshPtr_ is not defined.\n"
        << abort(FatalError);
    }
    return *blockMeshPtr_;
}

const Foam::block& Foam::cartesianMeshMap::getBlock() const
{
    if (!isBlockMeshOneBlock(getBlockMesh()))
    {
        FatalErrorInFunction
        << "blockMesh has more than one block.\n"
        << abort(FatalError);
    }
    return getBlockMesh()[0];
}

void Foam::cartesianMeshMap::set_h(const Foam::fvMesh& mesh)
{
    h_ = pow(mag(mesh.delta())()[0], -1);
}

Foam::scalar Foam::cartesianMeshMap::h() const
{
    return h_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cartesianMeshMap::cartesianMeshMap(const fvMesh& mesh)
    :
        mesh_(mesh),
        h_(0)
    {
        blockMesh::verboseOutput = false;
        CreateBlockMeshFromDict();
        isMeshCartesian(mesh);
        set_h(mesh),
        isMeshBlockMesh(mesh, getBlockMesh());
        isBlockMeshOneBlock(getBlockMesh());
    }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
const Foam::fvMesh& Foam::cartesianMeshMap::mesh() const
{
    return mesh_;
}


const Foam::labelVector& Foam::cartesianMeshMap::sizes() const
{
    return getBlock().sizes();
}


Foam::label Foam::cartesianMeshMap::operator()(label i, label j, label k) const
{
    getBlock().checkIndex(i, j, k);
    return getBlock().index(i, j, k);
}


Foam::labelVector Foam::cartesianMeshMap::operator()(label idx) const
{
    return getBlock().index(idx);
}


// ************************************************************************* //
