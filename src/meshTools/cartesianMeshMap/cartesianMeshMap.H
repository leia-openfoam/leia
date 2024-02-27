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

Class
    Foam::cartesianMeshMap

Description
    This class implements a mapping from global cell indices to the indices of
    a 3D orthogonal equidistant mesh.

SourceFiles
    cartesianMeshMapI.H
    cartesianMeshMap.C
    cartesianMeshMapIO.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshMap_H
#define cartesianMeshMap_H

#include "blockMesh.H"
#include "Time.H"
#include "surfaceFields.H"
#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Free Functions  * * * * * * * * * * * * * * * * * //
        bool isBlockMeshOneBlock(const blockMesh&);
        
        bool isMeshBlockMesh(const polyMesh& mesh, const blockMesh& blockMesh);

        bool isMeshCartesian(const fvMesh& mesh);


/*---------------------------------------------------------------------------*\
                         Class cartesianMeshMap Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshMap
{
    // Private Data

        //- Description of data_
        const fvMesh& mesh_;
        scalar h_;
        autoPtr<blockMesh> blockMeshPtr_ = nullptr;


    // Private Member Functions
        const blockMesh& CreateBlockMeshFromDict();

        const blockMesh& getBlockMesh() const;

        const block& getBlock() const;


        void set_h(const fvMesh& mesh);

public:

    // Constructors

        
        cartesianMeshMap(const fvMesh& mesh);

    // Member Functions

    // Access

        const fvMesh& mesh() const;

        scalar h() const;
        const labelVector& sizes() const;

        label operator()(label i, label j, label k) const;

        labelVector operator()(label idx) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
