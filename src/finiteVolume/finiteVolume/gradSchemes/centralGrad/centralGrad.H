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
    Foam::fv::centralGrad

Group
    grpFvGradSchemes

Description
    Basic second-order gradient scheme using the second order accurate central
    and equidistant meshes. At the boundaries the second order accurate 
    forward gradient scheme is used.

SourceFiles
    centralGrad.C

\*---------------------------------------------------------------------------*/

#ifndef centralGrad_H
#define centralGrad_H

#include "gradScheme.H"
#include "cartesianMeshMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class centralGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class centralGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

    cartesianMeshMap map_;


    // Private Member Functions

        //- No copy construct
        centralGrad(const centralGrad&) = delete;

        //- No copy assignment
        void operator=(const centralGrad&) = delete;

        

public:

    //- Runtime type information
    TypeName("centralGrad");


    // Constructors

        //- Construct from mesh
        centralGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            map_(mesh)
        {
            forAll(map_.sizes(), cmp)
            {
                if (map_.sizes()[cmp] == 2 )
                {
                    FatalErrorInFunction
                    << "A cartesian mesh with sizes[" << cmp << "] == 2 is not supported."
                    << abort(FatalError);
                }
            }
        }

        //- Construct from mesh and Istream
        centralGrad(const fvMesh& mesh, Istream& is)
        :
            centralGrad(mesh)
        {}


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //- for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;

        //- Correct the boundary values of the gradient using the patchField
        //- snGrad functions
        static void correctBoundaryConditions
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "centralGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
