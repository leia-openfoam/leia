/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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
    Foam::fv::noBvGrad

Group
    grpFvGradSchemes

Description
    Gradient scheme that decorates another gradient scheme. For the internal
    field, the decorated scheme is used to calculate the gradient. At the 
    boundary cells, a least-squares gradient calculation is used without the 
    influence of boundary surfaces. For level set fields this scheme can be 
    useful.

SourceFiles
    noBvGrad.C

\*---------------------------------------------------------------------------*/

#ifndef noBvGrad_H
#define noBvGrad_H

#include "gradScheme.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class noBvGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class noBvGrad
:
    public fv::gradScheme<Type>
{
private:
    // Private Data

        tmp<fv::gradScheme<Type>> tinnerGradScheme_;

    // Private Member Functions

        //- No copy construct
        noBvGrad(const noBvGrad&) = delete;

        //- No copy assignment
        void operator=(const noBvGrad&) = delete;

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > modified_leastSquaresGrad_calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;


public:

    //- Runtime type information
    TypeName("noBvGrad");


    // Constructors

        //- Construct from mesh
        noBvGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinnerGradScheme_(new gaussGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        noBvGrad(const fvMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh),
            tinnerGradScheme_(nullptr)
        {
            if (is.eof())
            {
                tinnerGradScheme_ =
                    tmp<fv::gradScheme<Type>>
                    (
                        new gaussGrad<Type>(mesh)
                    );
            }
            else
            {
                tinnerGradScheme_ =
                    tmp<fv::gradScheme<Type>>
                    (
                        fv::gradScheme<Type>::New(mesh, is)
                    );
            }
        }


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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noBvGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
