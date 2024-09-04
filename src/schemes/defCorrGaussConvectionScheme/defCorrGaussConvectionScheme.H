/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fv::defCorrGaussConvectionScheme

Group
    grpFvConvectionSchemes

Description
    Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
    defCorrGaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef defCorrGaussConvectionScheme_H
#define defCorrGaussConvectionScheme_H

#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class defCorrGaussConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class defCorrGaussConvectionScheme
:
    public fv::gaussConvectionScheme<Type>
{

public:

    //- Runtime type information
    TypeName("GaussDefCorr");


    // Constructors

        //- Construct from flux and interpolation scheme
        defCorrGaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<surfaceInterpolationScheme<Type>>& scheme
        )
        :
            gaussConvectionScheme<Type>(mesh, faceFlux, scheme)
        {}

        //- Construct from flux and Istream
        defCorrGaussConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            gaussConvectionScheme<Type>(mesh, faceFlux, is)
        {}


    // Member Functions

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "defCorrGaussConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
