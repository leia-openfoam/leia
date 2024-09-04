/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Tomislav Maric, TU Darmstadt 
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
    Foam::upwindSecondOrderDefCorr

Group
    grpFvLimitedSurfaceInterpolationSchemes

Description
    A second-order explicit defect correction for the upwind scheme. 

SourceFiles
    upwindSecondOrderDefCorr.C

\*---------------------------------------------------------------------------*/

#ifndef upwindSecondOrderDefCorr_H
#define upwindSecondOrderDefCorr_H

#include "upwind.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwindSecondOrderDefCorr Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindSecondOrderDefCorr
:
    public upwind<Type>
{

public:

    //- Runtime type information
    TypeName("upwindSecondOrderDefCorr");


    // Constructors

        //- Construct from faceFlux
        upwindSecondOrderDefCorr
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        upwindSecondOrderDefCorr
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            upwind<Type>(mesh, is)
        {}

        //- Construct from faceFlux and Istream
        upwindSecondOrderDefCorr
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            upwind<Type>(mesh, faceFlux)
        {}


    // Member Functions
    
        //- Return true if this scheme uses an explicit correction            
        virtual bool corrected() const                                        
        {                                                                     
            return true;                                                     
        }              

	    //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>& vf) const; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialize for Tensor, SymmTensor, SphericalTensor as not implemented,  
// because the defect correction requires fvc::grad, which is not available 
// for tensor fields, with the Gauss gradient scheme. 
template<>
tmp<GeometricField<Tensor<scalar>, fvsPatchField, surfaceMesh>>
upwindSecondOrderDefCorr<Tensor<scalar>>::correction
(
    const GeometricField<Tensor<scalar>, fvPatchField, volMesh>& vf
) const
{
    notImplemented("Defect requires fvc::grad(vf), unavailable for Tensor.");
    return tmp<GeometricField<Tensor<scalar>, fvsPatchField, surfaceMesh>>(nullptr);
}

template<>
tmp<GeometricField<SymmTensor<scalar>, fvsPatchField, surfaceMesh>>
upwindSecondOrderDefCorr<SymmTensor<scalar>>::correction
(
    const GeometricField<SymmTensor<scalar>, fvPatchField, volMesh>& vf
) const
{
    notImplemented("Defect requires fvc::grad(vf), unavailable for Tensor.");
    return tmp<GeometricField<SymmTensor<scalar>, fvsPatchField, surfaceMesh>>(nullptr);
}

template<>
tmp<GeometricField<SphericalTensor<scalar>, fvsPatchField, surfaceMesh>>
upwindSecondOrderDefCorr<SphericalTensor<scalar>>::correction
(
    const GeometricField<SphericalTensor<scalar>, fvPatchField, volMesh>& vf
) const
{
    notImplemented("Defect requires fvc::grad(vf), unavailable for Tensor.");
    return tmp<GeometricField<SphericalTensor<scalar>, fvsPatchField, surfaceMesh>>(nullptr);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "upwindSecondOrderDefCorr.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
