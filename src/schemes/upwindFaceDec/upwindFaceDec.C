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

\*---------------------------------------------------------------------------*/

#include "upwindFaceDec.H"
#include "error.H"
#include "fvcSnGrad.H"
#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<typename Type> 
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
upwindFaceDec<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Implements the defect correction, which is used by the  Gauss divergence
    // scheme (gaussConvectionScheme.C) by adding an explicit source to the
    // matrix in the form of
    // fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vf)); 
    // Therefore, the correction is the interpolation error estimate vfErr, 
    // and the error must have a negative sign as it is added as positive on
    // the right hand side of the linear system fvm (fvMatrix source). TM.

    using surfaceField = GeometricField<Type, fvsPatchField, surfaceMesh>;

    tmp<surfaceField> vfErrTmp
    (
        new surfaceField
        (
            IOobject
            (
                "vfErr",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(), 
            dimensioned<Type>("vfErr", vf.dimensions(), pTraits<Type>::zero) 
        )
    ); 
    surfaceField& vfErr = vfErrTmp.ref();

    const fvMesh& mesh = this->mesh();

    // 1) Face-centered upwind interpolation of cell centers C
    tmp<surfaceVectorField> CfUpwindTmp = upwind<vector>(mesh, this->faceFlux_).interpolate(mesh.C());
    const surfaceVectorField& CfUpwind  = CfUpwindTmp();

    // 2) Face-centered surface-normal gradient
    tmp<surfaceField> snGradVfTmp = fvc::snGrad(vf, snGradSchemeName_);
    const surfaceField& snGradVf  = snGradVfTmp();

    // 3) Compute connector, projection, and defect in one go
    const surfaceVectorField& Cf    = mesh.Cf();
    const surfaceVectorField& Sf    = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    // - connector vector for all faces
    surfaceVectorField dUf = CfUpwind - Cf;
    // - projection of connector onto face normal for all faces
    surfaceScalarField proj = (dUf & Sf) / magSf;

    // 4) Compute the defect for all faces: assign entire field including boundaries
    vfErr = - snGradVf * proj;

    return vfErrTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
