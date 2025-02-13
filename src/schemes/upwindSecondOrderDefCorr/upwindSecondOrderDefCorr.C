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

#include "upwindSecondOrderDefCorr.H"
#include "error.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<typename Type> 
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
upwindSecondOrderDefCorr<Type>::correction
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

    // Compute cell-centered gradient of vf.
    auto gradVfTmp = fvc::grad(vf); 
    const auto& gradVf = gradVfTmp();
    const auto& mesh = this->mesh();
    // Get owner-neighbour addressing.
    const auto& own = mesh.owner();
    const auto& nei = mesh.neighbour();
    // Get cell centers.
    const auto& C = mesh.C();
    // Get face centers.
    const auto& Cf = mesh.Cf();

    const auto& faceFlux = this->faceFlux_;

    // For internal faces
    forAll(own, faceI)
    {
        // If flux is positive, owner-cell is upwind. 
        if (faceFlux[faceI] > 0)
        {
            vfErr[faceI] = (gradVf[own[faceI]] & (Cf[faceI] - C[own[faceI]]));
        }
        else // If flux is negative, neighbor-cell is upwind. 
        {
            vfErr[faceI] = (gradVf[nei[faceI]] & (Cf[faceI] - C[nei[faceI]])); 
        }
    }

    // Computing vfErr on coupled boundaries.
    auto& vfErrBdryField = vfErr.boundaryFieldRef(); 
    // Boundary data for vfErr calculation.
    const auto& gradVfBdryField = gradVf.boundaryField(); 
    const auto& cfBdryField = Cf.boundaryField(); 
    const auto& cBdryField = C.boundaryField(); 
    const auto& faceFluxBdryField = faceFlux.boundaryField();
    const auto& patches = mesh.boundary(); 
    const auto& faceOwner = mesh.faceOwner();

    // For all boundary patches (faces).
    forAll(vfErrBdryField, patchI)
    {
        const fvPatch& patch = patches[patchI];

        // Face centers patch field
        const auto& cfPatch = cfBdryField[patchI];

        // Face-centered error patch field for calculation 
        auto& vfErrPatch = vfErrBdryField[patchI];

        // Face-centered volumetric flux patch field
        const auto& faceFluxPatchField = faceFluxBdryField[patchI];

        // Compute vfErr on all outflow patches
        forAll(faceFluxPatchField, faceI)
        {
            const label faceG = faceI + patch.start(); // Global label. 
            // If flux is positive, owner-cell is upwind. 
            if (faceFluxPatchField[faceI] > 0) // If flux is positive
            {
                vfErrPatch[faceI] = 
                (   
                    gradVf[faceOwner[faceG]] &  
                    (cfPatch[faceI] - C[faceOwner[faceG]])
                );
            }
        }

        if (isA<coupledFvPatch>(patch)) // coupled patch
        {
            // Get gradVf across coupled patch boundary 
            const auto& gradVfPatch = gradVfBdryField[patchI];
            auto gradVfNeiTmp = gradVfPatch.patchNeighbourField();
            const auto& gradVfNei = gradVfNeiTmp();

            // Get cell centers across coupled patch boundary 
            const auto& cPatch = cBdryField[patchI];
            auto cNeiTmp = cPatch.patchNeighbourField();
            const auto& cNei = cNeiTmp();

            // Compute vfErr on coupled patch.
            forAll(faceFluxPatchField, faceI)
            {
                if (faceFluxPatchField[faceI] < 0) // If flux is negative 
                {
                    // Coupled-patch neighbor is upwind
                    vfErrPatch[faceI] = 
                    (   
                        gradVfNei[faceI] &  
                        (cfPatch[faceI] - cNei[faceI])
                    );
                }
            }
        }
    }

    return vfErrTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
