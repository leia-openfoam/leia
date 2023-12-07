/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "noBvGrad.H"
#include "MeshObject.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class modified_leastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, modified_leastSquaresVectors>
{

    // Private Data

        //- Owner least-squares gradient vectors
        surfaceVectorField pVectors_;

        //- Neighbour least-squares gradient vectors
        surfaceVectorField nVectors_;


    // Private Member Functions

        //- Construct Least-squares gradient vectors
        void calcLeastSquaresVectors()
        {
            DebugInFunction << "Calculating least square gradient vectors" << nl;

            const fvMesh& mesh = mesh_;

            // Set local references to mesh data
            const labelUList& owner = mesh_.owner();
            const labelUList& neighbour = mesh_.neighbour();

            const volVectorField& C = mesh.C();

            // Set up temporary storage for the dd tensor (before inversion)
            symmTensorField dd(mesh_.nCells(), Zero);

            forAll(owner, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];

                const vector d(C[nei] - C[own]);
                const symmTensor wdd(sqr(d)/magSqr(d));
                dd[own] += wdd;
                dd[nei] += wdd;
            }


            // surfaceVectorField::Boundary& blsP =
            //     pVectors_.boundaryFieldRef();

            // forAll(blsP, patchi)
            // {
            //     const fvsPatchVectorField& patchLsP = blsP[patchi];

            //     const fvPatch& p = patchLsP.patch();
            //     const labelUList& faceCells = p.patch().faceCells();

            //     // Build the d-vectors
            //     const vectorField pd(p.delta());

            //     forAll(pd, patchFacei)
            //     {
            //         const vector& d = pd[patchFacei];

            //         dd[faceCells[patchFacei]] += sqr(d)/magSqr(d);
            //     }
            // }


            // Invert the dd tensor
            const symmTensorField invDd(inv(dd));


            // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
            forAll(owner, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];

                const vector d(C[nei] - C[own]);

                pVectors_[facei] = (invDd[own] & d)/magSqr(d);
                nVectors_[facei] = -(invDd[nei] & d)/magSqr(d);
            }

            // forAll(blsP, patchi)
            // {
            //     fvsPatchVectorField& patchLsP = blsP[patchi];

            //     const fvPatch& p = patchLsP.patch();
            //     const labelUList& faceCells = p.faceCells();

            //     // Build the d-vectors
            //     const vectorField pd(p.delta());

            //     forAll(pd, patchFacei)
            //     {
            //         const vector& d = pd[patchFacei];

            //         patchLsP[patchFacei] = (invDd[faceCells[patchFacei]] & d)/magSqr(d);
            //     }
            // }

            DebugInfo << "Finished calculating least square gradient vectors" << endl;
        }


public:

    // Declare name of the class and its debug switch
    TypeName("modified_leastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit modified_leastSquaresVectors(const fvMesh& mesh)
        :
            MeshObject<fvMesh, Foam::MoveableMeshObject, modified_leastSquaresVectors>(mesh),
            pVectors_
            (
                IOobject
                (
                    "LeastSquaresP",
                    mesh_.pointsInstance(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedVector(dimless/dimLength, Zero)
            ),
            nVectors_
            (
                IOobject
                (
                    "LeastSquaresN",
                    mesh_.pointsInstance(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedVector(dimless/dimLength, Zero)
            )
        {
            calcLeastSquaresVectors();
        }


    //- Destructor
    virtual ~modified_leastSquaresVectors(){}

    // Member functions

        //- Return const reference to owner least square vectors
        const surfaceVectorField& pVectors() const
        {
            return pVectors_;
        }

        //- Return const reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const
        {
            return nVectors_;
        }

        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints()
        {
            calcLeastSquaresVectors();
            return true;
        }
};
// }
// namespace Foam
// {
    defineTypeNameAndDebug(modified_leastSquaresVectors, 0);
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::noBvGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    // Delegate gradient calculation 
    tmp<GradFieldType> tgrad = tinnerGradScheme_().calcGrad(vsf, name);
    
    // GradFieldType& grad = tgrad.ref();
    // auto& gradbf = grad.boundaryFieldRef();

    // // Fix gradients in boundary cells

    // // fvc::reconstruct(fvc::snGrad(vsf)*mesh.magSf())

    // forAll(vsf.boundaryField(), patchi)
    // {
    //     const fvPatchField&     pfield = vsf.boundaryField()[patchi];
    //     fvPatchField::Internal& ifield = vsf.boundaryField()[patchi].internalField();

    //     if (!vsf.boundaryField()[patchi].coupled())
    //     {
    //         forAll(pfield, facei)
    //         {
    //             // grad(own(face)) = fvc::snGrad(oppositeFace(face)
    //             ifield.[facei] = 
    //         }


    //     }
    // }


    return tgrad;
}

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::noBvGrad<Type>::modified_leastSquaresGrad_calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;
 
    const fvMesh& mesh = vsf.mesh();
 
    tmp<GradFieldType> tlsGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>(vsf.dimensions()/dimLength, Zero),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GradFieldType& lsGrad = tlsGrad.ref();
 
    // Get reference to least square vectors
    const modified_leastSquaresVectors& lsv = modified_leastSquaresVectors::New(mesh);
 
    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();
 
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
 
    forAll(own, facei)
    {
        const label ownFacei = own[facei];
        const label neiFacei = nei[facei];
 
        const Type deltaVsf(vsf[neiFacei] - vsf[ownFacei]);
 
        lsGrad[ownFacei] += ownLs[facei]*deltaVsf;
        lsGrad[neiFacei] -= neiLs[facei]*deltaVsf;
    }
 
    // Boundary faces
    forAll(vsf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];
 
        const labelUList& faceCells =
            vsf.boundaryField()[patchi].patch().faceCells();
 
        if (vsf.boundaryField()[patchi].coupled())
        {
            const Field<Type> neiVsf
            (
                vsf.boundaryField()[patchi].patchNeighbourField()
            );
 
            forAll(neiVsf, patchFacei)
            {
                lsGrad[faceCells[patchFacei]] +=
                    patchOwnLs[patchFacei]
                   *(neiVsf[patchFacei] - vsf[faceCells[patchFacei]]);
            }
        }
        else
        {
            // const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
 
            // forAll(patchVsf, patchFacei)
            // {
            //     lsGrad[faceCells[patchFacei]] +=
            //          patchOwnLs[patchFacei]
            //         *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
            // }
        }
    }
 
 
    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);
 
    return tlsGrad;
}

// ************************************************************************* //
