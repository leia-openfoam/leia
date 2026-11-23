/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slQuadraticFitMeshTransfer

Description
    Transfer of psi with the quadratic fit of the reconstruction, for hexahedral
    meshes with 2:1 refinement (blastAMR, dynamicRefineFvMesh).

    store(psi) fits psi on the present mesh and writes four REGISTERED fields:
    the gradient g and the Hessian H of the fit at the cell centre, the cell
    centre, and the refinement level with the cell size. The mesh maps
    registered fields: a new child cell gets the values of its parent, a merged
    cell gets the volume-weighted mean of its children.

    apply(psi) then corrects psi cell by cell, with d = (new centre) - (mapped
    centre):
        child of a refined cell   psi = psi_parent + g.d + d.H.d/2
        merged cell               psi = mean(psi_children) - (h_c^2/8) tr(H)
        every other cell          unchanged, bit for bit
    h_c is the size of the children. Both forms are exact for a quadratic psi:
    the eight children of a cube sit at +-h_c/2 in every direction, so the mean
    of a quadratic over them is its centre value plus (h_c^2/8) tr(H).
    A child and a merged cell are recognised by the level: the mapped level
    against the present level. One level per mesh update is assumed.

    The fit is the reconstruction of levelSet.semiLagrangian.reconstruction; it
    must provide fitDerivatives. The object is built for one mesh.

    Selected with `meshTransfer quadraticFit`.

SourceFiles
    slQuadraticFitMeshTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef slQuadraticFitMeshTransfer_H
#define slQuadraticFitMeshTransfer_H

#include "slMeshTransfer.H"
#include "slReconstruction.H"

namespace Foam
{

class slQuadraticFitMeshTransfer
:
    public slMeshTransfer
{
    autoPtr<slReconstruction> fit_;

    //- Registered fields that the mesh maps.
    volVectorField grad_;
    volSymmTensorField hess_;
    volVectorField centre_;
    volScalarField level_;
    volScalarField size_;

    bool stored_;

    //- Present refinement level of every cell (registered cellLevel, else 0).
    labelList presentLevel() const;

public:

    TypeName("quadraticFit");

    explicit slQuadraticFitMeshTransfer(const fvMesh& mesh);

    virtual ~slQuadraticFitMeshTransfer() = default;

    virtual void store(const volScalarField& psi);

    virtual label apply(volScalarField& psi);
};

} // End namespace Foam

#endif
