/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Class
    Foam::slMeshTransfer

Description
    Transfer of the level set psi to a mesh that a dynamic mesh has refined,
    unrefined or redistributed. Run-time-selectable family; the base class is
    the inert model `none`: psi keeps the values that the mesh gave it
    (OpenFOAM copies the parent value to new child cells, which is first order,
    and takes the volume-weighted mean of merged cells).

    A model is used in three calls per time step, by the owner (the solver):
        store(psi)   before mesh.update(): keep what the transfer needs;
        apply(psi)   after a topology change: correct the mapped values;
    and the owner then constructs a new object, because a model may hold
    objects that are built for one mesh.

    Dictionary: levelSet { semiLagrangian { meshTransfer <type>; } },
    default none.

SourceFiles
    slMeshTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef slMeshTransfer_H
#define slMeshTransfer_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class slMeshTransfer
{
protected:

    const fvMesh& mesh_;

public:

    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slMeshTransfer,
        Mesh,
        (const fvMesh& mesh),
        (mesh)
    );

    explicit slMeshTransfer(const fvMesh& mesh);

    static autoPtr<slMeshTransfer> New(const fvMesh& mesh);

    virtual ~slMeshTransfer() = default;

    //- Keep the data of psi on the present mesh. Call before mesh.update().
    virtual void store(const volScalarField&)
    {}

    //- Correct the values that the mesh change gave to psi. Call once, after a
    //  topology change. Returns the number of cells that were corrected.
    virtual label apply(volScalarField&)
    {
        return 0;
    }
};

} // End namespace Foam

#endif
