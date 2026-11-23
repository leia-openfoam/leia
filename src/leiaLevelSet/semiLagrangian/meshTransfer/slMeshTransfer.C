/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slMeshTransfer.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(slMeshTransfer, 0);
    defineRunTimeSelectionTable(slMeshTransfer, Mesh);
    addToRunTimeSelectionTable(slMeshTransfer, slMeshTransfer, Mesh);
}

Foam::slMeshTransfer::slMeshTransfer(const fvMesh& mesh)
:
    mesh_(mesh)
{}

Foam::autoPtr<Foam::slMeshTransfer> Foam::slMeshTransfer::New(const fvMesh& mesh)
{
    const word modelType
    (
        mesh.solutionDict().subOrEmptyDict("levelSet")
            .subOrEmptyDict("semiLagrangian")
            .getOrDefault<word>("meshTransfer", "none")
    );

    auto* ctorPtr = MeshConstructorTable(modelType);
    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "slMeshTransfer", modelType, *MeshConstructorTablePtr_
        ) << exit(FatalError);
    }
    return autoPtr<slMeshTransfer>(ctorPtr(mesh));
}
