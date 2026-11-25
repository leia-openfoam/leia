/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slFootIntegrator.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

namespace Foam
{
    defineTypeNameAndDebug(slFootIntegrator, 0);
    defineRunTimeSelectionTable(slFootIntegrator, Mesh);
    addToRunTimeSelectionTable(slFootIntegrator, slFootIntegrator, Mesh);
}


Foam::slFootIntegrator::slFootIntegrator(const fvMesh& mesh)
:
    mesh_(mesh)
{}


Foam::autoPtr<Foam::slFootIntegrator>
Foam::slFootIntegrator::New(const fvMesh& mesh)
{
    const word modelType
    (
        mesh.solutionDict().subOrEmptyDict("levelSet")
            .subOrEmptyDict("semiLagrangian")
            .getOrDefault<word>("footIntegrator", "taylor")
    );

    auto* ctorPtr = MeshConstructorTable(modelType);
    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "levelSet.semiLagrangian.footIntegrator", modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalError);
    }
    return autoPtr<slFootIntegrator>(ctorPtr(mesh));
}


void Foam::slFootIntegrator::feet
(
    const volVectorField& uNew,
    const volVectorField& uOld,
    const scalar dt,
    pointField& feet
) const
{
    // The historical Taylor foot of pointValueScheme, moved here unchanged:
    //   x_d = x_c - u^{n+1} dt + 1/2 [ du/dt + (u^{n+1}.grad)u^{n+1} ] dt^2
    // (KEEP the 1/2 on the dt^2 term.) The gradient is built here and the loop
    // keeps its stack form: an autoPtr dereference inside the loop once changed
    // the compiler's FP contraction and moved the answer in the 9th digit.
    const volVectorField& C = mesh_.C();
    const volTensorField gradU(fvc::grad(uNew, "gradU"));
    forAll(C, c)
    {
        const vector& uN = uNew[c];
        const vector& uO = uOld[c];
        const vector accel = (uN - uO)/dt + (uN & gradU[c]);
        feet[c] = C[c] - uN*dt + 0.5*accel*dt*dt;
    }
}

// ************************************************************************* //
