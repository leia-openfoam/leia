/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slQuadraticFitMeshTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "labelIOList.H"
#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(slQuadraticFitMeshTransfer, 0);
    addToRunTimeSelectionTable(slMeshTransfer, slQuadraticFitMeshTransfer, Mesh);
}

namespace
{
template<class Type>
Foam::IOobject transferIO(const Foam::word& name, const Foam::fvMesh& mesh)
{
    // Registered (the mesh maps it), never written.
    return Foam::IOobject
    (
        name,
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::NO_READ,
        Foam::IOobject::NO_WRITE,
        true
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slQuadraticFitMeshTransfer::slQuadraticFitMeshTransfer(const fvMesh& mesh)
:
    slMeshTransfer(mesh),
    fit_(slReconstruction::New(mesh)),
    grad_
    (
        transferIO<vector>("slTransferGrad", mesh), mesh,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    hess_
    (
        transferIO<symmTensor>("slTransferHess", mesh), mesh,
        dimensionedSymmTensor(dimless, Zero),
        zeroGradientFvPatchSymmTensorField::typeName
    ),
    centre_
    (
        transferIO<vector>("slTransferCentre", mesh), mesh,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    level_
    (
        transferIO<scalar>("slTransferLevel", mesh), mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    size_
    (
        transferIO<scalar>("slTransferSize", mesh), mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    stored_(false)
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::slQuadraticFitMeshTransfer::presentLevel() const
{
    if (mesh_.foundObject<labelIOList>("cellLevel"))
    {
        const labelList& cl = mesh_.lookupObject<labelIOList>("cellLevel");
        if (cl.size() == mesh_.nCells()) { return cl; }
    }
    return labelList(mesh_.nCells(), Zero);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::slQuadraticFitMeshTransfer::store(const volScalarField& psi)
{
    fit_->update(psi);

    const labelList lev(presentLevel());
    const vectorField& C = mesh_.C();
    const scalarField& V = mesh_.V();
    vectorField& g = grad_.primitiveFieldRef();
    symmTensorField& H = hess_.primitiveFieldRef();
    vectorField& x = centre_.primitiveFieldRef();
    scalarField& L = level_.primitiveFieldRef();
    scalarField& h = size_.primitiveFieldRef();
    forAll(g, c)
    {
        fit_->fitDerivatives(c, g[c], H[c]);   // order 1 or 0: H or g is zero
        x[c] = C[c];
        L[c] = lev[c];
        h[c] = Foam::cbrt(V[c]);
    }
    stored_ = true;
}


Foam::label Foam::slQuadraticFitMeshTransfer::apply(volScalarField& psi)
{
    if (!stored_)
    {
        FatalErrorInFunction
            << "apply() without store() on the mesh before the change"
            << exit(FatalError);
    }
    const labelList lev(presentLevel());
    const vectorField& C = mesh_.C();
    const vectorField& g = grad_.primitiveField();
    const symmTensorField& H = hess_.primitiveField();
    const vectorField& x = centre_.primitiveField();
    const scalarField& L = level_.primitiveField();
    const scalarField& h = size_.primitiveField();
    scalarField& p = psi.primitiveFieldRef();

    label nChild = 0, nMerged = 0;
    forAll(p, c)
    {
        const label Lold = label(std::floor(L[c] + 0.5));
        if (lev[c] > Lold)
        {
            const vector d = C[c] - x[c];
            p[c] += (g[c] & d) + 0.5*(d & (H[c] & d));
            ++nChild;
        }
        else if (lev[c] < Lold)
        {
            p[c] -= 0.125*sqr(h[c])*tr(H[c]);
            ++nMerged;
        }
    }
    psi.correctBoundaryConditions();
    stored_ = false;

    reduce(nChild, sumOp<label>());
    reduce(nMerged, sumOp<label>());
    Info<< "slMeshTransfer quadraticFit: " << nChild << " child cells, "
        << nMerged << " merged cells corrected" << endl;
    return nChild + nMerged;
}
