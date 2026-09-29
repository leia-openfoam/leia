/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCurvature (davofCurvature.H)
\*---------------------------------------------------------------------------*/

#include "davofCurvature.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(davofCurvature, false);
    defineRunTimeSelectionTable(davofCurvature, Dictionary);
    addToRunTimeSelectionTable(davofCurvature, davofCurvature, Dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::davofCurvature::davofCurvature
(
    const word& name,
    const fvMesh& mesh,
    const davofState& state,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    state_(state),
    dict_(dict),
    kappaS_
    (
        IOobject
        (
            "kappa.davof." + name, mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh, dimensionedScalar(dimless/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    KS_
    (
        IOobject
        (
            "K.davof." + name, mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh, dimensionedScalar(dimless/dimArea, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa1S_
    (
        IOobject
        (
            "kappa1.davof." + name, mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh, dimensionedScalar(dimless/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa2S_
    (
        IOobject
        (
            "kappa2.davof." + name, mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        mesh, dimensionedScalar(dimless/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    filled_(mesh.nCells(), false),
    nFallback_(0),
    cpuSeconds_(0)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::davofCurvature> Foam::davofCurvature::New
(
    const word& name,
    const fvMesh& mesh,
    const davofState& state,
    const dictionary& dict
)
{
    const word type = dict.getOrDefault<word>("type", "none");

    auto* ctorPtr = DictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "davofCurvature",
            type,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Info<< "Selecting davof curvature model " << name << " : " << type << endl;

    return autoPtr<davofCurvature>(ctorPtr(name, mesh, state, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::davofCurvature::write() const
{
    kappaS_.write();
    KS_.write();
    kappa1S_.write();
    kappa2S_.write();
}


// ************************************************************************* //
