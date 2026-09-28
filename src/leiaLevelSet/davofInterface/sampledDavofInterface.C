/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: sampledDavofInterface (see the header)
\*---------------------------------------------------------------------------*/

#include "sampledDavofInterface.H"
#include "dictionary.H"
#include "volFields.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledDavofInterface, 0);

    // Selected as `davofInterface`: OpenFOAM's libgeometricVoF registers
    // `interface`, TwoPhaseFlow `reconstructedInterface`.
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledDavofInterface,
        word,
        davofInterface
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledDavofInterface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }
    prevTimeIndex_ = fvm.time().timeIndex();

    surfPtr_.clear();
    clearGeom();

    const volVectorField& m = fvm.lookupObject<volVectorField>(normalName_);
    const volVectorField& xp = fvm.lookupObject<volVectorField>(centreName_);
    const scalarField& V = fvm.V();

    boolList isInterface(fvm.nCells(), false);
    forAll(isInterface, c)
    {
        isInterface[c] =
            mag(m[c]) > wispTol_*Foam::pow(V[c], 2.0/3.0) && mag(m[c]) > 0;
    }
    surfPtr_.reset
    (
        new plicInterfaceSurface
        (
            fvm, m.primitiveField(), xp.primitiveField(), &isInterface
        )
    );
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledDavofInterface::sampledDavofInterface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    normalName_(dict.getOrDefault<word>("normal", "m.davof")),
    centreName_(dict.getOrDefault<word>("centre", "xPlane.davof")),
    wispTol_(dict.getOrDefault<scalar>("wispTol", 0)),
    surfPtr_(nullptr),
    prevTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledDavofInterface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());
    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledDavofInterface::expire()
{
    surfPtr_.clear();
    clearGeom();
    if (prevTimeIndex_ == -1)
    {
        return false;
    }
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledDavofInterface::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField> Foam::sampledDavofInterface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledDavofInterface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledDavofInterface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDavofInterface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledDavofInterface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledDavofInterface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledDavofInterface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledDavofInterface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDavofInterface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledDavofInterface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledDavofInterface::print(Ostream& os, int level) const
{
    os  << "sampledDavofInterface: " << name()
        << " (normal " << normalName_ << ", centre " << centreName_ << ")";
    if (level)
    {
        os << " faces:" << faces().size() << " points:" << points().size();
    }
}


// ************************************************************************* //
