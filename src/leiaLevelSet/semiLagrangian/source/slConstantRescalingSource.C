/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slConstantRescalingSource.H"
#include "addToRunTimeSelectionTable.H"
#include "slReconstruction.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
    defineTypeNameAndDebug(slConstantRescalingSource, 0);
    addToRunTimeSelectionTable(slSource, slConstantRescalingSource, Mesh);
}

Foam::slConstantRescalingSource::slConstantRescalingSource(const fvMesh& mesh)
:
    slSource(mesh),
    theta_(dict_.getOrDefault<scalar>("theta", 0.1))
{
    if (!(theta_ > 0))
    {
        FatalIOErrorInFunction(dict_)
            << "semiLagrangian.source constantRescaling needs theta > 0; got "
            << theta_ << exit(FatalIOError);
    }
}

void Foam::slConstantRescalingSource::apply
(
    volScalarField& psi,
    const volVectorField&,
    slReconstruction& geometryFit,
    const scalar
)
{
    if (!mesh_.foundObject<volScalarField>("NarrowBand"))
    {
        FatalErrorInFunction
            << "semiLagrangian.source constantRescaling needs the registered"
            << " field NarrowBand: set levelSet.narrowBand.type to signChange."
            << exit(FatalError);
    }
    const volScalarField& nb =
        mesh_.lookupObject<volScalarField>("NarrowBand");

    psi.correctBoundaryConditions();
    geometryFit.update(psi);

    scalar sumG = 0;
    label n = 0;
    vector g(Zero);
    symmTensor H(Zero);
    forAll(nb, celli)
    {
        if (nb[celli] > 0.5 && geometryFit.fitDerivatives(celli, g, H) >= 1)
        {
            sumG += mag(g);
            ++n;
        }
    }
    reduce(sumG, sumOp<scalar>());
    reduce(n, sumOp<label>());
    if (n == 0)
    {
        return;
    }
    const scalar G = sumG/scalar(n);

    const word refName("slConstantRescalingReference");
    if (!mesh_.foundObject<uniformDimensionedScalarField>(refName))
    {
        regIOobject::store
        (
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    refName,
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                dimensionedScalar(refName, dimless, G)
            )
        );
        Info<< "slConstantRescalingSource: reference G0 = " << G
            << " (mean |grad psi| of " << n << " interface cells), theta "
            << theta_ << endl;
        return;
    }
    const scalar G0 =
        mesh_.lookupObject<uniformDimensionedScalarField>(refName).value();

    const scalar dev = mag(G/G0 - 1);
    if (dev <= theta_)
    {
        return;
    }
    const scalar factor = G0/G;
    psi.primitiveFieldRef() *= factor;
    psi.correctBoundaryConditions();
    Info<< "slConstantRescalingSource: event at time index "
        << mesh_.time().timeIndex() << ", G " << G << ", G0 " << G0
        << ", |G/G0 - 1| " << dev << ", psi multiplied by " << factor
        << " in every cell" << endl;
}
