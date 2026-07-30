/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM and is distributed under the GPL v3 or later.
\*---------------------------------------------------------------------------*/

#include "integralConormalSurfaceTension.H"
#include "isoGeometricCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

defineTypeNameAndDebug(integralConormalSurfaceTension, false);
addToRunTimeSelectionTable
(
    surfaceTensionForce,
    integralConormalSurfaceTension,
    Mesh
);

integralConormalSurfaceTension::integralConormalSurfaceTension
(
    const fvMesh& mesh
)
:
    surfaceTensionForce(mesh),
    fvSolutionDict_(mesh_),
    levelSetDict_(fvSolutionDict_.subDict("levelSet")),
    surfTensionDict_(levelSetDict_.subDict("surfaceTensionForce")),
    field_
    (
        mesh_.lookupObject<volScalarField>
        (
            surfTensionDict_.getOrDefault<word>("field", "psi")
        )
    ),
    isoValue_(surfTensionDict_.getOrDefault<scalar>("isoValue", 0.0)),
    sgn_(surfTensionDict_.getOrDefault<scalar>("sign", 1.0))
{
    if (surfTensionDict_.found("delivery"))
    {
        const word legacyDelivery(surfTensionDict_.get<word>("delivery"));
        if (legacyDelivery == "directCell")
        {
            FatalIOErrorInFunction(surfTensionDict_)
                << "The directCell production bypass has been removed. "
                << "integralConormalSurfaceTension now always returns the "
                << "integrated scalar face force flux G_sigma,f. Remove the "
                << "delivery entry."
                << exit(FatalIOError);
        }
        WarningInFunction
            << "Ignoring obsolete delivery entry '" << legacyDelivery
            << "': the model now always returns G_sigma,f directly." << endl;
    }

    Info<< "integralConormalSurfaceTension: curvature-free shared-face "
        << "traction on {" << field_.name() << "=" << isoValue_ << "}; "
        << "pseudo-2D internal-face prototype; delivery=integratedFaceFlux"
        << endl;
}


tmp<surfaceScalarField>
integralConormalSurfaceTension::calcFaceSurfaceTensionForceFlux() const
{
    tmp<volVectorField> tfCell = conservativeCellSurfaceTensionForce();

    tmp<surfaceScalarField> tGSigma
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GSigmaConormalFromCell",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            linearInterpolate(tfCell()) & mesh_.Sf()
        )
    );

    Info<< "integralConormalSurfaceTension: max|Gsigma|="
        << gMax(mag(tGSigma().primitiveField())) << " N/m" << endl;

    return tGSigma;
}


tmp<volVectorField>
integralConormalSurfaceTension::conservativeCellSurfaceTensionForce() const
{
    surfaceVectorField traction
    (
        IOobject
        (
            "integralConormalTraction",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("integralConormalTraction", dimForce, Zero)
    );

    const label nCut = computeIsoConormalFaceTraction
    (
        mesh_, field_, isoValue_, sgn_, sigma_.value(), traction
    );

    tmp<volVectorField> tf = fvc::surfaceIntegrate(traction);
    tf.ref().rename("fSigmaIntegralConormal");

    const vectorField& fI = tf().primitiveField();
    const scalarField& V = mesh_.V();
    vector netForce(Zero);
    scalar maxForceDensity = 0.0;
    forAll(fI, c)
    {
        netForce += V[c]*fI[c];
        maxForceDensity = max(maxForceDensity, mag(fI[c]));
    }
    reduce(netForce, sumOp<vector>());
    reduce(maxForceDensity, maxOp<scalar>());

    Info<< "integralConormalSurfaceTension: cut internal faces="
        << returnReduce(nCut, sumOp<label>())
        << ", |fCell|_max=" << maxForceDensity << " N/m^3"
        << ", net force=" << netForce << " N" << endl;

    return tf;
}

} // End namespace Foam

// ************************************************************************* //
