/*---------------------------------------------------------------------------*\
    Production curvature for the Eulerian two-phase solver: the SL
    reconstruction's kappa, delivered exactly as in the semi-Lagrangian solver.

Why this exists
---------------
The first coupled SDPLS study used divGradPsiSnGradAlpha -- kappa from raw FV
derivatives of psi -- and measured a source-curvature feedback loop: the source
modifies psi near the interface from the COMPUTED velocity's strain, kappa
inherits every perturbation instantly, the force follows, the currents grow, the
source amplifies. The droplet drains to nothing by t = 0.1 at every resolution
while the sourceless arm's currents DECAY (so the CSF sign pairing is right; the
coupling path is the defect). The topological band cut-off does not change the
drain rate -- the loop closes at the interface, not in the far field.

The production model of the semi-Lagrangian line (reconstructedCurvature fed by
the quadratic LSQ fit with foot-point re-referencing) evaluates kappa from a
smooth local fit of psi rather than from its raw derivatives, so high-wavenumber
psi perturbations do not pass straight through to the force. That machinery is a
SPATIAL operation on psi -- nothing in it depends on how psi was advected -- so
this solver constructs slReconstruction standalone and reuses the SL face
delivery headers read-only. No SL solver file is modified.

Included from createFields.H BEFORE surfaceTensionForce::New: the
reconstructedCurvature constructor looks `kappa` up, so the fields must be
registered first (same ordering constraint the SL solver documents).

Everything is dict-driven and constructed ONLY when
levelSet.surfaceTensionForce.type == reconstructedCurvature, so every other
configuration pays nothing and behaves exactly as before.
\*---------------------------------------------------------------------------*/

// Is the production curvature selected? Read the same dict the force model
// will read, before constructing it.
const bool useReconstructedCurvature = []
(const fvMesh& m)
{
    const fvSolution& fvSol(m);
    return fvSol.subDict("levelSet")
        .subOrEmptyDict("surfaceTensionForce")
        .getOrDefault<word>("type", "none") == "reconstructedCurvature";
}(mesh);

autoPtr<volScalarField> kappaPtr;
autoPtr<surfaceScalarField> kappaInterfaceFacePtr;
autoPtr<slReconstruction> curvatureRecon;

if (useReconstructedCurvature)
{
    Info<< "Creating field kappa (reconstructed production curvature)\n" << endl;
    kappaPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("kappa", dimless/dimLength, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    kappaInterfaceFacePtr.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "kappaInterfaceFace",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("kappaInterfaceFace", dimless/dimLength, 0.0)
        )
    );

    // The quadratic LSQ fit machinery, standalone. Reads
    // levelSet.semiLagrangian.reconstruction (default
    // quadraticWeightedLeastSquares -- the production type) from fvSolution.
    curvatureRecon = slReconstruction::New(mesh);

    // Initialise kappa from the t=0 psi, so the very first Young-Laplace
    // balance and the force constructed below see a real curvature rather than
    // the zero the field was created with -- the same ordering the SL solver
    // documents for its own initialisation.
    curvatureRecon->update(psi);
    curvatureRecon->meanCurvatureClosestPoint(kappaPtr());
    computeStabilizedFootPointFaceCurvature
    (
        mesh, psi, alpha1, kappaPtr(), curvatureRecon(),
        kappaInterfaceFacePtr()
    );
}
