// After a topology change of the mesh (refinement, unrefinement, redistribution).
// psi holds psi^n with the values that the mesh gave it. Included by the solver
// inside `if (dynamicMesh && mesh.topoChanging())`.
{
    // 1. psi on the new mesh, second order.
    meshTransfer->apply(psi);

    // 2. Every object that is built for one mesh is constructed again. Each
    //    autoPtr is cleared FIRST: two objects must not register a field of the
    //    same name. The stencil object of OpenFOAM is deleted by the topology
    //    change, so the old reconstruction must not be used again.
    meshTransfer.clear();
    slAdv.clear();
    narrowBand.clear();
    phaseInd.clear();
    velocityModel.clear();

    velocityModel = Foam::velocityModel::New(mesh);
    // The prescribed flux and velocity, exact on the new mesh (not mapped values).
    // The step rescales them to t^{n+1} from U0 and phi0 as on a static mesh.
    velocityModel->setVolumetricFlux(phi);
    velocityModel->setVelocity(U);
    U0 == U;
    Uold == U;
    if (phi0Ptr.valid())
    {
        phi0Ptr.reset(new surfaceScalarField("phi0", phi));
    }
    if (traceFluxExtension)
    {
        velExtPtr.clear();
        velExtPtr = velocityExtension::New(mesh);
        velExtPtr->correct();
    }
    if (traceFromFlux)
    {
        UtracePtr.reset(new volVectorField("Utrace", fvc::reconstruct(phi)));
        UtraceOldPtr.reset(new volVectorField("UtraceOld", UtracePtr()));
    }

    narrowBand = Foam::narrowBand::New(mesh, psi);
    narrowBand->calc();
    phaseInd = phaseIndicator::New(mesh);
    slAdv = slAdvection::New(mesh);
    meshTransfer = slMeshTransfer::New(mesh);

    Info<< "Mesh changed: " << returnReduce(mesh.nCells(), sumOp<label>())
        << " cells; the semi-Lagrangian objects are built again" << nl << endl;
}
