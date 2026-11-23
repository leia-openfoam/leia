// The initial adaptation of a dynamic mesh, before the first time step. Included by
// the solver after errorCalculation.H, inside `if (dynamicMesh)`.
//
// The case starts on the coarse mesh with psi from leiaSetFields. The mesh refines one
// level per update and updates once per time index. So the time INDEX is advanced
// without a change of the time VALUE, the mesh is updated, and psi is set again from
// the exact profile on the new mesh, until an update changes nothing. The run then
// starts with the exact initial field on the adapted mesh. No function object runs
// here, and nothing is written.
{
    const label maxPasses = 8;
    label pass = 0;
    for (; pass < maxPasses; ++pass)
    {
        runTime.setTime(runTime.value(), runTime.timeIndex() + 1);
        refinementBand->calc(psi);
        meshTransfer->store(psi);
        mesh.update();
        if (!mesh.topoChanging()) { break; }

        #include "rebuildAfterMeshChange.H"

        // The exact initial field on the new mesh replaces the transferred one.
        const vectorField& Cc = mesh.C();
        scalarField& p = psi.primitiveFieldRef();
        forAll(p, c) { p[c] = referenceProfile->value(Cc[c]); }
        psi.correctBoundaryConditions();
        narrowBand->calc();
    }
    if (pass == maxPasses)
    {
        FatalErrorInFunction
            << "the initial adaptation did not settle in " << maxPasses
            << " mesh updates" << exit(FatalError);
    }
    phaseInd->calcPhaseIndicator(alpha, psi);
    Info<< "Initial adaptation: " << pass << " mesh update(s), "
        << returnReduce(mesh.nCells(), sumOp<label>()) << " cells" << nl << endl;
}
