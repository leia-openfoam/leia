// The exact reference fields psi0 and alpha0 on the PRESENT mesh, for the error
// report of a run with a dynamic mesh. Included by the solver on a report step.
// The phase indicator reads the registered narrow band, which belongs to psi.
// So psi is replaced by the reference for the evaluation of alpha0 and restored
// bit for bit afterwards, and the band is computed again for psi.
{
    const volScalarField psiSave("psiSave", psi);
    const vectorField& Cc = mesh.C();
    scalarField& p = psi.primitiveFieldRef();
    forAll(p, c)
    {
        p[c] = referenceProfile->value(Cc[c] - referenceShift);
    }
    psi.correctBoundaryConditions();
    narrowBand->calc();
    phaseInd->calcPhaseIndicator(alpha0, psi);
    psi0 == psi;
    psi == psiSave;
    narrowBand->calc();
}
