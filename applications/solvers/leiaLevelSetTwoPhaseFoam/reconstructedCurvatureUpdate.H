/*---------------------------------------------------------------------------*\
    Per-iteration update of the production curvature. Included after alphaEqn.H
    (psi and alpha are current) and before UEqn.H (the momentum predictor reads
    the surface-tension flux built from kappa). Mirrors the SL solver's
    slAlphaEqn.H sequence: fit update -> closest-point cell curvature ->
    stabilized foot-point face delivery. No-op unless the production curvature
    was selected in createFields.
\*---------------------------------------------------------------------------*/

if (useReconstructedCurvature)
{
    curvatureRecon->update(psi);
    curvatureRecon->meanCurvatureClosestPoint(kappaPtr());
    computeStabilizedFootPointFaceCurvature
    (
        mesh, psi, alpha1, kappaPtr(), curvatureRecon(),
        kappaInterfaceFacePtr()
    );
}
