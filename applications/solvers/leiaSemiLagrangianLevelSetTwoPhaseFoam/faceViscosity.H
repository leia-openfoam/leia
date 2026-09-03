/*---------------------------------------------------------------------------*\
    Face dynamic viscosity mu_f for the one-field DNS viscous term (UEqn.H).

    TWO ORTHOGONAL CHOICES, set in fvSolution levelSet.massFlux:

      viscosityFaceAlpha : algebraic | geometric      -- where alpha_f comes from
      viscosityFaceModel : linear | harmonic | normalBlend  -- how mu_f is formed

    They are independent, and an earlier single token conflated them, leaving the
    (algebraic, harmonic) combination untested. Note that (algebraic, linear) is
    exactly OpenFOAM's default face viscosity: because the linear model commutes
    with interpolation,

        < alpha_c mu1 + (1-alpha_c) mu2 >_f = <alpha_c>_f mu1 + (1-<alpha_c>_f) mu2,

    interpolating mu and interpolating alpha then mixing give the same field. That
    coincidence is what allowed the two axes to be collapsed unnoticed; it does NOT
    hold for the harmonic or blended models.

    normalBlend carries no free parameter. Across a layered medium the correct face
    value depends on orientation: NORMAL to the interface the phases are in SERIES
    (harmonic), TANGENT to it they conduct in PARALLEL (arithmetic). The weight is
    therefore the alignment of the interface normal with the face normal,

        w = |n_hat . Sf| / |Sf|,     mu_f = w mu_harm + (1 - w) mu_lin,

    a geometric quantity. Away from the interface alpha_f is 0 or 1 and every model
    returns the same value, so the blend only acts where n_hat is well defined.
\*---------------------------------------------------------------------------*/
{
    const surfaceScalarField mu1f(fvc::interpolate(rho1*mixture.nuModel1().nu()));
    const surfaceScalarField mu2f(fvc::interpolate(rho2*mixture.nuModel2().nu()));

    // --- axis 1: alpha_f -------------------------------------------------
    // ?: cannot mix a field and a tmp<field>, so build it explicitly.
    surfaceScalarField alphafVisc("alphafVisc", alphaf);
    if (!viscosityFaceAlphaGeometric)
    {
        alphafVisc = fvc::interpolate(alpha1);   // linear interpolation of the indicator
    }

    // --- axis 2: the mixing model ----------------------------------------
    const dimensionedScalar muSmall("muSmall", mu1f.dimensions(), VSMALL);
    const surfaceScalarField muLin
    (
        alphafVisc*mu1f + (scalar(1) - alphafVisc)*mu2f
    );

    if (viscosityFaceModelName == "linear")
    {
        muf = muLin;
    }
    else
    {
        const surfaceScalarField muHarm
        (
            scalar(1)
           /max
            (
                alphafVisc/max(mu1f, muSmall)
              + (scalar(1) - alphafVisc)/max(mu2f, muSmall),
                dimensionedScalar("s", inv(mu1f.dimensions()), VSMALL)
            )
        );

        if (viscosityFaceModelName == "harmonic")
        {
            muf = muHarm;
        }
        else    // normalBlend
        {
            // Interface normal at the face, from the phase indicator gradient. Only its
            // DIRECTION is used, and only where the interface is: elsewhere alphafVisc is
            // 0 or 1, muLin == muHarm, and the weight is irrelevant.
            const surfaceVectorField nHatf
            (
                fvc::interpolate(fvc::grad(alpha1))
              / max
                (
                    mag(fvc::interpolate(fvc::grad(alpha1))),
                    dimensionedScalar("g", dimless/dimLength, VSMALL)
                )
            );
            const surfaceScalarField w(mag(nHatf & mesh.Sf())/mesh.magSf());
            muf = w*muHarm + (scalar(1) - w)*muLin;
        }
    }
}
