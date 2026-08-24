{
    // Spurious-current amplitudes: L_inf and volume-averaged |U| (both global).
    const volScalarField magU(mag(U));
    const scalar maxMagU  = gMax(magU.primitiveField());
    const scalar meanMagU = magU.weightedAverage(mesh.V()).value();

    // Translation-aware disturbance metrics. For a translating equilibrium,
    // U_ref is the prescribed bulk velocity and U' = U-U_ref is the physically
    // relevant parasitic current. The cell loop keeps reductions explicit.
    const scalarField& V = mesh.V();
    scalar totalV = 0;
    scalar maxMagUPrime = 0;
    scalar integralMagUPrime = 0;
    scalar disturbanceKineticEnergy = 0;
    scalar phaseVolume = 0;
    vector phaseFirstMoment(vector::zero);
    const vectorField& cellCentres = mesh.C();
    forAll(U, c)
    {
        const vector UPrime = U[c] - dropletReferenceVelocity;
        const scalar magUPrime = mag(UPrime);
        totalV += V[c];
        maxMagUPrime = max(maxMagUPrime, magUPrime);
        integralMagUPrime += V[c]*magUPrime;
        disturbanceKineticEnergy += 0.5*V[c]*rho[c]*magSqr(UPrime);
        phaseVolume += V[c]*alpha1[c];
        phaseFirstMoment += V[c]*alpha1[c]*cellCentres[c];
    }
    reduce(totalV, sumOp<scalar>());
    reduce(maxMagUPrime, maxOp<scalar>());
    reduce(integralMagUPrime, sumOp<scalar>());
    reduce(disturbanceKineticEnergy, sumOp<scalar>());
    reduce(phaseVolume, sumOp<scalar>());
    reduce(phaseFirstMoment, sumOp<vector>());
    const scalar meanMagUPrime =
        (totalV > VSMALL) ? integralMagUPrime/totalV : 0;
    const scalar phaseVolumeChange =
        (dropletInitialPhaseVolume > VSMALL)
      ? phaseVolume/dropletInitialPhaseVolume - 1 : 0;
    const scalar phaseVolumeRelError = mag(phaseVolumeChange);
    const point phaseCentroid =
        (phaseVolume > VSMALL)
      ? point(phaseFirstMoment/phaseVolume) : point::zero;
    const point expectedCentroid =
        dropletInitialCenter + dropletReferenceVelocity*runTime.value();
    const scalar centroidError = mag(phaseCentroid - expectedCentroid);

    // Geometric zero-set error for the translating-circle gate. Locate each
    // internal owner-neighbour psi sign crossing by linear interpolation, then
    // compare its radius about center_0 + U_ref*t with the exact radius. This
    // follows the moving reference, unlike the legacy stationary error object.
    const labelUList& fown = mesh.owner();
    const labelUList& fnei = mesh.neighbour();
    const vectorField& centres = mesh.C();
    scalar zeroSetCount = 0;
    scalar zeroSetRadialS2 = 0;
    scalar zeroSetRadialLinf = 0;
    scalar m2SumCos = 0, m2SumSin = 0;
    scalar m2SumCos2 = 0, m2SumCosSin = 0, m2SumSin2 = 0;
    scalar m2SumDr = 0, m2SumDrCos = 0, m2SumDrSin = 0;
    forAll(fnei, f)
    {
        const label own = fown[f];
        const label nei = fnei[f];
        const scalar psiOwn = psi[own];
        const scalar psiNei = psi[nei];
        if (psiOwn*psiNei > 0 || mag(psiOwn - psiNei) < VSMALL)
        {
            continue;
        }
        const scalar lambda = psiOwn/(psiOwn - psiNei);
        const point crossing = centres[own] + lambda*(centres[nei] - centres[own]);
        const vector radial = crossing - expectedCentroid;
        const scalar crossingRadius = mag(radial);
        const scalar signedRadialError = crossingRadius - dropletR;
        const scalar radialError = mag(signedRadialError);
        zeroSetCount += 1;
        zeroSetRadialS2 += sqr(radialError);
        zeroSetRadialLinf = max(zeroSetRadialLinf, radialError);
        if (crossingRadius > VSMALL)
        {
            const scalar cosTheta = radial.x()/crossingRadius;
            const scalar sinTheta = radial.y()/crossingRadius;
            const scalar cosMode2 = sqr(cosTheta) - sqr(sinTheta);
            const scalar sinMode2 = 2*cosTheta*sinTheta;
            m2SumCos += cosMode2;
            m2SumSin += sinMode2;
            m2SumCos2 += sqr(cosMode2);
            m2SumCosSin += cosMode2*sinMode2;
            m2SumSin2 += sqr(sinMode2);
            m2SumDr += signedRadialError;
            m2SumDrCos += signedRadialError*cosMode2;
            m2SumDrSin += signedRadialError*sinMode2;
        }
    }
    reduce(zeroSetCount, sumOp<scalar>());
    reduce(zeroSetRadialS2, sumOp<scalar>());
    reduce(zeroSetRadialLinf, maxOp<scalar>());
    reduce(m2SumCos, sumOp<scalar>());
    reduce(m2SumSin, sumOp<scalar>());
    reduce(m2SumCos2, sumOp<scalar>());
    reduce(m2SumCosSin, sumOp<scalar>());
    reduce(m2SumSin2, sumOp<scalar>());
    reduce(m2SumDr, sumOp<scalar>());
    reduce(m2SumDrCos, sumOp<scalar>());
    reduce(m2SumDrSin, sumOp<scalar>());
    const scalar zeroSetRadialL2 =
        (zeroSetCount > 0) ? Foam::sqrt(zeroSetRadialS2/zeroSetCount) : 0;

    // Least-squares r-R = mean + a*cos(2 theta) + b*sin(2 theta).  Removing
    // the fitted mean makes the modal coefficients insensitive to the
    // nonuniform angular density of Cartesian face crossings.
    scalar m2MeanRadialOffset = 0;
    scalar m2CosCoefficient = 0;
    scalar m2SinCoefficient = 0;
    if (zeroSetCount > 3)
    {
        const scalar invCount = 1/zeroSetCount;
        const scalar cc = m2SumCos2 - sqr(m2SumCos)*invCount;
        const scalar ss = m2SumSin2 - sqr(m2SumSin)*invCount;
        const scalar cs = m2SumCosSin - m2SumCos*m2SumSin*invCount;
        const scalar dc = m2SumDrCos - m2SumDr*m2SumCos*invCount;
        const scalar ds = m2SumDrSin - m2SumDr*m2SumSin*invCount;
        const scalar determinant = cc*ss - sqr(cs);
        if (mag(determinant) > VSMALL)
        {
            m2CosCoefficient = (dc*ss - ds*cs)/determinant;
            m2SinCoefficient = (ds*cc - dc*cs)/determinant;
        }
        m2MeanRadialOffset =
            (m2SumDr - m2CosCoefficient*m2SumCos
             - m2SinCoefficient*m2SumSin)*invCount;
    }
    const scalar m2Amplitude =
        Foam::sqrt(sqr(m2CosCoefficient) + sqr(m2SinCoefficient));
    const scalar m2Phase =
        0.5*Foam::atan2(m2SinCoefficient, m2CosCoefficient);

    // Laplace jump: alpha-partitioned mean pressure, drop minus ambient,
    // volume-weighted so it is independent of the mesh and of any probe location.
    //
    // The two means are taken over the PURE phases (alpha1 > 1-eps and
    // alpha1 < eps), NOT over alpha1 >/< 1/2. Partitioning at 1/2 puts the whole
    // interface transition into the "drop" average, and those cells sit partway
    // across the pressure jump, so the reported jump is biased LOW by a few
    // percent -- an artefact of the probe, not of the solver.
    //
    // Measured on a frozen circle at N=128 with EXACT constant curvature
    // (kappa = 1000 exactly, max|U| ~ 3e-11 m/s, i.e. nothing physical moving):
    //   alpha1 >/< 1/2          : 70.4456 Pa  (-3.15% vs sigma/R = 72.74 Pa)
    //   alpha1 > 0.99 / < 0.01  : 72.7394 Pa  (-0.00%)
    //   pure phases             : 72.7400 Pa  (-0.00%)
    // The discrete Laplace jump is therefore EXACT to five significant figures;
    // the old partition was measuring the width of the alpha transition. It was
    // also non-monotone in solver quality -- at N=64 a computed curvature scored
    // -1.63% against -2.74% for an exact one, because the curvature error
    // partially cancelled the probe bias, and improving the curvature 27x made
    // the reported jump worse. Any earlier reported Laplace-jump convergence
    // (e.g. 4.2% -> 1.5% over a mesh ladder) is that artefact and not a solver
    // result.
    const scalar alphaPure = 1e-2;
    scalar vIn = 0, vOut = 0, pvIn = 0, pvOut = 0;
    forAll(alpha1, c)
    {
        if (alpha1[c] > 1 - alphaPure) { vIn  += V[c]; pvIn  += V[c]*p[c]; }
        else if (alpha1[c] < alphaPure) { vOut += V[c]; pvOut += V[c]*p[c]; }
        // cells inside the transition contribute to neither mean
    }
    reduce(vIn,   sumOp<scalar>());
    reduce(vOut,  sumOp<scalar>());
    reduce(pvIn,  sumOp<scalar>());
    reduce(pvOut, sumOp<scalar>());
    const scalar pIn      = (vIn  > SMALL) ? pvIn/vIn   : 0;
    const scalar pOut     = (vOut > SMALL) ? pvOut/vOut : 0;
    const scalar pLaplace = pIn - pOut;

    // --- curvature-error diagnostic in the snGrad(alpha) band --------------
    // kappa holds this step's curvature (slAlphaEqn). Compare with the exact 1/R
    // and record the min |grad psi| in the same band, to time-resolve whether the
    // curvature degrades before or after the spurious currents grow.
    const surfaceScalarField snA(fvc::snGrad(alpha1));
    const scalarField& snIi = snA.primitiveField();
    const scalar snMx = max(gMax(mag(snIi)), VSMALL);
    const volScalarField magGradPsi(mag(fvc::grad(psi)));
    boolList bnd(mesh.nCells(), false);
    forAll(snIi, f)
    {
        if (mag(snIi[f]) > 1e-8*snMx) { bnd[fown[f]] = true; bnd[fnei[f]] = true; }
    }
    scalar bV = 0, bs2 = 0, kLinf = 0;
    scalar gMinBand = GREAT, gMaxBand = -GREAT, gErrS2 = 0;
    // Popinet (JCP 228:5838) Fig. 10 quantity: the STANDARD DEVIATION of the
    // curvature along the interface, not its error against 1/R. The balanced-force
    // argument says a spatially UNIFORM curvature error is snGrad-free and the
    // pressure solve absorbs it exactly, so only the variation reaches momentum;
    // kErrL2Band mixes the two, and at numerical equilibrium Popinet's std dev
    // falls to round-off (1e-11) while the error against 1/R does NOT (his
    // equilibrium shape differs from the exact circle at O(h^2)). Accumulated as
    // raw moments in this same loop so it costs nothing.
    scalar sK = 0, sK2 = 0;
    forAll(bnd, c)
    {
        if (!bnd[c]) { continue; }
        const scalar e = mag(kappa[c] - kappaExactDrop);
        const scalar gMag = magGradPsi[c];
        bV += V[c]; bs2 += V[c]*e*e; kLinf = max(kLinf, e);
        sK += V[c]*kappa[c]; sK2 += V[c]*kappa[c]*kappa[c];
        gMinBand = min(gMinBand, gMag);
        gMaxBand = max(gMaxBand, gMag);
        gErrS2 += V[c]*sqr(gMag - 1);
    }
    reduce(bV, sumOp<scalar>());   reduce(bs2, sumOp<scalar>());
    reduce(kLinf, maxOp<scalar>()); reduce(gMinBand, minOp<scalar>());
    reduce(gMaxBand, maxOp<scalar>()); reduce(gErrS2, sumOp<scalar>());
    reduce(sK, sumOp<scalar>());   reduce(sK2, sumOp<scalar>());
    const scalar kErrL2 = (bV > 0) ? Foam::sqrt(bs2/bV) : 0;
    // mean and standard deviation of the CELL curvature over the band. The mean
    // minus the exact value is the ABSORBABLE part of the error; the standard
    // deviation is the part the balanced CSF force actually delivers. To
    // round-off, kErrL2^2 = (mean - exact)^2 + stdDev^2.
    const scalar kMeanBand = (bV > 0) ? sK/bV : 0;
    const scalar kStdBand =
        (bV > 0) ? Foam::sqrt(max(sK2/bV - kMeanBand*kMeanBand, scalar(0))) : 0;

    // The same statistic on the FACES the force is assembled on, |Sf|-weighted
    // over the CSF-active set -- the closest available analogue of "along the
    // interface", since the force consumes face values.
    scalar kStdFaces = 0, kMeanFaces = 0;
    {
        const surfaceScalarField kappaF(fvc::interpolate(kappa));
        const scalarField& kf = kappaF.primitiveField();
        const scalarField& aSf = mesh.magSf().primitiveField();
        scalar wSum = 0, wK = 0, wK2 = 0;
        forAll(snIi, f)
        {
            if (mag(snIi[f]) <= 1e-8*snMx) { continue; }
            wSum += aSf[f]; wK += aSf[f]*kf[f]; wK2 += aSf[f]*kf[f]*kf[f];
        }
        reduce(wSum, sumOp<scalar>()); reduce(wK, sumOp<scalar>());
        reduce(wK2, sumOp<scalar>());
        kMeanFaces = (wSum > 0) ? wK/wSum : 0;
        kStdFaces = (wSum > 0)
            ? Foam::sqrt(max(wK2/wSum - kMeanFaces*kMeanFaces, scalar(0))) : 0;
    }
    const scalar gMin = (bV > 0) ? gMinBand : 0;
    const scalar gMax = (bV > 0) ? gMaxBand : 0;
    const scalar gradPsiL2Error = (bV > 0) ? Foam::sqrt(gErrS2/bV) : 0;
    const scalar gradPsiRatio =
        (bV > 0) ? gMax/Foam::max(gMin, VSMALL) : 0;

    // --- WP0: band normal-profile mode spectrum (bandModeSpectrum.H) -------
    // A{2h,4h,8h} amplitudes along the LLS-plane normal; the leading
    // indicator candidates for the parasitic feedback (minGradPsiBand lags).
    const bandModeSpectrumNorms modeSpec = computeBandModeSpectrum
    (
        mesh, psi, mesh.lookupObject<volScalarField>("NarrowBand")
    );

    // --- WP8.1: where the delivered curvature varies (across the support
    // along the interface normal, or along the interface). Re-evaluates the
    // force model's face flux with the current fields and recovers kappa_f
    // from it, so the measure is model-agnostic.
    const capillaryDriverSplitNorms driverSplit = computeCapillaryDriverSplit
    (
        mesh, psi, alpha1, fSigma->faceSurfaceTensionForceFlux()(),
        slAdv->reconstruction()
    );

    // ---- solenoidality of the ADVECTING field (see createDropletMetricsFile.H)
    // div(phi) is zero to solver tolerance by construction; div(U_cell) from
    // fvc::reconstruct is not, and psi is transported with U_cell. Both are
    // measured on the same band (|psi| < 2h) so the ratio is meaningful, and
    // the last column converts div(U) into the per-step spurious normal
    // displacement it implies, in units of the cell size.
    scalar divUBandL2 = 0, divUBandLinf = 0;
    scalar divPhiBandL2 = 0, divPhiBandLinf = 0;
    scalar spuriousDispPerStepOverH = 0;
    {
        // fvc::div(volVectorField) would require a div(U) scheme entry in
        // every case; the flux form needs none (fvc::div of a
        // surfaceScalarField is surfaceIntegrate) and is the DIRECTLY
        // comparable quantity: the same discrete divergence operator applied
        // to the interpolated cell velocity and to the pressure-corrected
        // face flux.
        const volScalarField divU
        (
            fvc::div(fvc::interpolate(U) & mesh.Sf())
        );
        const volScalarField divPhi(fvc::div(phi));
        const scalarField& dU = divU.primitiveField();
        const scalarField& dP = divPhi.primitiveField();
        const scalarField& psiIn = psi.primitiveField();
        const scalarField& Vc = mesh.V();

        // Use the SAME band the other metrics use (the registered NarrowBand
        // field) rather than a |psi| threshold against a guessed cell size: on
        // a 2D mesh the cell volume is h^2*Lz, so no root of V is the cell
        // size, and a first attempt with pow(V, 1/nGeometricD) selected an
        // EMPTY band and silently reported zeros.
        const scalarField& nb =
            mesh.lookupObject<volScalarField>("NarrowBand").primitiveField();
        const scalar hCell = Foam::sqrt(gAverage(mesh.magSf()));

        scalar sU = 0, sP = 0, vSum = 0;
        forAll(psiIn, c)
        {
            if (nb[c] < 0.5) { continue; }
            sU += Vc[c]*sqr(dU[c]);
            sP += Vc[c]*sqr(dP[c]);
            vSum += Vc[c];
            divUBandLinf = max(divUBandLinf, mag(dU[c]));
            divPhiBandLinf = max(divPhiBandLinf, mag(dP[c]));
        }
        reduce(sU, sumOp<scalar>());
        reduce(sP, sumOp<scalar>());
        reduce(vSum, sumOp<scalar>());
        reduce(divUBandLinf, maxOp<scalar>());
        reduce(divPhiBandLinf, maxOp<scalar>());
        if (vSum > VSMALL)
        {
            divUBandL2 = Foam::sqrt(sU/vSum);
            divPhiBandL2 = Foam::sqrt(sP/vSum);
        }
        spuriousDispPerStepOverH =
            divUBandL2*hCell*runTime.deltaTValue()/max(hCell, VSMALL);
    }

    if (Pstream::master())
    {
        dropletFilePtr()
            << runTime.value() << ',' << maxMagU << ',' << meanMagU << ','
            << maxMagUPrime << ',' << meanMagUPrime << ','
            << disturbanceKineticEnergy << ',' << phaseVolume << ','
            << phaseVolumeChange << ',' << phaseVolumeRelError << ','
            << phaseCentroid.x() << ','
            << phaseCentroid.y() << ',' << phaseCentroid.z() << ','
            << centroidError << ',' << zeroSetRadialL2 << ','
            << zeroSetRadialLinf << ','
            << m2MeanRadialOffset << ',' << m2CosCoefficient << ','
            << m2SinCoefficient << ',' << m2Amplitude << ','
            << m2Phase << ','
            << pLaplace << ',' << kErrL2 << ',' << kLinf << ','
            << kMeanBand << ',' << kStdBand << ','
            << kMeanFaces << ',' << kStdFaces << ','
            << gMin << ',' << gMax << ',' << gradPsiL2Error << ','
            << gradPsiRatio << ',' << rhoMassResidualL1 << ','
            << rhoMassResidualLinf << ',' << rhoMassResidualRelL1 << ','
            << rhoLENTResetL1 << ','
            << modeSpec.aMax[0] << ',' << modeSpec.aL2[0] << ','
            << modeSpec.aMax[1] << ',' << modeSpec.aL2[1] << ','
            << modeSpec.aMax[2] << ',' << modeSpec.aL2[2] << ','
            << modeSpec.nExcluded[1] << ',' << modeSpec.nExcluded[2] << ','
            << driverSplit.acrossSupportL2 << ','
            << driverSplit.alongInterfaceL2 << ','
            << driverSplit.acrossSupportMax << ','
            << driverSplit.alongInterfaceMax << ','
            << driverSplit.nAcross << ',' << driverSplit.nAlong << ','
            << driverSplit.normalDisagreementMeanDeg << ','
            << driverSplit.normalDisagreementMaxDeg << ','
            << driverSplit.nPlaneNormalFallback << ','
            << divUBandL2 << ',' << divUBandLinf << ','
            << divPhiBandL2 << ',' << divPhiBandLinf << ','
            << spuriousDispPerStepOverH << nl;
    }
}
