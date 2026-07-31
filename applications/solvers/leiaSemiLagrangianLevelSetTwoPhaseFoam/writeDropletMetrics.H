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
    forAll(bnd, c)
    {
        if (!bnd[c]) { continue; }
        const scalar e = mag(kappa[c] - kappaExactDrop);
        const scalar gMag = magGradPsi[c];
        bV += V[c]; bs2 += V[c]*e*e; kLinf = max(kLinf, e);
        gMinBand = min(gMinBand, gMag);
        gMaxBand = max(gMaxBand, gMag);
        gErrS2 += V[c]*sqr(gMag - 1);
    }
    reduce(bV, sumOp<scalar>());   reduce(bs2, sumOp<scalar>());
    reduce(kLinf, maxOp<scalar>()); reduce(gMinBand, minOp<scalar>());
    reduce(gMaxBand, maxOp<scalar>()); reduce(gErrS2, sumOp<scalar>());
    const scalar kErrL2 = (bV > 0) ? Foam::sqrt(bs2/bV) : 0;
    const scalar gMin = (bV > 0) ? gMinBand : 0;
    const scalar gMax = (bV > 0) ? gMaxBand : 0;
    const scalar gradPsiL2Error = (bV > 0) ? Foam::sqrt(gErrS2/bV) : 0;
    const scalar gradPsiRatio =
        (bV > 0) ? gMax/Foam::max(gMin, VSMALL) : 0;

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
            << gMin << ',' << gMax << ',' << gradPsiL2Error << ','
            << gradPsiRatio << ',' << rhoMassResidualL1 << ','
            << rhoMassResidualLinf << ',' << rhoMassResidualRelL1 << ','
            << rhoLENTResetL1 << nl;
    }
}
