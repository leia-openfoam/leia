// Capillary-droplet / disturbance-current metrics CSV.
//
// One row per time step is appended to leiaSemiLagrangianLevelSetTwoPhaseFoam.csv
// at the case root by writeDropletMetrics.H:
//
//   TIME      simulated time [s]
//   maxMagU   max_c |U_c| (only a spurious-current metric for U_ref=0)  [m/s]
//   meanMagU  volume-averaged |U|                                       [m/s]
//   maxMagUPrime   max_c |U_c-U_ref|                                    [m/s]
//   meanMagUPrime  volume-averaged |U-U_ref|                            [m/s]
//   l2MagUPrime    volume-weighted L2 of |U-U_ref|                      [m/s]
//                  Report this, not maxMagUPrime alone: L_inf is a single-cell
//                  extremum set by the sub-cell interface position and cannot be
//                  expected to converge under refinement.
//   disturbanceKineticEnergy  integral 0.5*rho*|U-U_ref|^2 dV             [J]
//   phaseVolume/phaseVolumeChange/phaseVolumeRelError
//             Detrixhe-Aslam liquid volume and signed/absolute change  [m3,-,-]
//   centroid{X,Y,Z}/centroidError  liquid centroid and error relative to
//             center_0 + U_ref*t                                      [m,m]
//   zeroSetRadialL2/Linf  error of linearly interpolated psi=0 crossings
//             relative to the analytically translated circle             [m]
//   m2MeanRadialOffset/m2CosCoefficient/m2SinCoefficient/m2Amplitude
//             least-squares polar mode-2 decomposition of r-R             [m]
//   m2Phase  orientation phase, 0.5 atan2(m2Sin,m2Cos)                    [rad]
//   pLaplace  Laplace jump, <p> over the drop (alpha1>1/2) minus <p>    [Pa]
//             over the ambient (alpha1<1/2); mesh-robust (no fixed probe).
//   kErrL2Band    L2 error of the run-time curvature kappa vs 1/R in the
//                 snGrad(alpha) band                                     [1/m]
//   kErrLinfBand  Linf of the same                                       [1/m]
//   minGradPsiBand  min |grad psi| in the band (SDF drift diagnostic)     [-]
//   maxGradPsiBand  max |grad psi| in the band                            [-]
//   gradPsiL2ErrorBand  volume-weighted L2 norm of |grad psi|-1            [-]
//   gradPsiRatioBand  max/min |grad psi| (profile conditioning)            [-]
//   rhoMassResidualL1/Linf  residual of the auxiliary rhoLENT equation
//   rhoMassResidualRelL1  L1 residual scaled by |ddt(rho)|+|div(rhoPhi)|
//   rhoLENTResetL1  mean |rho_aux-rho_DA| before the end-step reset
//   rhoClipL1       volume-mean |clip| applied by massFlux.boundRho   [kg/m3]
//   rhoClipFraction fraction of cells clipped that step                  [-]
//
// The curvature/gradient columns are a diagnostic: they reveal whether the
// curvature error grows BEFORE (curvature-driven) or AFTER (coupling-driven) the
// spurious currents, and whether psi flattens (|grad psi| -> 0, which the |g|^3 in
// the curvature would amplify). The reductions are collective; only the master
// writes the row, so the output is race-free and lands at the case root.
// fvMesh IS-A fvSolution: bind a reference to the mesh's own registered
// fvSolution (as slAdvection does) -- constructing a value here is ambiguous.
const fvSolution& fvSolDrop = mesh;
const dictionary& levelSetDrop = fvSolDrop.subDict("levelSet");
const dictionary& implicitSurfaceDrop = levelSetDrop.subDict("implicitSurface");
const scalar dropletR =
    implicitSurfaceDrop.get<scalar>("radius");
const scalar kappaExactDrop = scalar(mesh.nGeometricD() - 1)/dropletR;   // 1/R (2D), 2/R (3D)
const vector dropletReferenceVelocity =
    levelSetDrop.getOrDefault<vector>("dropletReferenceVelocity", vector::zero);
const point dropletInitialCenter = implicitSurfaceDrop.get<point>("center");

scalar dropletInitialPhaseVolume = 0;
forAll(alpha1, c)
{
    dropletInitialPhaseVolume += mesh.V()[c]*alpha1[c];
}
reduce(dropletInitialPhaseVolume, sumOp<scalar>());

scalar capillaryResidualWorkAccum = 0;

autoPtr<OFstream> dropletFilePtr;
if (Pstream::master())
{
    dropletFilePtr.reset
    (
        new OFstream("leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
    );
    dropletFilePtr().precision(10);
    dropletFilePtr()
        << "TIME,maxMagU,meanMagU,maxMagUPrime,meanMagUPrime,l2MagUPrime,"
        << "disturbanceKineticEnergy,phaseVolume,phaseVolumeChange,"
        << "phaseVolumeRelError,"
        << "centroidX,centroidY,centroidZ,centroidError,"
        << "zeroSetRadialL2,zeroSetRadialLinf,"
        << "m2MeanRadialOffset,m2CosCoefficient,m2SinCoefficient,"
        << "m2Amplitude,m2Phase,"
        << "pLaplace,kErrL2Band,kErrLinfBand,"
        // Popinet Fig. 10: the STANDARD DEVIATION of the curvature along the
        // interface. kErrL2Band is an error against 1/R and so mixes the
        // absorbable uniform offset with the non-absorbable variation; the
        // balanced CSF force consumes ONLY the variation, because a uniform
        // kappa error is snGrad of a cell field and the pressure solve absorbs it
        // exactly. Popinet's std dev falls to 1e-11 (round-off) once the
        // interface reaches its NUMERICAL equilibrium, while the error against
        // 1/R does not, because that equilibrium shape differs from the exact
        // circle at O(h^2). Band = volume-weighted over the snGrad(alpha) band;
        // Faces = |Sf|-weighted over the CSF-active faces, the set the force is
        // actually assembled on. To round-off,
        //   kErrL2Band^2 = (kappaMeanBand - kappaExact)^2 + kappaStdDevBand^2.
        << "kappaMeanBand,kappaStdDevBand,"
        << "kappaMeanActiveFaces,kappaStdDevActiveFaces,"
        // Work done by the part of the capillary force the pressure projection
        // cannot absorb -- the only part that can drive flow. Zero for an exact
        // constant-curvature interface, so on the stationary droplet it is
        // entirely spurious energy input. Compare P/E against the measured
        // envelope growth and against the viscous decay 2*nu*k^2.
        << "capillaryResidualPower,capillaryResidualWork,"
        << "minGradPsiBand,maxGradPsiBand,gradPsiL2ErrorBand,"
        << "gradPsiRatioBand,rhoMassResidualL1,rhoMassResidualLinf,"
        << "rhoMassResidualRelL1,rhoLENTResetL1,rhoClipL1,rhoClipFraction,"
        << "massProjResidualBefore,massProjResidualAfter,"
        << "massProjGlobalDrift,massProjCorrFraction,"
        // WP0 band mode spectrum (bandModeSpectrum.H): normal-profile mode
        // amplitudes A(kh) = |psi_{+k} - 2 psi_0 + psi_{-k}|/4 along the
        // LLS-plane normal, k = 1 (2h), 2 (4h), 4 (8h); band max and
        // volume-weighted L2 per mode; A4h/A8h chain-truncation exclusions.
        << "A2hMaxBand,A2hL2Band,A4hMaxBand,A4hL2Band,"
        << "A8hMaxBand,A8hL2Band,A4hExcluded,A8hExcluded,"
        // WP8.1 (capillaryDriverSplit.H): relative variation rate [1/m] of the
        // delivered face curvature ACROSS the force support (separations along
        // the interface normal) and ALONG the interface, plus the pair counts.
        << "driverAcrossSupportL2,driverAlongInterfaceL2,"
        << "driverAcrossSupportMax,driverAlongInterfaceMax,"
        << "driverNAcross,driverNAlong,"
        // angle between the quadratic-fit normal (used for the split) and
        // the linear-plane normal, over band cells [deg], + fallback count
        << "driverNormalDisagreementMeanDeg,driverNormalDisagreementMaxDeg,"
        << "driverNPlaneNormalFallback" << ','
        // SOLENOIDALITY OF THE ADVECTING FIELD. psi is transported along
        // characteristics with the CELL velocity U (slAlphaEqn.H, identity path
        // when velocityExtension is none), while incompressibility is enforced
        // on the FACE flux phi: div(phi) = 0 to solver tolerance, div(U_cell)
        // is NOT. interFoam/Basilisk advect the colour function with the flux
        // that satisfies discrete continuity exactly; a reinit-free level set
        // has no bounded quantity to absorb the mismatch, so it accumulates in
        // |grad psi|. These columns measure it: the band L2/Linf of div(U),
        // the same for div(phi) as the reference floor, and the implied
        // per-step spurious normal displacement div(U)*h*dt in units of h.
        << "divUBandL2,divUBandLinf,divPhiBandL2,divPhiBandLinf,"
        << "spuriousDispPerStepOverH" << nl;
}
