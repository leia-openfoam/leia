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

autoPtr<OFstream> dropletFilePtr;
if (Pstream::master())
{
    dropletFilePtr.reset
    (
        new OFstream("leiaSemiLagrangianLevelSetTwoPhaseFoam.csv")
    );
    dropletFilePtr().precision(10);
    dropletFilePtr()
        << "TIME,maxMagU,meanMagU,maxMagUPrime,meanMagUPrime,"
        << "disturbanceKineticEnergy,phaseVolume,phaseVolumeChange,"
        << "phaseVolumeRelError,"
        << "centroidX,centroidY,centroidZ,centroidError,"
        << "zeroSetRadialL2,zeroSetRadialLinf,"
        << "m2MeanRadialOffset,m2CosCoefficient,m2SinCoefficient,"
        << "m2Amplitude,m2Phase,"
        << "pLaplace,kErrL2Band,kErrLinfBand,"
        << "minGradPsiBand,maxGradPsiBand,gradPsiL2ErrorBand,"
        << "gradPsiRatioBand,rhoMassResidualL1,rhoMassResidualLinf,"
        << "rhoMassResidualRelL1,rhoLENTResetL1,"
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
        << "driverNAcross,driverNAlong" << nl;
}
