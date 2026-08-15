// Registered mean-curvature field (1/length), filled each step from the
// semi-Lagrangian quadratic reconstruction (symbolic gradient + Hessian, evaluated
// at the CELL CENTRE with NO foot-point normal extension) and consumed by the
// reconstructedCurvature surfaceTensionForce model. Also filled once at t=0 during
// initialisation (see the solver .C) so the initial Young-Laplace balance uses it
// and the curvature on the clean signed-distance field can be inspected.
// Created BEFORE createFields.H so the force's constructor can look it up.
Info<< "Creating field kappa\n" << endl;
volScalarField kappa
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
    "zeroGradient"
);

// Face-native interface curvature.  connectedInterfaceCurvature fills this
// directly by extending one value per connected zero-surface element along its
// normal ray.  reconstructedCurvature consumes it without cell-to-face
// interpolation when faceInterpolation=connectedInterface.
Info<< "Creating field kappaInterfaceFace\n" << endl;
surfaceScalarField kappaInterfaceFace
(
    IOobject
    (
        "kappaInterfaceFace",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaInterfaceFace", dimless/dimLength, 0.0)
);

// Stabilized foot-point FACE curvature (stabilizedFootPointFaceCurvature.H):
// the arithmetic kappa_f re-referenced to the interface with the stabilized
// foot-point offset of each active face centre. Filled each step when
// curvatureExtension type == stabilizedFootPointFace; consumed by
// reconstructedCurvature via faceCurvatureSource registered /
// faceCurvatureField kappaStableFootFace.
Info<< "Creating field kappaStableFootFace\n" << endl;
surfaceScalarField kappaStableFootFace
(
    IOobject
    (
        "kappaStableFootFace",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaStableFootFace", dimless/dimLength, 0.0)
);

// Semi-Lagrangian level-set advection + symbolic curvature evaluator. Selects the
// reconstruction/correction from fvSolution levelSet.semiLagrangian.
autoPtr<slAdvection> slAdv = slAdvection::New(mesh);

// Optional harmonic (Laplace) smoothing of the mean curvature (see
// harmonicMeanCurvatureExtensionEqn.H): pin kappa at the interface cells and solve
// laplacian(kappa)=0 to replace the run-time oscillating band curvature with a
// smooth harmonic extension. Toggled by fvSolution:
//   levelSet { curvatureExtension { type harmonicLaplace; } }   // "none" disables
word curvatureExtensionType("none");
{
    const fvSolution& fvSol = mesh;   // fvMesh is-a fvSolution
    const dictionary& lsDict = fvSol.subDict("levelSet");
    curvatureExtensionType =
        lsDict.subOrEmptyDict("curvatureExtension").getOrDefault<word>("type", "none");
}
const bool harmonicCurvatureExtension = (curvatureExtensionType == "harmonicLaplace");
// Foot-point height-function curvature (footPointCurvature.H): kappa from a
// local parabola (height-function) fit of the interface feet of the linear
// LSQ planes, delivered normal-constant by uniform nearest-foot assignment --
// the unstructured-compatible height function (geometry-only, integral,
// normal-constant), the equilibrium-restoring option.
const bool footPointCurvatureExtension =
    (curvatureExtensionType == "footPointHeightFunction");
// Closest-point Newton foot on the cell's OWN quadratic (type
// closestPointNewton): symbolic kappa evaluated at the TRUE closest point of
// {P=0}, found by a trust-region KKT Newton (fit Hessian) instead of the
// gradient-projection flow of the original meanCurvature extension; per-cell
// fallback to the cell-centre value on guard failure.
const bool closestPointCurvatureExtension =
    (curvatureExtensionType == "closestPointNewton");
// Diagnostic null-space control: retain the CSF-support-weighted mean of the
// reconstructed curvature, but make it spatially constant.  This distinguishes
// curvature fluctuations from mean-curvature bias using the complete pressure
// and momentum solve.
const bool interfaceMeanCurvatureExtension =
    (curvatureExtensionType == "interfaceMean");
// Connected zero-surface reconstruction: one curvature per conforming 2-D
// interface segment, tangential Helmholtz regularisation on the segment chain,
// and normal-ray delivery directly to kappaInterfaceFace for the CSF force.
const bool connectedInterfaceCurvatureExtension =
    (curvatureExtensionType == "connectedInterface");
// type fvm: the selected surfaceTensionForce model computes its own curvature
// with FVM operators (e.g. divGradPsiSnGradAlpha) -- SKIP the symbolic
// quadratic-fit kappa fill entirely. Required for non-quadratic
// reconstructions (linearTaylor has no Hessian; the symbolic fill would
// FatalError with NotImplemented). kappa stays 0 (the kErr CSV columns are
// then meaningless -- ignore them for such runs).
const bool fvmCurvature = (curvatureExtensionType == "fvm");
// Stabilized foot-point FACE delivery: keep the raw contour-referenced cell
// fill (the default meanCurvatureNoExtension branch) and re-reference the
// interpolated FACE value to the interface with the stabilized foot-point
// offset (stabilizedFootPointFaceCurvature.H). The measured face-centered
// curvature gate: h^2.0 active-face L2 vs h^1.1 for the plain arithmetic
// delivery on the static exact-SDF circle.
const bool stabilizedFootPointFaceExtension =
    (curvatureExtensionType == "stabilizedFootPointFace");
// CUT-CELL variant: the same parallel-surface inverse, but applied ONCE PER
// CUT CELL (0 < alpha < 1) and assigned to that cell's active faces, so the
// curvature difference across the force support is identically zero rather
// than merely small (WP8.1). Dimension-agnostic: cut cells, the fit
// derivatives and the Gaussian term are defined identically in 2D and 3D.
const bool cutCellFootPointFaceExtension =
    (curvatureExtensionType == "cutCellFootPointFace");
// CELL-MEAN variant: the PER-FACE inversion (each face keeps its own foot
// distance), then averaged over the active faces of each cut cell. Structurally
// the same delivery as cutCellFootPointFace -- one value per cut cell -- but
// built from n averaged inversions instead of one concentrated inversion, which
// is what separates them numerically: WP8.2 measured that the gain, not the
// across-support structure, governs the coupled growth rate, and averaging
// lowers the gain where concentrating raises it.
const bool cellMeanFootPointFaceExtension =
    (curvatureExtensionType == "cellMeanFootPointFace");
// FOOT-POINT-EVALUATED variant: each adjacent cell finds the foot point OF THE
// FACE CENTRE on its own quadratic's zero set and evaluates the fit curvature
// AT that point -- the value belongs to the interface point itself, so no
// parallel-surface conversion (and no level-set-spacing assumption) enters.
// The face takes the LINEAR interpolation of the two per-cell values (both
// estimate the SAME interface point), one-sided where only one fit converges.
// Measured motivation (plan sec. 17): the conversion-based delivery drops to
// first order when the level sets are not normal offsets of the interface;
// foot evaluation holds second order there (0.187 vs 8.69 at N=512 on the
// quadratic-form ellipse) and loses where the fit error dominates instead.
const bool footPointEvaluatedFaceExtension =
    (curvatureExtensionType == "footPointEvaluatedFace");
// CELL-foot-evaluated variant: the same stabilized foot-point evaluation, but
// ONCE PER SUPPORT CELL at that cell's OWN centre, the resulting CELL field
// delivered to the faces by plain arithmetic interpolation. Keeps what the
// face-evaluated variant is for (no parallel-surface inverse, hence no
// psi = f(signed distance) hypothesis) and restores the structure the pressure
// projection can absorb: sigma kappaBar snGrad(alpha)|Sf| is snGrad of a CELL
// field, so a uniform curvature error is cancelled exactly and only the
// face-to-face VARIATION reaches the velocity. Evaluating per FACE resamples
// the same quadratic at a different point of its own zero set for every face
// -- a conic's curvature is not constant along itself -- which is why that
// delivery's relative face-to-face jump is O(h) where production's is O(h^2)
// (measured 112x apart at N=256 on the exact-SDF circle, WP8.1 driver split).
const bool cellFootPointEvaluatedFaceExtension =
    (curvatureExtensionType == "cellFootPointEvaluatedFace");
// SYMMETRIC FACE-MEAN variant: the per-face inversion followed by a smoothing
// pass over each face's own owner+neighbour active-face set, which is symmetric
// about that face and therefore annihilates linear fields -- second order is
// preserved where the cell-centred lumping of the two variants above loses an
// order (plan sec. 12: 1.02 and 1.03 against 1.98 on the varying-curvature gate).
// BAND RENORMALIZATION (fvSolution levelSet.psiRenormalization.type):
//   none            -- default, byte-identical to every previous run
//   footPointGradient -- psi <- psi/|grad psi at the foot point| on the band each
//                      step. Restores the property the parallel-surface inverse
//                      assumes; a strictly positive rescaling, so the zero set and
//                      every sign are preserved exactly (bandRenormalization.H).
// COST: one extra per-cell quadratic fit per step, because the fit must be
// refreshed on the freshly advected psi before the foot-point gradient is read,
// and the curvature fill later fits the renormalized psi. Roughly doubles the
// reconstruction cost of a step; reported rather than hidden.
word psiRenormalizationType("none");
{
    const fvSolution& fvSol = mesh;
    psiRenormalizationType =
        fvSol.subDict("levelSet").subOrEmptyDict("psiRenormalization")
             .getOrDefault<word>("type", "none");
}
const bool bandRenormalization = (psiRenormalizationType == "footPointGradient");
// Gate threshold on the dimensionless band non-parallelism |P H n|/beta * h
// (<= 0 fires every step), smoothing sweeps on the divisor, and the RELAXATION
// exponent omega: psi <- psi/beta_Gamma^omega. See bandRenormalization.H for why
// omega < 1 is the difference between removing the drift and kicking the field.
scalar psiRenormThreshold = 5e-3;
label  psiRenormSweeps = 2;
scalar psiRenormRelaxation = 1.0;
{
    const dictionary& d =
        static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("psiRenormalization");
    psiRenormThreshold = d.getOrDefault<scalar>("threshold", 5e-3);
    psiRenormSweeps = d.getOrDefault<label>("smoothingSweeps", 2);
    psiRenormRelaxation = d.getOrDefault<scalar>("relaxation", 1.0);
}
if (bandRenormalization)
{
    Info<< "  gate threshold (aT*h) = " << psiRenormThreshold
        << ", smoothing sweeps = " << psiRenormSweeps
        << ", relaxation = " << psiRenormRelaxation << nl;
}
if (psiRenormalizationType != "none" && !bandRenormalization)
{
    FatalErrorInFunction
        << "levelSet.psiRenormalization.type must be none or footPointGradient "
        << "(found " << psiRenormalizationType << ")." << exit(FatalError);
}
Info<< "psi renormalization: " << psiRenormalizationType << nl << endl;

const bool symmetricFaceMeanExtension =
    (curvatureExtensionType == "symmetricFaceMeanFootPointFace");
const scalar symmetricFaceMeanTheta =
    static_cast<const fvSolution&>(mesh).subDict("levelSet")
    .subOrEmptyDict("curvatureExtension")
    .getOrDefault<scalar>("faceSmoothing", 0.5);
if (stabilizedFootPointFaceExtension || cutCellFootPointFaceExtension
 || cellMeanFootPointFaceExtension || symmetricFaceMeanExtension
 || footPointEvaluatedFaceExtension || cellFootPointEvaluatedFaceExtension)
{
    // The parallel-curve inverse acts at the FACE: a cell-level
    // offsetCorrection would be applied twice -- silently wrong curvature.
    const word oc =
        static_cast<const fvSolution&>(mesh).subDict("levelSet")
        .subOrEmptyDict("semiLagrangian")
        .getOrDefault<word>("offsetCorrection", "none");
    if (oc != "none")
    {
        FatalErrorInFunction
            << "curvatureExtension stabilizedFootPointFace requires "
            << "semiLagrangian offsetCorrection none (found " << oc
            << "): the face delivery applies the parallel-curve inverse "
            << "itself and would double-correct."
            << exit(FatalError);
    }
}
Info<< "Curvature extension: " << curvatureExtensionType
    << (harmonicCurvatureExtension ? " (harmonic Laplace smoothing ON)" : "")
    << (footPointCurvatureExtension ? " (foot-point height function ON)" : "")
    << (closestPointCurvatureExtension ? " (closest-point Newton foot ON)" : "")
    << (interfaceMeanCurvatureExtension ? " (interface-mean diagnostic ON)" : "")
    << (connectedInterfaceCurvatureExtension ? " (connected face-native curvature ON)" : "")
    << (stabilizedFootPointFaceExtension ? " (stabilized foot-point FACE delivery ON)" : "")
    << (cutCellFootPointFaceExtension ? " (CUT-CELL foot-point delivery ON)" : "")
    << (cellMeanFootPointFaceExtension ? " (CELL-MEAN of per-face inversions ON)" : "")
    << (symmetricFaceMeanExtension ? " (SYMMETRIC face-mean smoothing ON)" : "")
    << (footPointEvaluatedFaceExtension ? " (FOOT-POINT-EVALUATED face delivery ON)" : "")
    << (cellFootPointEvaluatedFaceExtension ? " (CELL-foot-evaluated delivery ON)" : "")
    << (fvmCurvature ? " (symbolic kappa fill OFF -- force model computes its own)" : "")
    << nl << endl;

// Optional grid-scale psi filter (fvSolution levelSet.psiFilter): the damping
// term the reinit-free SL transport lacks. type biharmonicBand applies
// psi -= theta * L(L(psi)) in the narrow band + 1 ring each step, with
// L(psi) = psi - average(interpolate(psi)). L annihilates linear fields and
// maps quadratics to a constant, so L(L(.)) annihilates smooth curved signed
// distances to leading order (no systematic interface displacement); the 2h
// checkerboard mode -- the carrier of the parasitic feedback (Defect 2) --
// is damped at near-maximal rate. Default: none (off, byte-identical).
word psiFilterType("none");
scalar psiFilterTheta = 0;
{
    const fvSolution& fvSol = mesh;
    const dictionary& lsDict = fvSol.subDict("levelSet");
    const dictionary& pfDict = lsDict.subOrEmptyDict("psiFilter");
    psiFilterType = pfDict.getOrDefault<word>("type", "none");
    psiFilterTheta = pfDict.getOrDefault<scalar>("theta", 0.05);
}
const bool psiFilterOn = (psiFilterType == "biharmonicBand");
if (psiFilterOn)
{
    Info<< "psi filter: biharmonicBand, theta = " << psiFilterTheta
        << nl << endl;
}

// Optional curvature clamp (DIAGNOSTIC probe, not a fix): if fvSolution provides
//   levelSet { curvatureClamp { min <scalar>; max <scalar>; } }
// the cell-centred kappa is clamped to [min,max] after each evaluation. Purpose:
// discriminate whether the sign-indefinite grid-scale kappa excursions (the
// aliasing channel) are load-bearing for the parasitic runaway -- if clamping to
// a physical range (e.g. [0, 2/R]) stabilises, they are. Absent dict = off.
bool kappaClampOn = false;
scalar kappaClampMin = 0, kappaClampMax = 0;
{
    const fvSolution& fvSol = mesh;
    const dictionary& ccDict =
        fvSol.subDict("levelSet").subOrEmptyDict("curvatureClamp");
    if (ccDict.found("min") && ccDict.found("max"))
    {
        kappaClampOn = true;
        kappaClampMin = ccDict.get<scalar>("min");
        kappaClampMax = ccDict.get<scalar>("max");
        Info<< "Curvature clamp (diagnostic): kappa in ["
            << kappaClampMin << ", " << kappaClampMax << "] 1/m" << nl << endl;
    }
}
