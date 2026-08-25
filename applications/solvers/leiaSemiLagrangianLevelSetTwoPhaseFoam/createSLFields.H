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

// WHERE IN TIME THE CAPILLARY FORCE IS EVALUATED.
//
//   endStep  (default, and what every historical run did): slAlphaEqn advects
//            psi over the FULL step before UEqn, so the curvature -- and hence
//            the force -- is built from psi^{n+1}. Writing the linearised
//            capillary oscillator eta' = u, u' = -omega^2 eta, that is
//              eta^{n+1} = eta^n + dt u^n,  u^{n+1} = u^n - dt omega^2 eta^{n+1}
//            i.e. SYMPLECTIC EULER: det(M) = 1 exactly, so |lambda| = 1 while
//            the eigenvalues stay complex, which holds for omega dt < 2. It is
//            already neutrally stable and area preserving -- NOT the fully
//            explicit force-at-t^n scheme, whose |1 + i omega dt| > 1 always.
//
//   midpoint: Stormer-Verlet. Drift a HALF step to psi^{n+1/2}, build the force
//            there, solve momentum, then drift the remaining half with u^{n+1}.
//            This is Popinet's arrangement (JCP 228:5838): the volume fraction
//            is staggered at n +/- 1/2 against velocity at n, with the force
//            term (sigma kappa delta_s n)_{n+1/2} between the two velocity
//            levels.
//
// MEASURED EXPECTATION, so nobody reads more into this switch than it carries:
// for the LINEAR oscillator the two are SPECTRALLY IDENTICAL -- same det = 1 and
// the same trace 2 - omega^2 dt^2, hence the same eigenvalues, the same
// amplification and the same phase error per step. Centring the force cannot by
// itself change linear stability or the step limit. What it changes is the
// offset between the velocity and interface samples, i.e. the accuracy of the
// coupling, and it makes our arrangement match Popinet's. The stationary
// droplet's growth is attributed to the curvature error depending on where the
// interface sits RELATIVE TO THE MESH (the measured m = 8 mesh-locked mode),
// which makes the force non-conservative and voids the symplectic guarantee --
// a defect no time centring addresses. This switch exists to settle that by
// measurement rather than by argument.
//
// Cost: midpoint runs the interface pipeline twice per step (once at n+1/2 for
// the force, once at n+1 so alpha, rho and the metrics are at the step end), so
// expect roughly 2x the per-step cost. endStep is untouched and byte-identical.
word capillaryForceCentring("endStep");
{
    const fvSolution& fvSolCfc = mesh;          // fvMesh is-a fvSolution
    capillaryForceCentring = fvSolCfc.subDict("levelSet")
        .getOrDefault<word>("capillaryForceCentring", "endStep");
}
const bool midpointForceCentring = (capillaryForceCentring == "midpoint");
if (capillaryForceCentring != "endStep" && !midpointForceCentring)
{
    FatalErrorInFunction
        << "levelSet.capillaryForceCentring must be endStep or midpoint, got "
        << capillaryForceCentring << exit(FatalError);
}
Info<< "Capillary force centring: " << capillaryForceCentring
    << (midpointForceCentring ? " (Stormer-Verlet, force at n+1/2)"
                              : " (symplectic Euler, force at n+1)") << endl;

//- Set while slAlphaEqn.H is being re-entered for the SECOND half drift of the
//  Stormer-Verlet split, after momentum has produced u^{n+1}.
bool slSecondHalfDrift = false;

// FREEZE-KAPPA DIAGNOSTIC (env SL_FREEZE_KAPPA=1, research only, default off):
// fill the cell curvature on the FIRST step and never refill it. The force then
// acts with a frozen, near-exact kappa while psi, alpha and the flow all stay
// live. This cuts exactly one arrow of the feedback loop
//     kappa -> force -> velocity -> psi -> kappa
// (the last one). Bounded => the runaway REQUIRES the kappa feedback, closing
// the attribution chain; still growing => the transport/state side self-excites
// through a frozen force, and the amplifier is upstream of the curvature. The
// companion clamp experiment already showed AMPLITUDE limiting does not help
// (growth at +388 1/s with kappa clamped to [-2000, 4000] around the exact
// 1000), so this tests the PHASE: a frozen kappa cannot correlate with the
// oscillating interface, whatever its amplitude error.
const bool freezeKappa = !Foam::getEnv("SL_FREEZE_KAPPA").empty();
if (freezeKappa)
{
    Info<< "SL_FREEZE_KAPPA: cell curvature is filled on step 1 "
        << "and FROZEN afterwards (diagnostic)" << nl << endl;
}

// THE NON-GRADIENT CAPILLARY VELOCITY, kept so the WORK it does can be measured.
//
// The pressure projection IS a Helmholtz decomposition in the rAUf-weighted
// metric: it removes exactly the part of the capillary flux that lies in the
// range of snGrad, and what it leaves, R_f = phig - p_rghEqn.flux(), is the
// non-gradient remainder that fvc::reconstruct turns into velocity. In CSF form
// f = sigma*kappa*grad(H), so
//     curl(f) = sigma * grad(kappa) x grad(H),
// which vanishes if and only if kappa is CONSTANT ALONG THE INTERFACE. Hence
// "non-gradient content" and "curvature varying along Gamma" are the same
// statement -- the quantity kappaStdDevBand measures.
//
// On a stationary droplet the exact kappa is constant, so every bit of this
// residual is discretisation error, and the work it does on the flow is exactly
// the spurious energy input. pEqn.H already builds the field; this store keeps
// it alive past pEqn's scope so writeDropletMetrics.H can form
//     P = sum_c rho_c V_c (U_c . dU_c) / dt
// and compare the energy it injects per capillary period against the measured
// growth of the envelope. That comparison is what makes the "the force does net
// work because the curvature error is mesh-position dependent" claim falsifiable.
volVectorField USigmaResStore
(
    IOobject("USigmaResStore", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("USigmaResStore", dimVelocity, Zero)
);
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
// CELL-CENTRE parallel-surface inverse applied IN PLACE to the cell curvature,
// then delivered by ordinary interpolation (faceCurvatureSource model). NOT the
// cut-cell delivery: that computes the same inverse and then lumps one value
// onto every active face of the cell, discarding the cell-to-cell smoothness.
// Measured motivation (cellCentreInverseCurvature.H): the non-gradient content
// alpha_f*snGrad(kappa_c) converges at order +2.01 over the full CSF support,
// against +0.09 for every other cell curvature model, and is 3200x smaller at
// N=512 -- the only cell field that makes the unbalanceable part of the force
// converge under refinement.
const bool cellCentreInverseExtension =
    (curvatureExtensionType == "cellCentreInverse");
// GAUSSIAN-CURVATURE TERM of the parallel-surface inverse, on/off:
//   levelSet { curvatureExtension { type cellCentreInverse; gaussianCurvature no; } }
// Default YES = the K-aware inverse, byte-identical to every previous run.
//
// WHY THIS SWITCH EXISTS. The inverse is one dimension-seamless expression,
//     kappa^Gamma = (kappa - 2 K d)/(1 - d kappa + K d^2),
// in which K is IDENTICALLY ZERO in 2D and active only in 3D. Setting it to no
// runs the 3D case through the same algebra 2D uses, so it is the one available
// control that isolates a dimension-specific term of the delivery. It is
// motivated by measurement, not symmetry: with the K-aware form the coupled 3D
// ladder is stable at R/h = 10.0 and 12.7 (max|U| 1.00e-3, 1.06e-3) and
// destabilises at 15.8 (7.04e-2, growth rate 90.6 1/s against 40.5 and 52.3),
// while the WP0 normal-direction band spectrum stays pinned at 1.00x of its t=0
// value in every arm -- so whatever feeds the 3D growth is NOT the corrugation
// the filter damps, and K is one of only two dimension-specific candidates.
//
// EXPECTED COST if K is innocent: the static 3D sphere gate scored the K-aware
// inverse at h^1.95 and the K-less control at h^1.02, so switching K off should
// lose an order in the delivery. If the 3D instability disappears ANYWAY, the
// growth is being fed by K's fit noise rather than by the delivery's order, and
// the two can then be traded deliberately. If it persists at h^1.02, K is
// exonerated and the tangential-structure candidate is the remaining one.
bool cellCentreInverseGaussian = true;
if (cellCentreInverseExtension)
{
    const fvSolution& fvSolK = mesh;
    cellCentreInverseGaussian =
        fvSolK.subDict("levelSet").subOrEmptyDict("curvatureExtension")
            .getOrDefault<bool>("gaussianCurvature", true);
    Info<< "cellCentreInverse: Gaussian-curvature term "
        << (cellCentreInverseGaussian ? "ON (K-aware, 3D-active)"
                                      : "OFF (K == 0, the 2D algebra)") << endl;
}
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
    << (cellCentreInverseExtension ? " (CELL-CENTRE parallel-surface inverse ON)" : "")
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
