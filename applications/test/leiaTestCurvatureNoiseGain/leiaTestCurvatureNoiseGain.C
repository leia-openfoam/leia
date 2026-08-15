/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    leiaTestCurvatureNoiseGain

Description
    How strongly does a face-centred curvature delivery AMPLIFY a perturbation
    of psi?

    leiaTestMeanCurvature measures ACCURACY: the error of kappa_f against the
    exact value on an exact signed distance field. That is the wrong number for
    predicting the coupled run. In the stationary-droplet solve psi is never
    exact: every step it carries the semi-Lagrangian reconstruction error, and
    the delivery converts that error into a face force. What decides whether
    the parasitic-current loop grows is therefore not how accurate kappa_f is
    on a clean psi, but how much kappa_f MOVES when psi moves -- the gain

        G = || d kappa_f || / || d psi ||        [1/m^2]

    measured on the faces that actually carry the CSF force (|snGrad(alpha)|
    above the active threshold).

    MEASURED (this is what motivated the utility): over the three deliveries
    that have coupled stationary-droplet blow-up times at N=128, accuracy and
    stability are ANTI-correlated -- the more accurate the delivery on the
    static gate, the sooner the droplet blows up. The gain is the candidate
    explanation, and it is measured here rather than argued.

    Method:
      1. baseline: recon->update(psi), cell curvature from
         meanCurvatureNoExtension, then each delivery's kappa_f;
      2. perturbed: psi_e = psi + eps*h*r with r in [-1,1] a reproducible
         per-cell pseudo-random field (splitmix64 of the cell index and a
         seed), recon->update(psi_e), and the SAME three deliveries;
      3. alpha and the active-face set are held FIXED, so the measurement
         isolates d kappa_f / d psi and does not mix in the (separate)
         question of how the phase indicator responds.

    Averaged over NSEEDS realisations. Sweeping EPS checks that the response
    is linear -- if it is not, the delivery is behaving like a switch and the
    gain is not a meaningful single number.

    Reported per delivery (leiaTestCurvatureNoiseGain.csv):
      GAIN_L2      || d kappa_f ||_2 / || d psi ||_inf on active faces [1/m^2]
      GAIN_LINF    max | d kappa_f | / || d psi ||_inf                [1/m^2]
      GAIN_DIMLESS GAIN_L2 * h^2 -- the amplification RELATIVE to a plain
                   second difference of psi, which is the h^-2 that any
                   curvature operator must pay. 1 means "as sensitive as a
                   second difference"; 10 means ten times as sensitive.
      E_L2         the static accuracy of the same delivery, so the two
                   numbers can be read side by side.

    SERIAL (like leiaTestMeanCurvature: internal faces only).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvSolution.H"
#include "slReconstruction.H"

// The production face deliveries, measured exactly as the solver applies them.
#include "stabilizedFootPointFaceCurvature.H"
#include "capillaryDriverSplit.H"
#include "levelSetImplicitSurfaces.H"

using namespace Foam;

// Reproducible per-cell pseudo-random number in [-1, 1]: splitmix64 finaliser
// on (seed, cell index). Deterministic across runs and compilers, unlike
// rand(), so a reported gain can be reproduced exactly.
static inline scalar unitNoise(const label cellI, const label seed)
{
    uint64_t z = uint64_t(cellI) + 0x9E3779B97F4A7C15ull*uint64_t(seed + 1);
    z += 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30))*0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27))*0x94D049BB133111EBull;
    z = z ^ (z >> 31);
    // 53-bit mantissa -> [0,1) -> [-1,1]
    const scalar u = scalar(z >> 11)/scalar(uint64_t(1) << 53);
    return 2.0*u - 1.0;
}


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "alphaName", "word",
        "Volume-fraction field defining the CSF force support (default alpha.water)."
    );
    argList::addOption
    (
        "eps", "list",
        "Noise amplitudes as a fraction of h (default \"(0.001 0.01 0.1)\")."
    );
    argList::addOption
    (
        "nSeeds", "label",
        "Noise realisations averaged per amplitude (default 8)."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "leiaTestCurvatureNoiseGain is a SERIAL test: the norms cover "
            << "internal faces only, so a decomposed run would drop the "
            << "processor-seam part of the CSF force support."
            << exit(FatalError);
    }

    const word alphaName = args.getOrDefault<word>("alphaName", "alpha.water");
    const label nSeeds = args.getOrDefault<label>("nSeeds", 8);
    List<scalar> epsList;
    if (args.found("eps"))
    {
        epsList = args.get<List<scalar>>("eps");
    }
    else
    {
        epsList = List<scalar>({1e-3, 1e-2, 1e-1});
    }

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    const volScalarField alpha
    (
        IOobject(alphaName, runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    const fvSolution& fvSol(mesh);
    const dictionary& implSurf =
        fvSol.subDict("levelSet").subDict("implicitSurface");
    const label nd = mesh.nGeometricD();
    const bool is2D = (nd == 2);
    // Position-dependent exact curvature on the varying-curvature gates; the
    // radius-based constant elsewhere (implicitSphere::curvature() returns 1/R
    // in 3D too, so it must not be used there). Only the E_L2 accuracy column
    // depends on this -- the GAIN never references the exact value.
    const word surfType = implSurf.get<word>("type");
    const bool varyingKappa =
        (surfType == "signedDistanceEllipse" || surfType == "implicitEllipsoid");
    autoPtr<implicitSurface> exactSurf;
    scalar kappaExact = 0;
    if (surfType == "implicitEllipsoid")
    {
        // Quadratic-form psi (non-parallel foliation): the reference is still
        // the zero-set ellipse's curvature at the closest point -- see the
        // identical branch in leiaTestMeanCurvature.C.
        exactSurf.reset
        (
            new signedDistanceEllipse
            (
                implSurf.get<vector>("center"),
                implSurf.get<vector>("axes")
            )
        );
    }
    else if (varyingKappa)
    {
        exactSurf = implicitSurface::New(surfType, implSurf);
    }
    else
    {
        kappaExact = scalar(nd - 1)/implSurf.get<scalar>("radius");
    }
    auto kappaExactAt = [&](const point& x) -> scalar
    {
        return varyingKappa ? exactSurf->curvature(x) : kappaExact;
    };

    // The deliveries measured here apply the parallel-surface inverse
    // themselves; a reconstruction that already offset-corrects would be
    // corrected twice (same contract as leiaTestMeanCurvature).
    {
        const word offsetCorrection =
            fvSol.subDict("levelSet").subOrEmptyDict("semiLagrangian")
                 .getOrDefault<word>("offsetCorrection", "none");
        if (offsetCorrection != "none")
        {
            FatalErrorInFunction
                << "requires offsetCorrection none (found "
                << offsetCorrection << ")." << exit(FatalError);
        }
    }

    const boundBox bb(mesh.points());
    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar Lz = bb.max().z() - bb.min().z();
    const scalar h =
        is2D ? Foam::sqrt(Lx*Ly/mesh.nCells())
             : Foam::cbrt(Lx*Ly*Lz/mesh.nCells());

    autoPtr<slReconstruction> recon = slReconstruction::New(mesh);

    volScalarField kappa
    (
        IOobject("kappaNG", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0), "zeroGradient"
    );
    surfaceScalarField kappaFace
    (
        IOobject("kappaFaceNG", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar("k", dimless/dimLength, 0.0)
    );

    // The force support: exactly the faces the CSF term acts on.
    const surfaceScalarField snA(fvc::snGrad(alpha));
    const scalarField& snI = snA.primitiveField();
    const scalarField& magSfIn = mesh.magSf().primitiveField();
    const vectorField& CfIn = mesh.Cf().primitiveField();
    const label nIntFaces = mesh.nInternalFaces();
    const scalar snMax = max(gMax(mag(snI)()), VSMALL);
    boolList activeFace(nIntFaces, false);
    label nActive = 0;
    scalar sumAf = 0;
    forAll(activeFace, f)
    {
        if (mag(snI[f]) > 1e-8*snMax)
        {
            activeFace[f] = true;
            ++nActive;
            sumAf += magSfIn[f];
        }
    }
    Info<< "h = " << h << ", kappa_exact = " << kappaExact
        << ", active faces = " << nActive << nl << endl;

    if (nActive == 0)
    {
        FatalErrorInFunction << "no active faces." << exit(FatalError);
    }

    // The three deliveries with coupled stationary-droplet blow-up times, plus
    // the cell-mean variant: the SAME per-face inversions as perFaceInverse,
    // averaged over each cut cell's active faces. It has the same across-support
    // structure as cutCellInverse but is built by averaging n inversions rather
    // than concentrating on one, so its gain is the direct test of whether
    // averaging lowers the gain.
    const wordList models
    (
        {"arithmetic", "perFaceInverse", "cutCellInverse", "cellMeanInverse",
         "symFaceMean050", "footEvalFace", "cellFootEvalFace"}
    );

    // kappa_f for every delivery, from the CURRENT psi in the reconstruction.
    auto deliver = [&](const volScalarField& psiUse, List<scalarField>& kf)
    {
        recon->update(psiUse);
        recon->meanCurvatureNoExtension(kappa);
        kappa.correctBoundaryConditions();

        kf.setSize(models.size());

        // arithmetic: the classical CSF delivery, interpolate of the cell field
        kappaFace = fvc::interpolate(kappa);
        kf[0] = kappaFace.primitiveField();

        // per-face parallel-surface inverse (curvatureExtension
        // stabilizedFootPointFace)
        computeStabilizedFootPointFaceCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace
        );
        kf[1] = kappaFace.primitiveField();

        // one inverted curvature per cut cell (curvatureExtension
        // cutCellFootPointFace)
        computeCutCellFootPointFaceCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace
        );
        kf[2] = kappaFace.primitiveField();

        // per-face inversions, averaged over each cut cell's active faces
        // (curvatureExtension cellMeanFootPointFace)
        computeCellMeanFootPointFaceCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace
        );
        kf[3] = kappaFace.primitiveField();

        // per-face inversions smoothed SYMMETRICALLY about each face
        // (curvatureExtension symmetricFaceMeanFootPointFace, theta = 0.5)
        computeSymmetricFaceMeanCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace, 0.5
        );
        kf[4] = kappaFace.primitiveField();

        // fit curvature EVALUATED at the face centre's foot on each adjacent
        // cell's own zero set, linearly interpolated -- no parallel-surface
        // conversion (curvatureExtension footPointEvaluatedFace)
        computeFootPointEvaluatedFaceCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace
        );
        kf[5] = kappaFace.primitiveField();

        // the SAME stabilized foot-point evaluation, but once per support CELL
        // at that cell's own centre, delivered by plain interpolation
        // (curvatureExtension cellFootPointEvaluatedFace) -- the interpolate of
        // a cell field, which is what the pressure projection can absorb
        computeCellFootPointEvaluatedFaceCurvature
        (
            mesh, psiUse, alpha, kappa, recon(), kappaFace
        );
        kf[6] = kappaFace.primitiveField();
    };

    List<scalarField> kfBase;
    deliver(psi, kfBase);

    // Static accuracy of each delivery on the unperturbed field, so accuracy
    // and gain are reported from the same run and cannot drift apart.
    scalarList eL2(models.size(), scalar(0));
    forAll(models, m)
    {
        scalar s2 = 0;
        forAll(activeFace, f)
        {
            if (!activeFace[f]) { continue; }
            const scalar e = kfBase[m][f] - kappaExactAt(CfIn[f]);
            s2 += magSfIn[f]*e*e;
        }
        eL2[m] = Foam::sqrt(s2/sumAf);
    }

    OFstream os(runTime.path()/"leiaTestCurvatureNoiseGain.csv");
    os << "MODEL,N_ACTIVE_FACES,DELTA_X,KAPPA_EXACT,EPS,N_SEEDS,"
       << "E_L2,GAIN_L2,GAIN_LINF,GAIN_DIMLESS,"
       << "GAIN_DRIVER_ACROSS,GAIN_DRIVER_ALONG,"
       << "GAIN_DRIVER_ACROSS_DIMLESS,GAIN_DRIVER_ALONG_DIMLESS" << nl;

    volScalarField psiPerturbed("psiPerturbed", psi);

    for (const scalar eps : epsList)
    {
        const scalar a = eps*h;                       // ||d psi||_inf

        scalarList gL2(models.size(), scalar(0));
        scalarList gLi(models.size(), scalar(0));
        // The part of the response the pressure projection CANNOT remove: the
        // face-pair variation RATE of d kappa_f, binned across the force
        // support and along the interface exactly as the solver's WP8.1
        // diagnostic bins the delivered curvature itself.
        scalarList gDrAcross(models.size(), scalar(0));
        scalarList gDrAlong(models.size(), scalar(0));

        for (label seed = 0; seed < nSeeds; ++seed)
        {
            scalarField& pIn = psiPerturbed.primitiveFieldRef();
            const scalarField& p0 = psi.primitiveField();
            forAll(pIn, c) { pIn[c] = p0[c] + a*unitNoise(c, seed); }
            psiPerturbed.correctBoundaryConditions();

            List<scalarField> kfP;
            deliver(psiPerturbed, kfP);

            // Restore the BASELINE fit before the driver split, so the normals
            // that assign each face pair to a bin are a property of the
            // geometry and not of the particular noise realisation.
            recon->update(psi);

            forAll(models, m)
            {
                scalar s2 = 0, li = 0;
                scalarField dK(nIntFaces, 0.0);
                forAll(activeFace, f)
                {
                    if (!activeFace[f]) { continue; }
                    const scalar d = kfP[m][f] - kfBase[m][f];
                    dK[f] = d;
                    s2 += magSfIn[f]*d*d;
                    li = max(li, mag(d));
                }
                gL2[m] += Foam::sqrt(s2/sumAf)/a;
                gLi[m] = max(gLi[m], li/a);

                const capillaryDriverSplitNorms dr =
                    computeFacePairVariation
                    (
                        mesh, psi, dK, activeFace, recon(), a
                    );
                gDrAcross[m] += dr.acrossSupportL2;
                gDrAlong[m] += dr.alongInterfaceL2;
            }
        }

        Info<< "eps = " << eps << " (|d psi|_inf = " << a << " m)" << nl;
        forAll(models, m)
        {
            const scalar g2 = gL2[m]/scalar(nSeeds);
            const scalar dAc = gDrAcross[m]/scalar(nSeeds);
            const scalar dAl = gDrAlong[m]/scalar(nSeeds);
            os  << models[m] << ',' << nActive << ',' << h << ','
                << kappaExact << ',' << eps << ',' << nSeeds << ','
                << eL2[m] << ',' << g2 << ',' << gLi[m] << ','
                << g2*h*h << ',' << dAc << ',' << dAl << ','
                << dAc*h*h*h << ',' << dAl*h*h*h << nl;

            Info<< "  " << models[m]
                << ": E_L2 = " << eL2[m] << " 1/m"
                << ", gain = " << g2 << " 1/m^2"
                << " (x" << g2*h*h << " a plain second difference)"
                << ", driver gain across/along = "
                << dAc*h*h*h << " / " << dAl*h*h*h << " (dimensionless)" << nl;
        }
        Info<< endl;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
