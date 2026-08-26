/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sdplsRdiv.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "fvSolution.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sdplsRdiv, false);
addToRunTimeSelectionTable(sdplsSource, sdplsRdiv, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sdplsRdiv::sdplsRdiv(const dictionary& dict, const fvMesh& mesh)
:
    sdplsR(dict, mesh)
{
    // The conservative assembly cannot compose with a cut-off: a mollifier
    // multiplying a FLUX destroys the conservation this form exists to
    // provide. Refuse loudly rather than run something the dictionary does
    // not describe.
    const word mollType =
        dict.subOrEmptyDict("mollifier").getOrDefault<word>("type", "none");
    if (mollType != "none")
    {
        FatalErrorInFunction
            << "sdplsSource type `Rdiv` is the CONSERVATIVE (divergence-form) "
            << "source; a mollifier (`" << mollType << "`) multiplying its "
            << "fluxes would break exactly the conservation it provides." << nl
            << "Use mollifier type `none` with Rdiv, or type `R` if a cut-off "
            << "is wanted."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvScalarMatrix>
sdplsRdiv::fvmsdplsSource(const volScalarField& psi, const volVectorField& U)
{
    // Refresh R_ (the strain a = n.grad(v).n) and the lagged normal, exactly
    // as the algebraic source does -- same Picard structure, measured
    // converged at nDefCorr = 3.
    update(psi, U);

    const fvMesh& mesh = psi.mesh();

    // The lagged interface normal, from the same dedicated gradient the
    // strain uses (gradPsi model + gradPsiSdpls scheme).
    volVectorField const gradPsi(grad(psi));
    dimensioned<scalar> const eps =
        dimensioned<scalar>("eps", gradPsi.dimensions(), SMALL);
    volVectorField const nHat(gradPsi/(mag(gradPsi) + eps));

    // Normal velocity w = (n.v) n and its face flux. The flux is NAMED so
    // fvSchemes can give div(phiW,psi) its own scheme; the solver's own mass
    // flux phi is looked up rather than rebuilt from U, so the same flux that
    // advects psi in the base equation appears in the correction terms.
    volVectorField const w(("wSdpls"), (nHat & U)*nHat);
    surfaceScalarField phiW("phiW", fvc::flux(w));

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    // ---- OPT-IN RESEARCH SEAM: single-operator (composite-flux) assembly ----
    //
    // Default off, so an absent key reproduces every Rdiv run ever made.
    //
    // WHAT IT CHANGES AND WHY. eulerianAdvection assembles
    //     ddt(psi) + div(phi, psi) - Sp(div phi, psi) == S,
    // and fvMatrix::operator== subtracts S termwise, so the DEFAULT
    //     S = div(phiW,psi) - div(phi,psi) - Sp(div(phiW) - a, psi)
    // leaves the system matrix
    //     ddt + 2 div(phi, psi) - div(phiW, psi) + (diagonal terms),
    // a DIFFERENCE of two upwind convection operators. Subtracting an M-matrix
    // flips its off-diagonal signs: measured at t = 0, N = 32, 687 of 1984
    // internal faces (35%) carry sign(phi_f) != sign(phiW_f), the two operators
    // upwind from opposite cells there, and 799 of 3968 off-diagonal entries
    // (20%) of the assembled psi matrix are POSITIVE.
    //
    // Setting compositeFlux on returns instead
    //     S = div(phi,psi) - div(2 phi - phiW, psi) - Sp(div(phiW) - a, psi),
    // which leaves
    //     ddt + div(2 phi - phiW, psi) + (the same diagonal terms).
    // Identical at the continuum -- the effective advecting velocity is
    // 2U - w = U + U_t and the tangential part U_t transports nothing because
    // U_t . grad(psi) = 0 -- but assembled as ONE fvm::div, whose upwind
    // off-diagonals are <= 0 in every row on every mesh by construction. The
    // flux is NAMED "phiW" so it reads the same div(phiW,psi) fvSchemes entry
    // and the discretisation family is unchanged.
    //
    // getOrDefault<word> + an explicit Switch, NOT getOrDefault<Switch>(key,
    // "off"). MEASURED HERE: with a STRING-LITERAL default the template
    // parameter is Switch and the literal reaches the Switch parameter as
    // const char*, which converts to bool (non-null => TRUE) instead of being
    // parsed as the word "off". The A/B that found it: an arm whose dictionary
    // carried NO compositeFlux entry took the compositeFlux branch anyway --
    // nPosOffDiag = 0 in a run that must have had 799 -- and produced a file
    // byte-identical to the arm that asked for it. src/leiaLevelSet/
    // velocityModel/prescribedVelocityModels.C:185 has the same construct,
    // `getOrDefault<Switch>("oscillation", "off")`, whose comment states that
    // the safe default is off; it is only masked there because every rendered
    // case writes the key explicitly. Reading a word and constructing the
    // Switch from it has one meaning.
    const Switch compositeFlux
    (
        sourceDict_.getOrDefault<word>("compositeFlux", word("off"))
    );

    // S_matrix = div(phiW, psi) - div(phi, psi) - Sp(div(phiW) - a, psi).
    //
    // Derivation (see the header): substituting S = -div_s(psi v) into the
    // base equation and moving everything the base already assembles back to
    // the left leaves exactly these three terms on the right. At the continuum
    // they collapse to psi (a - div v): the algebraic source plus the same
    // compressibility correction the base equation carries. Discretely the
    // two div terms are CONSERVATIVE -- their global sum telescopes to
    // boundary fluxes -- which is the property the algebraic source lacks and
    // the coupled droplet measurements show it needs.
    // ---- PRODUCT-RULE RESIDUAL DIAGNOSTIC (added 2026-08-20) ----------------
    //
    // WHY THIS IS MEASURED HERE AND NOT INFERRED FROM A CONVERGENCE TABLE.
    // Rdiv and the algebraic source R are CONTINUUM-IDENTICAL for incompressible
    // flow, yet on the 2D reversed vortex with a prescribed divergence-free
    // velocity -- no capillary force, no curvature, no feedback of any kind --
    // R converges at +1.112 in the band gradient error while Rdiv DIVERGES at
    // -3.800 (0.391, 0.256, 1.56, 8.00, 20.5, 103, 669 over N = 32..256, every
    // run COMPLETING, relative volume error reaching 136% at N=256). The only
    // structural difference between them is that Rdiv routes part of the
    // operator through fvm::div FLUXES. So the defect must lie in that routing.
    //
    // The sign-split variant RdivStrictSp tested and FALSIFIED the first
    // hypothesis, that the unsigned Sp coefficient f = a - div(w) destroys
    // diagonal dominance: it demonstrably removed the M-matrix violation (512 of
    // 1024 cells carried f > 0 at N=32) and the divergence barely moved (-3.554
    // against -3.800, 635 against 669 at the finest rung).
    //
    // WHAT THIS RESIDUAL TESTS. At the continuum the product rule gives
    //     div(w psi) - psi div(w) - w . grad(psi) = 0
    // identically. Rdiv RELIES on that identity: the psi*div(w) content carried
    // implicitly inside fvm::div(phiW, psi) is supposed to cancel against the
    // -Sp(div(phiW), psi) term, leaving only the convective part. Discretely
    // there is NO guarantee that it does -- fvm::div(phiW, psi) interpolates psi
    // to faces with the div(phiW,psi) scheme while fvc::div(phiW) is a plain
    // flux divergence, and the two need not be consistent. If this residual
    // fails to vanish, and worse GROWS with refinement, then the cancellation
    // Rdiv is built on does not happen discretely, and that is the defect.
    //
    // The reference norm is printed alongside deliberately: the two terms nearly
    // cancel, so a residual is only meaningful relative to the magnitude of what
    // is being differenced.
    {
        volScalarField const prDiv("prDiv", fvc::div(phiW, psi));
        volScalarField const prPsiDivW("prPsiDivW", psi*fvc::div(phiW));
        volScalarField const prConv("prConv", w & gradPsi);
        volScalarField const prResid("prResid", prDiv - prPsiDivW - prConv);

        const scalarField& r = prResid.primitiveField();
        const scalarField& d = prDiv.primitiveField();
        const scalarField& V = mesh.V();
        const scalar Vtot = gSum(V);

        Info<< "sdplsRdiv product-rule residual"
            << " |div(w psi) - psi div w - w.grad psi|:"
            << " Linf = " << gMax(mag(r))
            << ", L2 = " << Foam::sqrt(gSum(sqr(r)*V)/Vtot)
            << " | reference |div(w psi)|: Linf = " << gMax(mag(d))
            << ", L2 = " << Foam::sqrt(gSum(sqr(d)*V)/Vtot)
            << endl;
    }
    // -------------------------------------------------------------------------

    // ---- OFF-DIAGONAL SIGN OF THE ASSEMBLED psi MATRIX (added 2026-08-26) ---
    //
    // WHY A SECOND M-MATRIX PROBE AFTER RdivStrictSp FALSIFIED THE FIRST ONE.
    // RdivStrictSp sign-split the ALGEBRAIC coefficient f = a - div(w) only. It
    // demonstrably acted (512 of 1024 cells carried f > 0 at N = 32) and the
    // divergence barely moved (-3.554 against -3.800), which rules out the
    // algebraic Sp term -- and nothing else. The TRANSPORT terms were never
    // examined, and they carry the same concern in a much larger form.
    //
    // WHAT THE SOLVER ACTUALLY ASSEMBLES. eulerianAdvection builds
    //     ddt(psi) + div(phi, psi) - Sp(div phi, psi) == S
    // and fvMatrix::operator== subtracts the right-hand matrix termwise, so
    // with S = div(phiW,psi) - div(phi,psi) - Sp(div(phiW) - a, psi) the system
    // matrix is
    //     ddt(psi) + 2 div(phi, psi) - div(phiW, psi) + (diagonal terms).
    // The continuum collapse is exact: the effective advecting velocity is
    // 2U - w = U + U_t with U_t = U - (U.n)n the TANGENTIAL velocity, and
    // U_t . grad(psi) = 0 identically because grad(psi) is parallel to n. So
    // the tangential half of that operator transports nothing -- at the
    // continuum.
    //
    // DISCRETELY IT IS NOT A NULL OPERATOR, AND IT IS NOT AN M-MATRIX.
    // gaussConvectionScheme sets, per internal face, lower = -w_f phi_f and
    // upper = lower + phi_f, so for upwind (w_f = 1 if phi_f > 0, else 0) BOTH
    // off-diagonals of a single fvm::div are <= 0 whatever the sign of the flux,
    // and negSumDiag() makes the operator an M-matrix. SUBTRACTING one such
    // operator flips its off-diagonal signs. On a face where phi_f and phiW_f
    // have the SAME sign the composite off-diagonal is still <= 0 as long as
    // 2|phi_f| >= |phiW_f|; on a face where they have OPPOSITE signs the two
    // operators upwind from opposite cells and the composite off-diagonal is
    // exactly +|phiW_f| -- positive, i.e. an anti-diffusive coupling. The
    // number of such faces is not small a priori: phiW = flux((n.U)n) has the
    // face-normal component of the INTERFACE NORMAL in it, and n turns through
    // 2 pi around the circle while Sf does not.
    //
    // THIS IS THE ONE STRUCTURAL DIFFERENCE BETWEEN Rdiv AND R. The algebraic
    // source assembles ddt + div(phi,psi) - Sp(div phi) - Sp(a): its only
    // off-diagonals come from a SINGLE fvm::div, so they are <= 0 in every cell,
    // on every mesh, at every step, by construction. R cannot have this defect
    // and Rdiv has it structurally.
    //
    // ddt and both Sp terms touch the DIAGONAL only, so the off-diagonals of
    // the composite below are exactly those of the full system matrix. The
    // sharp criterion is therefore whether the positive off-diagonal row sum
    // exceeds the ddt diagonal V/dt, because that is the cell whose row cannot
    // be made diagonally dominant by any time step this run is taking.
    //
    // Internal faces only: on a coupled patch fvm::div contributes through
    // internal/boundaryCoeffs, not upper/lower, so a parallel run undercounts.
    {
        // The composite MUST be the operator this call is about to assemble,
        // or the A/B against compositeFlux measures nothing: the whole claim is
        // that the two forms differ in their off-diagonal signs and in nothing
        // else.
        tmp<fvScalarMatrix> tcomp
        (
            compositeFlux
          ? tmp<fvScalarMatrix>
            (
                fvm::div
                (
                    surfaceScalarField("phiW", 2*phi - phiW),
                    psi
                )
            )
          : tmp<fvScalarMatrix>
            (
                2*fvm::div(phi, psi) - fvm::div(phiW, psi)
            )
        );
        const fvScalarMatrix& comp = tcomp();

        const scalarField& up = comp.upper();
        const scalarField& lo = comp.lower();
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& V = mesh.V();
        const scalar dt = mesh.time().deltaTValue();

        scalarField posOffSum(mesh.nCells(), 0.0);
        label nPosOff = 0;
        label nSignMismatch = 0;
        scalar mismatchFlux = 0;
        scalar totalFlux = 0;

        forAll(own, faceI)
        {
            if (up[faceI] > 0) { ++nPosOff; posOffSum[own[faceI]] += up[faceI]; }
            if (lo[faceI] > 0) { ++nPosOff; posOffSum[nei[faceI]] += lo[faceI]; }

            const scalar pf = phi[faceI];
            const scalar pw = phiW[faceI];
            totalFlux += mag(pf);
            if (pf*pw < 0)
            {
                ++nSignMismatch;
                mismatchFlux += mag(pw);
            }
        }

        label nRowsBroken = 0;
        forAll(posOffSum, cellI)
        {
            if (posOffSum[cellI] > V[cellI]/dt) ++nRowsBroken;
        }

        const label nIntFaces = own.size();
        Info<< "sdplsRdiv matrixSign"
            << " t=" << mesh.time().timeOutputValue()
            << " nIntFaces=" << returnReduce(nIntFaces, sumOp<label>())
            << " nSignMismatchFaces=" << returnReduce(nSignMismatch, sumOp<label>())
            << " signMismatchFrac="
            << scalar(returnReduce(nSignMismatch, sumOp<label>()))
              /max(scalar(returnReduce(nIntFaces, sumOp<label>())), SMALL)
            << " mismatchFluxFrac="
            << returnReduce(mismatchFlux, sumOp<scalar>())
              /max(returnReduce(totalFlux, sumOp<scalar>()), SMALL)
            << " nPosOffDiag=" << returnReduce(nPosOff, sumOp<label>())
            << " posOffDiagFrac="
            << scalar(returnReduce(nPosOff, sumOp<label>()))
              /max(scalar(2*returnReduce(nIntFaces, sumOp<label>())), SMALL)
            << " nRowsNotDominant=" << returnReduce(nRowsBroken, sumOp<label>())
            << "/" << returnReduce(mesh.nCells(), sumOp<label>())
            << " maxPosOffSumTimesDtOverV="
            << gMax(posOffSum*dt/V)
            << endl;
    }
    // -------------------------------------------------------------------------

    // ---- EXACT-COEFFICIENT DIAGNOSTIC (added 2026-08-26) --------------------
    //
    // WHAT THE PRODUCT-RULE RESIDUAL ABOVE CANNOT SEE.
    // That residual differences div(w psi), psi div(w) and w.grad(psi) built
    // from the SAME discrete coefficient div(w) = fvc::div(phiW). A coefficient
    // that is noisy but used CONSISTENTLY cancels out of it, so it isolates
    // ASSEMBLY inconsistency (suspect ii) and says NOTHING about the accuracy
    // of the coefficient itself (suspect i): fvc::div(phiW) contains
    // div(nHat) -- SECOND derivatives of psi -- whose noise grows as h shrinks
    // and could be feeding the inconsistency rather than acting independently.
    // The two suspects imply opposite repairs, so they must be measured apart.
    //
    // THE CLOSED-FORM REFERENCE. On cases/2Dvortex at t = 0 the level set is
    // the EXACT signed distance to a circle and the velocity is the analytic
    // shear2D field, so div(w) is known in closed form and needs no finite
    // differencing:
    //
    //   psi(x) = |x - c| - R        (leiaSetFields writes exactly this at the
    //                                cell centres; implicitSphere::value)
    //   n      = grad(psi)          = (x - c)/r = e_r,  r = |x - c|,  |n| = 1
    //   w      = (n.U) n = s n,     s = U.n
    //
    //   div(w) = div(s n) = n.grad(s) + s div(n)
    //   n.grad(s) = n_i d_i (U_j n_j) = n_i n_j d_i U_j + U_j (n_i d_i n_j)
    //
    // The second term is U.((n.grad) n), and it VANISHES IDENTICALLY for a
    // signed-distance field: n is constant along its own rays (for the circle
    // n = e_r depends on the polar angle only, so d(e_r)/dr = 0). The first
    // term is exactly the normal strain a = n.grad(U).n that R_ approximates.
    // div(n) of a signed distance to a sphere is (d-1)/r -- 1/r for a circle in
    // 2D, 2/r for a sphere in 3D -- hence
    //
    //   div(w)_exact = n.grad(U).n + (U.n) (d-1)/r.                       (*)
    //
    // Both terms are invariant under n -> -n, so the sign convention of psi is
    // irrelevant here.
    //
    // With shear2D scaled by the reversal factor f(t) = cos(pi t/tau):
    //   u  =  sin(2 pi y) sin^2(pi x) f
    //   v  = -sin(2 pi x) sin^2(pi y) f
    //   d_x u =  pi sin(2 pi x) sin(2 pi y) f      [d/dx sin^2(pi x) = pi sin(2 pi x)]
    //   d_y u =  2 pi cos(2 pi y) sin^2(pi x) f
    //   d_x v = -2 pi cos(2 pi x) sin^2(pi y) f
    //   d_y v = -pi sin(2 pi x) sin(2 pi y) f
    // so d_x u + d_y v = 0 exactly: the reference velocity is solenoidal to the
    // last bit, which is what makes Rdiv and R continuum-identical here.
    //   a_exact = n_x^2 d_x u + n_x n_y d_x v + n_y n_x d_y u + n_y^2 d_y v
    // in OpenFOAM's convention grad(U)_ij = d_i U_j, the same contraction
    // sdplsSource::R() performs discretely as (grad_U & n) & n.
    //
    // THE THREE ERRORS REPORTED, and why the split is the whole point:
    //   divW_err    = fvc::div(fvc::flux(w_h))  - (*)   TOTAL coefficient error
    //   divW_opErr  = fvc::div(fvc::flux(w_ex)) - (*)   the flux-divergence
    //                 OPERATOR alone, fed an exact, smooth w -- no discrete
    //                 normal anywhere in it
    //   divW_normalErr = the difference of the two = the contribution of the
    //                 discrete nHat, i.e. of the second derivatives of psi
    // If divW_opErr converges while divW_err grows, the amplifier is the
    // fitted normal and not the divergence operator, and suspect (i) is live.
    //
    // AND THE CONTROLLED SUBSTITUTION THAT SEPARATES (i) FROM (ii): the product
    // rule residual is recomputed with w_ex in place of w_h, everything else
    // (the discrete psi, the div(phiW,psi) scheme, the gradPsi model) held
    // fixed. A residual that still GROWS on a perfectly smooth coefficient is
    // assembly inconsistency that no better coefficient can repair; a residual
    // that collapses when the coefficient is smoothed is coefficient noise
    // feeding the inconsistency. This substitution is only possible inside this
    // assembly, which is why the diagnostic lives here and not in a standalone
    // test application.
    //
    // NORMS ARE TAKEN OVER A FIXED PHYSICAL BAND |psi_exact| <= W, NOT over a
    // multiple of h. The same argument the mollifier block in the case
    // fvSolution records: a norm whose domain moves with the mesh has no
    // convergence order, because the order of a sequence of norms over
    // different regions is not the order of anything.
    //
    // VALIDITY. The reference is exact only while psi is still the initial
    // circle, i.e. at the FIRST assembly of the first time step. psiRefDefect
    // (max |psi - psi_exact| over the band) is printed on the same line so the
    // instant at which the reference stops describing the field is visible in
    // the log rather than assumed. Take the first printed line per run.
    //
    // Guarded, not switched: it runs only where the closed form is the truth
    // (implicitSphere + shear2D, serial), and prints one warning and skips
    // otherwise. Parallel is excluded on purpose -- a coupled patch field holds
    // the NEIGHBOUR CELL value, not the face value, so an analytic normal
    // written at the face centre would corrupt fvc::flux exactly as documented
    // in velocityModel::setVelocity.
    {
        const fvSolution& fvSolutionDict(mesh);
        const dictionary& levelSetDict = fvSolutionDict.subDict("levelSet");
        const dictionary& velocityDict = fvSolutionDict.subDict("velocityModel");

        const word surfaceType =
            levelSetDict.subOrEmptyDict("implicitSurface")
                .getOrDefault<word>("type", "none");
        const word velocityType =
            velocityDict.getOrDefault<word>("type", "none");

        const bool referenceAvailable =
            (surfaceType == "implicitSphere")
         && (velocityType == "shear2D")
         && !Pstream::parRun();

        static bool skipReported = false;

        if (!referenceAvailable)
        {
            if (!skipReported)
            {
                skipReported = true;
                Info<< "sdplsRdiv exactRef: SKIPPED (implicitSurface '"
                    << surfaceType << "', velocityModel '" << velocityType
                    << "', parRun " << Switch(Pstream::parRun())
                    << "); the closed form div(w) = n.grad(U).n + (U.n)(d-1)/r"
                    << " is derived for implicitSphere + shear2D in serial."
                    << endl;
            }
        }
        else
        {
            const dictionary& surfaceDict =
                levelSetDict.subDict("implicitSurface");
            const vector centre = surfaceDict.get<vector>("center");
            const scalar radius = surfaceDict.get<scalar>("radius");

            // Rebuilt, not looked up: the velocity model object is not visible
            // from here. The reconstruction is VERIFIED against the U field
            // below, which checks the model type, the oscillation switch and
            // tau in one assertion instead of trusting three dictionary reads.
            // getOrDefault<word> + explicit Switch: see the note at the
            // compositeFlux read above -- a string-literal default to
            // getOrDefault<Switch> is silently TRUE.
            const Switch oscillating
            (
                velocityDict.getOrDefault<word>("oscillation", word("on"))
            );
            const scalar tau =
                velocityDict.getOrDefault<scalar>
                ("tau", mesh.time().endTime().value());
            const scalar tNow = mesh.time().timeOutputValue();
            const scalar fOsc =
                oscillating ? Foam::cos(M_PI*tNow/tau) : 1.0;

            const label nGeomD = mesh.nGeometricD();
            const vectorField& Cc = mesh.C().primitiveField();

            // n_exact = (x - c)/|x - c|, at cell centres and at boundary face
            // centres (no coupled patches here: serial only, guarded above).
            volVectorField nEx
            (
                IOobject
                (
                    "nExactSdpls",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedVector(dimless, Zero)
            );

            {
                vectorField& nExI = nEx.primitiveFieldRef();
                forAll(nExI, cellI)
                {
                    const vector dvec = Cc[cellI] - centre;
                    nExI[cellI] = dvec/max(mag(dvec), SMALL);
                }

                const auto& CfBf = mesh.Cf().boundaryField();
                auto& nExBf = nEx.boundaryFieldRef();
                forAll(nExBf, patchI)
                {
                    const vectorField& Cfp = CfBf[patchI];
                    vectorField& nExP = nExBf[patchI];
                    forAll(nExP, faceI)
                    {
                        const vector dvec = Cfp[faceI] - centre;
                        nExP[faceI] = dvec/max(mag(dvec), SMALL);
                    }
                }
            }

            // The analytic reference (*), the analytic velocity used to verify
            // the reconstruction, and the exact signed distance defining the
            // measurement band -- one loop, one set of trigonometric calls.
            volScalarField divWexact
            (
                IOobject
                (
                    "divWexactSdpls",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 0)
            );

            scalarField& divWexactI = divWexact.primitiveFieldRef();
            vectorField Uexact(mesh.nCells(), Zero);
            scalarField psiExact(mesh.nCells(), 0);

            forAll(divWexactI, cellI)
            {
                const scalar x = Cc[cellI].x();
                const scalar y = Cc[cellI].y();

                const scalar sx  = Foam::sin(M_PI*x);
                const scalar sy  = Foam::sin(M_PI*y);
                const scalar s2x = Foam::sin(2*M_PI*x);
                const scalar s2y = Foam::sin(2*M_PI*y);
                const scalar c2x = Foam::cos(2*M_PI*x);
                const scalar c2y = Foam::cos(2*M_PI*y);

                const scalar u =  s2y*sx*sx*fOsc;
                const scalar v = -s2x*sy*sy*fOsc;
                Uexact[cellI] = vector(u, v, 0);

                const scalar dudx =  M_PI*s2x*s2y*fOsc;
                const scalar dudy =  2*M_PI*c2y*sx*sx*fOsc;
                const scalar dvdx = -2*M_PI*c2x*sy*sy*fOsc;
                const scalar dvdy = -M_PI*s2x*s2y*fOsc;

                const vector dvec = Cc[cellI] - centre;
                const scalar r = max(mag(dvec), SMALL);
                const vector n = dvec/r;
                psiExact[cellI] = r - radius;

                const scalar aExact =
                    n.x()*n.x()*dudx + n.x()*n.y()*dvdx
                  + n.y()*n.x()*dudy + n.y()*n.y()*dvdy;

                divWexactI[cellI] =
                    aExact + (u*n.x() + v*n.y())*scalar(nGeomD - 1)/r;
            }

            // The reconstruction must BE the field the solver is advecting
            // with, or every number below describes a different problem. This
            // catches a changed velocity model, a changed oscillation switch
            // and a changed tau at once.
            const scalar Uscale = max(gMax(mag(U.primitiveField())), SMALL);
            const scalar UrefDefect = gMax(mag(U.primitiveField() - Uexact));
            if (UrefDefect > 1e-10*Uscale)
            {
                FatalErrorInFunction
                    << "The analytic shear2D reconstruction does not reproduce"
                    << " the solver's U field: max|U - U_exact| = "
                    << UrefDefect << " (max|U| = " << Uscale << ")." << nl
                    << "The closed-form div(w) reference would then be measured"
                    << " against a velocity the run does not use."
                    << exit(FatalError);
            }

            // w built from the EXACT normal and the same U. The field NAMES
            // must stay "wSdpls"/"phiW": fvc::flux looks up
            // interpolate(wSdpls) and fvc::div(phiW, psi) looks up the
            // div(phiW,psi) scheme, so renaming them would silently change the
            // discretisation being measured. Neither field is registered
            // (registerObject = false in this constructor), so the duplicate
            // names cost nothing.
            volVectorField const wEx("wSdpls", (nEx & U)*nEx);
            surfaceScalarField phiWex("phiW", fvc::flux(wEx));

            volScalarField const divWh("divWh", fvc::div(phiW));
            volScalarField const divWopEx("divWopEx", fvc::div(phiWex));

            // Product-rule residual with the discrete coefficient (the same
            // quantity the block above prints globally) and with the exact one,
            // differing in NOTHING but w.
            volScalarField const prResidH
            (
                "prResidH",
                fvc::div(phiW, psi) - psi*divWh - (w & gradPsi)
            );
            volScalarField const prResidEx
            (
                "prResidEx",
                fvc::div(phiWex, psi) - psi*divWopEx - (wEx & gradPsi)
            );
            volScalarField const prRef("prRef", fvc::div(phiW, psi));

            // Fixed physical band, from the EXACT distance: identical annulus
            // on every mesh of the ladder, and independent of the state of the
            // solution.
            const scalar bandWidth =
                sourceDict_.getOrDefault<scalar>("exactRefBandWidth", 0.05);

            const scalarField& V = mesh.V();
            scalarField bandMask(mesh.nCells(), 0.0);
            label nBandCells = 0;
            forAll(bandMask, cellI)
            {
                if (mag(psiExact[cellI]) <= bandWidth)
                {
                    bandMask[cellI] = 1.0;
                    ++nBandCells;
                }
            }
            const scalar bandVol = max(gSum(bandMask*V), VSMALL);

            auto bandL2 = [&](const scalarField& e)
            {
                return Foam::sqrt(gSum(bandMask*sqr(e)*V)/bandVol);
            };
            auto bandLinf = [&](const scalarField& e)
            {
                return gMax(bandMask*mag(e));
            };

            const scalarField eTot
            (
                divWh.primitiveField() - divWexact.primitiveField()
            );
            const scalarField eOp
            (
                divWopEx.primitiveField() - divWexact.primitiveField()
            );
            const scalarField eNormal
            (
                divWh.primitiveField() - divWopEx.primitiveField()
            );

            // Same "unstructured FVM length approximation" the solver's error
            // CSV writes as DELTA_X, so the two tables share an abscissa.
            const scalar deltaX =
                gMin(pow(mesh.deltaCoeffs().internalField(), -1)());

            Info<< "sdplsRdiv exactRef"
                << " t=" << tNow
                << " h=" << deltaX
                << " bandWidth=" << bandWidth
                << " nBandCells=" << nBandCells
                << " psiRefDefect=" << bandLinf(psi.primitiveField() - psiExact)
                << " divW_err_Linf=" << bandLinf(eTot)
                << " divW_err_L2=" << bandL2(eTot)
                << " divW_opErr_Linf=" << bandLinf(eOp)
                << " divW_opErr_L2=" << bandL2(eOp)
                << " divW_normalErr_Linf=" << bandLinf(eNormal)
                << " divW_normalErr_L2=" << bandL2(eNormal)
                << " divW_exact_Linf=" << bandLinf(divWexact.primitiveField())
                << " divW_exact_L2=" << bandL2(divWexact.primitiveField())
                << " prResid_h_Linf=" << bandLinf(prResidH.primitiveField())
                << " prResid_h_L2=" << bandL2(prResidH.primitiveField())
                << " prResid_exactW_Linf="
                << bandLinf(prResidEx.primitiveField())
                << " prResid_exactW_L2=" << bandL2(prResidEx.primitiveField())
                << " prRef_L2=" << bandL2(prRef.primitiveField())
                << endl;

            // ---- WHERE THE RESIDUAL LIVES ----------------------------------
            //
            // The global norms the block above this one prints turned out NOT
            // to describe the band. MEASURED on the short probe
            // (config/sdplsRdivProductRule2Dvortex.yaml, END_TIME 0.1, one full
            // reversal, N = 32/64/128/256), at the END of the run:
            //   GLOBAL residual L2  1.515e-02  8.802e-03  2.260e-02  6.852e-02
            //   BAND   residual L2  8.890e-03  2.118e-03  5.068e-04  1.267e-04
            // The global norm turns and grows (this is the -0.789 order already
            // recorded in the study config); the fixed-annulus band norm falls
            // at h^2.04 through the same runs. The same short study's own error
            // table agrees: gradientErrorBand +2.003, gradientErrorBandHalf
            // +1.523, volumeErrorHalf +1.926, while the GLOBAL max gradient
            // error gradientErrorMax runs at -1.347. So the identity fails
            // somewhere OUTSIDE the interface neighbourhood and a global norm
            // cannot say where.
            //
            // Four disjoint regions are therefore reported, plus the position
            // of the global maximum:
            //   inner  psi_exact < -W : inside the initial circle. Contains the
            //          CENTRE r = 0, where the signed distance has a kink, nHat
            //          is undefined and div(n) = (d-1)/r is singular. A
            //          non-converging global Linf is expected to sit here and
            //          it must be identified, not attributed to the scheme.
            //   band   |psi_exact| <= W (as above)
            //   outer  psi_exact > W
            //   hband  |psi| <= 2h of the CURRENT psi: the moving interface
            //          neighbourhood, which is where gradientErrorBand -- the
            //          metric Rdiv diverges in -- is measured. Its domain moves
            //          with h, so it carries no convergence order of its own;
            //          it is reported because it is the region the failure is
            //          reported in. || |grad psi| - 1 || is printed over the
            //          same cells so the residual and the failure metric can be
            //          read off one line at one instant.
            const scalarField& prH = prResidH.primitiveField();
            const scalarField& psiI = psi.primitiveField();

            scalarField innerMask(mesh.nCells(), 0.0);
            scalarField outerMask(mesh.nCells(), 0.0);
            scalarField hbandMask(mesh.nCells(), 0.0);
            label nHband = 0;

            forAll(psiExact, cellI)
            {
                if (psiExact[cellI] < -bandWidth) innerMask[cellI] = 1.0;
                if (psiExact[cellI] >  bandWidth) outerMask[cellI] = 1.0;
                if (mag(psiI[cellI]) <= 2*deltaX)
                {
                    hbandMask[cellI] = 1.0;
                    ++nHband;
                }
            }

            auto regionL2 = [&](const scalarField& m, const scalarField& e)
            {
                return Foam::sqrt(gSum(m*sqr(e)*V)/max(gSum(m*V), VSMALL));
            };
            auto regionLinf = [&](const scalarField& m, const scalarField& e)
            {
                return gMax(m*mag(e));
            };

            const scalarField gradErr
            (
                mag(mag(gradPsi.primitiveField()) - 1.0)
            );

            // Serial-only block (guarded above), so a plain loop is the global
            // argmax.
            label iMax = 0;
            scalar vMax = 0;
            forAll(prH, cellI)
            {
                if (mag(prH[cellI]) > vMax)
                {
                    vMax = mag(prH[cellI]);
                    iMax = cellI;
                }
            }
            const scalar rMax = mag(Cc[iMax] - centre);

            Info<< "sdplsRdiv exactRefLoc"
                << " t=" << tNow
                << " h=" << deltaX
                << " inner_prResid_Linf=" << regionLinf(innerMask, prH)
                << " inner_prResid_L2=" << regionL2(innerMask, prH)
                << " band_prResid_Linf=" << bandLinf(prH)
                << " band_prResid_L2=" << bandL2(prH)
                << " outer_prResid_Linf=" << regionLinf(outerMask, prH)
                << " outer_prResid_L2=" << regionL2(outerMask, prH)
                << " hband_nCells=" << nHband
                << " hband_prResid_Linf=" << regionLinf(hbandMask, prH)
                << " hband_prResid_L2=" << regionL2(hbandMask, prH)
                << " hband_gradErr_Linf=" << regionLinf(hbandMask, gradErr)
                << " hband_gradErr_L2=" << regionL2(hbandMask, gradErr)
                << " argmax_val=" << vMax
                << " argmax_x=" << Cc[iMax].x()
                << " argmax_y=" << Cc[iMax].y()
                << " argmax_r=" << rMax
                << " argmax_rOverH=" << rMax/deltaX
                << " argmax_psi=" << psiI[iMax]
                << " argmax_psiExact=" << psiExact[iMax]
                << endl;
        }
    }
    // -------------------------------------------------------------------------

    if (compositeFlux)
    {
        surfaceScalarField const phiComposite("phiW", 2*phi - phiW);

        return
        (
            fvm::div(phi, psi)
          - fvm::div(phiComposite, psi)
          - fvm::Sp(fvc::div(phiW) - strain(), psi)
        );
    }

    return
    (
        fvm::div(phiW, psi)
      - fvm::div(phi, psi)
      - fvm::Sp(fvc::div(phiW) - strain(), psi)
    );
}

} // End namespace Foam

// ************************************************************************* //
