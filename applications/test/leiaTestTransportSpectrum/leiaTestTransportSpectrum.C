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

Application
    leiaTestTransportSpectrum

Description
    FROZEN-VELOCITY POWER ITERATION of the semi-Lagrangian transport operator.

    WHY THIS EXISTS. The amplification bound Lambda_c is a ONE-STEP, per-cell
    sensitivity: it bounds |psi^{n+1}_c| by Lambda_c times the stencil maximum.
    It does NOT decide repeated stability, and there is an exact counterexample:
    the anchored quadratic evaluated at x_i - C h on a uniform 1D grid is
    Lax-Wendroff,

        psi_i^{n+1} = (C^2+C)/2 psi_{i-1} + (1-C^2) psi_i + (C^2-C)/2 psi_{i+1},

    whose row norm is Lambda = 1 + C - C^2 > 1 while its Fourier modulus obeys
    |G|^2 = 1 - 4C^2(1-C^2) sin^4(theta/2) <= 1. So Lambda > 1 is compatible
    with L2 stability, and the correlation between Lambda_max = 1.2608 on pMesh
    (which fails) and 1.05 on hexahedra (which do not) is CORRELATION, not a
    mechanism.

    WHAT THIS MEASURES INSTEAD. With the velocity FROZEN and the quasi-monotone
    clip OFF, the transport stage is an exactly LINEAR operator B: the departure
    feet are fixed, and the per-cell fit rank is decided from stencil POSITIONS
    alone (quadraticPivotTol), never from the field. This application applies B
    repeatedly to a seed field, renormalising each time, and reports

      * the per-iteration growth factor  ||B psi_k|| / ||psi_k||, whose limit is
        the spectral radius rho(B) -- norm independent, so L1 and L2 agree in the
        limit and both are reported;
      * the finite-time power norms ||B^n psi_0|| / ||psi_0|| WITHOUT
        renormalisation, which lower-bound ||B^n||. For a NON-NORMAL operator
        these can greatly exceed rho^n, and transient growth is what destroys a
        level set in a finite horizon, so they are reported separately.

    rho(B) > 1 proves asymptotic growth of this operator on this mesh.
    rho(B) <= 1 does NOT prove robust stability for a non-normal B; read the
    transient column as well.

    SEEDS. The dominant mode is unknown a priori, so the seed must not privilege
    one: `random` fills psi with a reproducible pseudo-random field (identical on
    every rank count, because it is keyed on the cell CENTRE, not the cell
    index -- a index-keyed seed would make the experiment decomposition
    dependent). `checkerboard` uses sign(sin) of a high-frequency function of the
    centre, the closest mesh-agnostic analogue of the 2h mode that the reported
    failure grows. Both are recorded.

Usage
    leiaTestTransportSpectrum [-nIter 200] [-nTransient 40] [-seed random]
        [-writeModes]

    Run in a case whose 0/ holds psi (container and boundary conditions; the
    VALUES are overwritten by the seed) and U. The velocity is read once and
    held frozen. fvSolution/levelSet/semiLagrangian selects the reconstruction
    and the corrector exactly as the solver does.

    Writes leiaTestTransportSpectrum.csv and returns non-zero when rho > 1 + tol.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "slAdvection.H"
#include "slReconstruction.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A reproducible pseudo-random value keyed on a POSITION, so the seed field is
//  identical for every decomposition. An index-keyed generator would give a
//  different field on 1 and on 32 ranks and the comparison would be void.
static inline scalar positionHash(const point& p, const scalar L)
{
    // Three incommensurate frequencies, then the fractional part: cheap, and it
    // decorrelates neighbouring cells at the mesh scale.
    const scalar s =
        std::sin(12.9898*p.x()/L + 78.233*p.y()/L + 37.719*p.z()/L)*43758.5453;
    return 2.0*(s - std::floor(s)) - 1.0;      // in [-1, 1)
}


int main(int argc, char *argv[])
{
    argList::addOption("nIter", "label", "power iterations (default 200)");
    argList::addOption
    (
        "nTransient", "label",
        "un-renormalised iterations for the power norms (default 40)"
    );
    argList::addOption("seed", "word", "random | checkerboard (default random)");
    argList::addOption
    (
        "trace", "word",
        "projectedFlux (default, PRODUCTION) | cellCentred"
    );
    argList::addOption("tol", "scalar", "rho pass tolerance (default 1e-6)");
    argList::addOption
    (
        "mode", "word",
        "spectrum (default, linear power iteration) | growth"
        " (nonlinear finite-time amplification about the physical state)"
    );
    argList::addOption
    (
        "amp", "scalar",
        "growth mode: perturbation amplitude as a MULTIPLE of the mean cell size"
        " (default 1; the review asks for below and above O(h))"
    );
    argList::addBoolOption("writeModes", "write the converged mode as psiMode");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const label nIter = args.getOrDefault<label>("nIter", 200);
    const label nTrans = args.getOrDefault<label>("nTransient", 40);
    const word seedType = args.getOrDefault<word>("seed", "random");
    const word traceType = args.getOrDefault<word>("trace", "projectedFlux");
    if (traceType != "projectedFlux" && traceType != "cellCentred")
    {
        FatalErrorInFunction << "trace must be projectedFlux or cellCentred, got "
            << traceType << exit(FatalError);
    }
    const scalar tol = args.getOrDefault<scalar>("tol", 1e-6);
    const word mode = args.getOrDefault<word>("mode", "spectrum");
    const scalar ampH = args.getOrDefault<scalar>("amp", 1.0);
    if (mode != "spectrum" && mode != "growth")
    {
        FatalErrorInFunction << "mode must be spectrum or growth, got "
            << mode << exit(FatalError);
    }

    if (seedType != "random" && seedType != "checkerboard")
    {
        FatalErrorInFunction << "seed must be random or checkerboard, got "
            << seedType << exit(FatalError);
    }

    // ---- the field container: values are overwritten, the BCs are kept -----
    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    // ---- the FROZEN velocity ----------------------------------------------
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    const volVectorField Ufrozen(U);

    // THE TRACE VELOCITY IS PART OF THE OPERATOR, so it must be the production one.
    // The two-phase solver runs SL_TRACE_VELOCITY projectedFlux: it traces the foot with
    // fvc::reconstruct(phi), NOT with the cell velocity, and STATUS records that the
    // projectedFlux win IS the reconstruct operator. On a uniform stream the two agree on
    // hexahedra and differ on polyhedra, which is exactly where the defect lives -- so
    // measuring the cell-velocity operator would measure the wrong B.
    // phi is registered because slReconstruction looks it up by name when
    // stencilBoundaryFaces is inflowOnly.
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        fvc::flux(Ufrozen)
    );
    volVectorField Utrace
    (
        IOobject("Utrace", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        Ufrozen
    );
    if (traceType == "projectedFlux")
    {
        Utrace == fvc::reconstruct(phi);
    }
    // The velocity is FROZEN, so the old level equals the new one and both are built once.
    const volVectorField UtraceOld(Utrace);

    autoPtr<slAdvection> slAdv = slAdvection::New(mesh);

    const scalar dt = runTime.deltaTValue();
    const scalar Umax = gMax(mag(Ufrozen.primitiveField()));
    const scalar L = Foam::max(gMax(mesh.C().component(0)()), SMALL);

    Info<< nl << "leiaTestTransportSpectrum" << nl
        << "  cells            = " << returnReduce(mesh.nCells(), sumOp<label>()) << nl
        << "  deltaT           = " << dt << nl
        << "  max|U|           = " << Umax << nl
        << "  |d| = max|U| dt  = " << Umax*dt << nl
        << "  seed             = " << seedType << nl
        << "  trace velocity   = " << traceType
        << (traceType == "projectedFlux" ? "  (PRODUCTION)" : "  (control)") << nl
        << "  max|Utrace|      = " << gMax(mag(Utrace.primitiveField())) << nl
        << "  nIter/nTransient = " << nIter << " / " << nTrans << nl << endl;

    // The operator is linear ONLY with no bound acting: every bound is a
    // nonlinear, data-dependent clip. A spectral radius is not defined for a
    // nonlinear map, so the SPECTRUM mode refuses to run with one -- and the
    // external review of 2026-09-09 is explicit that this may not be worked
    // around: "The linear diagnostic Lambda cannot simply be assigned the value
    // one after nonlinear clipping", and "For the proposed nonlinear schemes,
    // test the mapping itself rather than reusing a frozen linear-weight
    // argument."
    //
    // THE GROWTH MODE IS THAT TEST OF THE MAPPING. It runs with any bound.
    const bool bounded = !slAdv->corrector().bound().inert();
    if (mode == "spectrum" && bounded)
    {
        FatalErrorInFunction
            << "valueBound " << slAdv->corrector().bound().type()
            << " is active. The transport stage is then NONLINEAR and a spectral"
            << " radius is not defined for it. Either set"
            << " levelSet/semiLagrangian/valueBound none, or run -mode growth,"
            << " which measures the nonlinear map's finite-time amplification"
            << " about the physical state instead."
            << exit(FatalError);
    }

    // ---- the seed ----------------------------------------------------------
    // EVERY collective is evaluated ONCE, by every rank, OUTSIDE every loop. gAverage is
    // a reduction; calling it inside forAll(f, c) makes each rank enter it nCells times,
    // and the ranks do not hold equal cell counts, so they deadlock. That hung the
    // checkerboard arm for its whole 50-minute timeout while the random arm on the SAME
    // mesh finished in one minute. Same class as a collective inside Pstream::master().
    const scalar hMean = Foam::pow(gAverage(mesh.V()), 1.0/3.0);
    const scalar kWave = constant::mathematical::pi/Foam::max(hMean, SMALL);

    auto fillSeed = [&](volScalarField& f)
    {
        forAll(f, c)
        {
            const point& p = mesh.C()[c];
            if (seedType == "random")
            {
                f[c] = positionHash(p, L);
            }
            else
            {
                // The mesh-agnostic analogue of the 2h mode: a high-frequency
                // sign pattern in space. sign(sin(...)) keeps the amplitude at
                // one so the seed norm is mesh independent.
                const scalar sv = std::sin(kWave*p.x())*std::sin(kWave*p.y())
                                *std::sin(kWave*p.z());
                f[c] = (sv >= 0 ? 1.0 : -1.0);
            }
        }
        f.correctBoundaryConditions();
    };

    auto l2 = [&](const volScalarField& f)
    {
        return Foam::sqrt(gSum(sqr(f.primitiveField())*mesh.V())/gSum(mesh.V()));
    };
    auto l1 = [&](const volScalarField& f)
    {
        return gSum(mag(f.primitiveField())*mesh.V())/gSum(mesh.V());
    };

    OFstream* csv = nullptr;
    if (Pstream::master())
    {
        csv = new OFstream(runTime.path()/"leiaTestTransportSpectrum.csv");
        *csv << "PHASE,ITER,L2,L1,GROWTH_L2,GROWTH_L1" << endl;
    }

    // ======================= MODE: growth ================================= //
    // WHY THIS MODE EXISTS. Every value bound is a nonlinear, data-dependent
    // clip, so the transport stage stops being linear the moment one is active:
    // a spectral radius is undefined, the power iteration has nothing to
    // converge to, and the linear amplification bound Lambda may NOT be
    // reassigned. The external review of 2026-09-09 states the requirement
    // directly -- "For the proposed nonlinear schemes, test the mapping itself
    // rather than reusing a frozen linear-weight argument."
    //
    // WHAT IT MEASURES. The amplification of a PERTURBATION about the physically
    // relevant state, which is what stability means here. Two frozen-velocity
    // trajectories run side by side: one from the case's own psi (the exact
    // signed distance field that leiaSetFields wrote) and one from that field
    // plus a perturbation. The read-out is
    //
    //     g_N = ||psi_pert^N - psi_base^N||_2 / ||psi_pert^0 - psi_base^0||_2 ,
    //     per-step growth = g_N^(1/N) .
    //
    // A random seed is NOT used as the state here, deliberately: it is not a
    // distance function, so a bound premised on the Lipschitz property of one
    // would fire everywhere and the measurement would describe a field the
    // method never transports.
    //
    // ITS SELF-VALIDATION IS PRE-REGISTERED. With valueBound none the difference
    // of two linear trajectories evolves under the SAME linear operator, so the
    // per-step growth must converge to rho(B) -- the already measured 1.00441 on
    // production hexahedra and 1.01028 on production pMesh. If it does not, the
    // instrument is wrong and no verdict about any bound may be read from it.
    if (mode == "growth")
    {
        volScalarField psiBase
        (
            IOobject("psiBase", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            psi
        );
        volScalarField psiPert
        (
            IOobject("psiPert", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            psi
        );

        // The perturbation amplitude is scaled by the mean cell size, so "below
        // and above O(h)" is what -amp selects and the number transfers across
        // meshes. hMean is a collective, evaluated once by every rank above.
        const scalar eps = ampH*hMean;
        volScalarField dseed
        (
            IOobject("dseed", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            psi
        );
        fillSeed(dseed);
        psiPert.primitiveFieldRef() += eps*dseed.primitiveField();
        psiPert.correctBoundaryConditions();

        volScalarField diff
        (
            IOobject("diff", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            psiPert - psiBase
        );
        const scalar e0 = l2(diff);

        Info<< "MODE growth: nonlinear finite-time amplification" << nl
            << "  valueBound       = "
            << slAdv->corrector().bound().type()
            << (bounded ? "  (NONLINEAR)" : "  (linear -> must reproduce rho)") << nl
            << "  base state       = the case's psi (exact signed distance)" << nl
            << "  perturbation     = " << seedType << ", amp = " << ampH
            << " h = " << eps << nl
            << "  ||d psi^0||_2    = " << e0 << nl
            << "  steps            = " << nIter << nl
            << "  iter  ||dpsi||_2  g_N  per-step g_N^(1/N)" << endl;

        if (e0 < SMALL)
        {
            FatalErrorInFunction << "the perturbation is zero" << exit(FatalError);
        }

        scalar gN = 1.0, perStep = 1.0, perStepPrev = 0.0;
        label nDoneG = 0;
        for (label it = 1; it <= nIter; ++it)
        {
            // BOTH trajectories advance under the SAME frozen velocity. Any
            // difference between them is the operator's action on the
            // perturbation, which is the whole measurement.
            slAdv->advect(psiBase, Utrace, UtraceOld);
            slAdv->advect(psiPert, Utrace, UtraceOld);

            diff.primitiveFieldRef() =
                psiPert.primitiveField() - psiBase.primitiveField();
            diff.correctBoundaryConditions();

            // Collectives on EVERY rank, outside every guard.
            const scalar e = l2(diff);
            nDoneG = it;
            if (!std::isfinite(e))
            {
                Info<< "  iter " << it << ": NON-FINITE -- the map overflowed"
                    << endl;
                break;
            }
            gN = e/e0;
            perStep = Foam::pow(gN, 1.0/scalar(it));
            if (it <= 10 || it % 20 == 0 || it == nIter)
            {
                Info<< "  " << it << "  " << e << "  " << gN << "  " << perStep
                    << "   dperStep = "
                    << (it > 1 ? perStep - perStepPrev : 0.0) << endl;
            }
            if (Pstream::master())
            {
                *csv << "growth," << it << ',' << e << ',' << gN << ','
                     << perStep << ',' << gN << endl;
            }
            perStepPrev = perStep;
        }

        const scalar lastIncrement = perStep - perStepPrev;
        const bool growsG = (perStep > 1.0 + tol);

        Info<< nl << "RESULT" << nl
            << "  mode              = growth (the NONLINEAR map itself)" << nl
            << "  valueBound        = "
            << slAdv->corrector().bound().type() << nl
            << "  perturbation      = " << seedType
            << ", amp = " << ampH << " h" << nl
            << "  steps             = " << nDoneG << nl
            << "  total growth g_N  = " << gN << nl
            << "  per-step growth   = " << perStep << nl
            << "  per-step - 1      = " << perStep - 1.0 << nl
            << "  last increment    = " << lastIncrement
            << (Foam::mag(lastIncrement) > 0.1*Foam::mag(perStep - 1.0)
                 ? "   *** NOT CONVERGED: raise -nIter ***" : "")
            << nl << endl;

        Info<< (growsG
              ? "GROWS: the map amplifies a perturbation of the distance field"
              : "BOUNDED: per-step growth <= 1 + tol over this horizon")
            << endl;

        if (Pstream::master()) { delete csv; }
        Info<< "End\n" << endl;
        return growsG ? 1 : 0;
    }

    // ---- PHASE 1: un-renormalised power norms (transient growth) ----------
    // ||B^n psi_0|| / ||psi_0|| lower-bounds ||B^n||. For a non-normal B this
    // can far exceed rho^n, and a level set only has to survive a finite
    // horizon, so this column is reported in its own right.
    Info<< "PHASE 1: finite-time power norms, no renormalisation" << nl
        << "  iter  L2  L1  L2/L2_0" << endl;
    fillSeed(psi);
    const scalar l20 = l2(psi);
    const scalar l10 = l1(psi);
    scalar worstAmp = 1.0;
    label worstAt = 0;
    for (label it = 1; it <= nTrans; ++it)
    {
        slAdv->advect(psi, Utrace, UtraceOld);

        const scalar a = l2(psi), b = l1(psi);
        if (!std::isfinite(a))
        {
            Info<< "  iter " << it << ": NON-FINITE -- the operator overflowed"
                << endl;
            break;
        }
        if (a/l20 > worstAmp) { worstAmp = a/l20; worstAt = it; }
        if (it <= 10 || it % 10 == 0 || it == nTrans)
        {
            Info<< "  " << it << "  " << a
                << "  " << b << "  " << a/l20 << endl;
        }
        if (Pstream::master())
        {
            *csv << "transient," << it << ',' << a << ',' << b << ','
                 << a/l20 << ',' << b/l10 << endl;
        }
    }
    Info<< "  largest power norm ||B^n psi_0||/||psi_0|| = " << worstAmp
        << " at n = " << worstAt << nl << endl;

    // ---- PHASE 2: renormalised power iteration -> rho(B) -------------------
    Info<< "PHASE 2: renormalised power iteration, growth -> rho(B)" << nl
        << "  iter  growth_L2  growth_L1" << endl;
    fillSeed(psi);
    scalar gL2 = 1.0, gL1 = 1.0;
    scalar gL2prev = 0.0;
    label nDone = 0;
    for (label it = 1; it <= nIter; ++it)
    {
        const scalar n2 = l2(psi), n1 = l1(psi);
        if (n2 < SMALL) { Info<< "  seed annihilated at iter " << it << endl; break; }
        psi.primitiveFieldRef() /= n2;
        psi.correctBoundaryConditions();

        slAdv->advect(psi, Utrace, UtraceOld);

        // EVERY collective is evaluated by EVERY rank, OUTSIDE any master guard.
        // l2() and l1() call gSum. Calling them inside Pstream::master() blocks
        // rank 0 in the reduction while the others never enter it -- the exact
        // deadlock CLAUDE.md records for semiImplicitCapillaryForce's diagnostics
        // (four cluster arms alive but silent for 76 minutes). It hung this app at
        // power iteration 1 on 8 ranks.
        const scalar a2 = l2(psi);
        const scalar a1 = l1(psi);
        gL2 = a2;
        gL1 = a1/Foam::max(n1/n2, SMALL);
        nDone = it;
        if (!std::isfinite(gL2))
        {
            Info<< "  iter " << it << ": NON-FINITE" << endl;
            break;
        }
        if (it <= 10 || it % 20 == 0 || it == nIter)
        {
            Info<< "  " << it << "  " << gL2
                << "  " << gL1
                << "   dgrowth = " << (it > 1 ? gL2 - gL2prev : 0.0) << endl;
        }
        if (Pstream::master())
        {
            *csv << "power," << it << ',' << a2 << ',' << a1 << ','
                 << gL2 << ',' << gL1 << endl;
        }
        gL2prev = gL2;
    }

    if (args.found("writeModes"))
    {
        volScalarField mode
        (
            IOobject("psiMode", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            psi
        );
        mode.write();
    }

    const scalar rho = gL2;
    // The power iteration converges only when the increment DECAYS. On the first run
    // of this probe the factor climbed monotonically with GROWING increments over 40
    // iterations, so 40 was far too few and rho <= 1 was not a conclusion. Report the
    // last increment so an unconverged estimate cannot be read as a verdict.
    const scalar dLast = mag(gL2 - gL2prev);
    Info<< nl << "RESULT" << nl
        << "  rho(B) estimate (L2)      = " << rho << nl
        << "  rho(B) estimate (L1)      = " << gL1 << nl
        << "  power iterations used     = " << nDone << nl
        << "  last increment |dgrowth|  = " << dLast
        << (dLast > 1e-8 ? "   NOT CONVERGED -- rho is a rising lower bound" : "   converged")
        << nl
        << "  largest transient growth  = " << worstAmp << " at n = " << worstAt
        << nl
        << "  growth per unit time      = " << (rho > 0 ? Foam::log(rho)/dt : 0)
        << " 1/s" << nl << endl;

    const bool grows = (rho > 1.0 + tol);
    Info<< (grows ? "GROWS: rho > 1 + tol -- this operator amplifies on this mesh"
                  : "BOUNDED: rho <= 1 + tol (non-normal transients may still grow)")
        << nl;

    if (Pstream::master()) { delete csv; }

    Info<< nl << "End" << nl << endl;
    return grows ? 1 : 0;
}

// ************************************************************************* //
