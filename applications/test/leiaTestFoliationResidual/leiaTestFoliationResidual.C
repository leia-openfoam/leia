/*---------------------------------------------------------------------------*\
Application
    leiaTestFoliationResidual

Description
    Is the loss of the PARALLEL FOLIATION what makes the delivered curvature go
    bad in a coupled run?

    The parallel-surface inverse used by every foot-point curvature delivery is
    exact if and only if the level sets of psi are the parallel offsets of the
    interface -- equivalently beta := |grad psi| is constant ALONG each level set.
    It does NOT need |grad psi| = 1 (plan sec. 11.1). When the foliation is not
    parallel, the residual of the inverse is the integral along the normal of

        D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2 = -beta Lap_Gamma(1/beta)

    so the delivered interface curvature is biased by d*D + O(d^2), FIRST order in
    the offset d (plan sec. 11.2).

    This utility evaluates D on a SAVED psi field, using the volumetric identity

        Lap_Gamma u = Lap u - n.(grad grad u).n - kappa (n.grad u),  n = grad psi/beta

    with u = ln beta, so no third-derivative fit of psi is needed in the solver --
    the wide finite-volume stencil is acceptable here BECAUSE this is a diagnostic
    and not a force term. Noise costs accuracy in the reported number; it cannot
    destabilise anything.

    Reported per time, band-restricted to |psi| < bandWidth*h:
        beta         L2/Linf of (beta - 1)          -- how far from a signed distance
        aT           L2/Linf of |P grad(ln beta)|   -- the non-parallelism VECTOR
                                                      (equals |P H n|/beta; zero iff
                                                       the level sets are parallel)
        D            L2/Linf                        [1/m^2]
        biasHalfH    L2 of |D|*h/2                  [1/m] the bias a face at the
                                                      typical offset h/2 inherits
        biasPsi      L2 of |D * psi/beta|           [1/m] using each cell's own
                                                      offset

    THE DECISION RULE, stated before looking at the data: in the cell-mean N=128
    coupled run the band curvature error grows 70 -> 187 -> 1100 1/m over
    t = 0.05 -> 0.08 -> 0.10 while min(beta) in the band falls 0.988 -> 0.858 ->
    0.724. If biasHalfH is of the same ORDER as that error and grows with it, the
    foliation residual explains the growth and the next lever is psi (keep the
    level sets parallel). If biasHalfH stays one or more orders BELOW the measured
    error, the correlation is spurious and the lever is elsewhere.

    Runs over every time directory holding psi. Parallel-safe.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "bandWidth", "scalar",
        "Band half-width in cells for the norms (default 2)."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const scalar bandWidth = args.getOrDefault<scalar>("bandWidth", 2.0);

    // Representative cell size (uniform hex cases).
    const boundBox bb(mesh.points());
    const label nd = mesh.nGeometricD();
    const scalar Lx = bb.max().x() - bb.min().x();
    const scalar Ly = bb.max().y() - bb.min().y();
    const scalar Lz = bb.max().z() - bb.min().z();
    scalar nCellsG = mesh.nCells();
    reduce(nCellsG, sumOp<scalar>());
    const scalar h =
        (nd == 2) ? Foam::sqrt(Lx*Ly/nCellsG) : Foam::cbrt(Lx*Ly*Lz/nCellsG);

    autoPtr<OFstream> os;
    if (Pstream::master())
    {
        // globalPath, not path(): in a decomposed run the latter is the
        // processor0 directory and the CSV would hide inside it.
        os.reset
        (
            new OFstream(runTime.globalPath()/"leiaTestFoliationResidual.csv")
        );
        os() << "TIME,DELTA_X,N_BAND,BETA_MIN,BETA_MAX,BETA_ERR_L2,"
             << "AT_L2,AT_LINF,D_L2,D_LINF,BIAS_HALFH_L2,BIAS_PSI_L2" << nl;
    }

    forAll(timeDirs, ti)
    {
        runTime.setTime(timeDirs[ti], ti);

        IOobject psiIO
        (
            "psi", runTime.timeName(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        );
        if (!psiIO.typeHeaderOk<volScalarField>(true)) { continue; }

        Info<< "Time = " << runTime.timeName() << endl;
        const volScalarField psi(psiIO, mesh);

        const volVectorField g(fvc::grad(psi));
        volScalarField beta("beta", mag(g));
        beta.max(SMALL);
        const volVectorField n(g/beta);

        volScalarField lnBeta("lnBeta", Foam::log(beta));
        lnBeta.correctBoundaryConditions();

        const volVectorField gLn(fvc::grad(lnBeta));
        const volTensorField HLn(fvc::grad(gLn));

        // kappa = div(grad psi/|grad psi|) from the field gradient and Hessian
        // rather than fvc::div(n): the same algebraic convention the solver uses,
        // and it needs no divScheme entry in the case's fvSchemes (the saved
        // droplet cases do not carry one for this expression).
        const volTensorField Hpsi(fvc::grad(g));
        const volScalarField kappa
        (
            (tr(Hpsi)*magSqr(g) - (g & (Hpsi & g)))/pow3(beta)
        );

        // tangential part of grad(ln beta): the non-parallelism vector a
        const volVectorField aT(gLn - (n & gLn)*n);

        // D = Lap_Gamma(ln beta) - |grad_Gamma ln beta|^2
        const volScalarField D
        (
            (tr(HLn) - (n & (HLn & n)) - kappa*(n & gLn)) - magSqr(aT)
        );

        // band: |psi| < bandWidth*h
        const scalarField& psiIn = psi.primitiveField();
        const scalarField& V = mesh.V();
        const scalarField& bIn = beta.primitiveField();
        const scalarField& DIn = D.primitiveField();
        const scalarField aTm(mag(aT)().primitiveField());

        scalar sumV = 0, bE = 0, aL = 0, dL = 0, bhL = 0, bpL = 0;
        scalar aI = 0, dI = 0, bMin = GREAT, bMax = -GREAT;
        label nBand = 0;

        forAll(psiIn, c)
        {
            if (mag(psiIn[c]) > bandWidth*h) { continue; }
            ++nBand;
            sumV += V[c];
            bE += V[c]*sqr(bIn[c] - 1.0);
            aL += V[c]*sqr(aTm[c]);
            dL += V[c]*sqr(DIn[c]);
            bhL += V[c]*sqr(mag(DIn[c])*0.5*h);
            bpL += V[c]*sqr(DIn[c]*psiIn[c]/bIn[c]);
            aI = max(aI, aTm[c]);
            dI = max(dI, mag(DIn[c]));
            bMin = min(bMin, bIn[c]);
            bMax = max(bMax, bIn[c]);
        }
        reduce(nBand, sumOp<label>());
        reduce(sumV, sumOp<scalar>());
        reduce(bE, sumOp<scalar>()); reduce(aL, sumOp<scalar>());
        reduce(dL, sumOp<scalar>()); reduce(bhL, sumOp<scalar>());
        reduce(bpL, sumOp<scalar>());
        reduce(aI, maxOp<scalar>()); reduce(dI, maxOp<scalar>());
        reduce(bMin, minOp<scalar>()); reduce(bMax, maxOp<scalar>());

        if (sumV <= 0) { continue; }
        const scalar rb = Foam::sqrt(bE/sumV);
        const scalar ra = Foam::sqrt(aL/sumV);
        const scalar rd = Foam::sqrt(dL/sumV);
        const scalar rbh = Foam::sqrt(bhL/sumV);
        const scalar rbp = Foam::sqrt(bpL/sumV);

        Info<< "  band cells " << nBand
            << ", beta in [" << bMin << ", " << bMax << "]" << nl
            << "  |beta-1| L2 = " << rb
            << ", |a_T| L2 = " << ra << " (max " << aI << ") 1/m" << nl
            << "  D  L2 = " << rd << " (max " << dI << ") 1/m^2" << nl
            << "  bias at h/2 = " << rbh << " 1/m, bias at psi/beta = "
            << rbp << " 1/m" << endl;

        if (Pstream::master())
        {
            os() << runTime.timeName() << ',' << h << ',' << nBand << ','
                 << bMin << ',' << bMax << ',' << rb << ','
                 << ra << ',' << aI << ',' << rd << ',' << dI << ','
                 << rbh << ',' << rbp << nl;
        }
    }

    Info<< "End\n" << endl;
    return 0;
}
