/*---------------------------------------------------------------------------*\
Application
    leiaTestDeparturePoint

Description
    Direct verification of the departure-point (x_d) kernel of the
    semi-Lagrangian advection -- the SHIPPED code path, not a replica.

    The prescribed velocity is a time-modulated solid-body rotation about the
    domain centre c,

        u(x, t) = cos(omega t) * Omega * e_z ^ (x - c),

    whose backward characteristic is EXACT: the trajectory arriving at x_c at
    t^{n+1} departs from a rotation of x_c about c by the angle

        dphi = (Omega/omega) * [ sin(omega t^{n+1}) - sin(omega t^n) ].

    The feet the scheme ACTUALLY uses are extracted by advecting the linear
    fields psiX = x and psiY = y through slAdvection::advect: the quadratic
    least-squares fit reproduces linear fields exactly, so after one step
    psiX(x_c) = x_d.x and psiY(x_c) = y_d exactly (up to the fit's roundoff).

    Two velocity-supply variants, matching the two production callers:

      exact : Unew = u(., t^{n+1}), Uold = u(., t^n) -- what the prescribed-
              velocity kinematic solvers supply after ++runTime.
      ab2   : Unew = u^n + (dt/dt0)(u^n - u^{n-1}), Uold = u(., t^n) -- the
              coupled solver's Adams-Bashforth substitution (slAlphaEqn.H).

    Swept: omega*dt in {0, 0.105, 0.26, 0.52, 1.04} (0.52 = the capillary
    operating point omega_max*dt at CAPILLARY_DT_COEFF = 0.010861 -- the
    stiffest 2h capillary mode is resolved by ~12 steps at EVERY N), and a
    dt-halving ladder at fixed omega for the per-step order (expected 3:
    the kernel truncates the Taylor trace after the material-acceleration
    term).

    Error norms are restricted to cells farther than 4 h from any boundary
    (the CPC fit stencils of near-wall cells are one-sided).

    Writes departurePointErrors.csv into the case directory:
      VARIANT,OMEGA_DT,DT,CFL,N_INTERIOR,E_FOOT_L2,E_FOOT_LINF,E_REL_LINF

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "slAdvection.H"
#include "OFstream.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The advection object of the production solvers (reads
    // fvSolution.levelSet.semiLagrangian of the case).
    autoPtr<slAdvection> slAdv = slAdvection::New(mesh);

    const vectorField& C = mesh.C().primitiveField();
    const boundBox bb(mesh.points());
    const vector centre = bb.centre();

    // Interior mask: > 4 h from every boundary (fit stencils are one-sided
    // near walls; the kernel itself is stencil-free but the psi evaluation
    // at the foot is not).
    const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());
    boolList interior(mesh.nCells(), false);
    forAll(C, c)
    {
        const vector d = C[c] - centre;
        const scalar Lx = 0.5*(bb.max().x() - bb.min().x());
        const scalar Ly = 0.5*(bb.max().y() - bb.min().y());
        interior[c] =
            (mag(d.x()) < Lx - 4.0*h) && (mag(d.y()) < Ly - 4.0*h);
    }
    label nInt = 0;
    forAll(interior, c) { if (interior[c]) ++nInt; }
    reduce(nInt, sumOp<label>());

    // Rotation rate chosen so the advective CFL matches the coupled runs'
    // order of magnitude: |u|max = OmegaRot * Rmax, CFL = |u|max*dt/h.
    // The foot displacement stays well inside the fit stencil.
    const scalar Rmax = 0.5*Foam::sqrt
    (
        sqr(bb.max().x() - bb.min().x()) + sqr(bb.max().y() - bb.min().y())
    );

    auto uAt = [&](const scalar t, const scalar omega, const scalar OmegaRot)
    {
        auto tfld = tmp<volVectorField>::New
        (
            IOobject("uRot", runTime.timeName(), mesh, IOobject::NO_READ,
                     IOobject::NO_WRITE),
            mesh,
            dimensionedVector(dimVelocity, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        vectorField& u = tfld.ref().primitiveFieldRef();
        const scalar f = (omega > VSMALL) ? Foam::cos(omega*t) : 1.0;
        forAll(u, c)
        {
            const vector d = C[c] - centre;
            u[c] = f*OmegaRot*vector(-d.y(), d.x(), 0);
        }
        tfld.ref().correctBoundaryConditions();
        return tfld;
    };

    // psi = x and psi = y as separate scalar fields with zeroGradient walls.
    auto makeCoord = [&](const direction cmpt, const word& nm)
    {
        auto tfld = tmp<volScalarField>::New
        (
            IOobject(nm, runTime.timeName(), mesh, IOobject::NO_READ,
                     IOobject::NO_WRITE),
            mesh,
            dimensionedScalar(dimLength, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        scalarField& s = tfld.ref().primitiveFieldRef();
        forAll(s, c) { s[c] = C[c][cmpt]; }
        tfld.ref().correctBoundaryConditions();
        return tfld;
    };

    OFstream csv(runTime.globalPath()/"departurePointErrors.csv");
    csv << "VARIANT,OMEGA_DT,DT,CFL,N_INTERIOR,E_FOOT_L2,E_FOOT_LINF,"
        << "E_REL_LINF\n";

    // Base step: the coupled N=128 operating point's dt magnitude relative
    // to this mesh is irrelevant -- what matters is (omega*dt, CFL); pick dt
    // so CFL = 0.5 at ladder level 0 and halve it twice.
    const scalar dt0base = 0.5*h/(1.0*Rmax);   // OmegaRot = 1 rad/s below

    for (const scalar omegaDt :
         {0.0, 0.105, 0.26, 0.52, 1.04})
    {
        for (label lvl = 0; lvl < 3; ++lvl)
        {
            const scalar dt = dt0base/Foam::pow(2.0, lvl);
            const scalar OmegaRot = 1.0;             // rad/s
            const scalar omega = omegaDt/dt;         // omega*dt fixed per row
            const scalar CFL = OmegaRot*Rmax*dt/h;

            // Arrival at t^{n+1} = dt (t^n = 0, t^{n-1} = -dt).
            const scalar tNew = dt, tOld = 0.0, tPrev = -dt;

            // Exact backward rotation angle over [t^n, t^{n+1}].
            const scalar dphi =
                (omega > VSMALL)
              ? (OmegaRot/omega)*(Foam::sin(omega*tNew) - Foam::sin(omega*tOld))
              : OmegaRot*dt;
            const scalar cphi = Foam::cos(dphi), sphi = Foam::sin(dphi);

            for (const word& variant : {word("exact"), word("ab2")})
            {
                tmp<volVectorField> tUold = uAt(tOld, omega, OmegaRot);
                tmp<volVectorField> tUnew;
                if (variant == "exact")
                {
                    tUnew = uAt(tNew, omega, OmegaRot);
                }
                else
                {
                    // AB2 substitution at fixed step: 2 u^n - u^{n-1}.
                    tmp<volVectorField> tUprev = uAt(tPrev, omega, OmegaRot);
                    tUnew = tmp<volVectorField>::New
                    (
                        "UextStar",
                        2.0*tUold() - tUprev()
                    );
                }

                // runTime.deltaT drives the kernel's dt.
                runTime.setDeltaT(dt, false);

                tmp<volScalarField> tPsiX = makeCoord(0, "psiX");
                tmp<volScalarField> tPsiY = makeCoord(1, "psiY");
                slAdv->advect(tPsiX.ref(), tUnew(), tUold());
                slAdv->advect(tPsiY.ref(), tUnew(), tUold());

                const scalarField& fx = tPsiX().primitiveField();
                const scalarField& fy = tPsiY().primitiveField();

                scalar s2 = 0, sInf = 0;
                forAll(C, c)
                {
                    if (!interior[c]) { continue; }
                    const vector d = C[c] - centre;
                    // Backward rotation by dphi (departure lies BEHIND the
                    // rotation): R(-dphi) d.
                    const vector xdEx = centre + vector
                    (
                        cphi*d.x() + sphi*d.y(),
                       -sphi*d.x() + cphi*d.y(),
                        d.z()
                    );
                    const scalar e = Foam::sqrt
                    (
                        sqr(fx[c] - xdEx.x()) + sqr(fy[c] - xdEx.y())
                    );
                    s2 += e*e; sInf = max(sInf, e);
                }
                reduce(s2, sumOp<scalar>());
                reduce(sInf, maxOp<scalar>());
                const scalar eL2 = Foam::sqrt(s2/max(nInt, 1));
                const scalar disp = max(OmegaRot*Rmax*dt, VSMALL);

                csv << variant << ',' << omegaDt << ',' << dt << ','
                    << CFL << ',' << nInt << ','
                    << eL2 << ',' << sInf << ',' << sInf/disp << nl;

                Info<< variant << " omegaDt=" << omegaDt << " dt=" << dt
                    << " CFL=" << CFL
                    << " E_L2=" << eL2 << " E_Linf=" << sInf
                    << " rel=" << sInf/disp << endl;
            }
        }
    }

    Info<< "\nWrote " << (runTime.globalPath()/"departurePointErrors.csv")
        << nl << "End" << endl;

    return 0;
}

// ************************************************************************* //
