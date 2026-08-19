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
    leiaTestSdplsSource

Description
    Analytic regression test for the SDPLS source term
    (src/leiaLevelSet/sdplsSource/), on the AFFINE exact solution.

    WHY THIS EXISTS. Every SDPLS number ever published from this repository was
    produced with a source whose EXPLICIT branch carried the wrong sign: the
    assembled matrix solved

        Dpsi/Dt = Sp*psi^{n+1} - Sc      instead of      + Sc,

    because fvMatrix::operator== subtracts the right-hand matrix termwise, so a
    contribution written into source_ enters the transport equation with a minus
    (OpenFOAM's own fvm::Su encodes this by doing `source_ -= V*su`). The
    consequence was that `explicit` DOUBLED the gradient drift instead of
    cancelling it, and `strictNegativeSpLinearImplicit` did so on its f_nl > 0
    branch only. `simpleLinearImplicit` (Sc = 0) was never affected, which is
    exactly why the defect survived: the one variant that worked was never the
    one swept over a mesh ladder.

    A single assertion catches the whole family: EVERY LINEARIZATION OF THE
    SAME CONTINUUM TERM must apply the SAME physical source at psi = psi^n.
    That is Layer 1 below, and it covers `explicit`, `simpleLinearImplicit` and
    `strictNegativeSpLinearImplicit`. `exponential` is deliberately OUTSIDE
    that family -- it delivers the exact amplification factor exp(f dt) rather
    than the tangent coefficient f -- so it is measured against its own closed
    form instead, and is exempt from the cross-discretization agreement check.

    THE SETUP is affine, so every discrete operator involved is exact and the
    test has no truncation error to tolerate:

      psi = c x                  =>  grad psi = (c,0,0) exactly, |grad psi| = c
      U   = alpha (x, -y, 0)     =>  grad U constant, div U = 0
      n   = e_x                  =>  R = n.(grad U).n = alpha exactly

    Boundary values are assigned analytically (calculated patches), so the
    gradients are exact on boundary cells too and every cell is a valid sample.

    LAYER 1 -- the applied-source contract (the sign gate).
    For a matrix B returned by fvmsdplsSource, the physical source the equation
    applies is (B.diag*psi - B.source)/V. This must equal f_nl*psi for EVERY
    discretization, where f_nl = R = alpha for `R`, and beta - |grad psi| for
    `beta`. Under the old defect the `explicit` and (for alpha > 0)
    `strictNegativeSpLinearImplicit` rows returned -f_nl*psi.

    LAYER 2 -- composition with ddt.
    One implicit-Euler step of Dpsi/Dt = f_nl psi assembled as
    `fvm::ddt(psi) == fvmsdplsSource(psi, U)`. The system is diagonal, so
    psi^{n+1} = source/diag exactly, and each variant has a closed discrete form:

        Sc-carried (explicit branch):  psi^{n+1} = psi^n (1 + f dt)
        Sp-carried (implicit branch):  psi^{n+1} = psi^n / (1 - f dt)
        exponential (exact factor):    psi^{n+1} = psi^n e^{f dt}

    the first two first-order consistent with the exact psi^n e^{f dt} and the
    third equal to it. Running BOTH signs of alpha exercises both branches of
    strictNegativeSpLinearImplicit, whose split is by the sign of f_nl. The
    case fvSchemes pins `ddt(psi) Euler`, which is also the pairing under which
    `exponential` is exact.

    LAYER 3 -- the gradient law, which is what SDPLS is actually for.
    On the interface the exact statement is D|grad psi|/Dt = |grad psi| (F - a).
    With F = a (the `R` source) the interfacial gradient magnitude is preserved;
    with no source it decays as e^{-a t}. The affine model reproduces both in
    closed form, so the test checks the SIGN OF THE EFFECT directly: one step
    must move |grad psi| toward its analytic target, not away from it. This is
    the assertion that fails loudest under a sign reversal, because a reversed
    source moves the gradient in the wrong direction at first order in dt.

    Exits non-zero on any mismatch (usable as a repo regression test).

Usage
    Run in the dedicated meshed case:
        cd cases/sdplsSourceUnit && blockMesh && leiaTestSdplsSource
    or simply
        bash cases/sdplsSourceUnit/Allrun.sh

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "calculatedFvPatchFields.H"

#include "sdplsSource.H"
#include "sdplsR.H"
#include "sdplsBeta.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * Test harness  * * * * * * * * * * * * * * //

static label nFail = 0;
static label nPass = 0;

static void check
(
    const string& name,
    const scalar got,
    const scalar want,
    const scalar tol = 1e-11
)
{
    const scalar scale = max(mag(want), scalar(1));
    if (!std::isfinite(got) || mag(got - want) > tol*scale)
    {
        Info<< "FAIL " << name << ": got " << got << ", want " << want
            << " (|diff| = " << mag(got - want) << ')' << nl;
        ++nFail;
    }
    else
    {
        ++nPass;
        Info<< "ok   " << name << " (" << got << ')' << nl;
    }
}


//- Max over cells of |a - b|, reduced over ranks.
static scalar maxDiff(const scalarField& a, const scalarField& b)
{
    scalar m = 0;
    forAll(a, c)
    {
        m = max(m, mag(a[c] - b[c]));
    }
    return returnReduce(m, maxOp<scalar>());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The affine fields are exact only if the boundary values are the analytic
    // ones; a decomposed run would additionally need the halo values, and this
    // test is cheap enough to keep serial.
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "leiaTestSdplsSource is a SERIAL analytic test." << exit(FatalError);
    }

    const scalar dt = runTime.deltaTValue();
    const scalarField& V = mesh.V().field();

    // ------------------------------------------------------------------ //
    // Affine exact fields. psi = cSlope*x  ->  |grad psi| = cSlope.
    // ------------------------------------------------------------------ //
    const scalar cSlope = 1.0;          // signed-distance scaling of the plane

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("psi", dimLength, 0),
        calculatedFvPatchScalarField::typeName
    );

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector("U", dimVelocity, Zero),
        calculatedFvPatchVectorField::typeName
    );

    // alpha is set per case below; fill() rewrites both fields.
    auto fill = [&](const scalar alpha)
    {
        forAll(mesh.C(), c)
        {
            const point& x = mesh.C()[c];
            psi[c] = cSlope*x.x();
            U[c] = vector(alpha*x.x(), -alpha*x.y(), 0);
        }
        forAll(psi.boundaryField(), patchi)
        {
            fvPatchScalarField& pp = psi.boundaryFieldRef()[patchi];
            fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];
            const vectorField& Cf = mesh.boundary()[patchi].Cf();
            forAll(pp, i)
            {
                pp[i] = cSlope*Cf[i].x();
                pU[i] = vector(alpha*Cf[i].x(), -alpha*Cf[i].y(), 0);
            }
        }
        // Refresh the old-time level so each case starts from its own state.
        psi.oldTime() == psi;
    };

    // Builds the source-model dictionary; `type` is informational (the model
    // class is chosen by construction, not by the RTS table, so the test is
    // independent of the case's fvSolution).
    auto srcDict = [](const word& disc, const scalar beta) -> dictionary
    {
        dictionary d;
        d.add("discretization", disc);
        d.add("beta", beta);
        return d;
    };

    // `explicit` MUST stay first: it is the reference of the
    // cross-discretization agreement check below. `exponential` is appended
    // last and is measured against its OWN contract -- it is the only member
    // that does not apply f_nl*psi at psi = psi^n (see exactFactor below).
    const List<word> discs
    {
        "explicit",
        "simpleLinearImplicit",
        "strictNegativeSpLinearImplicit",
        "exponential"
    };

    // ------------------------------------------------------------------ //
    // Cases: both signs of the normal strain, so that the sign-split
    // strictNegativeSpLinearImplicit is exercised on BOTH of its branches.
    //   alpha > 0 : the gradient is compressed, |grad psi| decays without a
    //               source; f_nl = R = alpha > 0 -> the Sc (explicit) branch.
    //   alpha < 0 : the gradient is stretched; f_nl < 0 -> the Sp branch.
    // ------------------------------------------------------------------ //
    const List<scalar> alphas{ 0.7, -0.7, 0.0 };

    const scalar beta = 1.3;

    for (const scalar alpha : alphas)
    {
        fill(alpha);

        const word tag = "alpha=" + name(alpha);
        Info<< nl << "==== " << tag << " (exact R = " << alpha << ") ====" << nl;

        // The exact affine values every model is measured against.
        const scalar fR = alpha;                 // R      : f_nl = n.(grad U).n
        const scalar fB = beta - cSlope;         // beta   : f_nl = beta - |grad psi|

        // Layer 1 + 2, per model, per discretization.
        struct ModelSpec { word label; bool isBeta; scalar f; };
        const List<ModelSpec> models
        {
            { "R",    false, fR },
            { "beta", true,  fB }
        };

        for (const ModelSpec& m : models)
        {
            // Reference: the applied source of the FIRST discretization; every
            // other must reproduce it. Filled on the first pass.
            scalarField appliedRef;

            for (const word& disc : discs)
            {
                const string id = tag + " " + m.label + "/" + disc;

                // One source object alive at a time: sdplsSource registers the
                // diagnostic fields "R" and "SDPLS_nonLinearPart" in the mesh
                // registry, so two coexisting instances would collide.
                autoPtr<sdplsSource> src;
                const dictionary d = srcDict(disc, beta);
                if (m.isBeta)
                {
                    src.reset(new sdplsBeta(d, mesh));
                }
                else
                {
                    src.reset(new sdplsR(d, mesh));
                }

                const fvScalarMatrix B(src->fvmsdplsSource(psi, U));

                // ---- Layer 1: the applied-source contract ---------------- //
                // For a matrix on the right of `==`, the physical source is
                // (diag*psi - source)/V. This is the invariant the sign defect
                // violated, and it must hold for every linearization.
                scalarField applied(mesh.nCells());
                forAll(applied, c)
                {
                    applied[c] =
                        (B.diag()[c]*psi[c] - B.source()[c])/V[c];
                }

                // `exponential` is NOT a linearization of the continuum
                // coefficient: it delivers the EXACT one-step amplification
                // factor exp(f dt), so the source it applies at psi = psi^n is
                // psi (e^{f dt} - 1)/dt, not f psi. The two differ at
                // O(f^2 dt): 2.5e-3 at f = 0.7, dt = 0.01, seven orders above
                // the 1e-10 tolerance, so it gets its own expectation here and
                // is exempt from the agreement check below. It is still gated
                // on the SIGN (the defect this test exists for) and on the
                // consistency layer, which it passes with zero deviation.
                const bool exactFactor = (disc == "exponential");
                const scalar appliedCoeff =
                    exactFactor
                  ? (Foam::exp(m.f*dt) - 1.0)/dt
                  : m.f;

                scalarField expected(mesh.nCells());
                forAll(expected, c)
                {
                    expected[c] = appliedCoeff*psi[c];
                }

                check
                (
                    id
                  + (exactFactor
                        ? " : applied source == psi (e^{f dt} - 1)/dt"
                        : " : applied source == f_nl*psi"),
                    maxDiff(applied, expected),
                    0.0,
                    1e-10
                );

                // ---- all discretizations agree ------------------------- //
                if (appliedRef.empty())
                {
                    appliedRef = applied;
                }
                else if (!exactFactor)
                {
                    check
                    (
                        id + " : agrees with " + discs[0],
                        maxDiff(applied, appliedRef),
                        0.0,
                        1e-10
                    );
                }

                // ---- Layer 2: one implicit-Euler step ------------------- //
                // Diagonal system (no div term), so psi^{n+1} = source/diag.
                fvScalarMatrix eqn
                (
                    fvm::ddt(psi) == src->fvmsdplsSource(psi, U)
                );

                // Which branch carries the term is the documented contract of
                // each linearization; f_nl is spatially constant here.
                const bool implicitBranch =
                    (disc == "simpleLinearImplicit")
                 || (disc == "strictNegativeSpLinearImplicit" && m.f < 0);

                scalarField psiNew(mesh.nCells());
                scalarField psiWant(mesh.nCells());
                forAll(psiNew, c)
                {
                    psiNew[c] = eqn.source()[c]/eqn.diag()[c];

                    // `exponential` is neither branch: Sp = 0 and
                    // Sc = psi^n (e^{f dt} - 1)/dt, so the Euler step
                    // psi^n + dt Sc collapses to the EXACT factor. This is the
                    // assertion that would fail if the class ever reverted to
                    // a truncated factor or read the outer ITERATE instead of
                    // psi.oldTime().
                    psiWant[c] =
                        exactFactor
                      ? psi[c]*Foam::exp(m.f*dt)
                      : (implicitBranch
                            ? psi[c]/(1.0 - m.f*dt)
                            : psi[c]*(1.0 + m.f*dt));
                }

                check
                (
                    id + " : one Euler step matches the closed form",
                    maxDiff(psiNew, psiWant),
                    0.0,
                    1e-10
                );

                // Both branches must be first-order consistent with the exact
                // pointwise solution psi e^{f dt} -- a reversed sign shows up
                // here as an O(dt) error, not an O(dt^2) one.
                scalarField psiExact(mesh.nCells());
                forAll(psiExact, c)
                {
                    psiExact[c] = psi[c]*Foam::exp(m.f*dt);
                }
                const scalar consistency = maxDiff(psiNew, psiExact);
                const scalar firstOrderBound =
                    max(mag(m.f)*dt, scalar(1e-30))*mag(m.f)*dt;
                if (consistency > 5*firstOrderBound + 1e-12)
                {
                    Info<< "FAIL " << id
                        << " : deviation from psi*exp(f dt) is " << consistency
                        << ", exceeds the O((f dt)^2) bound "
                        << 5*firstOrderBound << nl;
                    ++nFail;
                }
                else
                {
                    ++nPass;
                    Info<< "ok   " << id
                        << " : consistent with psi*exp(f dt) at O((f dt)^2)"
                        << nl;
                }
            }
        }

        // ------------------------------------------------------------------ //
        // Layer 3: the gradient law -- what SDPLS is actually for.
        //
        //   no source :  D|grad psi|/Dt = -a |grad psi|
        //   F = a     :  D|grad psi|/Dt = 0
        //
        // The source multiplies psi by (1 + a dt) (explicit branch), so the
        // affine slope becomes c(1 + a dt) while pure advection would compress
        // it to c(1 - a dt). A REVERSED source gives c(1 - a dt) -- i.e. it
        // moves the gradient the wrong way at first order, doubling the drift
        // instead of cancelling it. Assert the direction explicitly.
        // ------------------------------------------------------------------ //
        if (mag(alpha) > SMALL)
        {
            fill(alpha);

            autoPtr<sdplsSource> src;
            src.reset(new sdplsR(srcDict("explicit", beta), mesh));

            fvScalarMatrix eqn
            (
                fvm::ddt(psi) == src->fvmsdplsSource(psi, U)
            );

            // Slope of the updated affine field: psi^{n+1} = c(1 + a dt) x.
            scalar slopeNew = 0;
            forAll(mesh.C(), c)
            {
                if (mag(mesh.C()[c].x()) > SMALL)
                {
                    slopeNew =
                        (eqn.source()[c]/eqn.diag()[c])/mesh.C()[c].x();
                    break;
                }
            }

            check
            (
                tag + " R/explicit : affine slope -> c(1 + a dt)",
                slopeNew,
                cSlope*(1.0 + alpha*dt),
                1e-10
            );

            // The load-bearing direction test. Under pure advection the slope
            // would go to c(1 - a dt); the source must land on the OPPOSITE
            // side of c.
            const scalar advectedSlope = cSlope*(1.0 - alpha*dt);
            const bool correctSide =
                (slopeNew - cSlope)*(advectedSlope - cSlope) < 0;
            if (!correctSide)
            {
                Info<< "FAIL " << tag
                    << " R/explicit : source moves |grad psi| the SAME way as"
                    << " advection (slope " << slopeNew << ", advection-only "
                    << advectedSlope << ", initial " << cSlope
                    << ") -- this is the reversed-sign signature" << nl;
                ++nFail;
            }
            else
            {
                ++nPass;
                Info<< "ok   " << tag
                    << " R/explicit : source opposes the advective drift" << nl;
            }
        }
    }

    // ------------------------------------------------------------------ //
    // noSource must be exactly inert.
    // ------------------------------------------------------------------ //
    {
        fill(0.7);
        autoPtr<sdplsSource> src;
        src.reset(new sdplsSource(srcDict("explicit", beta), mesh));

        const fvScalarMatrix B(src->fvmsdplsSource(psi, U));
        scalar m = 0;
        forAll(psi, c)
        {
            m = max(m, mag(B.diag()[c]*psi[c] - B.source()[c]));
        }
        check("noSource : applies nothing", returnReduce(m, maxOp<scalar>()), 0.0);
    }

    Info<< nl << "==== " << nPass << " passed, " << nFail << " failed ===="
        << nl << endl;

    if (nFail)
    {
        Info<< "leiaTestSdplsSource FAILED" << nl << endl;
        return 1;
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
