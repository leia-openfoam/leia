/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestParallelSurfaceInverse

Description
    Mesh-free analytic unit test of the dimension-seamless parallel-surface
    curvature inverse (stabilizedFootPointFaceCurvature.H):

        kappa^G = (kappa_d - 2 K_d d) / (1 - kappa_d d + K_d d^2)

    Given EXACT offset-surface inputs (kappa_d, K_d) at signed distance d
    (sign of psi: positive outside), the inverse must return the interface
    mean-curvature sum exactly (to round-off) on:
      * sphere  R:            kappa_d = 2/(R+d),   K_d = 1/(R+d)^2  -> 2/R
      * cylinder R (K = 0):   kappa_d = 1/(R+d),   K_d = 0          -> 1/R
        (also the K = 0 reduction identity kappa_d/(1 - d kappa_d))
      * torus (R0, r) at poloidal angle theta -- including the SADDLE side
        (theta = pi, K < 0): offset torus has kappa1_d = 1/(r+d),
        kappa2_d = cos(theta)/(R0 + (r+d) cos(theta)); the inverse must
        return 1/r + cos(theta)/(R0 + r cos(theta)).
    Plus: the guard branch (|1 - kappa_d d + K_d d^2| <= 1/2 returns the
    input unchanged) and the implicit-surface Gaussian-curvature expression
    K = (g.(cof(H).g))/|g|^4 on the analytic sphere Hessian.

    Exits non-zero on any mismatch (usable as a repo regression test).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "stabilizedFootPointFaceCurvature.H"

using namespace Foam;

static label nFail = 0;

void check(const word& name, const scalar got, const scalar want,
           const scalar tol = 1e-12)
{
    const scalar scale = max(mag(want), scalar(1));
    if (mag(got - want) > tol*scale)
    {
        Info<< "FAIL " << name << ": got " << got << ", want " << want
            << " (|diff| = " << mag(got - want) << ")" << nl;
        ++nFail;
    }
    else
    {
        Info<< "ok   " << name << " (" << got << ")" << nl;
    }
}


int main(int argc, char *argv[])
{
    // --- sphere: kappa = 2/R, both offset signs -----------------------------
    {
        const scalar R = 0.15;
        for (const scalar d : {0.3*R, -0.3*R, 0.05*R, -0.05*R})
        {
            const scalar s = R + d;
            const scalar kd = 2.0/s;
            const scalar Kd = 1.0/(s*s);
            check("sphere R=0.15 d=" + name(d),
                  parallelSurfaceInverse(kd, d, Kd), 2.0/R);
        }
    }

    // --- cylinder: K = 0 exactly, kappa = 1/R -------------------------------
    {
        const scalar R = 1e-3;
        for (const scalar d : {0.2*R, -0.2*R})
        {
            const scalar kd = 1.0/(R + d);
            check("cylinder R=1e-3 d=" + name(d),
                  parallelSurfaceInverse(kd, d, scalar(0)), 1.0/R);

            // K = 0 reduction identity: must equal the 2D scalar inverse
            // bit-for-bit.
            const scalar twoD = kd/(1.0 - d*kd);
            if (parallelSurfaceInverse(kd, d, scalar(0)) != twoD)
            {
                Info<< "FAIL cylinder K=0 reduction is not bit-identical" << nl;
                ++nFail;
            }
        }
    }

    // --- torus (R0, r): general kappa1 != kappa2, saddle K < 0 --------------
    {
        const scalar R0 = 1.0, r = 0.3;
        for (const scalar theta : {0.0, 2.0*constant::mathematical::pi/3.0,
                                   constant::mathematical::pi})
        {
            const scalar ct = Foam::cos(theta);
            const scalar kExact = 1.0/r + ct/(R0 + r*ct);
            for (const scalar d : {0.1*r, -0.1*r})
            {
                const scalar k1d = 1.0/(r + d);
                const scalar k2d = ct/(R0 + (r + d)*ct);
                const scalar kd = k1d + k2d;
                const scalar Kd = k1d*k2d;   // K < 0 for theta > pi/2
                check("torus theta=" + name(theta) + " d=" + name(d)
                      + " (K_d=" + name(Kd) + ")",
                      parallelSurfaceInverse(kd, d, Kd), kExact);
            }
        }
    }

    // --- guard branch: past half the local curvature radius -> input kept ---
    {
        const scalar R = 0.15;
        const scalar d = 0.6*R;                 // sphere denom (R/(R+d))^2 < 1/2
        const scalar s = R + d;
        const scalar kd = 2.0/s, Kd = 1.0/(s*s);
        const scalar denom = 1.0 - kd*d + Kd*d*d;
        if (mag(denom) > 0.5)
        {
            Info<< "FAIL guard test setup: denom = " << denom
                << " expected <= 1/2" << nl;
            ++nFail;
        }
        else if (parallelSurfaceInverse(kd, d, Kd) != kd)
        {
            Info<< "FAIL guard branch did not return the input unchanged" << nl;
            ++nFail;
        }
        else
        {
            Info<< "ok   guard branch (denom = " << denom << ")" << nl;
        }
    }

    // --- implicit-surface Gaussian curvature on the analytic sphere ---------
    {
        // psi = |x| - R at radius rr, point on the z-axis: g = e_z,
        // H = (I - zz)/rr -> K = 1/rr^2.
        const scalar rr = 0.2;
        const vector g(0, 0, 1);
        const symmTensor H(1.0/rr, 0, 0,
                                   1.0/rr, 0,
                                           0);
        const scalar K = (g & (cof(H) & g))/pow4(mag(g));
        check("Gaussian K on sphere Hessian", K, 1.0/(rr*rr));

        // pseudo-2D fit (zero z-components everywhere): K must be EXACTLY +0.
        const vector g2(0.6, 0.8, 0);
        const symmTensor H2(3.0, 1.0, 0,
                                 2.0, 0,
                                      0);
        const scalar K2 = (g2 & (cof(H2) & g2))/pow4(mag(g2));
        if (K2 != 0.0)
        {
            Info<< "FAIL pseudo-2D Gaussian curvature is not exactly zero: "
                << K2 << nl;
            ++nFail;
        }
        else
        {
            Info<< "ok   pseudo-2D K == +0 exactly" << nl;
        }
    }

    if (nFail)
    {
        Info<< nl << nFail << " CHECK(S) FAILED" << nl << endl;
        return 1;
    }
    Info<< nl << "All parallel-surface inverse checks passed." << nl << endl;
    return 0;
}

// ************************************************************************* //
