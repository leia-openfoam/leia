/*---------------------------------------------------------------------------*\
Verification gate for signedDistanceEllipsoid (no mesh, pure arithmetic).

Exits non-zero on any failure, so it doubles as a repo test:
 1. DEGENERATE SPHERE: axes (R,R,R) must reproduce implicitSphere's signed
    distance to round-off and curvature 2/R at random points.
 2. SDF PROPERTIES on a triaxial ellipsoid: the closest point lies ON the
    surface (algebraic residual ~ 0), |grad| = 1, and Eikonal consistency
    |value(x + eps g) - value(x) - eps| < tol.
 3. MID-PLANE CROSS-CHECK: in the z = c_z plane the closest point of the 3D
    triaxial ellipsoid with axes (a, b, c) does NOT generally coincide with
    the 2D ellipse's (the 3D surface curves away in z), BUT for the special
    prolate case a > b = c the x-y plane section IS the distance-realising
    set, so signedDistanceEllipse(a, b) must agree there.
 4. CURVATURE CONVENTION: on the sphere limit curvature = 2/R (div convention);
    at the ends of the axes of a triaxial ellipsoid the exact total curvature
    kappa1 + kappa2 has the closed form  a*(1/b^2 + 1/c^2)  (at the a-end).
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "levelSetImplicitSurfaces.H"
#include <random>

using namespace Foam;

int main(int argc, char *argv[])
{
    label nFail = 0;
    auto expect = [&](const bool ok, const word& what, const scalar got,
                      const scalar want)
    {
        if (!ok)
        {
            ++nFail;
            Info<< "FAIL " << what << ": got " << got
                << ", want " << want << nl;
        }
    };

    std::mt19937 rng(42);
    std::uniform_real_distribution<scalar> u(-1.0, 1.0);

    // ---- 1. degenerate sphere ------------------------------------------
    {
        const vector c(0.35, 0.35, 0.35);
        const scalar R = 0.15;
        signedDistanceEllipsoid e(c, vector(R, R, R));
        implicitSphere s(c, R);
        scalar worstV = 0, worstK = 0;
        for (label i = 0; i < 2000; ++i)
        {
            const vector x = c + 0.4*vector(u(rng), u(rng), u(rng));
            worstV = max(worstV, mag(e.value(x) - s.value(x)));
            worstK = max(worstK, mag(e.curvature(x) - 2.0/R));
        }
        expect(worstV < 1e-10, "sphere-degenerate |value diff|", worstV, 0);
        expect(worstK < 1e-8*(2.0/R), "sphere-degenerate curvature", worstK, 0);
    }

    // ---- 2. SDF properties, triaxial ------------------------------------
    {
        const vector c(0.35, 0.35, 0.35);
        const vector ax(0.21, 0.15, 0.107);
        signedDistanceEllipsoid e(c, ax);
        implicitEllipsoid f(c, ax);
        scalar worstAlg = 0, worstG = 0, worstEik = 0;
        const scalar eps = 1e-6;
        for (label i = 0; i < 2000; ++i)
        {
            const vector x = c + 0.35*vector(u(rng), u(rng), u(rng));
            // closest point on surface: reconstruct via x - value*grad
            const scalar d = e.value(x);
            const vector g = e.grad(x);
            const vector q = x - d*g;
            worstAlg = max(worstAlg, mag(f.value(q)));
            worstG = max(worstG, mag(mag(g) - 1.0));
            worstEik = max(worstEik,
                mag(e.value(x + eps*g) - d - eps)/eps);
        }
        expect(worstAlg < 1e-9, "closest point on surface (algebraic)",
               worstAlg, 0);
        expect(worstG < 1e-12, "|grad| = 1", worstG, 0);
        expect(worstEik < 1e-5, "Eikonal step consistency", worstEik, 0);
    }

    // ---- 3. prolate mid-plane vs 2D ellipse ------------------------------
    {
        const vector c(0.0, 0.0, 0.0);
        const scalar a = 0.21, b = 0.13;
        signedDistanceEllipsoid e3(c, vector(a, b, b));   // prolate: b = c
        signedDistanceEllipse   e2(c, vector(a, b, 1.0));
        scalar worst = 0;
        for (label i = 0; i < 500; ++i)
        {
            const vector x(0.3*u(rng), 0.3*u(rng), 0.0);
            worst = max(worst, mag(e3.value(x) - e2.value(x)));
        }
        expect(worst < 1e-9, "prolate mid-plane vs signedDistanceEllipse",
               worst, 0);
    }

    // ---- 4. curvature closed forms at the axis ends ----------------------
    {
        const vector c(0.35, 0.35, 0.35);
        const scalar a = 0.21, b = 0.15, cc = 0.107;
        signedDistanceEllipsoid e(c, vector(a, b, cc));
        // evaluate slightly OUTSIDE each axis end; closest point = the end
        struct End { vector dir; scalar want; };
        const End ends[3] =
        {
            {vector(1,0,0), a*(1.0/sqr(b) + 1.0/sqr(cc))},
            {vector(0,1,0), b*(1.0/sqr(a) + 1.0/sqr(cc))},
            {vector(0,0,1), cc*(1.0/sqr(a) + 1.0/sqr(b))},
        };
        const scalar ax_[3] = {a, b, cc};
        for (label j = 0; j < 3; ++j)
        {
            const vector x = c + (ax_[j] + 0.01)*ends[j].dir;
            const scalar got = e.curvature(x);
            expect(mag(got - ends[j].want) < 1e-6*ends[j].want,
                   "axis-end total curvature", got, ends[j].want);
        }
    }

    if (nFail)
    {
        Info<< nFail << " FAILURES" << nl;
        return 1;
    }
    Info<< "signedDistanceEllipsoid: ALL CHECKS PASS" << nl;
    return 0;
}
