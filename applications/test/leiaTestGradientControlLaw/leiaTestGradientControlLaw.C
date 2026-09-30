/*---------------------------------------------------------------------------*\
Unit gate for libleiaGradientControl (no mesh, pure arithmetic;
docs/plan-halo-limited-gradient-control.md, D1). Exits non-zero on any failure.

 1. Every law: G(q = 1) = 0 exactly; G has the sign opposite to q - 1;
    dG/dq equals a central difference of G.
 2. The q- and z-pairs have the same slope at q = 1 for the same mu.
 3. Soft wall: tanh(gamma) = 0.9 at |q - 1| = deltaS with gamma = artanh(0.9);
    the response table of the dossier (3 to 10 percent); the rate scales with
    sqrt(sigma^2 + epsD^2); an even exponent is refused.
 4. Strain weights: none = 0, full = 1, omega(1) = 0 and
    omega(1 +- deltaS) = 1 - exp(-beta); with the full weight the strain
    cancels exactly on the interface law q' = q (F - a).
 5. RK4 integration of q' = q G (strain cancelled) against the closed forms:
    the logistic law (linearQ), the z-logistic law z = q^2 (linearZ), the
    implicit cubic solution (cubicQ); the other laws relax q monotonically
    towards 1 without crossing it.
 6. Selection: every law name selects its class; an unknown name and a
    missing coefficient are refused.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "gradientControlLaw.H"
#include "IOstreams.H"
#include "IStringStream.H"
#include <cmath>
#include <functional>

using namespace Foam;

namespace
{

label nFail = 0;
label nPass = 0;

void expect(const bool ok, const std::string& what, const scalar got, const scalar want)
{
    if (ok)
    {
        ++nPass;
    }
    else
    {
        ++nFail;
        Info<< "FAIL " << what.c_str() << ": got " << got << ", want " << want << nl;
    }
}

dictionary lawDict(const std::string& text)
{
    IStringStream is(text);
    return dictionary(is);
}

// The reference coefficients (the default.parameter values; mu = 2 to test
// that mu enters). Every entry is present, as in the case templates.
std::string coeffs(const word& type, const word& weight = "none")
{
    return "type " + type + "; mu 2; eps 0.02; c 1; cKappa 1.25; deltaS 0.08;"
        " p 5; gamma 1.4722194895832204; epsD 0;"   // artanh(0.9) to 17 digits
        " strainWeight { type " + weight + "; beta 1; m 2; deltaS 0.08; }";
}

// One RK4 step of q' = q (F(q, sigma, a) - a).
scalar rk4(const gradientControlLaw& law, scalar q, const scalar dt,
           const scalar sigma, const scalar a)
{
    auto f = [&](const scalar x) { return x*(law.F(x*x, sigma, a) - a); };
    const scalar k1 = f(q);
    const scalar k2 = f(q + 0.5*dt*k1);
    const scalar k3 = f(q + 0.5*dt*k2);
    const scalar k4 = f(q + dt*k3);
    return q + dt/6*(k1 + 2*k2 + 2*k3 + k4);
}

scalar integrate(const gradientControlLaw& law, scalar q, const scalar T,
                 const label n, const scalar sigma = 0, const scalar a = 0)
{
    const scalar dt = T/n;
    for (label i = 0; i < n; ++i)
    {
        q = rk4(law, q, dt, sigma, a);
    }
    return q;
}

// The implicit solution of e' = -mu (1 + e) e^3 (cubicQ, q = 1 + e):
// H(e) = -1/(2 e^2) + 1/e + ln|e/(1 + e)|, H(e(T)) - H(e0) = -mu T.
scalar Hcubic(const scalar e)
{
    return -0.5/(e*e) + 1/e + std::log(std::fabs(e/(1 + e)));
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noBanner();
    argList args(argc, argv, false, false, false);

    FatalError.throwing(true);
    FatalIOError.throwing(true);

    const wordList laws
    ({
        "linearQ", "linearZ", "cubicQ", "cubicZ", "twoThirdsZReg",
        "saturatedLinearZ", "softWall"
    });

    // ---- 1. G(1) = 0, the sign, and dG/dq ------------------------------------
    for (const word& name : laws)
    {
        auto law = gradientControlLaw::New(lawDict(coeffs(name)));
        expect(law->type() == name, "selected type " + name, 0, 0);
        expect(law->active() && law->hasRate(), name + " active", 0, 1);
        for (const scalar sigma : {0.0, 0.7, 3.0})
        {
            const scalar g1 = law->rate(1, sigma);
            expect(g1 == 0, name + " G(q = 1) == 0 exactly", g1, 0);
        }
        scalar worstSign = 1, worstDer = 0;
        for (label i = 1; i < 60; ++i)
        {
            const scalar q = 0.05*i;              // 0.05 ... 2.95
            if (std::fabs(q - 1) < 1e-9) continue;
            const scalar G = law->rate(q*q, 1.0);
            worstSign = std::min(worstSign, -G*(q - 1));
            const scalar h = 1e-6;
            const scalar fd =
                (law->rate(sqr(q + h), 1.0) - law->rate(sqr(q - h), 1.0))/(2*h);
            const scalar d = law->dRateDq(q*q, 1.0);
            worstDer = std::max(worstDer, std::fabs(d - fd)/std::max(std::fabs(fd), 1e-3));
        }
        // softWall: tanh(gamma s^5) underflows to 0 only below |q - 1| ~ 1e-62;
        // every grid point is far from that, so the sign test is strict.
        expect(worstSign > 0, name + " G has the sign of 1 - q (strict)", worstSign, 0);
        expect(worstDer < 1e-6, name + " dG/dq against a central difference", worstDer, 0);
    }

    // ---- 2. the q/z pairs: equal slopes at q = 1 --------------------------------
    {
        auto lq = gradientControlLaw::New(lawDict(coeffs("linearQ")));
        auto lz = gradientControlLaw::New(lawDict(coeffs("linearZ")));
        auto cq = gradientControlLaw::New(lawDict(coeffs("cubicQ")));
        auto cz = gradientControlLaw::New(lawDict(coeffs("cubicZ")));
        expect(lq->dRateDq(1, 0) == -2, "linearQ slope -mu at q = 1", lq->dRateDq(1, 0), -2);
        expect(lz->dRateDq(1, 0) == -2, "linearZ slope -mu at q = 1", lz->dRateDq(1, 0), -2);
        expect(cq->dRateDq(1, 0) == 0, "cubicQ slope 0 at q = 1", cq->dRateDq(1, 0), 0);
        expect(cz->dRateDq(1, 0) == 0, "cubicZ slope 0 at q = 1", cz->dRateDq(1, 0), 0);
        // second order: G_z - G_q = -mu e^2/2 for the linear pair
        const scalar e = 1e-3, q = 1 + e;
        const scalar dG = lz->rate(q*q, 0) - lq->rate(q*q, 0);
        expect(std::fabs(dG + 2*0.5*e*e) < 1e-15, "linear pair differs by -mu e^2/2", dG, -e*e);
    }

    // ---- 3. the soft wall ------------------------------------------------------
    {
        auto sw = gradientControlLaw::New(lawDict(coeffs("softWall")));
        expect(sw->needsStrainRate(), "softWall needs the strain rate", 0, 1);
        const scalar cKappa = 1.25;
        // Psi(s) = -G(1 + s, sigma = 1)/cKappa with epsD = 0
        auto Psi = [&](const scalar s) { return -sw->rate(sqr(1 + s), 1.0)/cKappa; };
        expect(std::fabs(Psi(0.08) - 0.9) < 1e-12, "tanh(gamma) = 0.9 at |q - 1| = deltaS", Psi(0.08), 0.9);
        // robust design point: deltaA - eps_q = 0.10 - 0.02 = deltaS -> Theta* = 0.9
        expect(std::fabs(Psi(0.10 - 0.02) - 0.9) < 1e-12, "Theta* = 0.9 at deltaA - eps_q", Psi(0.08), 0.9);
        const scalar s[7] = {0.03, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10};
        const scalar table[7] = {0.011, 0.139, 0.336, 0.638, 0.900, 0.990, 0.9997};
        for (label i = 0; i < 7; ++i)
        {
            expect(std::fabs(Psi(s[i]) - table[i]) <= 5e-4, "dossier table, q = 1 + " + std::to_string(s[i]), Psi(s[i]), table[i]);
            expect(std::fabs(-Psi(-s[i]) - table[i]) <= 5e-4, "dossier table, q = 1 - " + std::to_string(s[i]), -Psi(-s[i]), table[i]);
        }
        // the rate scales with sqrt(sigma^2 + epsD^2)
        std::string c = coeffs("softWall");
        c.replace(c.find("epsD 0;"), 7, "epsD 0.3;");
        auto swe = gradientControlLaw::New(lawDict(c));
        const scalar r = swe->rate(sqr(1.05), 0.4)/swe->rate(sqr(1.05), 0.0);
        expect(std::fabs(r - 0.5/0.3) < 1e-13, "softWall rate ~ sqrt(sigma^2 + epsD^2)", r, 0.5/0.3);
        // an even exponent is refused
        for (const std::string& bad : {std::string("p 4;"), std::string("p 0;")})
        {
            std::string cb = coeffs("softWall");
            cb.replace(cb.find("p 5;"), 4, bad);
            bool refused = false;
            try { auto x = gradientControlLaw::New(lawDict(cb)); }
            catch (const Foam::error&) { refused = true; }
            expect(refused, "softWall refuses " + bad, refused, 1);
        }
    }

    // ---- 4. strain weights -------------------------------------------------------
    {
        auto none = gradientControlLaw::New(lawDict(coeffs("none", "none")));
        auto full = gradientControlLaw::New(lawDict(coeffs("none", "full")));
        auto omeg = gradientControlLaw::New(lawDict(coeffs("none", "omega")));
        expect(!none->active(), "law none + weight none is inactive", none->active(), 0);
        expect(full->active() && full->needsNormalStrain(), "weight full is active", full->active(), 1);
        scalar worst0 = 0, worst1 = 0;
        for (label i = 1; i < 60; ++i)
        {
            const scalar q = 0.05*i, a = 0.37;
            worst0 = std::max(worst0, std::fabs(none->F(q*q, 1, a)));
            worst1 = std::max(worst1, std::fabs(full->F(q*q, 1, a) - a));
        }
        expect(worst0 == 0, "weight none: F = 0", worst0, 0);
        expect(worst1 == 0, "weight full + law none: F = a exactly", worst1, 0);
        const scalar w1 = omeg->weight().w(1);
        expect(w1 == 0, "omega(q = 1) = 0 exactly", w1, 0);
        for (const scalar q : {1.08, 0.92})
        {
            const scalar w = omeg->weight().w(q*q);
            expect(std::fabs(w - (1 - std::exp(-1.0))) < 1e-12, "omega(1 +- deltaS) = 1 - exp(-beta)", w, 1 - std::exp(-1.0));
        }
        expect(omeg->weight().w(sqr(1.24)) > 1 - 1e-15, "omega -> 1 at |q - 1| = 3 deltaS", omeg->weight().w(sqr(1.24)), 1);
        // the SDPLS cancellation: q' = q (F - a) with w = 1 does not see a
        auto lzFull = gradientControlLaw::New(lawDict(coeffs("linearZ", "full")));
        auto lzNone = gradientControlLaw::New(lawDict(coeffs("linearZ", "none")));
        const scalar qa = integrate(*lzFull, 1.4, 2.0, 400, 1.0, 0.8);
        const scalar qb = integrate(*lzNone, 1.4, 2.0, 400, 1.0, 0.0);
        expect(std::fabs(qa - qb) < 1e-14, "weight full cancels the strain on q' = q (F - a)", qa, qb);
    }

    // ---- 5. RK4 against the closed forms ---------------------------------------
    {
        const scalar mu = 2, T = 3.0/mu;
        auto lq = gradientControlLaw::New(lawDict(coeffs("linearQ")));
        auto lz = gradientControlLaw::New(lawDict(coeffs("linearZ")));
        auto cq = gradientControlLaw::New(lawDict(coeffs("cubicQ")));
        for (const scalar q0 : {0.5, 1.6})
        {
            const scalar exact = q0/(q0 + (1 - q0)*std::exp(-mu*T));
            const scalar e1 = std::fabs(integrate(*lq, q0, T, 150) - exact);
            const scalar e2 = std::fabs(integrate(*lq, q0, T, 300) - exact);
            expect(e2 < 1e-9, "linearQ: logistic q(T), q0 = " + std::to_string(q0), e2, 0);
            expect(std::fabs(std::log2(e1/e2) - 4) < 0.3, "linearQ: RK4 order 4", std::log2(e1/e2), 4);
            const scalar z0 = q0*q0;
            const scalar zex = z0/(z0 + (1 - z0)*std::exp(-mu*T));
            const scalar qz = integrate(*lz, q0, T, 300);
            expect(std::fabs(qz*qz - zex) < 1e-9, "linearZ: z-logistic z(T), q0 = " + std::to_string(q0), qz*qz, zex);
        }
        for (const scalar q0 : {0.8, 1.25})
        {
            const scalar e0 = q0 - 1;
            const scalar qT = integrate(*cq, q0, T, 3000);
            const scalar res = Hcubic(qT - 1) - Hcubic(e0) + mu*T;
            expect(std::fabs(res) < 1e-9*mu*T, "cubicQ: implicit solution, q0 = " + std::to_string(q0), res, 0);
        }
        for (const word& name : laws)
        {
            auto law = gradientControlLaw::New(lawDict(coeffs(name)));
            for (const scalar q0 : {0.7, 1.3})
            {
                scalar q = q0, worst = 0;
                bool monotone = true;
                for (label i = 0; i < 2000; ++i)
                {
                    const scalar qn = rk4(*law, q, 1e-3, 1.0, 0.0);
                    if (std::fabs(qn - 1) > std::fabs(q - 1) || (qn - 1)*(q0 - 1) < 0)
                    {
                        monotone = false;
                        worst = qn;
                    }
                    q = qn;
                }
                expect(monotone, name + " relaxes q monotonically without crossing 1, q0 = " + std::to_string(q0), worst, 1);
            }
        }
    }

    // ---- 6. selection refusals ---------------------------------------------------
    {
        bool refused = false;
        try { auto x = gradientControlLaw::New(lawDict("type noSuchLaw;")); }
        catch (const Foam::error&) { refused = true; }
        expect(refused, "an unknown law name is refused", refused, 1);
        refused = false;
        try { auto x = gradientControlLaw::New(lawDict("type linearZ;")); }
        catch (const Foam::error&) { refused = true; }
        expect(refused, "a missing mu is refused", refused, 1);
        refused = false;
        try { auto x = gradientControlLaw::New(lawDict("type linearZ; mu 0;")); }
        catch (const Foam::error&) { refused = true; }
        expect(refused, "mu = 0 is refused", refused, 1);
        refused = false;
        try { auto x = gradientControlLaw::New(lawDict("type none; strainWeight { type omega; beta 1; m 2; }")); }
        catch (const Foam::error&) { refused = true; }
        expect(refused, "omega without deltaS is refused", refused, 1);
    }

    Info<< nl << "==== " << nPass << " passed, " << nFail << " failed ====" << endl;
    return nFail ? 1 : 0;
}
