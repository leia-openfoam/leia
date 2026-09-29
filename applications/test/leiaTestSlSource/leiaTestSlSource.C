/*---------------------------------------------------------------------------*\
Unit gate for slGradientControlSource, the semi-Lagrangian gradient-control
source (docs/plan-halo-limited-gradient-control.md, D3). Run it in
cases/slSourceUnit (serial and on 4 ranks: bash cases/slSourceUnit/Allrun.sh).
Exits non-zero on any failure.

Fields: the planar level set psi = g0 x and the affine velocity
u = (alpha x, -alpha y, 0), so that q = |grad psi| = g0, n = e_x,
a = n.D.n = alpha and sigma = |D|_F = sqrt(2)|alpha| are exact for the
leastSquares gradients. One apply() with dt must give:

 1. psi = psi^0 exp(dt F) in the band |psi|/q <= 3h (6 x 40 = 240 cells),
    F from the law at the exact values, to round-off;
 2. psi = psi^0 bit for bit outside the band;
 3. no sign change anywhere;
 4. slSourceF = F in the band and 0 outside; slSourceClamp = 1 exactly where
    |dt F| > 30 (and psi = psi^0 exp(-+30) there);
 5. the band cell size h = 0.05 in every cell, processor neighbours included;
 6. the processor-patch values of psi are current after apply().

On 4 ranks the band crosses the processor boundary x = 0, so the gradient,
the band and the counts must still be exact.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "calculatedFvPatchFields.H"
#include "slSource.H"
#include "slGradientControlSource.H"
#include "slReconstruction.H"

using namespace Foam;

static label nFail = 0;
static label nPass = 0;

static void check(const string& name, const scalar got, const scalar want, const scalar tol)
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


static dictionary lawDict(const word& type, const word& weight, const scalar mu)
{
    dictionary l;
    l.add("type", type);
    l.add("mu", mu);
    l.add("eps", 0.02);
    l.add("c", 1.0);
    l.add("cKappa", 1.25);
    l.add("deltaS", 0.08);
    l.add("p", label(5));
    l.add("gamma", 1.4722194895832204);
    l.add("epsD", 0.0);
    dictionary w;
    w.add("type", weight);
    w.add("beta", 1.0);
    w.add("m", label(2));
    w.add("deltaS", 0.08);
    l.add("strainWeight", w);
    return l;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar dt = runTime.deltaTValue();
    const scalar hExact = 0.05;
    const scalar bandHalfWidth = 3*hExact;

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(dimLength, 0),
        calculatedFvPatchScalarField::typeName
    );
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        calculatedFvPatchVectorField::typeName
    );

    auto fill = [&](const scalar g0, const scalar alpha)
    {
        forAll(mesh.C(), c)
        {
            const point& x = mesh.C()[c];
            psi[c] = g0*x.x();
            U[c] = vector(alpha*x.x(), -alpha*x.y(), 0);
        }
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled()) continue;
            fvPatchScalarField& pp = psi.boundaryFieldRef()[patchi];
            fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];
            const vectorField& Cf = mesh.boundary()[patchi].Cf();
            forAll(pp, i)
            {
                pp[i] = g0*Cf[i].x();
                pU[i] = vector(alpha*Cf[i].x(), -alpha*Cf[i].y(), 0);
            }
        }
        // Processor patches: the halo values the leastSquares gradient reads.
        psi.correctBoundaryConditions();
        U.correctBoundaryConditions();
    };

    // The test replaces the law of levelSet.semiLagrangian.source in memory.
    fvSolution& fvs = const_cast<fvSolution&>(static_cast<const fvSolution&>(mesh));
    dictionary& srcDict =
        fvs.subDict("levelSet").subDict("semiLagrangian").subDict("source");

    autoPtr<slReconstruction> fit(slReconstruction::New(mesh));

    // ---- the inert model ------------------------------------------------------
    {
        srcDict.set("type", word("none"));
        autoPtr<slSource> src(slSource::New(mesh));
        check("type none is inactive", src->active(), 0, 0);
        srcDict.set("type", word("gradientControl"));
    }

    struct Case { word law; word weight; scalar mu; scalar g0; scalar alpha; };
    const List<Case> cases
    {
        {"linearZ",  "none",  2,    1.4,  0.0},   // q > 1: F < 0
        {"linearZ",  "none",  2,    0.7,  0.0},   // q < 1: F > 0
        {"none",     "full",  2,    1.0,  0.8},   // F = a = alpha
        {"softWall", "none",  2,    1.12, 0.5},   // F = -1.25 sqrt(2) 0.5 tanh(...)
        {"linearZ",  "omega", 2,    1.2,  0.3},   // F = w(q) a + G(q)
        {"linearQ",  "none",  1000, 2.0,  0.0}    // dt F = -50: clamped at -30
    };

    for (const Case& k : cases)
    {
        const string tag =
            k.law + "/" + k.weight + " g0=" + name(k.g0) + " alpha=" + name(k.alpha);
        Info<< nl << "==== " << tag << " ====" << nl;
        srcDict.set("law", lawDict(k.law, k.weight, k.mu));
        fill(k.g0, k.alpha);
        const scalarField psi0(psi.primitiveField());

        autoPtr<slSource> src(slSource::New(mesh));
        const auto& gc = refCast<const slGradientControlSource>(src());

        const scalar F = gc.law().F(sqr(k.g0), Foam::sqrt(2.0)*mag(k.alpha), k.alpha);
        const scalar x = dt*F;
        const bool clamped = mag(x) > slGradientControlSource::maxExponent;
        const scalar factor =
            Foam::exp(clamped ? sign(x)*slGradientControlSource::maxExponent : x);

        src->apply(psi, U, fit(), dt);

        const volScalarField& Ff = mesh.lookupObject<volScalarField>("slSourceF");
        const volScalarField& Cl = mesh.lookupObject<volScalarField>("slSourceClamp");

        label nBand = 0, nFlip = 0, nClamp = 0;
        scalar worstIn = 0, worstOut = 0, worstF = 0, worstH = 0;
        forAll(psi0, c)
        {
            const bool inBand = mag(mesh.C()[c].x()) <= bandHalfWidth;
            if (inBand)
            {
                ++nBand;
                const scalar want = psi0[c]*factor;
                worstIn = max(worstIn, mag(psi[c] - want)/max(mag(want), VSMALL));
                worstF = max(worstF, mag(Ff[c] - F)/max(mag(F), scalar(1)));
            }
            else
            {
                worstOut = max(worstOut, mag(psi[c] - psi0[c]));
                worstF = max(worstF, mag(Ff[c]));
            }
            if (psi0[c]*psi[c] < 0 || (psi0[c] == 0 && psi[c] != 0))
            {
                ++nFlip;
            }
            if (Cl[c] > 0.5)
            {
                ++nClamp;
            }
            worstH = max(worstH, mag(gc.cellSize()[c] - hExact));
        }
        reduce(nBand, sumOp<label>());
        reduce(nFlip, sumOp<label>());
        reduce(nClamp, sumOp<label>());
        reduce(worstIn, maxOp<scalar>());
        reduce(worstOut, maxOp<scalar>());
        reduce(worstF, maxOp<scalar>());
        reduce(worstH, maxOp<scalar>());

        check(tag + " : band cells", nBand, 240, 0);
        check(tag + " : psi = psi0 exp(dt F) in the band", worstIn, 0, 1e-12);
        check(tag + " : psi unchanged outside the band (bit for bit)", worstOut, 0, 0);
        check(tag + " : no sign change", nFlip, 0, 0);
        check(tag + " : slSourceF = F in the band, 0 outside", worstF, 0, 1e-12);
        check(tag + " : clamped cells", nClamp, clamped ? 240 : 0, 0);
        check(tag + " : cell size h = 0.05 everywhere", worstH, 0, 1e-14);

        // The processor-patch values of psi after apply() must be the current
        // neighbour values: a fresh exchange on a copy may not change them.
        volScalarField psiFresh("psiFresh", psi);
        psiFresh.correctBoundaryConditions();
        scalar worstHalo = 0;
        forAll(psi.boundaryField(), patchi)
        {
            if (!psi.boundaryField()[patchi].coupled()) continue;
            worstHalo = max
            (
                worstHalo,
                max(mag(psi.boundaryField()[patchi] - psiFresh.boundaryField()[patchi]))
            );
        }
        reduce(worstHalo, maxOp<scalar>());
        check(tag + " : processor-patch values are current", worstHalo, 0, 0);
    }

    Info<< nl << "==== " << nPass << " passed, " << nFail << " failed ====" << nl << endl;
    if (nFail)
    {
        Info<< "leiaTestSlSource FAILED" << nl << endl;
        return 1;
    }
    Info<< "End\n" << endl;
    return 0;
}
