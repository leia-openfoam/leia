/*---------------------------------------------------------------------------*\
Unit gate for the halo-limited velocity extension
(docs/plan-halo-limited-gradient-control.md, D4). Run it in cases/haloLimitedUnit
(serial and on 4 ranks: bash cases/haloLimitedUnit/Allrun.sh). Exits non-zero
on any failure.

 1. Travel and weight: S(0) = 0, S'(0) = 1, S odd, |S| <= R up to |d| = 1e6 R
    (|S| < R in exact arithmetic; 1 + 1e-24 rounds to 1);
    c(0) = 1; 1 - w = (beta/(2m)) (d/R)^(2m) to leading order.
 2. Direction: d = psi/q and e = grad(psi)/q for q >= 0.5; Q2 > 0 for q -> 0.
 3. Planar signed distance psi = x - x0 with an affine velocity: the local
    quadratic models are exact, so the face correction is
    -w_f S_R(d_f) (dU/dx).S_f and Uext = U - w S_R(d) dU/dx in every cell, to
    round-off; physical patches keep phi; |Y - x|/R < 1 everywhere.
 4. Uniform velocity: phiExt == phi and Uext == U bit for bit.
 5. The divergence-free quadratic velocity (alpha x + beta x^2,
    -(alpha + 2 beta x) y) of the dossier with psi = x: the correction equals
    w_f [u(Y_f) - u(x_f)].S_f with the exact u at the exact sample point, on
    the faces whose cells carry the full quadratic model.
 6. Coupled faces carry exactly opposite corrections on the two sides.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "calculatedFvPatchFields.H"
#include "haloLimited.H"
#include "syncTools.H"
#include <functional>

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


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // ---- fields the extension looks up (psi, U, phi) ---------------------------
    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimLength, 0), calculatedFvPatchScalarField::typeName
    );
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimVelocity, Zero), calculatedFvPatchVectorField::typeName
    );
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimVelocity*dimArea, 0)
    );

    const scalar x0 = 0.01;
    auto setFields = [&](const std::function<vector(const point&)>& u)
    {
        forAll(mesh.C(), c)
        {
            psi[c] = mesh.C()[c].x() - x0;
            U[c] = u(mesh.C()[c]);
        }
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled()) continue;
            const vectorField& Cf = mesh.boundary()[patchi].Cf();
            forAll(Cf, i)
            {
                psi.boundaryFieldRef()[patchi][i] = Cf[i].x() - x0;
                U.boundaryFieldRef()[patchi][i] = u(Cf[i]);
            }
        }
        psi.correctBoundaryConditions();
        U.correctBoundaryConditions();
        // The exact face flux of u (not an interpolation), so that the test
        // measures the correction alone.
        forAll(phi, f)
        {
            phi[f] = u(mesh.Cf()[f]) & mesh.Sf()[f];
        }
        forAll(phi.boundaryField(), patchi)
        {
            fvsPatchScalarField& pp = phi.boundaryFieldRef()[patchi];
            const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
            const vectorField& Sp = mesh.Sf().boundaryField()[patchi];
            forAll(pp, i)
            {
                pp[i] = u(Cf[i]) & Sp[i];
            }
        }
    };
    setFields([](const point&) { return vector::zero; });

    autoPtr<velocityExtension> ext(velocityExtension::New(mesh));
    const auto& hl = refCast<const velocityExtensions::haloLimited>(ext());
    const extensionTravel& travel = hl.travelStrategy();
    const extensionWeight& weight = hl.weightStrategy();
    const extensionDirection& direction = hl.directionStrategy();

    // ---- 1. travel and weight ----------------------------------------------------
    {
        const scalar R = 0.05;
        check("S(0) = 0", travel.S(0, R), 0, 0);
        check("c(0) = 1", travel.fraction(0, R), 1, 0);
        const scalar eps = 1e-7*R;
        check("S'(0) = 1", (travel.S(eps, R) - travel.S(-eps, R))/(2*eps), 1, 1e-9);
        scalar worstOdd = 0, worstCap = -GREAT, worstT = 0;
        for (scalar t = 1e-3; t < 1e6; t *= 1.7)
        {
            worstOdd = max(worstOdd, mag(travel.S(t*R, R) + travel.S(-t*R, R)));
            const scalar r = mag(travel.S(t*R, R))/R;
            if (r > worstCap)
            {
                worstCap = r;
                worstT = t;
            }
        }
        Info<< "travel cap: max |S|/R - 1 = " << (worstCap - 1) << " at |d|/R " << worstT << nl;
        check("S odd", worstOdd, 0, 0);
        check("|S| <= R up to |d| = 1e6 R", worstCap <= 1 ? 0 : worstCap, 0, 0);
        // 1 - w = (beta/(2m)) t^(2m) + ..., beta = 1, m = 2
        const scalar t = 1e-2;
        const scalar ratio = (1 - weight.w(travel.fraction(t*R, R)))/Foam::pow(t, 4);
        check("1 - w = (beta/(2m)) (d/R)^(2m)", ratio, 0.25, 1e-3);
    }

    // ---- 2. direction --------------------------------------------------------------
    {
        scalar d; vector e;
        direction.evaluate(0.3, vector(0.8, 0, 0), d, e);
        check("direction: d = psi/q for q >= 0.5", d, 0.3/0.8, 1e-15);
        check("direction: e = g/q for q >= 0.5", mag(e - vector(1, 0, 0)), 0, 1e-15);
        scalar minQ2 = GREAT;
        for (scalar q = 0; q < 0.6; q += 1e-3)
        {
            minQ2 = min(minQ2, extensionDirections::levelSet::Q2(q*q));
        }
        check("direction: Q2 > 0.15 for every q", minQ2 > 0.15 ? 0 : 1, 0, 0);
    }

    // The per-face sample geometry, recomputed here from the model's own
    // strategies with the planar, exact d_f = x_f - x0 and e = e_x.
    const scalarField& hc = hl.cellSize().primitiveField();
    auto faceR = [&](const label f) -> scalar
    {
        return hl.radiusCells()*min(hc[mesh.owner()[f]], hc[mesh.neighbour()[f]]);
    };

    // ---- 3. planar level set, affine velocity --------------------------------------
    {
        const scalar a = 0.7, b = 0.3, c = -0.4;
        auto u = [&](const point& x) { return vector(a*x.x() + b*x.y(), c*x.x() - a*x.y(), 0); };
        setFields(u);
        ext->correct();
        const vector dUdx(a, c, 0);
        const surfaceScalarField corr(ext->phi() - phi);
        scalar worst = 0, worstB = 0;
        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            const scalar d = mesh.Cf()[f].x() - x0;
            const scalar R = faceR(f);
            const scalar want =
                -weight.w(travel.fraction(d, R))*travel.S(d, R)*(dUdx & mesh.Sf()[f]);
            worst = max(worst, mag(corr[f] - want)/mesh.magSf()[f]);
        }
        forAll(corr.boundaryField(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const scalarField& cp = corr.boundaryField()[patchi];
            if (!p.coupled())
            {
                if (cp.size()) worstB = max(worstB, max(mag(cp)));
                continue;
            }
            // coupled faces: the same formula, R from min(h_own, h_nbr)
            const scalarField& hN = hl.cellSize().boundaryField()[patchi];
            const labelUList& fc = p.faceCells();
            forAll(cp, i)
            {
                const scalar d = p.Cf()[i].x() - x0;
                const scalar R = hl.radiusCells()*min(hc[fc[i]], hN[i]);
                const scalar want =
                    -weight.w(travel.fraction(d, R))*travel.S(d, R)*(dUdx & p.Sf()[i]);
                worst = max(worst, mag(cp[i] - want)/p.magSf()[i]);
            }
        }
        reduce(worst, maxOp<scalar>());
        reduce(worstB, maxOp<scalar>());
        check("affine: face correction = -w S (dU/dx).Sf", worst, 0, 1e-12);
        check("affine: physical patches keep phi (bit for bit)", worstB, 0, 0);

        scalar worstU = 0;
        forAll(U, ci)
        {
            const scalar d = mesh.C()[ci].x() - x0;
            const scalar R = hl.radiusCells()*hc[ci];
            const vector want = U[ci] - weight.w(travel.fraction(d, R))*travel.S(d, R)*dUdx;
            worstU = max(worstU, mag(ext->Uext()[ci] - want));
        }
        reduce(worstU, maxOp<scalar>());
        check("affine: Uext = U - w S dU/dx", worstU, 0, 1e-12);

        const volScalarField& reach = mesh.lookupObject<volScalarField>("hlReach");
        check("reach |Y - x|/R < 1 in every cell", gMax(reach.primitiveField()) < 1 ? 0 : 1, 0, 0);

        // 6. coupled faces: the correction field carries exactly opposite values
        //    on the two sides, and phiExt adds no asymmetry to phi.
        auto antisymmetry = [&](const surfaceScalarField& fld, label& nCoupled) -> scalar
        {
            scalarField bc(mesh.nBoundaryFaces(), Zero);
            forAll(fld.boundaryField(), patchi)
            {
                const label b0 = mesh.boundary()[patchi].start() - mesh.nInternalFaces();
                forAll(fld.boundaryField()[patchi], i)
                {
                    bc[b0 + i] = fld.boundaryField()[patchi][i];
                }
            }
            scalarField bcNbr(bc);
            syncTools::swapBoundaryFaceList(mesh, bcNbr);
            scalar worstA = 0;
            nCoupled = 0;
            forAll(mesh.boundary(), patchi)
            {
                if (!mesh.boundary()[patchi].coupled()) continue;
                const label b0 = mesh.boundary()[patchi].start() - mesh.nInternalFaces();
                forAll(mesh.boundary()[patchi], i)
                {
                    worstA = max(worstA, mag(bc[b0 + i] + bcNbr[b0 + i]));
                    ++nCoupled;
                }
            }
            reduce(worstA, maxOp<scalar>());
            reduce(nCoupled, sumOp<label>());
            return worstA;
        };
        label nCoupled = 0;
        const scalar antiCorr =
            antisymmetry(mesh.lookupObject<surfaceScalarField>("hlCorrection"), nCoupled);
        const scalar antiPhi = antisymmetry(phi, nCoupled);
        const scalar antiExt = antisymmetry(ext->phi(), nCoupled);
        Info<< "coupled faces: " << nCoupled << nl;
        check("coupled faces: hlCorrection exactly opposite on the two sides", antiCorr, 0, 0);
        check("coupled faces: phiExt as antisymmetric as phi", antiExt, antiPhi, 0);
    }

    // ---- 4. uniform velocity -------------------------------------------------------
    {
        setFields([](const point&) { return vector(0.3, -0.2, 0); });
        ext->correct();
        scalar dphi = max(mag(ext->phi().primitiveField() - phi.primitiveField()));
        scalar dU = max(mag(ext->Uext().primitiveField() - U.primitiveField()));
        forAll(phi.boundaryField(), patchi)
        {
            if (phi.boundaryField()[patchi].size())
            {
                dphi = max(dphi, max(mag(ext->phi().boundaryField()[patchi] - phi.boundaryField()[patchi])));
            }
        }
        reduce(dphi, maxOp<scalar>());
        reduce(dU, maxOp<scalar>());
        check("uniform U: phiExt == phi (bit for bit)", dphi, 0, 0);
        check("uniform U: Uext == U (bit for bit)", dU, 0, 0);
    }

    // ---- 5. the dossier's quadratic velocity ---------------------------------------
    {
        const scalar alpha = 0.6, beta = 0.9;
        auto u = [&](const point& x)
        {
            return vector(alpha*x.x() + beta*sqr(x.x()), -(alpha + 2*beta*x.x())*x.y(), 0);
        };
        setFields(u);
        ext->correct();
        const surfaceScalarField corr(ext->phi() - phi);
        // interior faces only: both cells at least two cells from every wall,
        // where the stencil carries the full quadratic model
        const scalar inner = 1 - 2.5*0.05;
        scalar worst = 0;
        label nFaces = 0;
        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            const point& xf = mesh.Cf()[f];
            const point& xP = mesh.C()[mesh.owner()[f]];
            const point& xN = mesh.C()[mesh.neighbour()[f]];
            if (max(max(mag(xP.x()), mag(xP.y())), max(mag(xN.x()), mag(xN.y()))) > inner) continue;
            const scalar d = xf.x() - x0;
            const scalar R = faceR(f);
            const point Y = xf - travel.S(d, R)*vector(1, 0, 0);
            const scalar want = weight.w(travel.fraction(d, R))*((u(Y) - u(xf)) & mesh.Sf()[f]);
            worst = max(worst, mag(corr[f] - want)/mesh.magSf()[f]);
            ++nFaces;
        }
        reduce(worst, maxOp<scalar>());
        reduce(nFaces, sumOp<label>());
        Info<< "quadratic velocity: " << nFaces << " interior faces" << nl;
        check("quadratic u: correction = w [u(Y) - u(x)].Sf with the exact u", worst, 0, 1e-11);
    }

    Info<< nl << "==== " << nPass << " passed, " << nFail << " failed ====" << nl << endl;
    if (nFail)
    {
        Info<< "leiaTestHaloLimited FAILED" << nl << endl;
        return 1;
    }
    Info<< "End\n" << endl;
    return 0;
}
