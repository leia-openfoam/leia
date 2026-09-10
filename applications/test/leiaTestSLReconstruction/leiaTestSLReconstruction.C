/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestSLReconstruction

Description
    Unit tests for the semi-Lagrangian reconstructions (slReconstruction):

      (a) centre reproduction  -- evaluate(c, x_c) == psiOld[c] to machine
          precision, all reconstructions;
      (b) polynomial exactness -- with psi a global quadratic, quadraticTaylor
          and quadraticWeightedLeastSquares reproduce it at an interior off-centre
          point to ~1e-10; with psi a global linear field, linearTaylor and
          linearWeightedLeastSquares reproduce it exactly;
      (c) constant-velocity foot -- for uniform U the Taylor acceleration term
          (du/dt + (u.grad)u) vanishes, so the foot reduces to x_d = x_c - u*dt.

    Runs in a meshed case (e.g. cases/2Dvortex after blockMesh + leiaSetFields).
    Writes PASS/FAIL lines and leiaTestSLReconstruction.csv; returns non-zero
    if any assertion fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "linearTaylorReconstruction.H"
#include "linearWeightedLeastSquaresReconstruction.H"
#include "quadraticTaylorReconstruction.H"
#include "quadraticWeightedLeastSquaresReconstruction.H"
#include "uncachedQuadraticWeightedLeastSquaresReconstruction.H"
#include "signedDistanceQuadraticWeightedLeastSquaresReconstruction.H"
#include "bandQuadraticWeightedLeastSquaresReconstruction.H"
#include "defectCorrectedIDWReconstruction.H"
#include "slReconstruction.H"
#include "slValueBound.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field psi (container + BCs; values are overwritten)\n"
        << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const volVectorField& C = mesh.C();
    const Vector<label> gd = mesh.geometricD();

    // Analytic global polynomials. Quadratic:
    //   f(x) = a + g.x + 1/2 x^T Q x
    const scalar a0 = 0.37;
    const vector g0(0.7, -0.4, 0.0);
    const symmTensor Q0(1.3, 0.5, 0.0, -0.9, 0.0, 0.0);  // xx,xy,xz,yy,yz,zz
    auto fQuad = [&](const point& x) -> scalar
    {
        return a0 + (g0 & x) + 0.5*(x & (Q0 & x));
    };
    auto fLin = [&](const point& x) -> scalar
    {
        return a0 + (g0 & x);
    };

    // Cell length scale = max face-neighbour centre distance (robust on 2D
    // meshes of arbitrary thickness, unlike V^(1/nd)); used to keep the
    // exactness check on genuinely interior cells whose cell-point-cell stencil
    // does not touch a physical boundary (where the zeroGradient BC pollutes
    // the LSQ fit of an imposed global polynomial).
    const scalar h = gMax(1.0/mesh.deltaCoeffs().primitiveField());
    const boundBox bb(mesh.points());

    scalar worstCentre = 0;   // over all reconstructions, both fields
    scalar errLinInterp = 0;  // linearTaylor on the linear field
    scalar errLinWLSLin = 0;  // linearWeightedLeastSquares on the linear field (exact)
    scalar errLinWLSQuad = 0; // linearWeightedLeastSquares on the quadratic field (O(h^2))
    scalar errNestQuad = 0;   // quadraticTaylor (was nestedLSQ) on the quadratic field
    scalar errQuadQuad = 0;   // quadraticWeightedLeastSquares on the quadratic field
    scalar errUncachedQuad = 0;   // uncachedQuadraticWeightedLeastSquares on quadratic
    scalar errUncachedVsCached = 0; // uncached vs cached quad, cell-for-cell (parity gate)
    scalar errDefectLin = 0;  // defectCorrectedIDW on the linear field (linear-exact)

    // ---- helper: set psi to an analytic field -------------------------- //
    auto setPsi = [&](std::function<scalar(const point&)> f)
    {
        scalarField v(mesh.nCells());
        forAll(v, c) { v[c] = f(C[c]); }
        psi.primitiveFieldRef() = v;
        psi.correctBoundaryConditions();
    };

    // ---- helper: centre reproduction + interior off-centre exactness --- //
    auto testExactness = [&]
    (
        slReconstruction& R,
        std::function<scalar(const point&)> f
    ) -> scalar
    {
        R.update(psi);
        scalar wCentre = 0, wInterp = 0;
        forAll(C, c)
        {
            wCentre = Foam::max(wCentre, Foam::mag(R.evaluate(c, C[c]) - psi[c]));

            // interior only (skip stencils truncated by the domain box)
            const point& xc = C[c];
            const bool interior =
                (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
             && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
            if (!interior) { continue; }

            const scalar r = R.stencilRadius(c);
            vector dir(1.0, 1.0, 0.0);
            if (gd[2] == 1) { dir.z() = 1.0; }
            dir /= Foam::mag(dir);
            const point xe = xc + 0.25*r*dir;
            wInterp = Foam::max(wInterp, Foam::mag(R.evaluate(c, xe) - f(xe)));
        }
        reduce(wCentre, maxOp<scalar>());
        reduce(wInterp, maxOp<scalar>());
        worstCentre = Foam::max(worstCentre, wCentre);
        return wInterp;
    };

    // ================= (a)+(b) reconstruction tests ==================== //
    {
        linearTaylorReconstruction lin(mesh);
        linearWeightedLeastSquaresReconstruction lwls(mesh);
        quadraticTaylorReconstruction    nst(mesh);
        quadraticWeightedLeastSquaresReconstruction quad(mesh);
        uncachedQuadraticWeightedLeastSquaresReconstruction uq(mesh);
        defectCorrectedIDWReconstruction defect(mesh);

        // Linear field: linearTaylor + linearWeightedLeastSquares + defectCorrectedIDW
        // must be exact. (defectCorrectedIDW is 2nd-order ACCURATE, not quadratic-exact,
        // so it is asserted on the linear field only; on quadratics it is O(h^2), not 0.
        // linearWeightedLeastSquares fits a linear polynomial to VALUES -> exact on a
        // globally linear field; on the quadratic field it is only O(h^2), asserted
        // with a loose bound below.)
        setPsi(fLin);
        errLinInterp = testExactness(lin, fLin);
        errLinWLSLin = testExactness(lwls, fLin);
        errDefectLin = testExactness(defect, fLin);

        // Quadratic field: quadraticTaylor + quadraticWeightedLeastSquares must be exact.
        // uncachedQuadraticWeightedLeastSquares is the SAME weighted LS fit, so it
        // must be quadratic-exact too AND match the cached quad cell-for-cell.
        // linearWeightedLeastSquares is NOT exact on a quadratic (loose O(h^2) bound).
        setPsi(fQuad);
        errLinWLSQuad = testExactness(lwls, fQuad);
        errNestQuad = testExactness(nst, fQuad);
        errQuadQuad = testExactness(quad, fQuad);
        errUncachedQuad = testExactness(uq, fQuad);
        forAll(C, c)   // parity: uncached vs cached quad at an interior off-centre point
        {
            const point& xc = C[c];
            const bool interior =
                (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
             && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
            if (!interior) { continue; }
            const scalar r = quad.stencilRadius(c);
            vector dir(1.0, 1.0, 0.0);
            if (gd[2] == 1) { dir.z() = 1.0; }
            dir /= Foam::mag(dir);
            const point xe = xc + 0.25*r*dir;
            errUncachedVsCached =
                Foam::max(errUncachedVsCached,
                          Foam::mag(uq.evaluate(c, xe) - quad.evaluate(c, xe)));
        }
        reduce(errUncachedVsCached, maxOp<scalar>());
        // (also refresh centre-reproduction for lin on the quadratic field)
        lin.update(psi);
        scalar wc = 0;
        forAll(C, c) { wc = Foam::max(wc, Foam::mag(lin.evaluate(c, C[c]) - psi[c])); }
        reduce(wc, maxOp<scalar>());
        worstCentre = Foam::max(worstCentre, wc);
    }

    // ===== (b') bandQuadraticWeightedLeastSquares: exact in the band, == full quadratic ===== //
    scalar errBandExact = 0;    // band-cell reconstruction vs the analytic field
    scalar errBandVsFull = 0;   // band-cell reconstruction vs full quadraticWeightedLeastSquares
    {
        // A zero-crossing quadratic (fQuad minus its value at the box centre) so
        // a band |psi| <= (nLayersBand+bandGuard)*h exists mid-domain; it is
        // still a global quadratic, so the band fit must be exact there.
        const point ctr = 0.5*(bb.min() + bb.max());
        const scalar shift = fQuad(ctr);
        auto fBand = [&](const point& x) -> scalar { return fQuad(x) - shift; };
        setPsi(fBand);

        bandQuadraticWeightedLeastSquaresReconstruction band(mesh);
        quadraticWeightedLeastSquaresReconstruction     quad2(mesh);
        band.update(psi);
        quad2.update(psi);

        const scalar bandW = 4.0*h;   // safely inside the fitted band
        forAll(C, c)
        {
            const point& xc = C[c];
            const bool interior =
                (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
             && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
            if (!interior || Foam::mag(psi[c]) > bandW) { continue; }

            const scalar r = band.stencilRadius(c);
            vector dir(1.0, 1.0, 0.0);
            if (gd[2] == 1) { dir.z() = 1.0; }
            dir /= Foam::mag(dir);
            const point xe = xc + 0.25*r*dir;
            errBandExact =
                Foam::max(errBandExact, Foam::mag(band.evaluate(c, xe) - fBand(xe)));
            errBandVsFull =
                Foam::max(errBandVsFull,
                          Foam::mag(band.evaluate(c, xe) - quad2.evaluate(c, xe)));
        }
        reduce(errBandExact, maxOp<scalar>());
        reduce(errBandVsFull, maxOp<scalar>());

        scalar wc = 0;
        forAll(C, c) { wc = Foam::max(wc, Foam::mag(band.evaluate(c, C[c]) - psi[c])); }
        reduce(wc, maxOp<scalar>());
        worstCentre = Foam::max(worstCentre, wc);
    }

    // === (b'') signedDistanceQuadraticWeightedLeastSquares: reproject a clean SDF === //
    // On an EXACT sphere signed-distance field the value fit has |grad P| ~ 1, so the
    // reprojection d = P0/|grad P| must reproduce the true distance (to ~O(h^2)) near
    // the interface. U is absent in this harness, so the normal-strain rescaling is
    // skipped -- this isolates the reprojection itself. Cells near the SDF's central
    // kink (excluded via the near-interface band) are not tested.
    scalar errSDF = 0;
    {
        // From the mesh, not from unit-box literals: hardcoded coordinates put
        // the sphere outside an SI-scale case, the band filter then matched no
        // cell, and errSDF = 0 was reported as a PASS having tested nothing.
        const point x0 = 0.5*(bb.min() + bb.max());
        const vector extSdf = bb.max() - bb.min();
        const scalar Rsdf = 0.2*Foam::min(extSdf.x(), extSdf.y());
        auto fSDF = [&](const point& x) -> scalar { return Foam::mag(x - x0) - Rsdf; };
        setPsi(fSDF);
        label nSDFtested = 0;

        signedDistanceQuadraticWeightedLeastSquaresReconstruction sd(mesh);
        sd.update(psi);

        const scalar bandW = 4.0*h;   // near the interface, away from the centre kink
        forAll(C, c)
        {
            const point& xc = C[c];
            const bool interior =
                (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
             && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
            if (!interior || Foam::mag(psi[c]) > bandW) { continue; }

            const scalar r = sd.stencilRadius(c);
            vector dir(1.0, 1.0, 0.0);
            if (gd[2] == 1) { dir.z() = 1.0; }
            dir /= Foam::mag(dir);
            const point xe = xc + 0.25*r*dir;
            errSDF = Foam::max(errSDF, Foam::mag(sd.evaluate(c, xe) - fSDF(xe)));
            ++nSDFtested;
        }
        reduce(errSDF, maxOp<scalar>());
        reduce(nSDFtested, sumOp<label>());
        // A test that examined no cell has not passed.
        if (nSDFtested == 0)
        {
            errSDF = GREAT;
            Info<< "  (b'') WARNING: no cell met the SDF band filter" << endl;
        }
    }

    // ======= (d) the distance-cone bound: admissibility and contraction ======= //
    // The bound is
    //     l_c = max_j ( psi_j - L|x_d - x_j| ),  u_c = min_j ( psi_j + L|x_d - x_j| )
    // over the arrival cell's stencil, with L = 1 for a signed distance field.
    //
    // WHAT THIS GATE DECIDES. The falsified monotone clip is exact only on
    // monotone data, and a signed distance field is never monotone in a closed
    // domain: it has genuine interior extrema on its medial axis. The claim for
    // the cone bound is different and much stronger -- it is exact on the class
    // of fields the method transports, for ANY interface shape and ANY mesh:
    //
    //   1. For L-Lipschitz data the true departure value lies inside [l_c, u_c],
    //      so the interval is never EMPTY (admissibility) and never excludes the
    //      truth (containment).
    //   2. Clipping onto an interval that contains the true value moves the
    //      candidate toward it, so the guarded error cannot exceed the raw fit's
    //      error, cell by cell (contraction). This is why the bound cannot
    //      degrade the order on compatible data.
    //
    // If any of these fails on an EXACT distance field, the bound is wrong and
    // no coupled run is worth starting.
    scalar coneLo1D = 0, coneHi1D = 0, coneMonoErr1D = 0;   // the review's 1D case
    label nConeInadmissible = 0;      // empty intervals on exact data: must be 0
    scalar coneContain = 0;           // how far the truth fell OUTSIDE [l_c, u_c]
    scalar coneWorsened = 0;          // worst error INCREASE caused by clipping
    scalar coneTightPlane = 0;        // |l_c - psi_true| on a plane: must be ~0
    label nConeApexTested = 0;        // apex cells where the monotone clip errs
    scalar coneApexMonoErr = 0;       // the monotone clip's error there
    scalar coneApexConeErr = 0;       // the cone bound's error at the same cells
    {
        // ---- (d1) the review's verified 1D counterexample, pure arithmetic.
        // psi(x) = |x - 0.25| - 2 sampled at x = -1, 0, 1. The true value at the
        // departure point x_d = 0.25 is -2, which is BELOW the sample minimum
        // -1.75, so a clip to the sample range commits an error of 0.25. The cone
        // interval is [-2, -1.5] and admits the true value at its lower end.
        {
            const scalar xs[3] = {-1.0, 0.0, 1.0};
            const scalar xd = 0.25;
            scalar lo = -GREAT, hi = GREAT, smin = GREAT, smax = -GREAT;
            for (label i = 0; i < 3; ++i)
            {
                const scalar pj = Foam::mag(xs[i] - 0.25) - 2.0;
                const scalar r = Foam::mag(xd - xs[i]);
                lo = Foam::max(lo, pj - r);
                hi = Foam::min(hi, pj + r);
                smin = Foam::min(smin, pj);
                smax = Foam::max(smax, pj);
            }
            const scalar exact = -2.0;
            coneLo1D = lo;
            coneHi1D = hi;
            coneMonoErr1D = Foam::mag(Foam::min(Foam::max(exact, smin), smax) - exact);
        }

        // ---- (d2)/(d3)/(d4) on the mesh, through the public accessor.
        uncachedQuadraticWeightedLeastSquaresReconstruction R(mesh);

        // GEOMETRY FROM THE MESH, never from unit-box literals. A hardcoded
        // x0 = (0.5, 0.5) and R = 0.25 lie entirely OUTSIDE an SI-scale case
        // (the Popinet box is 5 mm x 2.5 mm), so every filter matches nothing
        // and the test reports a VACUOUS pass -- the same defect class as a
        // mesher that adds no cells and still exits 0.
        //
        // The sphere centre is offset by HALF A CELL off the box centre on
        // purpose. The distance cone's apex sitting exactly on a cell centre is
        // a measure-zero special case in which the apex value is itself a
        // stencil value and the monotone clip happens to be right. The generic
        // position -- apex BETWEEN cell centres -- is the one that matters, and
        // it is where the monotone clip is guaranteed wrong.
        const point boxCtr = 0.5*(bb.min() + bb.max());
        const vector ext = bb.max() - bb.min();
        const scalar Rs = 0.2*Foam::min(ext.x(), ext.y());
        const point x0 = boxCtr + 0.5*h*vector(1, 1, (gd[2] == 1 ? 1 : 0));
        vector nrm(0.6, 0.8, 0.0);            // |nrm| == 1 exactly
        auto fPlane = [&](const point& x) -> scalar { return nrm & (x - x0); };
        auto fSphere = [&](const point& x) -> scalar
        {
            return Foam::mag(x - x0) - Rs;
        };

        for (label field = 0; field < 2; ++field)
        {
            std::function<scalar(const point&)> f =
                (field == 0)
              ? std::function<scalar(const point&)>(fPlane)
              : std::function<scalar(const point&)>(fSphere);

            setPsi(f);
            R.update(psi);

            forAll(C, c)
            {
                const point& xc = C[c];
                const bool interior =
                    (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
                 && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
                if (!interior) { continue; }

                // A departure displacement inside the stencil hull, as the SL
                // trace produces.
                const scalar r = R.stencilRadius(c);
                vector dir(1.0, 1.0, 0.0);
                if (gd[2] == 1) { dir.z() = 1.0; }
                dir /= Foam::mag(dir);
                const point xd = xc - 0.25*r*dir;

                scalar lo, hi;
                const bool ok = R.stencilConeRange(c, xd, 1.0, false, lo, hi);
                if (!ok)
                {
                    ++nConeInadmissible;
                    continue;
                }

                const scalar exact = f(xd);
                // Containment: the truth must be inside. A POSITIVE value here is
                // a bound that would clip a correct value.
                coneContain = Foam::max
                (
                    coneContain,
                    Foam::max(lo - exact, exact - hi)
                );

                // Contraction: clipping must not increase the pointwise error.
                const scalar Hc = R.evaluateRaw(c, xd);
                const scalar Hb = Foam::min(Foam::max(Hc, lo), hi);
                coneWorsened = Foam::max
                (
                    coneWorsened,
                    Foam::mag(Hb - exact) - Foam::mag(Hc - exact)
                );

                if (field == 0)
                {
                    // On a plane the bound is TIGHT: a stencil point lying
                    // upstream along the gradient gives psi_j + |x_d - x_j|
                    // exactly equal to psi(x_d), so l_c reproduces the truth.
                    coneTightPlane =
                        Foam::max(coneTightPlane, Foam::mag(lo - exact));
                }
            }
        }

        // ---- THE CASE THAT KILLED THE MONOTONE CLIP: the distance cone's apex.
        // psi = |x - x0| - R has a genuine interior MINIMUM at x0. When the apex
        // lies BETWEEN cell centres, the true value at a departure point near it
        // is below every sampled value, so a clip to the stencil range is
        // guaranteed to be wrong -- and the exemption written to protect that
        // cell is what falsified the monotone clip at gate G4.
        //
        // SUB-CELL ALIGNMENT IS SWEPT, not chosen. The apex sitting exactly on a
        // cell centre is a measure-zero special case in which the apex value IS
        // a stencil value and the monotone clip happens to be right; a single
        // hardcoded offset can land there by accident, and did (an offset of half
        // a cell from the box centre put the apex exactly on a centre, because
        // the box centre falls on a FACE). The review of 2026-09-09 asks for
        // "a translating disk/sphere at several subcell alignments" for exactly
        // this reason.
        //
        // The foot is aimed AT the apex, which is the direction that puts the
        // departure point below the sampled minimum.
        for (label a = 0; a < 8; ++a)
        {
            const scalar frac = a/8.0;          // 0, 1/8, ... 7/8 of a cell
            const point xa =
                boxCtr + frac*h*vector(1, 0.7, (gd[2] == 1 ? 0.4 : 0));
            auto fApex = [&](const point& x) -> scalar
            {
                return Foam::mag(x - xa) - Rs;
            };
            setPsi(fApex);
            R.update(psi);

            forAll(C, c)
            {
                const point& xc = C[c];
                // Only cells whose stencil can reach the apex.
                if (Foam::mag(xc - xa) > 2.0*h) { continue; }
                const bool interior =
                    (xc.x() - bb.min().x() > 4*h) && (bb.max().x() - xc.x() > 4*h)
                 && (xc.y() - bb.min().y() > 4*h) && (bb.max().y() - xc.y() > 4*h);
                if (!interior) { continue; }

                const vector toApex = xa - xc;
                const scalar dm = Foam::mag(toApex);
                if (dm < SMALL) { continue; }    // the apex IS this cell centre

                // A departure point at the apex itself: the sharpest case, and
                // the one a translating droplet reaches every step.
                const point xd = xa;
                const scalar exact = fApex(xd);  // == -Rs

                scalar slo, shi;
                R.stencilRange(c, slo, shi);
                if (exact >= slo && exact <= shi) { continue; }   // no violation

                scalar lo, hi;
                const bool ok = R.stencilConeRange(c, xd, 1.0, false, lo, hi);
                ++nConeApexTested;
                if (!ok)
                {
                    ++nConeInadmissible;
                    continue;
                }
                coneApexMonoErr = Foam::max
                (
                    coneApexMonoErr,
                    Foam::mag(Foam::min(Foam::max(exact, slo), shi) - exact)
                );
                coneApexConeErr = Foam::max
                (
                    coneApexConeErr,
                    Foam::max(lo - exact, exact - hi)
                );
            }
        }

        reduce(nConeInadmissible, sumOp<label>());
        reduce(nConeApexTested, sumOp<label>());
        reduce(coneContain, maxOp<scalar>());
        reduce(coneWorsened, maxOp<scalar>());
        reduce(coneTightPlane, maxOp<scalar>());
        reduce(coneApexMonoErr, maxOp<scalar>());
        reduce(coneApexConeErr, maxOp<scalar>());
    }

    // ================= (c) constant-velocity foot ====================== //
    volVectorField U
    (
        IOobject("Utest", runTime.timeName(), mesh, IOobject::NO_READ,
                 IOobject::NO_WRITE),
        mesh,
        dimensionedVector("Utest", dimVelocity, vector(0.8, -0.3, 0.0)),
        "zeroGradient"
    );
    const volTensorField gradU(fvc::grad(U, "gradU"));
    scalar maxGradU = gMax(mag(gradU)().primitiveField());
    // Foot acceleration term for a uniform field: du/dt = 0 (u^{n+1}=u^n) and
    // (u.grad)u = u & gradU ~ 0 -> the dt^2 displacement term vanishes.
    scalar maxAccel = 0;
    forAll(C, c)
    {
        const vector accel = (U[c] & gradU[c]);   // du/dt term is exactly 0
        maxAccel = Foam::max(maxAccel, Foam::mag(accel));
    }
    reduce(maxAccel, maxOp<scalar>());

    // ========================= verdict ================================= //
    const scalar tol = 1e-9;
    const bool passCentre = worstCentre < tol;
    const bool passLin    = errLinInterp < 1e-8;
    const bool passLinWLS = errLinWLSLin < 1e-8;    // linearWLSQ exact on a linear field
    // linearWLSQ on a quadratic is only O(h^2) (a linear fit): loose bound, like the SDF gate.
    const bool passLinWLSQuad = errLinWLSQuad < 5e-2;
    const bool passNest   = errNestQuad < 1e-7;
    const bool passQuad   = errQuadQuad < 1e-5;
    const bool passUncached = errUncachedQuad < 1e-5;        // same fit -> quadratic-exact
    const bool passParity   = errUncachedVsCached < 1e-5;    // matches cached quad
    const bool passDefect = errDefectLin < 1e-8;   // defectCorrectedIDW: linear-exact
    const bool passFoot   = maxAccel < 1e-8;
    const bool passBandExact = errBandExact < 1e-7;
    const bool passBandVsFull = errBandVsFull < 1e-5;
    // SDF reprojection is only O(h^2)-accurate (SDF is non-polynomial): assert a loose
    // bound that a correct reprojection meets easily but a broken one (O(1)) fails.
    const bool passSDF = errSDF < 5e-2;

    // (d) the distance-cone bound. Every one of these is a hard assertion on
    // EXACT distance data, where the mathematics leaves no room: an empty
    // interval, a truth outside the interval, or an error made worse by clipping
    // would each mean the bound is wrong, not merely loose.
    const bool passCone1D =
        Foam::mag(coneLo1D + 2.0) < 1e-12          // l_c == -2
     && Foam::mag(coneHi1D + 1.5) < 1e-12          // u_c == -1.5
     && Foam::mag(coneMonoErr1D - 0.25) < 1e-12;   // the monotone clip errs 0.25
    const bool passConeAdmissible = (nConeInadmissible == 0);
    const bool passConeContain = coneContain < 1e-12;
    const bool passConeContract = coneWorsened < 1e-12;
    // Tightness on a plane is limited by the stencil's angular coverage, not by
    // round-off: l_c is exact only if some stencil point lies ON the gradient ray
    // through x_d. On a hexahedral stencil the nearest point is off that ray by
    // up to half a cell, so the gap is O(h). A loose bound that a correct
    // implementation meets easily and a broken one (O(1)) fails.
    const bool passConeTight = coneTightPlane < 2.0*h;
    // The apex must actually be REACHED, or the sharpest claim is untested.
    const bool passConeApex = (nConeApexTested > 0) && (coneApexConeErr < 1e-12);

    const bool allPass =
        passCentre && passLin && passLinWLS && passLinWLSQuad
        && passNest && passQuad && passDefect && passFoot
        && passBandExact && passBandVsFull && passUncached && passParity && passSDF
        && passCone1D && passConeAdmissible && passConeContain
        && passConeContract && passConeTight && passConeApex;

    Info<< nl << "=== leiaTestSLReconstruction ===" << nl
        << "  (a) centre reproduction  max|e| = " << worstCentre
        << "  [" << (passCentre ? "PASS" : "FAIL") << "]" << nl
        << "  (b) linearTaylor  (linear field)    = " << errLinInterp
        << "  [" << (passLin ? "PASS" : "FAIL") << "]" << nl
        << "  (b) linearWeightedLeastSquares (linear field, exact) = " << errLinWLSLin
        << "  [" << (passLinWLS ? "PASS" : "FAIL") << "]" << nl
        << "  (b) linearWeightedLeastSquares (quadratic field, O(h^2)) = " << errLinWLSQuad
        << "  [" << (passLinWLSQuad ? "PASS" : "FAIL") << "]" << nl
        << "  (b) quadraticTaylor (quadratic field) = " << errNestQuad
        << "  [" << (passNest ? "PASS" : "FAIL") << "]" << nl
        << "  (b) quadraticWeightedLeastSquares (quadratic field) = " << errQuadQuad
        << "  [" << (passQuad ? "PASS" : "FAIL") << "]" << nl
        << "  (b) uncachedQuadraticWeightedLeastSquares (quad fld) = " << errUncachedQuad
        << "  [" << (passUncached ? "PASS" : "FAIL") << "]" << nl
        << "  (b) uncached == cached quad (parity)  = " << errUncachedVsCached
        << "  [" << (passParity ? "PASS" : "FAIL") << "]" << nl
        << "  (b) defectCorrectedIDW (linear fld) = " << errDefectLin
        << "  [" << (passDefect ? "PASS" : "FAIL") << "]" << nl
        << "  (b')bandQuadraticWeightedLeastSquares exact in band = " << errBandExact
        << "  [" << (passBandExact ? "PASS" : "FAIL") << "]" << nl
        << "  (b')bandQuadraticWeightedLeastSquares == full quad  = " << errBandVsFull
        << "  [" << (passBandVsFull ? "PASS" : "FAIL") << "]" << nl
        << "  (b'')signedDistanceQuadratic SDF reproject = " << errSDF
        << "  [" << (passSDF ? "PASS" : "FAIL") << "]" << nl
        << "  (c) const-velocity foot accel |max| = " << maxAccel
        << " (gradU max " << maxGradU << ")"
        << "  [" << (passFoot ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound, review 1D case [" << coneLo1D << ", " << coneHi1D
        << "] vs monotone error " << coneMonoErr1D
        << "  [" << (passCone1D ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound, empty intervals on exact data = " << nConeInadmissible
        << "  [" << (passConeAdmissible ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound, truth outside interval = " << coneContain
        << "  [" << (passConeContain ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound, worst error INCREASE from clipping = " << coneWorsened
        << "  [" << (passConeContract ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound, plane tightness |l_c - exact| = " << coneTightPlane
        << " (tol " << 2.0*h << ")"
        << "  [" << (passConeTight ? "PASS" : "FAIL") << "]" << nl
        << "  (d) cone bound at " << nConeApexTested
        << " apex cells: monotone clip errs " << coneApexMonoErr
        << ", cone errs " << coneApexConeErr
        << "  [" << (passConeApex ? "PASS" : "FAIL") << "]" << nl
        << "  RESULT: " << (allPass ? "PASS" : "FAIL") << nl << endl;

    if (Pstream::master())
    {
        OFstream os("leiaTestSLReconstruction.csv");
        os << "TEST,VALUE,TOL,PASS\n";
        os << "centreReproduction," << worstCentre << "," << tol << ","
           << (passCentre ? 1 : 0) << "\n";
        os << "linearTaylorLinear," << errLinInterp << ",1e-8,"
           << (passLin ? 1 : 0) << "\n";
        os << "linearWeightedLeastSquaresLinear," << errLinWLSLin << ",1e-8,"
           << (passLinWLS ? 1 : 0) << "\n";
        os << "linearWeightedLeastSquaresQuadratic," << errLinWLSQuad << ",5e-2,"
           << (passLinWLSQuad ? 1 : 0) << "\n";
        os << "quadraticTaylorQuadratic," << errNestQuad << ",1e-7,"
           << (passNest ? 1 : 0) << "\n";
        os << "quadraticWeightedLeastSquaresQuadratic," << errQuadQuad << ",1e-5,"
           << (passQuad ? 1 : 0) << "\n";
        os << "uncachedQuadraticQuadratic," << errUncachedQuad << ",1e-5,"
           << (passUncached ? 1 : 0) << "\n";
        os << "uncachedVsCachedParity," << errUncachedVsCached << ",1e-5,"
           << (passParity ? 1 : 0) << "\n";
        os << "defectCorrectedIDWLinear," << errDefectLin << ",1e-8,"
           << (passDefect ? 1 : 0) << "\n";
        os << "bandQuadraticExactInBand," << errBandExact << ",1e-7,"
           << (passBandExact ? 1 : 0) << "\n";
        os << "bandQuadraticVsFull," << errBandVsFull << ",1e-5,"
           << (passBandVsFull ? 1 : 0) << "\n";
        os << "constVelocityFootAccel," << maxAccel << ",1e-8,"
           << (passFoot ? 1 : 0) << "\n";
        os << "signedDistanceSDFReproject," << errSDF << ",5e-2,"
           << (passSDF ? 1 : 0) << "\n";
        os << "coneBoundReview1Dcase," << coneMonoErr1D << ",0.25,"
           << (passCone1D ? 1 : 0) << "\n";
        os << "coneBoundEmptyIntervals," << nConeInadmissible << ",0,"
           << (passConeAdmissible ? 1 : 0) << "\n";
        os << "coneBoundTruthOutside," << coneContain << ",1e-12,"
           << (passConeContain ? 1 : 0) << "\n";
        os << "coneBoundErrorIncrease," << coneWorsened << ",1e-12,"
           << (passConeContract ? 1 : 0) << "\n";
        os << "coneBoundPlaneTightness," << coneTightPlane << "," << 2.0*h << ","
           << (passConeTight ? 1 : 0) << "\n";
        os << "coneBoundApexCells," << nConeApexTested << ",>0,"
           << (passConeApex ? 1 : 0) << "\n";
        os << "coneBoundApexMonotoneError," << coneApexMonoErr << ",NA,1\n";
    }

    Info<< "End\n" << endl;
    return allPass ? 0 : 1;
}


// ************************************************************************* //
