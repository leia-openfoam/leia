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
      (b) polynomial exactness -- with psi a global quadratic, nestedLSQ and
          quadraticWLSQ reproduce it at an interior off-centre point to ~1e-10;
          with psi a global linear field, linearTaylor reproduces it;
      (c) constant-velocity foot -- for uniform U the Taylor acceleration term
          (du/dt + (u.grad)u) vanishes, so the foot reduces to x_d = x_c - u*dt.

    Runs in a meshed case (e.g. cases/2Dvortex after blockMesh + leiaSetFields).
    Writes PASS/FAIL lines and leiaTestSLReconstruction.csv; returns non-zero
    if any assertion fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "linearTaylorReconstruction.H"
#include "nestedLSQReconstruction.H"
#include "quadraticWLSQReconstruction.H"
#include "bandQuadraticWLSQReconstruction.H"
#include "slReconstruction.H"

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
    scalar errNestQuad = 0;   // nestedLSQ on the quadratic field
    scalar errQuadQuad = 0;   // quadraticWLSQ on the quadratic field

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
        nestedLSQReconstruction    nst(mesh);
        quadraticWLSQReconstruction quad(mesh);

        // Linear field: linearTaylor must be exact.
        setPsi(fLin);
        errLinInterp = testExactness(lin, fLin);

        // Quadratic field: nestedLSQ + quadraticWLSQ must be exact.
        setPsi(fQuad);
        errNestQuad = testExactness(nst, fQuad);
        errQuadQuad = testExactness(quad, fQuad);
        // (also refresh centre-reproduction for lin on the quadratic field)
        lin.update(psi);
        scalar wc = 0;
        forAll(C, c) { wc = Foam::max(wc, Foam::mag(lin.evaluate(c, C[c]) - psi[c])); }
        reduce(wc, maxOp<scalar>());
        worstCentre = Foam::max(worstCentre, wc);
    }

    // ===== (b') bandQuadraticWLSQ: exact in the band, == full quadratic ===== //
    scalar errBandExact = 0;    // band-cell reconstruction vs the analytic field
    scalar errBandVsFull = 0;   // band-cell reconstruction vs full quadraticWLSQ
    {
        // A zero-crossing quadratic (fQuad minus its value at the box centre) so
        // a band |psi| <= (nLayersBand+bandGuard)*h exists mid-domain; it is
        // still a global quadratic, so the band fit must be exact there.
        const point ctr = 0.5*(bb.min() + bb.max());
        const scalar shift = fQuad(ctr);
        auto fBand = [&](const point& x) -> scalar { return fQuad(x) - shift; };
        setPsi(fBand);

        bandQuadraticWLSQReconstruction band(mesh);
        quadraticWLSQReconstruction     quad2(mesh);
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
    const bool passNest   = errNestQuad < 1e-7;
    const bool passQuad   = errQuadQuad < 1e-7;
    const bool passFoot   = maxAccel < 1e-8;
    const bool passBandExact = errBandExact < 1e-7;
    const bool passBandVsFull = errBandVsFull < 1e-9;
    const bool allPass =
        passCentre && passLin && passNest && passQuad && passFoot
        && passBandExact && passBandVsFull;

    Info<< nl << "=== leiaTestSLReconstruction ===" << nl
        << "  (a) centre reproduction  max|e| = " << worstCentre
        << "  [" << (passCentre ? "PASS" : "FAIL") << "]" << nl
        << "  (b) linearTaylor  (linear field)    = " << errLinInterp
        << "  [" << (passLin ? "PASS" : "FAIL") << "]" << nl
        << "  (b) nestedLSQ     (quadratic field) = " << errNestQuad
        << "  [" << (passNest ? "PASS" : "FAIL") << "]" << nl
        << "  (b) quadraticWLSQ (quadratic field) = " << errQuadQuad
        << "  [" << (passQuad ? "PASS" : "FAIL") << "]" << nl
        << "  (b')bandQuadraticWLSQ exact in band = " << errBandExact
        << "  [" << (passBandExact ? "PASS" : "FAIL") << "]" << nl
        << "  (b')bandQuadraticWLSQ == full quad  = " << errBandVsFull
        << "  [" << (passBandVsFull ? "PASS" : "FAIL") << "]" << nl
        << "  (c) const-velocity foot accel |max| = " << maxAccel
        << " (gradU max " << maxGradU << ")"
        << "  [" << (passFoot ? "PASS" : "FAIL") << "]" << nl
        << "  RESULT: " << (allPass ? "PASS" : "FAIL") << nl << endl;

    if (Pstream::master())
    {
        OFstream os("leiaTestSLReconstruction.csv");
        os << "TEST,VALUE,TOL,PASS\n";
        os << "centreReproduction," << worstCentre << "," << tol << ","
           << (passCentre ? 1 : 0) << "\n";
        os << "linearTaylorLinear," << errLinInterp << ",1e-8,"
           << (passLin ? 1 : 0) << "\n";
        os << "nestedLSQQuadratic," << errNestQuad << ",1e-7,"
           << (passNest ? 1 : 0) << "\n";
        os << "quadraticWLSQQuadratic," << errQuadQuad << ",1e-7,"
           << (passQuad ? 1 : 0) << "\n";
        os << "bandQuadraticExactInBand," << errBandExact << ",1e-7,"
           << (passBandExact ? 1 : 0) << "\n";
        os << "bandQuadraticVsFull," << errBandVsFull << ",1e-9,"
           << (passBandVsFull ? 1 : 0) << "\n";
        os << "constVelocityFootAccel," << maxAccel << ",1e-8,"
           << (passFoot ? 1 : 0) << "\n";
    }

    Info<< "End\n" << endl;
    return allPass ? 0 : 1;
}


// ************************************************************************* //
