/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestGhostCellReconstruction

Description
    Exact-solution gate of the ghost-cell reconstruction
    (ghostCellQuadraticWeightedLeastSquares, C-Cp/L2+g) on a mesh with 2:1
    refinement transitions. It runs in a meshed case (refineHexMesh has written
    polyMesh/cellLevel) that holds a field psi.

    PRE-REGISTERED (2026-11-12, before the first build).
      E1. A LINEAR field is returned to round-off (1e-10 of the field range) at
          an off-centre point of EVERY cell whose stencil has no boundary face,
          the cells at a transition included. Reason: the mean of the children
          of a linear field is its value at the parent centre, and the fit of a
          coarse cell returns a linear field exactly at a ghost position.
      E2. A QUADRATIC field is returned to round-off in every such cell that has
          no cell of another level in its stencil.
      E3. In the cells at a transition the quadratic field is NOT returned
          exactly: the mean of the children differs from the centre value by
          h_fine^2/8 times the Laplacian. The error there is reported and is not
          a failure. It must be larger than round-off; if it is not, the ghost
          path did not run and the gate FAILS.
      The baseline uncachedQuadraticWeightedLeastSquares is run on the same
      fields and reported: it returns the quadratic everywhere.
    The application exits non-zero if E1, E2 or E3 fails.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "slReconstruction.H"
#include "labelIOList.H"
#include "syncTools.H"
#include <functional>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    labelIOList cellLevel
    (
        IOobject
        (
            "cellLevel",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh.nCells(), Zero)
    );

    const volVectorField& C = mesh.C();
    const boundBox bb(mesh.points());

    // a cell is "at a transition" if a point neighbour has another level;
    // processor neighbours are covered by the levels on the points
    labelList pMin(mesh.nPoints(), labelMax), pMax(mesh.nPoints(), labelMin);
    const labelListList& cp = mesh.cellPoints();
    forAll(cp, c)
    {
        forAll(cp[c], i)
        {
            pMin[cp[c][i]] = min(pMin[cp[c][i]], cellLevel[c]);
            pMax[cp[c][i]] = max(pMax[cp[c][i]], cellLevel[c]);
        }
    }
    syncTools::syncPointList(mesh, pMin, minEqOp<label>(), labelMax);
    syncTools::syncPointList(mesh, pMax, maxEqOp<label>(), labelMin);
    boolList transition(mesh.nCells(), false);
    forAll(cp, c)
    {
        forAll(cp[c], i)
        {
            if (pMin[cp[c][i]] != pMax[cp[c][i]]) { transition[c] = true; }
        }
    }

    const scalar a0 = 0.37;
    const vector g0(0.7, -0.4, 0.55);
    const symmTensor Q0(1.3, 0.5, -0.3, -0.9, 0.4, 0.8);
    auto fLin = [&](const point& x) { return a0 + (g0 & x); };
    auto fQuad = [&](const point& x) { return a0 + (g0 & x) + 0.5*(x & (Q0 & x)); };

    // errors at an off-centre point: [0] cells away from a transition,
    // [1] cells at a transition; cells within 4 coarse cells of the box are skipped
    auto test = [&]
    (
        const word& name,
        std::function<scalar(const point&)> f,
        scalar& eAway,
        scalar& eAt
    )
    {
        forAll(C, c) { psi.primitiveFieldRef()[c] = f(C[c]); }
        psi.correctBoundaryConditions();
        autoPtr<slReconstruction> R(slReconstruction::New(mesh, name));
        R->update(psi);
        const scalar range = gMax(psi.primitiveField()) - gMin(psi.primitiveField());
        eAway = 0; eAt = 0;
        forAll(C, c)
        {
            const scalar r = R->stencilRadius(c);
            const point& xc = C[c];
            bool interior = true;
            for (direction k = 0; k < 3; ++k)
            {
                if (xc[k] - bb.min()[k] < 1.5*r || bb.max()[k] - xc[k] < 1.5*r)
                {
                    interior = false;
                }
            }
            if (!interior) { continue; }
            const point x = xc + 0.2*r*vector(0.5, -0.3, 0.4);
            const scalar e = mag(R->evaluateRaw(c, x) - f(x))/range;
            if (transition[c]) { eAt = max(eAt, e); } else { eAway = max(eAway, e); }
        }
        reduce(eAway, maxOp<scalar>());
        reduce(eAt, maxOp<scalar>());
    };

    label nT = 0;
    forAll(transition, c) { if (transition[c]) { ++nT; } }
    reduce(nT, sumOp<label>());
    Info<< "cells at a transition: " << nT << " of "
        << returnReduce(mesh.nCells(), sumOp<label>()) << nl << endl;

    const word ghost("ghostCellQuadraticWeightedLeastSquares");
    const word base("uncachedQuadraticWeightedLeastSquares");
    scalar gLinA, gLinT, gQuadA, gQuadT, bLinA, bLinT, bQuadA, bQuadT;
    test(ghost, fLin, gLinA, gLinT);
    test(ghost, fQuad, gQuadA, gQuadT);
    test(base, fLin, bLinA, bLinT);
    test(base, fQuad, bQuadA, bQuadT);

    Info<< "relative error at an off-centre point (away from / at a transition)" << nl
        << "  ghost     linear     " << gLinA << "  " << gLinT << nl
        << "  ghost     quadratic  " << gQuadA << "  " << gQuadT << nl
        << "  baseline  linear     " << bLinA << "  " << bLinT << nl
        << "  baseline  quadratic  " << bQuadA << "  " << bQuadT << nl << endl;

    const scalar tol = 1e-10;
    const bool E1 = gLinA < tol && gLinT < tol;
    const bool E2 = gQuadA < tol;
    const bool E3 = nT == 0 || gQuadT > 100*tol;
    Info<< (E1 ? "PASS" : "FAIL") << "  E1 linear field exact in every cell" << nl
        << (E2 ? "PASS" : "FAIL") << "  E2 quadratic field exact away from a transition" << nl
        << (E3 ? "PASS" : "FAIL") << "  E3 the ghost path acts at a transition" << nl
        << endl;

    return (E1 && E2 && E3) ? 0 : 1;
}
