/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.

Application
    leiaTestMeshTransfer

Description
    Unit gate of the mesh transfer of psi on a dynamic mesh (slMeshTransfer,
    STATUS 11.44, gate U1 of docs/plan-sl-amr-preprint.md).

    A sphere moves through the mesh in steps. In every pass the application
      1. sets psi to an EXACT test field on the present mesh,
      2. computes the refinement band of the sphere (slRefinementBand),
      3. calls store(psi), mesh.update(), apply(psi), as the solver does,
      4. compares psi with the exact field on the new mesh.
    Unchanged cells keep their value. A child cell and a merged cell get the
    value of the transfer. The error is reported for the cells within
    `radius + 0.15` of the sphere centre (away from the walls, where the patch
    values of a zeroGradient field are not exact data of the fit) and for all
    cells.

    Test fields (option -field):
      linear     a + b.r
      harmonic   a + b.r + 0.5 r.A.r with tr(A) = 0. The mean of the children of
                 a parent equals the centre value, so the restricted data of the
                 fit are exact, and the transfer must be exact to round-off.
      quadratic  the same with tr(A) = 3. The restricted data carry the term
                 (h^2/8) tr(A); the result measures that inconsistency.
      sphere     the signed distance of the sphere (not a polynomial: the error
                 is the third-order remainder).

    The application exits with 1 if the error of `linear` or `harmonic` near
    the sphere is above 1e-10 in any pass, or if no cell changed its level.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "slMeshTransfer.H"
#include "slRefinementBand.H"
#include "slReconstruction.H"
#include "labelIOList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption("field", "word", "linear | harmonic | quadratic | sphere");
    argList::addOption("passes", "label", "number of mesh updates (default 14)");
    argList::addOption("shift", "scalar", "sphere shift per pass (default 0.04)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    const word fieldType(args.getOrDefault<word>("field", "harmonic"));
    const label nPasses = args.getOrDefault<label>("passes", 14);
    const scalar shift = args.getOrDefault<scalar>("shift", 0.04);
    const scalar radius = 0.15;
    const point x0(1.0, 0.5, 0.5);
    const vector b(0.2, -0.1, 0.4);
    symmTensor A(1.0, 0.5, -0.3, -2.0, 0.7, 1.0);
    if (fieldType == "quadratic") { A += symmTensor(1, 0, 0, 1, 0, 1); }
    if (fieldType == "linear") { A = Zero; }

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField psiBand("psiBand", psi);

    point xc(0.5, 0.5, 0.5);
    auto exact = [&](const point& x) -> scalar
    {
        if (fieldType == "sphere") { return mag(x - xc) - radius; }
        const vector r(x - x0);
        return 0.3 + (b & r) + 0.5*(r & (A & r));
    };

    slRefinementBand band(mesh);
    autoPtr<slMeshTransfer> transfer(slMeshTransfer::New(mesh));

    scalar worstNear = 0;
    label nChanged = 0;
    for (label pass = 0; pass < nPasses; ++pass)
    {
        runTime.setTime(runTime.value(), runTime.timeIndex() + 1);
        xc = point(0.5 + shift*max(label(0), pass - 4), 0.5, 0.5);

        {
            const vectorField& C = mesh.C();
            scalarField& p = psi.primitiveFieldRef();
            scalarField& q = psiBand.primitiveFieldRef();
            forAll(p, c)
            {
                p[c] = exact(C[c]);
                q[c] = mag(C[c] - xc) - radius;
            }
            psi.correctBoundaryConditions();
            psiBand.correctBoundaryConditions();
        }

        band.calc(psiBand);
        transfer->store(psi);
        mesh.update();
        if (!mesh.topoChanging())
        {
            Info<< "pass " << pass << ": no mesh change" << endl;
            continue;
        }
        const label n = transfer->apply(psi);
        nChanged += n;
        transfer.clear();
        transfer = slMeshTransfer::New(mesh);

        const vectorField& C = mesh.C();
        const scalarField& V = mesh.V();
        const scalarField& p = psi.primitiveField();
        scalar eMaxNear = 0, eMaxAll = 0, eL1Near = 0, vNear = 0;
        forAll(p, c)
        {
            const scalar e = mag(p[c] - exact(C[c]));
            eMaxAll = max(eMaxAll, e);
            if (mag(C[c] - xc) < radius + 0.15)
            {
                eMaxNear = max(eMaxNear, e);
                eL1Near += e*V[c];
                vNear += V[c];
            }
        }
        reduce(eMaxNear, maxOp<scalar>());
        reduce(eMaxAll, maxOp<scalar>());
        reduce(eL1Near, sumOp<scalar>());
        reduce(vNear, sumOp<scalar>());
        worstNear = max(worstNear, eMaxNear);

        // FIT CHECK on the new mesh: the gradient and the value of the fit of the
        // EXACT field, per refinement level, for the cells near the sphere.
        if (fieldType != "sphere")
        {
            volScalarField ex("psiExact", psi);
            forAll(C, c) { ex.primitiveFieldRef()[c] = exact(C[c]); }
            ex.correctBoundaryConditions();
            autoPtr<slReconstruction> rec(slReconstruction::New(mesh));
            rec->update(ex);
            const labelIOList& lev = mesh.lookupObject<labelIOList>("cellLevel");
            const labelListList& cc = mesh.cellCells();
            scalarList gErr(4, Zero), vErr(4, Zero);
            labelList nBad(4, Zero), nLev(4, Zero);
            label printed = 0;
            forAll(C, c)
            {
                if (mag(C[c] - xc) >= radius + 0.15) { continue; }
                const label L = min(lev[c], label(3));
                vector g; symmTensor H;
                rec->fitDerivatives(c, g, H);
                const vector r(C[c] - x0);
                const scalar eg = mag(g - (b + (A & r)));
                const point xs(C[c] + 0.3*Foam::cbrt(V[c])*vector(1, -1, 1));
                const scalar ev = mag(rec->evaluateRaw(c, xs) - exact(xs));
                ++nLev[L];
                gErr[L] = max(gErr[L], eg);
                vErr[L] = max(vErr[L], ev);
                if (eg > 1e-9)
                {
                    ++nBad[L];
                    if (printed < 6)
                    {
                        ++printed;
                        label lo = 99, hi = -1;
                        for (const label n : cc[c]) { lo = min(lo, lev[n]); hi = max(hi, lev[n]); }
                        Pout<< "  bad fit: cell " << c << " centre " << C[c]
                            << " level " << lev[c] << " face-neighbour levels "
                            << lo << ".." << hi << " faces " << mesh.cells()[c].size()
                            << " |dg| " << eg << " H " << H << endl;
                    }
                }
            }
            for (label L = 0; L < 4; ++L)
            {
                reduce(nLev[L], sumOp<label>()); reduce(nBad[L], sumOp<label>());
                reduce(gErr[L], maxOp<scalar>()); reduce(vErr[L], maxOp<scalar>());
                if (nLev[L])
                {
                    Info<< "  fit check level " << L << ": cells " << nLev[L]
                        << ", gradient error above 1e-9 in " << nBad[L]
                        << ", max gradient error " << gErr[L]
                        << ", max value error " << vErr[L] << endl;
                }
            }
        }
        Info<< "pass " << pass << ": cells "
            << returnReduce(mesh.nCells(), sumOp<label>())
            << ", changed " << n
            << ", error near the sphere: max " << eMaxNear
            << ", L1 " << eL1Near/max(vNear, VSMALL)
            << "; all cells: max " << eMaxAll << endl;
    }

    Info<< nl << "leiaTestMeshTransfer field " << fieldType
        << ": changed cells " << nChanged
        << ", worst error near the sphere " << worstNear << endl;

    const bool exactField = (fieldType == "linear" || fieldType == "harmonic");
    if (nChanged == 0 || (exactField && worstNear > 1e-10))
    {
        Info<< "FAIL" << nl << endl;
        return 1;
    }
    Info<< (exactField ? "PASS" : "REPORT (no criterion for this field)") << nl << endl;
    return 0;
}
