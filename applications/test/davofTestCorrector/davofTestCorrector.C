/*---------------------------------------------------------------------------*\
    leia :: davof :: davofTestCorrector
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.

Application
    davofTestCorrector

Description
    The unit check of the iterative face-fraction corrector of the davof
    reconstruction (leia davof-iter, 2026-11-09; davof.H "The iterative
    face-fraction corrector") on a static DAVOF state (the fields alpha and
    alphaf of the time directory, as davofInitFields writes them): the
    corrector is applied PASSES times to the state, each pass one
    reconstruction from the current face fractions, one evaluation of the
    face fractions the planes induce and their replacement (correctorIterations
    1, correctorTolerance 0), and the max and L1 change of alpha_f of every
    pass is printed with the interface-cell count, the faces set, the cells
    without a plane and the identity defect.

    What the passes must show:
    - normalModel exact on the exact sphere state: pass 1 changes alpha_f by
      the realizability defect of a plane against the sphere (O(h^2));
      every later pass changes it by round-off only, because the normal is
      the exact one whatever the face fractions and CCS positions the plane
      by alpha alone: the clip of the same plane returns the same fraction.
      This is the identity check of the clip and of the mean rule.
    - normalModel ringArea: the passes contract towards a fixed point of
      (reconstruct, induce); the factors per pass are the static version of
      prediction P2.
    - a planar state (davofInitFields detrixheAslam): pass 1 changes alpha_f
      by round-off (the planes reproduce the transported fractions of a
      plane), the static version of the plane gate P1.

    Extension of 2026-11-25 (the static test of the corrector with the ring
    linear fit, branch davof-iter-linfit; the output above is unchanged):
    - pass 0 is the reconstruction on the state as read; after every pass p
      (0 .. PASSES) one line PCL is printed with the interface cells and the
      cells holding a plane, the max and L1 change of alpha_f (p >= 1), the
      max and L1 of |n^p - n^(p-1)| over the cells that hold a plane in both
      passes (p >= 1), the errors of the normal against the exact sphere
      (system/initAlphaFieldDict: type sphere, origin, radius) with the
      definitions of the reconstructionError function object (over the
      cells that hold a plane: L1 and Linf of |n - n_e|, L1 and Linf of
      1 - n.n_e, L1 of the distance of the plane centre to the sphere; n_e
      the exact normal at the plane centre), the identity defect and, for
      the model ringLinear, the cells with a negative weight, the cells with
      y = 0 and the max |W (q - qbar).y|;
    - -dump TAG writes for every pass the file cellErr_TAG_p<p>.csv with, per
      cell that holds a plane, the cell label, the unit normal, |n - n_e|,
      1 - n.n_e and the centre error, so that the norms on a common cell set
      of two runs can be formed outside;
    - -plicRDF evaluates the plicRDF reconstruction scheme on the same alpha
      (two calls of reconstruct(true), as the gate's solver makes: the
      second is the reference row) with the same line and dump (tag plicRDF);
    - -writeTime T writes alpha.water and alphaf.water after the passes into
      the time directory T, so that a second call can start from them with
      -time T;
    - -timing R times, in the same process, R repetitions of (a) one
      reconstruct(true) of the plicRDF scheme on this alpha, (b) one
      reconstruct(true) of the davof scheme, (c) one corrector pass of the
      davof scheme (correctFaceFractions: two reconstructions and one
      induced-fraction evaluation) on the state as read, alpha_f restored
      between the repetitions (the restore is not timed); the mean and the
      minimum wall time of each are printed.

Usage
    davofTestCorrector [-normalModel exact|ringArea|ringLinear] [-passes N]
                       [-time T] [-dump TAG] [-plicRDF] [-writeTime T]
                       [-timing R]
    in a rendered case of the DAVOF advection studies (system/fvSolution with
    the alpha.water solver dictionary, 0/alpha.water, 0/alphaf.water, 0/U).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "reconstructionSchemes.H"
#include "davof.H"
#include "plicRDF.H"
#include "ringLinearNormal.H"
#include "implicitFunction.H"

#include <chrono>
#include <fstream>

using namespace Foam;

namespace
{

// The state of the evaluation: the exact sphere and the normals of the
// previous pass
struct Evaluator
{
    autoPtr<implicitFunction> func;
    point origin;
    scalar radius;
    vectorField prevN;
    bool havePrev = false;

    // One PCL line of the planes (normal, centre) of a scheme; the changes
    // of alpha_f (dAfMax, dAfL1 < 0: none), the identity (< 0: none) and the
    // ringLinear diagnostics (negW < 0: none) are given
    void line
    (
        const word& tag,
        const label p,
        const fvMesh& mesh,
        const volVectorField& normal,
        const volVectorField& centre,
        const label nInterface,
        const scalar dAfMax,
        const scalar dAfL1,
        const scalar identity,
        const label negW,
        const label y0,
        const scalar maxWy,
        const bool dump
    )
    {
        const label nCells = mesh.nCells();
        if (prevN.size() != nCells) prevN.setSize(nCells, vector::zero);

        scalar sumD = 0, maxD = 0, sumCos = 0, maxCos = 0, sumC = 0;
        scalar sumDn = 0, maxDn = 0;
        label nPlane = 0, nCommon = 0;
        vectorField curN(nCells, vector::zero);

        std::ofstream os;
        if (dump)
        {
            const fileName fn
            (
                mesh.time().path()/("cellErr_" + tag + "_p" + Foam::name(p) + ".csv")
            );
            os.open(fn.c_str());
            os << "cell,nx,ny,nz,dn,cosErr,dc\n";
            os.precision(17);
        }

        forAll(normal, celli)
        {
            if (mag(normal[celli]) == 0) continue;
            const vector n = normal[celli]/mag(normal[celli]);
            vector ne = func->grad(centre[celli]);
            ne /= mag(ne);
            const scalar d = mag(n - ne);
            const scalar cosErr = 1 - (n & ne);
            const scalar dc = mag(radius - mag(centre[celli] - origin));
            ++nPlane;
            sumD += d; maxD = max(maxD, d);
            sumCos += cosErr; maxCos = max(maxCos, cosErr);
            sumC += dc;
            curN[celli] = n;
            if (havePrev && mag(prevN[celli]) > 0)
            {
                const scalar dn = mag(n - prevN[celli]);
                sumDn += dn; maxDn = max(maxDn, dn);
                ++nCommon;
            }
            if (dump)
            {
                os << celli << ',' << n.x() << ',' << n.y() << ',' << n.z()
                   << ',' << d << ',' << cosErr << ',' << dc << '\n';
            }
        }
        reduce(sumD, sumOp<scalar>()); reduce(maxD, maxOp<scalar>());
        reduce(sumCos, sumOp<scalar>()); reduce(maxCos, maxOp<scalar>());
        reduce(sumC, sumOp<scalar>());
        reduce(sumDn, sumOp<scalar>()); reduce(maxDn, maxOp<scalar>());
        reduce(nPlane, sumOp<label>()); reduce(nCommon, sumOp<label>());
        const scalar np = max(nPlane, label(1));

        Info<< "PCL tag " << tag << " pass " << p
            << " interfaceCells " << nInterface << " planeCells " << nPlane
            << " dAfMax " << dAfMax << " dAfL1 " << dAfL1
            << " dNCommon " << nCommon << " dNMax " << (nCommon ? maxDn : -1)
            << " dNL1 " << (nCommon ? sumDn/nCommon : -1)
            << " L1dn " << sumD/np << " Linfdn " << maxD
            << " L1cos " << sumCos/np << " Linfcos " << maxCos
            << " L1centre " << sumC/np << " identity " << identity
            << " negW " << negW << " y0 " << y0 << " maxWy " << maxWy << endl;

        prevN = curN;
        havePrev = true;
    }
};

}


int main(int argc, char *argv[])
{
    argList::addOption("normalModel", "word", "the normal model (default: the dictionary's)");
    argList::addOption("passes", "label", "the corrector passes (default 6)");
    argList::addOption("time", "scalar", "the time directory to read (default: the start time)");
    argList::addOption("dump", "word", "tag: write cellErr_<tag>_p<p>.csv for every pass");
    argList::addBoolOption("plicRDF", "evaluate the plicRDF reconstruction on the state (reference row)");
    argList::addOption("writeTime", "scalar", "write alpha.water and alphaf.water into this time directory after the passes");
    argList::addOption("timing", "label", "time R repetitions of plicRDF reconstruct, davof reconstruct and one corrector pass");
    #include "setRootCase.H"
    #include "createTime.H"
    if (args.found("time"))
    {
        runTime.setTime(args.get<scalar>("time"), 0);
    }
    #include "createMesh.H"

    volScalarField alpha1
    (
        IOobject("alpha.water", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        fvc::flux(U)
    );

    // the alpha solver dictionary with the corrector forced on (one iteration per
    // pass, no tolerance stop) and the model of the command line
    dictionary dict(mesh.solverDict(alpha1.name()));
    dict.set("faceFractionCorrector", word("mean"));
    dict.set("correctorIterations", label(1));
    dict.set("correctorTolerance", scalar(0));
    dict.set("correctorLog", true);
    if (args.found("normalModel"))
    {
        dict.set("normalModel", args.get<word>("normalModel"));
        dict.set("applyTo", word("all"));
    }
    // the keys of the ring linear fit: one pass, the correction on, the diagnostics printed
    dict.set("ringLinearPasses", label(1));
    dict.set("ringLinearCorrection", true);
    dict.set("ringLinearDiag", true);
    const label passes = args.getOrDefault<label>("passes", 6);

    autoPtr<reconstructionSchemes> surfPtr
    (
        reconstructionSchemes::New(alpha1, phi, U, dict)
    );
    if (!isA<reconstruction::davof>(*surfPtr))
    {
        FatalErrorInFunction
            << "reconstructionScheme " << surfPtr->type() << " is not davof"
            << exit(FatalError);
    }
    reconstruction::davof& dv = refCast<reconstruction::davof>(*surfPtr);

    // the exact sphere of the state (system/initAlphaFieldDict)
    Evaluator ev;
    {
        IOdictionary initDict
        (
            IOobject
            (
                "initAlphaFieldDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        ev.origin = initDict.get<point>("origin");
        ev.radius = initDict.get<scalar>("radius");
        dictionary fd;
        fd.add("type", word("sphere"));
        fd.add("origin", ev.origin);
        fd.add("radius", ev.radius);
        ev.func = implicitFunction::New("sphere", fd);
    }
    const word dumpTag = args.getOrDefault<word>("dump", word::null);
    const bool dump = !dumpTag.empty();

    Info<< nl << "davofTestCorrector: normalModel " << dv.normalModel().type()
        << ", " << passes << " passes on the state of t = " << runTime.timeName()
        << nl << endl;
    Info<< "exact sphere: origin " << ev.origin << " radius " << ev.radius << endl;

    // the plicRDF reference scheme (its own interface fields; the registry keeps the
    // davof scheme's objects, as in the plicRDF normal model)
    autoPtr<reconstruction::plicRDF> plicPtr;
    if (args.found("plicRDF") || args.found("timing"))
    {
        dictionary pd(dict);
        pd.set("reconstructionScheme", word("plicRDF"));
        pd.set("positioning", word("vofCutCell"));
        plicPtr.reset(new reconstruction::plicRDF(alpha1, phi, U, pd));
    }

    // the timings first, on the state as read
    if (args.found("timing"))
    {
        const label R = args.get<label>("timing");
        using clock = std::chrono::steady_clock;
        auto since = [](const clock::time_point& t0)
        {
            return std::chrono::duration<double>(clock::now() - t0).count();
        };
        scalarList tA(R), tB(R), tC(R);
        const surfaceScalarField afSaved(dv.alphaf());
        plicPtr().reconstruct(true);   // the first call starts from the gradient normals: not timed
        for (label r = 0; r < R; ++r)
        {
            auto t0 = clock::now();
            plicPtr().reconstruct(true);
            tA[r] = since(t0);
            dv.alphaf() = afSaved;
            t0 = clock::now();
            dv.reconstruct(true);
            tB[r] = since(t0);
            dv.alphaf() = afSaved;
            t0 = clock::now();
            dv.correctFaceFractions();
            tC[r] = since(t0);
            dv.alphaf() = afSaved;
        }
        auto mean = [](const scalarList& v) { scalar s = 0; for (scalar x : v) s += x; return s/v.size(); };
        auto minv = [](const scalarList& v) { scalar m = GREAT; for (scalar x : v) m = min(m, x); return m; };
        Info<< "TIMING R " << R << " normalModel " << dv.normalModel().type()
            << " (a) plicRDF reconstruct mean " << mean(tA) << " min " << minv(tA)
            << " (b) davof reconstruct mean " << mean(tB) << " min " << minv(tB)
            << " (c) davof corrector pass mean " << mean(tC) << " min " << minv(tC)
            << " cells " << mesh.nCells() << endl;
        dv.alphaf() = afSaved;
    }

    // the plicRDF reference row: two calls (the second is the reference)
    if (args.found("plicRDF"))
    {
        for (label c = 1; c <= 2; ++c)
        {
            plicPtr().reconstruct(true);
            Evaluator evp;
            evp.func = implicitFunction::New("sphere", [&]{ dictionary fd; fd.add("type", word("sphere")); fd.add("origin", ev.origin); fd.add("radius", ev.radius); return fd; }());
            evp.origin = ev.origin; evp.radius = ev.radius;
            evp.line
            (
                "plicRDF", c, mesh, plicPtr().normal(), plicPtr().centre(),
                plicPtr().interfaceLabels().size(), -1, -1, -1, -1, -1, -1,
                dump && c == 2
            );
        }
    }

    // pass 0: the reconstruction on the state as read
    auto ringDiag = [&](label& negW, label& y0, scalar& maxWy)
    {
        negW = -1; y0 = -1; maxWy = -1;
        if (isA<davofNormalModels::ringLinear>(dv.normalModel()))
        {
            const davofNormalModels::ringLinear& rl =
                refCast<const davofNormalModels::ringLinear>(dv.normalModel());
            negW = rl.nNegWeight(); y0 = rl.nYZero(); maxWy = rl.maxWy();
        }
    };
    {
        dv.reconstruct(true);
        label negW, y0; scalar maxWy;
        ringDiag(negW, y0, maxWy);
        ev.line
        (
            (dump ? dumpTag : word(dv.normalModel().type())), 0, mesh, dv.normal(), dv.centre(),
            dv.interfaceLabels().size(), -1, -1, -1, negW, y0, maxWy, dump
        );
    }

    Info<< "pass  maxChange  L1Change  iterationsUsed  identity" << endl;
    for (label p = 1; p <= passes; ++p)
    {
        dv.correctFaceFractions();
        const scalarList& mx = dv.correctorMaxChange();
        const scalarList& l1 = dv.correctorL1Change();
        Info<< "PASS " << p << " " << (mx.size() ? mx[0] : -1)
            << " " << (l1.size() ? l1[0] : -1)
            << " " << dv.correctorIterationsUsed()
            << " " << dv.correctorIdentity() << endl;
        label negW, y0; scalar maxWy;
        ringDiag(negW, y0, maxWy);
        ev.line
        (
            (dump ? dumpTag : word(dv.normalModel().type())), p, mesh, dv.normal(), dv.centre(),
            dv.interfaceLabels().size(), (mx.size() ? mx[0] : -1),
            (l1.size() ? l1[0] : -1), dv.correctorIdentity(), negW, y0, maxWy, dump
        );
    }

    if (args.found("writeTime"))
    {
        const scalar tw = args.get<scalar>("writeTime");
        runTime.setTime(tw, 0);
        volScalarField aOut
        (
            IOobject("alpha.water", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            alpha1
        );
        aOut.write();
        surfaceScalarField afOut
        (
            IOobject(dv.alphaf().name(), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            dv.alphaf()
        );
        afOut.write();
        // U for the next call
        volVectorField uOut
        (
            IOobject("U", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            U
        );
        uOut.write();
        Info<< "written alpha.water, " << dv.alphaf().name() << " and U into the time directory "
            << runTime.timeName() << endl;
    }

    Info<< nl << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
