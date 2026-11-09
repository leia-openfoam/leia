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

Usage
    davofTestCorrector [-normalModel exact|ringArea] [-passes N] [-time T]
    in a rendered case of the DAVOF advection studies (system/fvSolution with
    the alpha.water solver dictionary, 0/alpha.water, 0/alphaf.water, 0/U).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "reconstructionSchemes.H"
#include "davof.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption("normalModel", "word", "the normal model (default: the dictionary's)");
    argList::addOption("passes", "label", "the corrector passes (default 6)");
    argList::addOption("time", "scalar", "the time directory to read (default: the start time)");
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

    Info<< nl << "davofTestCorrector: normalModel " << dv.normalModel().type()
        << ", " << passes << " passes on the state of t = " << runTime.timeName()
        << nl << endl;
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
    }

    Info<< nl << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
