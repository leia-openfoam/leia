/*---------------------------------------------------------------------------*\
    leia :: functionObjects :: contactLineError (see contactLineError.H)
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, GPL-3.0-or-later.

    M. Fricke, T. Maric, D. Bothe, Contact line advection using the
    geometrical Volume-of-Fluid method, Journal of Computational Physics 407
    (2020) 109221, doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.
\*---------------------------------------------------------------------------*/

#include "contactLineError.H"
#include "davofAdvectionScheme.H"
#include "reconstructionSchemes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "globalIndex.H"
#include <cmath>

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(contactLineError, 0);
    addToRunTimeSelectionTable(functionObject, contactLineError, dictionary);
}
}

void Foam::functionObjects::contactLineError::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "contactLineError");
    writeCommented(os, "tRec");
    writeTabbed(os, "nCandL");
    writeTabbed(os, "xL");
    writeTabbed(os, "thetaL");
    writeTabbed(os, "cellL");
    writeTabbed(os, "nCandR");
    writeTabbed(os, "xR");
    writeTabbed(os, "thetaR");
    writeTabbed(os, "cellR");
    writeTabbed(os, "volume");
    writeTabbed(os, "L1toInitial");
    writeTabbed(os, "minAlpha");
    writeTabbed(os, "maxAlphaMinus1");
    os << endl;
}

Foam::functionObjects::contactLineError::contactLineError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    fieldName_(),
    wallPatch_(),
    alpha0_(),
    haveAlpha0_(false)
{
    nMissing_[0] = 0;
    nMissing_[1] = 0;
    read(dict);
    writeFileHeader(file());
}

bool Foam::functionObjects::contactLineError::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    fieldName_ = dict.get<word>("field");
    wallPatch_ = dict.get<word>("wallPatch");
    return true;
}

bool Foam::functionObjects::contactLineError::execute()
{
    const volScalarField& alpha = lookupObject<volScalarField>(fieldName_);
    const scalarField& a = alpha.primitiveField();

    if (!haveAlpha0_)
    {
        // The first call is at the start time, before the first step
        alpha0_ = a;
        haveAlpha0_ = true;
        return true;
    }

    davofAdvectionScheme& adv =
        mesh_.lookupObjectRef<davofAdvectionScheme>("davofAdvectionScheme");
    reconstructionSchemes& surf = adv.surf();

    const label patchi = mesh_.boundaryMesh().findPatchID(wallPatch_);
    if (patchi < 0)
    {
        FatalErrorInFunction
            << "wall patch " << wallPatch_ << " not found" << exit(FatalError);
    }
    const polyPatch& wall = mesh_.boundaryMesh()[patchi];

    // The interface cells of the reconstruction
    bitSet isInterface(mesh_.nCells());
    for (const label celli : surf.interfaceLabels())
    {
        isInterface.set(celli);
    }

    const vectorField& n = surf.normal();
    const vectorField& c = surf.centre();
    const pointField& pts = mesh_.points();

    // Per side: nCandidates, x, theta, global cell (as scalars)
    // side 0: left (n_Sigma,x < 0), side 1: right
    scalar cand[2] = {0, 0};
    scalar xs[2] = {0, 0};
    scalar ths[2] = {0, 0};
    scalar cells[2] = {-1, -1};

    const labelUList& own = mesh_.faceOwner();
    forAll(wall, pfi)
    {
        const label facei = wall.start() + pfi;
        const label celli = own[facei];
        if (!isInterface.test(celli)) continue;
        const scalar magn = mag(n[celli]);
        if (magn == 0) continue;

        const vector nOut(-n[celli]/magn);
        if (nOut.x() == 0) continue;

        // The wall face: x-range and centre
        const face& f = mesh_.faces()[facei];
        scalar xlo = GREAT, xhi = -GREAT;
        for (const label pi : f)
        {
            xlo = min(xlo, pts[pi].x());
            xhi = max(xhi, pts[pi].x());
        }
        const point& xf = mesh_.faceCentres()[facei];

        // The plane n.(p - c) = 0 on the wall plane z = 0, y = y_f
        const scalar x =
            c[celli].x()
          + (n[celli].z()*c[celli].z() - n[celli].y()*(xf.y() - c[celli].y()))
           /n[celli].x();
        if (x < xlo || x > xhi) continue;

        const scalar theta =
            std::acos(min(1.0, max(-1.0, nOut.z())));
        const label side = (nOut.x() < 0) ? 0 : 1;

        const bool better =
            (cand[side] == 0)
         || (side == 0 && x < xs[side])
         || (side == 1 && x > xs[side]);
        cand[side] += 1;
        if (better)
        {
            xs[side] = x;
            ths[side] = theta;
            cells[side] = scalar(celli);
        }
    }

    // Gather over the ranks: the outermost candidate per side
    const globalIndex gi(mesh_.nCells());
    List<List<scalar>> all(Pstream::nProcs());
    {
        List<scalar> mine(8);
        for (int s = 0; s < 2; ++s)
        {
            mine[4*s + 0] = cand[s];
            mine[4*s + 1] = xs[s];
            mine[4*s + 2] = ths[s];
            mine[4*s + 3] = (cells[s] >= 0) ? scalar(gi.toGlobal(label(cells[s]))) : -1;
        }
        all[Pstream::myProcNo()] = mine;
    }
    Pstream::gatherList(all);
    Pstream::scatterList(all);

    scalar gCand[2] = {0, 0};
    scalar gx[2] = {0, 0};
    scalar gth[2] = {0, 0};
    scalar gcell[2] = {-1, -1};
    for (int s = 0; s < 2; ++s)
    {
        bool have = false;
        forAll(all, proci)
        {
            const List<scalar>& m = all[proci];
            gCand[s] += m[4*s];
            if (m[4*s] > 0)
            {
                const bool better =
                    !have
                 || (s == 0 && m[4*s + 1] < gx[s])
                 || (s == 1 && m[4*s + 1] > gx[s]);
                if (better)
                {
                    have = true;
                    gx[s] = m[4*s + 1];
                    gth[s] = m[4*s + 2];
                    gcell[s] = m[4*s + 3];
                }
            }
        }
        if (!have) ++nMissing_[s];
    }

    // The liquid volume, the L1 distance from the initial alpha, extremes
    const scalarField& V = mesh_.V();
    scalar vol = 0, l1 = 0, amin = GREAT, amax = -GREAT;
    forAll(a, celli)
    {
        vol += a[celli]*V[celli];
        l1 += mag(a[celli] - alpha0_[celli])*V[celli];
        amin = min(amin, a[celli]);
        amax = max(amax, a[celli]);
    }
    reduce(vol, sumOp<scalar>());
    reduce(l1, sumOp<scalar>());
    reduce(amin, minOp<scalar>());
    reduce(amax, maxOp<scalar>());

    if (Pstream::master())
    {
        const scalar tRec = time_.value() - time_.deltaTValue();
        Ostream& os = file();
        os.precision(17);
        os << tRec;
        for (int s = 0; s < 2; ++s)
        {
            os << tab << gCand[s] << tab << gx[s] << tab << gth[s] << tab
               << gcell[s];
        }
        os << tab << vol << tab << l1 << tab << amin << tab << amax - 1
           << endl;
    }

    return true;
}

bool Foam::functionObjects::contactLineError::write()
{
    return true;
}
