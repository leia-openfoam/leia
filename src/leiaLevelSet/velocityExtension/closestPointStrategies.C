/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "closestPointStrategies.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(closestPointStep, 0);
    defineRunTimeSelectionTable(closestPointStep, dictionary);
    defineTypeNameAndDebug(closestPointFill, 0);
    defineRunTimeSelectionTable(closestPointFill, dictionary);

namespace closestPointSteps
{
    defineTypeNameAndDebug(distance, 0);
    addToRunTimeSelectionTable(closestPointStep, distance, dictionary);
    defineTypeNameAndDebug(newton, 0);
    addToRunTimeSelectionTable(closestPointStep, newton, dictionary);
}
namespace closestPointFills
{
    defineTypeNameAndDebug(steadyUpwind, 0);
    addToRunTimeSelectionTable(closestPointFill, steadyUpwind, dictionary);
    defineTypeNameAndDebug(footPointTaylor, 0);
    addToRunTimeSelectionTable(closestPointFill, footPointTaylor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::closestPointStep>
Foam::closestPointStep::New(const dictionary& dict)
{
    const word type = dict.getOrDefault<word>("cpStep", "distance");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "closestPointStep", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    Info<< "closestPoint: descent step " << type << endl;
    return autoPtr<closestPointStep>(ctorPtr(dict));
}


Foam::autoPtr<Foam::closestPointFill>
Foam::closestPointFill::New(const dictionary& dict)
{
    const word type = dict.getOrDefault<word>("cpFill", "steadyUpwind");

    auto* ctorPtr = dictionaryConstructorTable(type);
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict, "closestPointFill", type, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }
    Info<< "closestPoint: fill of failed descents " << type << endl;
    return autoPtr<closestPointFill>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::closestPointSteps::newton::newton(const dictionary& dict)
:
    closestPointStep(dict),
    qFloor_(dict.getOrDefault<scalar>("cpNewtonFloor", 1e-3))
{}


Foam::closestPointFills::footPointTaylor::footPointTaylor
(
    const dictionary& dict
)
:
    closestPointFill(dict),
    maxSweeps_(dict.getOrDefault<label>("cpFillMaxSweeps", 64))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace
{
    // Deterministic order of two candidate foot points at equal distance
    // (the same in serial and in parallel: coordinates only).
    bool lexLess(const Foam::point& a, const Foam::point& b)
    {
        if (a.x() != b.x()) return a.x() < b.x();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.z() < b.z();
    }
}


bool Foam::closestPointFills::footPointTaylor::fill
(
    const fvMesh& mesh,
    volVectorField& Uext,
    const boolList& fixed,
    const boolList& failed,
    const pointField& footP,
    const vectorField& footN,
    const tensorField& footG
) const
{
    const label nC = mesh.nCells();
    const vectorField& C = mesh.C().primitiveField();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const cellList& cells = mesh.cells();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // The foot-point data: p, n, u, G per cell; have = data present.
    boolList have(fixed);
    pointField P(footP);
    vectorField N(footN);
    vectorField Uf(Uext.primitiveField());
    tensorField G(footG);

    labelList todo;
    {
        DynamicList<label> t(nC/8);
        forAll(failed, c)
        {
            if (failed[c]) t.append(c);
        }
        todo.transfer(t);
    }

    // Exchange fields: zeroGradient on physical patches; coupled patches get
    // their constraint (processor) type and deliver the neighbour values.
    auto makeV = [&](const word& nm)
    {
        return volVectorField
        (
            IOobject(nm, fileName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, false),
            mesh, dimensionedVector(dimless, Zero), "zeroGradient"
        );
    };
    volScalarField vH
    (
        IOobject("cpFillHave", fileName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar(dimless, 0), "zeroGradient"
    );
    volVectorField vP(makeV("cpFillP"));
    volVectorField vN(makeV("cpFillN"));
    volVectorField vU(makeV("cpFillU"));
    volTensorField vG
    (
        IOobject("cpFillG", fileName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedTensor(dimless, Zero), "zeroGradient"
    );

    label nSweep = 0;
    label nFilledTot = 0;
    for (; nSweep < maxSweeps_; ++nSweep)
    {
        forAll(have, c)
        {
            vH[c] = have[c] ? 1.0 : 0.0;
        }
        vP.primitiveFieldRef() = P;
        vN.primitiveFieldRef() = N;
        vU.primitiveFieldRef() = Uf;
        vG.primitiveFieldRef() = G;
        vH.correctBoundaryConditions();
        vP.correctBoundaryConditions();
        vN.correctBoundaryConditions();
        vU.correctBoundaryConditions();
        vG.correctBoundaryConditions();

        // Neighbour values on the coupled patches (one evaluation per sweep)
        List<scalarField> nbH(pbm.size());
        List<vectorField> nbP(pbm.size()), nbN(pbm.size()), nbU(pbm.size());
        List<tensorField> nbG(pbm.size());
        forAll(pbm, pI)
        {
            if (mesh.boundary()[pI].coupled())
            {
                nbH[pI] = vH.boundaryField()[pI].patchNeighbourField();
                nbP[pI] = vP.boundaryField()[pI].patchNeighbourField();
                nbN[pI] = vN.boundaryField()[pI].patchNeighbourField();
                nbU[pI] = vU.boundaryField()[pI].patchNeighbourField();
                nbG[pI] = vG.boundaryField()[pI].patchNeighbourField();
            }
        }

        // Jacobi: read the previous sweep (have, P, ...), write new copies.
        boolList have1(have);
        pointField P1(P);
        vectorField N1(N), U1(Uf);
        tensorField G1(G);
        label nNew = 0;

        for (const label c : todo)
        {
            if (have[c]) continue;

            scalar best = GREAT;
            bool found = false;
            point bp(Zero);
            vector bn(Zero), bu(Zero);
            tensor bg(Zero);

            auto consider = [&]
            (
                const point& p, const vector& n, const vector& u, const tensor& g
            )
            {
                const scalar d = magSqr(C[c] - p);
                if (!found || d < best || (d == best && lexLess(p, bp)))
                {
                    found = true; best = d; bp = p; bn = n; bu = u; bg = g;
                }
            };

            for (const label f : cells[c])
            {
                if (mesh.isInternalFace(f))
                {
                    const label nb = (own[f] == c) ? nei[f] : own[f];
                    if (have[nb]) consider(P[nb], N[nb], Uf[nb], G[nb]);
                }
                else
                {
                    const label pI = pbm.whichPatch(f);
                    if (pI >= 0 && mesh.boundary()[pI].coupled())
                    {
                        const label i = f - pbm[pI].start();
                        if (nbH[pI][i] > 0.5)
                        {
                            consider(nbP[pI][i], nbN[pI][i], nbU[pI][i], nbG[pI][i]);
                        }
                    }
                }
            }

            if (found)
            {
                have1[c] = true; P1[c] = bp; N1[c] = bn; U1[c] = bu; G1[c] = bg;
                ++nNew;
            }
        }

        have.transfer(have1);
        P.transfer(P1); N.transfer(N1); Uf.transfer(U1); G.transfer(G1);

        reduce(nNew, sumOp<label>());
        nFilledTot += nNew;
        if (nNew == 0) break;
    }

    // Second-order Taylor value at the projected foot point
    label nUnreached = 0;
    vectorField& Ue = Uext.primitiveFieldRef();
    for (const label c : todo)
    {
        if (!have[c])
        {
            ++nUnreached;
            continue;   // keeps the base velocity
        }
        const vector& n = N[c];
        const vector r = C[c] - P[c];
        const point Pi = P[c] + (r - (r & n)*n);
        Ue[c] = Uf[c] + ((Pi - P[c]) & G[c]);
    }
    reduce(nUnreached, sumOp<label>());
    Uext.correctBoundaryConditions();

    Info<< "closestPoint footPointTaylor: " << nFilledTot << " cells filled in "
        << nSweep << " sweeps, " << nUnreached << " unreached" << endl;

    return true;
}

// ************************************************************************* //
