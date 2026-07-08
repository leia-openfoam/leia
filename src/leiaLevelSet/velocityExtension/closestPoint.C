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

#include "closestPoint.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "DynamicList.H"
#include "zoneDistribute.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(closestPoint, 0);
    addToRunTimeSelectionTable(velocityExtension, closestPoint, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::closestPoint::closestPoint(const fvMesh& mesh)
:
    steadyUpwind(mesh),
    cpTol_(velExtDict_.getOrDefault<scalar>("cpTol", 0.1)),
    cpMaxDist_(velExtDict_.getOrDefault<scalar>("cpMaxDist", 1.5)),
    cpHaloReach_(velExtDict_.getOrDefault<scalar>("cpHaloReach", 1.5))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::closestPoint::extend()
{
    interpolationCellPoint<scalar> psiInterp(psi_);
    interpolationCellPoint<vector> Uinterp(U_);
    const auto& C = mesh_.C();

    // Parallel finishing step (zoneDistribute; Scheufler & Roenby JCP 383
    // (2019), production-tested in geometricVoF/plicRDF): fetch ONE
    // cell-point-cell halo layer of remote (C, psi, grad psi, U, grad U)
    // around the band. A descent step that leaves the processor-local mesh
    // (findCell < 0) continues on first-order Taylor data anchored at the
    // nearest halo cell centre -- the same order as the local cell-point
    // interpolation -- instead of failing. Only excursions deeper than the
    // one-layer halo, physical-boundary overshoots and cyclic transforms
    // (whose untransformed centres misrepresent distances, so the reach cap
    // rejects them) still fall back to the pinned steadyUpwind fill.
    const bool haveHalo = Pstream::parRun();
    zoneDistribute* zd = nullptr;
    const labelListList* stencilPtr = nullptr;
    Map<vector> mapC;
    Map<scalar> mapPsi;
    Map<vector> mapGPsi;
    Map<vector> mapU;
    Map<tensor> mapGU;

    if (haveHalo)
    {
        boolList zone(mesh_.nCells(), false);
        forAll(band_, c)
        {
            zone[c] = band_[c] > 0.5;
        }

        zd = &zoneDistribute::New(mesh_);
        zd->setUpCommforZone(zone, true);

        // The same gradient flavours the family already uses: grad(psi)
        // resolves the fvSchemes 'grad(psi)' entry (as computeNormals);
        // grad(U) resolves the gradSchemes default.
        tmp<volVectorField> tgPsi = fvc::grad(psi_);
        tmp<volTensorField> tgU = fvc::grad(U_);

        mapC = zd->getDatafromOtherProc(zone, C);
        mapPsi = zd->getDatafromOtherProc(zone, psi_);
        mapGPsi = zd->getDatafromOtherProc(zone, tgPsi());
        mapU = zd->getDatafromOtherProc(zone, U_);
        mapGU = zd->getDatafromOtherProc(zone, tgU());
        stencilPtr = &zd->getStencil();
    }

    // Const views: Map::operator[] const FatalErrors on a missing key
    // instead of silently inserting.
    const Map<vector>& ccC = mapC;
    const Map<scalar>& ccPsi = mapPsi;
    const Map<vector>& ccGPsi = mapGPsi;
    const Map<vector>& ccU = mapU;
    const Map<tensor>& ccGU = mapGU;

    // Nearest REMOTE halo cell (global index) to point p, searched in the
    // stencil of the last local cell the walk visited. -1 if none lies
    // within cpHaloReach_*cellSize: one halo layer cannot represent points
    // beyond the immediate neighbour cells -- deeper walks must fall back.
    const auto nearestHalo = [&](const point& p, const label cLast) -> label
    {
        if (!haveHalo || cLast < 0)
        {
            return -1;
        }
        const labelList& st = (*stencilPtr)[cLast];
        const globalIndex& gn = zd->globalNumbering();
        label best = -1;
        scalar bestD = cpHaloReach_*cellSize_[cLast];
        for (const label gbl : st)
        {
            if (gn.isLocal(gbl))
            {
                continue;   // local containment already failed for p
            }
            const auto fnd = ccC.cfind(gbl);
            if (!fnd.good())
            {
                continue;
            }
            const scalar d = Foam::mag(p - *fnd);
            if (d < bestD)
            {
                bestD = d;
                best = gbl;
            }
        }
        return best;
    };

    // Every band cell gets the interface velocity at ITS closest interface
    // point: Newton descent along the LOCAL normal (re-evaluated at the cell
    // currently containing the iterate -- across the band the frozen starting
    // normal would be too crude), with acceptance guards. Successful cells
    // become Dirichlet data; failures (skeleton region, path cap, halo-reach
    // exceedance) are filled afterwards by the steadyUpwind solve pinned at
    // every successful cell.
    DynamicList<label> fixedCells(mesh_.nCells()/4);
    DynamicList<vector> fixedVals(mesh_.nCells()/4);
    label nFailed = 0;
    label nHalo = 0;

    forAll(band_, c)
    {
        if (band_[c] < 0.5)
        {
            continue;
        }

        const scalar hc = cellSize_[c];
        const scalar maxDist = cpMaxDist_*nLayers_*hc;
        // Taylor-vs-interpolation operator jumps at a processor crossing are
        // O(h^2); allow them without tripping the divergence guard.
        const scalar crossSlack = 1e-3*hc;

        point x = C[c];
        label cx = c;        // containing local cell (-1: halo/Taylor mode)
        label jg = -1;       // halo data cell (global index) when cx == -1
        label cLast = c;     // last local cell: stencil access + reach scale
        scalar apsi = Foam::mag(psi_[c]);
        bool ok = true;
        bool usedHalo = false;

        // Parallel: cap each step at ~one cell so the walk TRAVERSES the
        // cells on its way (cLast tracks to the true exit cell, whose halo
        // stencil covers the first remote layer at the crossing point). A
        // full Newton jump from an outer band cell would anchor the halo
        // search at the distant STARTING cell and miss. The capped walk
        // needs a larger iteration budget (~nLayers steps to reach the
        // interface + Newton polish). Serial keeps the exact original
        // full-jump iteration (bit-identical regression).
        const label nIter = haveHalo ? nDescent_ + nLayers_ + 2 : nDescent_;

        for (label it = 0; it < nIter; ++it)
        {
            scalar psix;
            vector nx;
            if (cx >= 0)
            {
                psix = psiInterp.interpolate(x, cx);
                nx = nHat_[cx];
            }
            else
            {
                const vector& Cj = ccC[jg];
                const vector& gj = ccGPsi[jg];
                psix = ccPsi[jg] + (gj & (x - Cj));
                nx = gj/Foam::max(Foam::mag(gj), SMALL);
            }

            vector dx = -psix*nx;
            if (haveHalo)
            {
                const scalar maxStep =
                    0.9*cellSize_[cLast >= 0 ? cLast : c];
                const scalar m = Foam::mag(dx);
                if (m > maxStep)
                {
                    dx *= maxStep/m;
                }
            }
            const point xNew = x + dx;

            if (Foam::mag(xNew - C[c]) > maxDist)
            {
                ok = false;   // walked too far (skeleton / bad normal)
                break;
            }

            // Re-localise: prefer the local mesh (a walk may re-enter it),
            // then the one-layer halo, then give up (deep excursion,
            // physical-boundary overshoot, cyclic transform).
            const label cn = mesh_.findCell(xNew);
            label jn = -1;
            if (cn < 0)
            {
                jn = nearestHalo(xNew, cLast);
                if (jn < 0)
                {
                    ok = false;
                    break;
                }
            }

            scalar apsiNew;
            if (cn >= 0)
            {
                apsiNew = Foam::mag(psiInterp.interpolate(xNew, cn));
            }
            else
            {
                apsiNew = Foam::mag
                (
                    ccPsi[jn] + (ccGPsi[jn] & (xNew - ccC[jn]))
                );
            }

            const bool pureLocal = (cx >= 0 && cn >= 0);
            if (apsiNew > apsi + (pureLocal ? 0.0 : crossSlack))
            {
                ok = false;   // diverging descent
                break;
            }

            x = xNew;
            cx = cn;
            jg = jn;
            apsi = apsiNew;
            if (cn >= 0)
            {
                cLast = cn;
            }
            else
            {
                usedHalo = true;
            }
        }

        if (ok && apsi < cpTol_*hc)
        {
            if (cx >= 0)
            {
                Uext_[c] = Uinterp.interpolate(x, cx);
            }
            else
            {
                // First-order Taylor of U at the foot point, anchored at the
                // halo data cell: U_j + (x - C_j) . (grad U)_j -- the same
                // order as the local linear (cell-point) interpolation.
                Uext_[c] = ccU[jg] + ((x - ccC[jg]) & ccGU[jg]);
            }
            if (usedHalo)
            {
                ++nHalo;
            }
            fixedCells.append(c);
            fixedVals.append(Uext_[c]);
        }
        else
        {
            ++nFailed;
        }
    }

    label nFixed = fixedCells.size();
    reduce(nFixed, sumOp<label>());
    reduce(nHalo, sumOp<label>());
    reduce(nFailed, sumOp<label>());
    Info<< "closestPoint: " << nFixed << " foot-pointed band cells ("
        << nHalo << " via halo), " << nFailed << " fallback cells" << endl;

    // Fill failures (and keep the far field consistent) with the steady
    // upwind transport pinned at ALL successful closest-point cells.
    if (nFailed > 0)
    {
        solveSteady(Uext_, fixedCells, fixedVals);
    }
}

// ************************************************************************* //
