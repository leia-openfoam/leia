/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Tomislav Maric & Julian Reitzel, TU Darmstadt
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

#include "neighbourNarrowBand.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvPatch.H"

namespace Foam
{
    defineTypeNameAndDebug(neighbourNarrowBand, false);
    addToRunTimeSelectionTable(narrowBand, neighbourNarrowBand, Dictionary);

    neighbourNarrowBand::neighbourNarrowBand(const dictionary& dict, const volScalarField& psi)
        :
            signChangeNarrowBand(dict, psi),
            ntimes_(dict.get<int>("n"))
    {}

    void neighbourNarrowBand::calc()
    {
        signChangeNarrowBand::calc();
        if (field().empty())
        {
            FatalErrorInFunction
            << "Initial NarrowBandField is empty! Can not find neighbours." << exit(FatalError);
        }

        // ONE LAYER AT A TIME, AND EACH LAYER IS FINISHED BEFORE THE NEXT BEGINS.
        //
        // The previous implementation dilated ntimes_ using only
        // mesh().cellCells(), which is LOCAL topology. A cell one layer across a
        // processor boundary from a marked cell was never reached, so the band
        // was torn at every processor boundary and the result depended on the
        // DECOMPOSITION -- the same case on 4 and on 16 ranks marked different
        // cells. signChangeNarrowBand already exchanges its SEED over coupled
        // patches; the dilation on top of it did not.
        //
        // Each pass is therefore: (1) mark locally from the current marker,
        // (2) exchange the marker over coupled patches and mark from the
        // neighbouring side, (3) commit. Doing the exchange once at the end
        // would not do: layer k+1 on this rank may only be reachable from a
        // cell that layer k put on the OTHER side of the boundary, so the
        // information has to cross once per layer.
        const labelListList& neighbourCells = mesh().cellCells();
        const auto& faceOwner = mesh().faceOwner();

        for (int layer = 0; layer < ntimes_; ++layer)
        {
            volScalarField& NarrowBand = field();

            // (1) LOCAL: mark neighbours of everything currently in the band.
            // Collected first and applied after, so a cell marked during this
            // pass cannot seed further growth within the same pass -- otherwise
            // one pass would advance several layers on the local side while the
            // coupled side advances exactly one, which is the same
            // decomposition dependence in a subtler form.
            labelList toMark;
            forAll(NarrowBand, cellI)
            {
                if (NarrowBand[cellI] <= 0.0) continue;
                const labelList& nbrs = neighbourCells[cellI];
                forAll(nbrs, j)
                {
                    if (NarrowBand[nbrs[j]] == 0.0) toMark.append(nbrs[j]);
                }
            }

            // (2) COUPLED: the marker on the other side of a processor patch.
            // Evaluating the boundary field performs the halo exchange; a cell
            // whose coupled neighbour is in the band joins this layer.
            NarrowBand.correctBoundaryConditions();
            const auto& nbBdry = NarrowBand.boundaryField();
            const auto& patches = mesh().boundary();
            forAll(nbBdry, patchI)
            {
                const fvPatch& patch = patches[patchI];
                if (!isA<coupledFvPatch>(patch)) continue;

                const auto nbNeiTmp = nbBdry[patchI].patchNeighbourField();
                const auto& nbNei = nbNeiTmp();
                forAll(nbNei, faceI)
                {
                    const label faceJ = faceI + patch.start();
                    const label ownCell = faceOwner[faceJ];
                    if (nbNei[faceI] > 0.0 && NarrowBand[ownCell] == 0.0)
                    {
                        toMark.append(ownCell);
                    }
                }
            }

            // (3) COMMIT this layer.
            forAll(toMark, i)
            {
                NarrowBand[toMark[i]] = 1.0;
            }
            NarrowBand.correctBoundaryConditions();
        }
    }

} // End namespace Foam

// ************************************************************************* //
