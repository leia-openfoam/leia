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

#include "distanceNarrowBand.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

namespace Foam
{
    defineTypeNameAndDebug(distanceNarrowBand, false);
    addToRunTimeSelectionTable(narrowBand, distanceNarrowBand, Dictionary);

    distanceNarrowBand::distanceNarrowBand(const dictionary& dict, const volScalarField& psi)
        :
            emptyNarrowBand(dict, psi),
            ncells_(dict.getOrDefault<scalar>("ncells", 5))
    {
        // Info << "NarrowBand type: " << TypeNameString << " with ncells: " << ncells_ << nl << endl;
    }

    void distanceNarrowBand::calc()
    {
        field() = dimensionedScalar(field().dimensions(), 0.);

        // gMin is a global reduction, so deltaX_min is already identical on
        // every rank -- the band width does not depend on the decomposition.
        const auto deltaX = pow(mesh().deltaCoeffs(),-1).cref();
        const auto deltaX_min = gMin(deltaX);

        // mag(). Without it the test reads `psi/dx < ncells`, which every cell
        // with psi < 0 satisfies -- and psi < 0 is the ENTIRE interior phase.
        // The "narrow band" was therefore the whole drop plus a one-sided shell
        // of ncells on the outside: not narrow, not a band, and not symmetric
        // about the interface. This is a purely local criterion, so it was
        // wrong identically on every rank, which is why it never showed up as a
        // parallel inconsistency.
        //
        // Note the criterion still presupposes that |psi| is a DISTANCE, i.e.
        // that psi is a signed distance function. Where |grad psi| leaves 1 the
        // marked width is not ncells cells -- see interfaceBandMollifier for a
        // band that makes no such assumption.
        const volScalarField& psiRef = psi();
        forAll(field(), cellID)
        {
            if (mag(psiRef[cellID])/deltaX_min < ncells_)
            {
                field()[cellID] = 1.0;
            }
        }
    }

} // End namespace Foam

// ************************************************************************* //
