/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt
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

#include "mollifier1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mollifier1, false);
    addToRunTimeSelectionTable(mollifier, mollifier1, Dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mollifier1::mollifier1(const dictionary& dict)
    :
        mollifier(dict),
        w1_(dict.get<scalar>("w1")),
        w2_(dict.getOrDefault<scalar>("w2", 1e-3))
{
    // w2 is "the distance at which G = 1e-3" (eq. 24 of arXiv:2208.01269), so
    // it must lie OUTSIDE the flat region. The default of 1e-3 is smaller than
    // any sensible w1 -- cases/3Dshear set w1 = 0.05 and no w2 at all, which
    // silently gave a decay length |w2 - w1| = 0.049 that happened to be
    // reasonable and was in no way what the dictionary said. Refuse instead:
    // the mollifier width is a modelling choice and must be stated.
    if (w2_ <= w1_)
    {
        FatalErrorInFunction
            << "SDPLS mollifier1 requires w2 > w1, got w1 = " << w1_
            << ", w2 = " << w2_ << "." << nl
            << "w1 is the half-width of the region where the cut-off is unity;"
            << " w2 is the distance at which it has decayed to 1e-3." << nl
            << "Both are PHYSICAL lengths, not multiples of the cell size."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * *  Member functions  * * * * * * * * * * * * * * //

double Foam::mollifier1::mollify(double x) const
{
    double s = log(1000.0)/pow((w2_-w1_),2); // help variable  
    if (x >= 0) {
        if (x < w1_) {
            return 1.0;
        } 
        else {
            return exp(- s*pow(x-w1_, 2));
        }
    } 
    else {
        return mollify(-x);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
