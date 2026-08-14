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

// Equation (24) of Fricke, Marić, Vučković, Roisman & Bothe, "A locally
// signed-distance preserving level set method (SDPLS) for moving interfaces",
// arXiv:2208.01269 (\cite fricke_locally_2022):
//
//              / 1                                        if 0 <= x <= w1,
//     G(x) =  |  exp( -ln(10^3) (x-w1)^2 / (w2-w1)^2 )     if w1 < x,
//              \ G(-x)                                     if x < 0,
//
// so G(w2) = 10^-3 exactly: at x = w2 the exponent evaluates to -ln(10^3).
// G reaches the transported field through eq. (23) of the same reference, where
// it multiplies the interfacial normal strain n.grad(v).n.
//
// The branch structure is the equation's and is kept literal so the code can be
// read against the paper. `s` is the Gaussian rate ln(10^3)/(w2-w1)^2; the
// constructor guarantees w2 > w1, so it is finite and positive.
double Foam::mollifier1::mollify(double x) const
{
    const double s = log(1000.0)/pow((w2_-w1_), 2);

    if (x >= 0) {
        if (x < w1_) {
            // Plateau: G == 1 on |x| <= w1. Being flat here (not merely
            // continuous) is what stops the cut-off from displacing the
            // discrete interface -- see the Description in mollifier1.H.
            return 1.0;
        }
        else {
            return exp(- s*pow(x-w1_, 2));
        }
    }
    else {
        // G is even; the negative side mirrors the positive one.
        return mollify(-x);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
