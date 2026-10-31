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

Description
    The pointwise distance estimates of D. Bothe (2026-10), shared by the
    rescaling redistancer and the projection-blend velocity extension:

      order 1  r    = psi/q,  q = |grad psi|
      order 2  dhat = 2 psi/(q + sqrt((q^2 - 2 psi d_nn psi)^+)),
               d_nn psi = n.H.n, n = grad psi/q, H the Hessian of psi

    The factor is positive, so the sign of psi is kept. The caller supplies
    grad psi and, for order 2, d_nn psi; a cell with q below gradFloor gets
    no estimate (the function returns false).

\*---------------------------------------------------------------------------*/

#ifndef rescaledDistance_H
#define rescaledDistance_H

#include "scalar.H"
#include "word.H"

namespace Foam
{
namespace rescaledDistance
{

//- The positive part: max(s, 0), or the smooth (s + sqrt(s^2 + eps^2))/2
inline scalar positivePart(const scalar s, const bool smooth, const scalar eps)
{
    return smooth ? 0.5*(s + Foam::sqrt(s*s + eps*eps)) : Foam::max(s, scalar(0));
}

//- The distance estimate of psi; false if q < gradFloor (dist unchanged).
//  negDisc is set when the discriminant of order 2 is negative.
inline bool estimate
(
    const label order,
    const scalar psi,
    const scalar q,
    const scalar dnn,
    const scalar gradFloor,
    const bool smooth,
    const scalar smoothEps,
    scalar& dist,
    bool& negDisc
)
{
    negDisc = false;
    if (q < gradFloor)
    {
        return false;
    }
    if (order == 1)
    {
        dist = psi/q;
        return true;
    }
    const scalar disc = q*q - 2*psi*dnn;
    negDisc = (disc < 0);
    dist = 2*psi/(q + Foam::sqrt(positivePart(disc, smooth, smoothEps)));
    return true;
}

} // End namespace rescaledDistance
} // End namespace Foam

#endif

// ************************************************************************* //
