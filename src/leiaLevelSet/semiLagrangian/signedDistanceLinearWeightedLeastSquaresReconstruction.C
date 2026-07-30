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

#include "signedDistanceLinearWeightedLeastSquaresReconstruction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        signedDistanceLinearWeightedLeastSquaresReconstruction,
        0
    );
    addToRunTimeSelectionTable
    (
        slReconstruction,
        signedDistanceLinearWeightedLeastSquaresReconstruction,
        Mesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::signedDistanceLinearWeightedLeastSquaresReconstruction::
signedDistanceLinearWeightedLeastSquaresReconstruction
(
    const fvMesh& mesh
)
:
    linearWeightedLeastSquaresReconstruction(mesh),
    sdfGradMin_(slDict_.getOrDefault<scalar>("sdfGradMin", 0.1))
{
    Info<< "signedDistanceLinearWeightedLeastSquares: signed-distance "
        << "renormalisation psi = P/max(|g|, " << sdfGradMin_ << ")" << endl;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::signedDistanceLinearWeightedLeastSquaresReconstruction::evaluateRaw
(
    const label c,
    const point& x
) const
{
    const scalar psiC = (*psiOldPtr_)[c];
    const label nc = ncoeff_[c];
    if (nc == 0) { return psiC; }        // constant fallback: no gradient info

    const vector d = x - mesh_.C()[c];
    scalar b[3];
    basis(d, nc, b);
    const scalar* cf = &coeffsFlat_[c*ncoeffFull_];
    scalar val = psiC;
    scalar g2 = 0;
    for (label k = 0; k < nc; ++k)
    {
        val += cf[k]*b[k];
        g2 += cf[k]*cf[k];               // fitted gradient = the coefficients
    }
    const scalar gm = Foam::sqrt(g2);

    // Signed-distance renormalisation psi = P(x)/max(|g|, sdfGradMin) -- the
    // plane's local signed distance (zero set unchanged, |grad| pinned to 1) --
    // applied ONLY in INTERFACE-BAND donor cells (a psi sign change within the
    // stencil): there the plane is a valid single-interface model and P/|g| is
    // its true signed distance. Elsewhere the plain fit value is returned.
    // Why band-restricted (measured, two failed designs):
    //  - a GRADIENT trust window switches the division on/off where val is
    //    LARGE -> O(1/gMin) value jumps ringing the SDF skeleton -> coupled run
    //    dead in ~10 steps;
    //  - ALWAYS dividing amplifies exactly where the plane is an INVALID
    //    interface model -- under-resolved filaments / the skeleton, where the
    //    stencil spans both interface sides so |g| -> 0 while val != 0 ->
    //    |grad psi| ran to 8e14 on the stretched 2D vortex.
    // The band-edge switch is benign because it cuts where val ~ 0: the two
    // branches differ by O(|val||1 - g|/g) -> 0 near the zero set.
    const List<scalar>& s = stencilPsi_[c];
    bool band = false;
    forAll(s, i)
    {
        if (s[i]*psiC < 0) { band = true; break; }
    }
    if (band)
    {
        return val/max(gm, sdfGradMin_);
    }
    return val;
}

// ************************************************************************* //
