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

Class
    Foam::deferredCorrector

Description
    Deferred-correction (fixed-point) semi-Lagrangian assignment: over the FIXED
    departure feet, iterate the reconstruction, rebuilding it from the CURRENT
    iterate each pass, so any lower-order gradient used in the lift is recomputed
    from the reconstructed field (self-consistent) rather than being a fixed
    external input:

        psi^{(0)} = psi^n
        psi^{(k+1)}(x_c) = (1-w) psi^{(k)}(x_c)
                           + w * reconstruct[psi^{(k)}]( x_d(c) )

    with under-relaxation w = defCorrRelax in (0,1] and a range-relative early-out
    defCorrTol. nDefCorr passes (hard-capped). Selected with
    `correction deferredCorrection;`.

    NOTE (finding on this level-set transport): this stabilises a correction ONLY
    when the reconstruction operator is contracting. For the value-inconsistent
    gradient lift (defectCorrectedIDW) the operator has spectral radius > 1, so the
    iteration DIVERGES (more passes -> faster blow-up); use `direct` with a
    value-fitting reconstruction (quadraticWeightedLeastSquares) there. Retained as
    the seam for developing a contracting deferred correction.

    Controls (semiLagrangian sub-dict): nDefCorr (default 3), nDefCorrMax (10),
    defCorrRelax (1.0), defCorrTol (0 = run all passes).

SourceFiles
    deferredCorrector.C

\*---------------------------------------------------------------------------*/

#ifndef deferredCorrector_H
#define deferredCorrector_H

#include "slCorrector.H"

namespace Foam
{

class deferredCorrector
:
    public slCorrector
{
    //- Deferred-correction passes (hard-capped by maxIters_).
    label nIters_;
    label maxIters_;

    //- Under-relaxation of the fixed-point update, in (0,1].
    scalar relax_;

    //- Range-relative early-out on the max nodal update (0 = run all passes).
    scalar tol_;

public:

    TypeName("deferredCorrection");

    deferredCorrector(const fvMesh& mesh, const dictionary& dict);

    virtual ~deferredCorrector() = default;

    virtual void correct
    (
        volScalarField& psi,
        const pointField& feet,
        slReconstruction& recon
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
