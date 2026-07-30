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
    Foam::slScheme

Description
    Abstract base for the semi-Lagrangian advection SCHEME -- how psi^n is
    turned into psi^{n+1} once the reconstruction and correction strategy are
    fixed. Runtime-selected from
    `levelSet { semiLagrangian { scheme ...; } }` (default pointValue -- every
    existing case dictionary stays valid).

      - pointValue : psi^{n+1}(x_c) = psi^n(x_d), the Lagrange interpolation of
                     the reconstructed departure value (non-conservative).
      - fluxForm   : a conservative finite-volume semi-Lagrangian -- a
                     swept-face flux update that conserves the discrete integral
                     of psi for a divergence-free flux (CFL <= 1).

    The scheme is handed the reconstruction and the correction strategy owned by
    slAdvection, so both variants reuse the SAME per-cell quadratic fit.

SourceFiles
    slScheme.C

\*---------------------------------------------------------------------------*/

#ifndef slScheme_H
#define slScheme_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "slReconstruction.H"
#include "slCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class slScheme Declaration
\*---------------------------------------------------------------------------*/

class slScheme
{
protected:

    const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("slScheme");

    declareRunTimeSelectionTable
    (
        autoPtr,
        slScheme,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );

    // Constructors

        explicit slScheme(const fvMesh& mesh);

    // Selectors

        //- Select from levelSet.semiLagrangian.scheme (default pointValue).
        static autoPtr<slScheme> New(const fvMesh& mesh);

    //- Destructor
    virtual ~slScheme() = default;

    // Member Functions

        //- Advect one step: psi holds psi^n on entry, psi^{n+1} on exit.
        //  recon/corrector are owned by slAdvection (fluxForm ignores the
        //  corrector). Unew = u^{n+1}, Uold = u^n at cell centres.
        virtual void advance
        (
            volScalarField& psi,
            const volVectorField& Unew,
            const volVectorField& Uold,
            slReconstruction& recon,
            slCorrector& corrector
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
