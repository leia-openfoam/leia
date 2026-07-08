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
    Foam::steadyUpwind

Description
    Velocity extension by the STEADY normal-transport equation

        sign(psi) (n . grad) Uext = 0,

    solved in ONE implicit linear solve per defect-correction pass (no
    pseudo-time truncation): the characteristics of w = sign(psi) n all
    emanate from the interface seed layer, so the upwind convection matrix is
    an M-matrix permutable to block-lower-triangular form -- a single solve
    has infinite reach (the implicit twin of fast marching / fast sweeping;
    Zhao et al. JCP 127 (1996), Peng et al. JCP 155 (1999), Aslam, Luo & Zhao
    SISC 36 (2014)).

    Assembled exactly like the solver's psi-equation (advective form of a
    non-solenoidal transport):

        fvm::div(phiW, Uext) - fvm::Sp(fvc::div(phiW), Uext)
        + fvm::Sp(lambda, Uext) == 0

    with the divergence scheme taken from fvSchemes (div(velExtPhiW,Uext)),
    case-controlled via the UEXT_DIV token. MEASURED: keep it Gauss upwind --
    the linearUpwind deferred correction (converged by the `nDefCorrExt`
    defect-correction loop; the implicit part remains the triangularizable
    upwind matrix) DIVERGES on this steady equation (no fvm::ddt damping).
    lambda is a tiny zero-inflow guard for discrete-skeleton cells (local
    maxima of |psi| with no upwind neighbour). Seeds are pinned via
    fvMatrix::setValues.

SourceFiles
    steadyUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef steadyUpwind_H
#define steadyUpwind_H

#include "interfaceExtension.H"

namespace Foam
{

class steadyUpwind
:
    public interfaceExtension
{
protected:

    //- Defect-correction passes converging the deferred (higher-order minus
    //  upwind) part of the divergence scheme. 1 = pure implicit (upwind).
    label nDefCorrExt_;

    //- Zero-inflow (discrete skeleton) regularization coefficient; the guard
    //  is lambda = zeroInflowGuard_/cellSize [1/length].
    scalar zeroInflowGuard_;

    //- Under-relaxation of the deferred correction between passes (the steady
    //  fixed-point iteration rings without damping). 1 = no relaxation.
    scalar defCorrRelax_;

    //- Solve the steady extension for the given field with the given seed
    //  Dirichlet data and an optional explicit source (used by the Aslam
    //  extrapolation cascade). Shared by steadyUpwind and derived models.
    void solveSteady
    (
        volVectorField& q,
        const labelUList& fixedCells,
        const UList<vector>& fixedVals,
        const volVectorField* source = nullptr
    );

    virtual void extend() override;

public:

    TypeName("steadyUpwind");

    steadyUpwind(const fvMesh& mesh);

    virtual ~steadyUpwind() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
