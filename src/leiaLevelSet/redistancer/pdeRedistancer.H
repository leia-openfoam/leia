/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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
    Foam::pdeRedistancer

Description
    Level-Set Redistancer using the reinitialization equation of Sussman, 
    Smereka, and Osher.

    Reference:
    \verbatim
    Sussman, M., Smereka, P., & Osher, S. (1994). A Level Set Approach for 
    Computing Solutions to Incompressible Two-Phase Flow. Journal of 
    Computational Physics, 114(1), 146–159.
    \endverbatim


SourceFiles
    pdeRedistancer.C

\*---------------------------------------------------------------------------*/

#ifndef pdeRedistancer_H
#define pdeRedistancer_H

#include "redistancer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pdeRedistancer Declaration
\*---------------------------------------------------------------------------*/

class pdeRedistancer
:
    public redistancer
{
private:
    // Data Members

    scalar deltaT_;

    //- Mesh-relative pseudo-time step: deltaT = deltaTCoeff*h_min when
    //  set (> 0). The explicit reinitialization equation is CFL-limited by
    //  the MESH (|grad psi| ~ 1 characteristic speed): a fixed deltaT that
    //  is stable on one mesh diverges on a finer one -- measured on the
    //  static circle gate (deltaT 0.1 at h = 1/64: band error 2.7e-2 ->
    //  7.5e6). deltaTCoeff <= 0 (default) keeps the legacy fixed deltaT.
    scalar deltaTCoeff_;

    uint ninterations_;

    //- Frozen-band (bulk-only) variant: freeze the sign-change band +
    //  anchorLayers guard layers as Dirichlet data (their TRANSPORTED values
    //  -- nothing is reconstructed, so no chord-plane displacement) and run
    //  the pseudo-time fill on the bulk only. Zero interface displacement by
    //  construction. Default false = classical whole-field reinit (existing
    //  behaviour, bitwise).
    Switch freezeBand_;

    //- |grad psi| discretization of the reinit Hamiltonian:
    //    central (default) -- the fvSchemes grad(psi); MEASURED to oscillate
    //                         and diverge (bulk grad-err 0.56 -> 59.7 over
    //                         300 pseudo-steps; runaway flips bulk signs).
    //    godunov           -- Rouy-Tourin upwind: axis-wise max of the
    //                         one-sided differences toward the interface.
    //                         Monotone and causal (information flows from
    //                         the anchors outward; no sign flips) -- the
    //                         classical cure.
    const word hamiltonian_;

    //- Rouy-Tourin upwind |grad psi| of psi_p (see hamiltonian_).
    tmp<volScalarField> magGradGodunov
    (
        const volScalarField& psi_p,
        const volScalarField& signPsi
    ) const;

    bool write_;

    void write(const volScalarField& psi) const;

public:

    // Static Data Members

    TypeName ("PDE");

    // Constructors

    //- Construct from components
    pdeRedistancer(const fvMesh& mesh);

    //- Destructor
    virtual ~pdeRedistancer() = default;

    // Member Functions
    
    virtual bool doRedistance(volScalarField& psi);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
