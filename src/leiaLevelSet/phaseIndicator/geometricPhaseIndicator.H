/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Tomislav Maric, TU Darmstadt
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
    Foam::geometricPhaseIndicator

Description
    Computes the phase indicator as a volume fraction from a field \f$\psi_c\f$. 
    
    In each cell \f$\Omega_c\f$, a piecewise-planer approximation \f$\psi^l_c\f$ is
    computed using linear least squares. This implicit plane is then used to
    truncate the cell and compute the phase-specific volume \f$\Omega_c^-\f$,
    resulting in the volume fraction 

    \f$\f$ \alpha_c = |\Omega_c^-| / |\Omega_c| \f$\f$

SourceFiles
    geometricPhaseIndicator.C

\*---------------------------------------------------------------------------*/

#ifndef geometricPhaseIndicator_H
#define geometricPhaseIndicator_H

#include "phaseIndicator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class geometricPhaseIndicator Declaration
\*---------------------------------------------------------------------------*/

class geometricPhaseIndicator
:
    public phaseIndicator
{
    
protected:

    // Narrow-band marker
    const volScalarField& narrowBand_;

    // Linear Least-Squares Approximation of \Psi
    // \psi(x,y,z)= nc_x x  + nc_y y + nc_z z + dc
    tmp<volVectorField> ncTmp_;
    tmp<volScalarField> dcTmp_;

public:

    // Static Data Members
    
    TypeName ("geometric");

    // Constructors

    //- Construct from components
    geometricPhaseIndicator(const fvMesh& mesh);

    //- Destructor
    virtual ~geometricPhaseIndicator() = default;

    // Member Functions
    virtual void calcPhaseIndicator
    (
        volScalarField& alpha, 
        const volScalarField& psi 
    ); 
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
