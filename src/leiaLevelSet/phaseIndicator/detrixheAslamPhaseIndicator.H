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
    Foam::detrixheAslamPhaseIndicator

Description
    Computes the phase indicator (volume fraction) from the level-set field
    \f$\psi_c\f$ as a run-time-selectable alternative to
    geometricPhaseIndicator.

    It uses the *same* linear least-squares reconstruction of \f$\psi\f$ in
    each narrow-band cell (\c leastSquaresPlaneCoeffs), normalised to a
    signed-distance plane \f$\phi(\mathbf{x}) = \mathbf{n}\cdot\mathbf{x}+d\f$,
    \f$|\mathbf{n}|=1\f$. Instead of clipping the cell against the plane
    (which is sensitive to geometric tolerances and has known issues in 2D),
    the cell is decomposed into tetrahedra using the OpenFOAM cell centroid,
    the face centroids and the cell corner points, and the fraction of each
    tetrahedron occupied by \f$\Omega_c^- = \{\phi<0\}\f$ is evaluated
    analytically with the closed-form expressions of

        M. Detrixhe and T. D. Aslam, "From level set to volume of fluid and
        back again at second-order accuracy", Int. J. Numer. Meth. Fluids
        2016; 80:231-255 (Eq. 7).

    The cell volume fraction is the volume-weighted average of the
    per-tetrahedron fractions. This decomposition is identical in 3D and
    pseudo-2D (one-cell-thick) meshes, so it is robust for 2D cases.

SourceFiles
    detrixheAslamPhaseIndicator.C

\*---------------------------------------------------------------------------*/

#ifndef detrixheAslamPhaseIndicator_H
#define detrixheAslamPhaseIndicator_H

#include "phaseIndicator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class detrixheAslamPhaseIndicator Declaration
\*---------------------------------------------------------------------------*/

class detrixheAslamPhaseIndicator
:
    public phaseIndicator
{

protected:

    //- Narrow-band marker
    const volScalarField& narrowBand_;

public:

    // Static Data Members

    TypeName ("detrixheAslam");

    // Constructors

    //- Construct from components
    detrixheAslamPhaseIndicator(const fvMesh& mesh);

    //- Destructor
    virtual ~detrixheAslamPhaseIndicator() = default;

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
