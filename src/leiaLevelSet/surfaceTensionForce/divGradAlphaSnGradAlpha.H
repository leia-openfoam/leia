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
    Foam::divGradAlphaSnGradAlpha

Description
    CSF surface-tension force with the curvature computed from the VOLUME
    FRACTION field itself (interFoam-style):

    \f[
        \kappa = -\nabla\cdot\frac{\nabla\alpha}{|\nabla\alpha| + \delta_N},
        \qquad
        f_\sigma = \sigma\,\kappa_f\,\nabla_f\alpha .
    \f]

    Because alpha comes from the geometric phase indicator, this kappa is a
    functional of the INTERFACE GEOMETRY ONLY: it cannot read the level-set
    normal profile, so the profile-mode feedback channel of the reconstruction
    curvature is closed by construction. Companion/falsifier model for the
    parasitic-current instability analysis (accuracy is interFoam-class, i.e.
    lower than the reconstruction curvature on a clean interface).

SourceFiles
    divGradAlphaSnGradAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef divGradAlphaSnGradAlpha_H
#define divGradAlphaSnGradAlpha_H

#include "surfaceTensionForce.H"
#include "fvSolution.H"
#include "volFieldsFwd.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class divGradAlphaSnGradAlpha Declaration
\*---------------------------------------------------------------------------*/

class divGradAlphaSnGradAlpha
:
    public surfaceTensionForce
{
    protected:

        const fvSolution& fvSolutionDict_;
        const dictionary& levelSetDict_;
        const dictionary& surfTensionDict_;

        // Volume fraction (geometric phase indicator)
        const volScalarField& alpha_;

        // Number of cell->face->cell averaging passes applied to alpha BEFORE
        // the curvature evaluation (smoothed-CSF, Lafaurie et al.). The
        // geometric indicator is a ~one-cell-sharp step, on which div(grad
        // alpha/|grad alpha|) picks up O(1/h) support-edge artifacts; smoothing
        // widens the transition to 2-3 cells like an algebraic-VoF alpha. The
        // FORCE still uses the sharp snGrad(alpha); only kappa sees the
        // smoothed field, which remains a functional of interface geometry.
        const label nSmooth_;

public:

    // Static Data Members
    TypeName ("divGradAlphaSnGradAlpha");

    // Constructors
    divGradAlphaSnGradAlpha(const fvMesh& mesh);

    //- Destructor
    virtual ~divGradAlphaSnGradAlpha() = default;

    // Member Functions
    virtual tmp<surfaceScalarField> calcFaceSurfaceTensionForceFlux() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
