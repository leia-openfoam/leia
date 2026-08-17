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
    Foam::curvaturePressurePotential

Description
    Capillary force as the discrete gradient of a capillary PRESSURE POTENTIAL
    built from the COMPUTED curvature field:

        G_sigma,f = snGrad(sigma * kappa * alpha)_f * |S_f|

    This is the "remainder-off" arm of the identity (user's decomposition,
    2026-08-17):

        sigma*kappa*grad(alpha) = grad(sigma*kappa*alpha) - sigma*alpha*grad(kappa)

    Keeping only the first term makes the delivered force an EXACT discrete
    gradient for ANY kappa field, so the pressure equation reproduces it
    identically (p_rgh = sigma*kappa*alpha solves it, and rAUf cancels between
    phig and p_rghEqn.flux()) and the velocity contribution is zero to solver
    tolerance. It therefore CANNOT relax a deformed interface -- the discarded
    term is the entire shape-relaxation driver, non-gradient exactly where
    grad(alpha) x grad(kappa) != 0. It is a DIAGNOSTIC, not a production model.

    WHAT IT SEPARATES, that constantCurvature{SurfaceTension,PressurePotential}
    cannot. Those two use a CONSTANT curvature, for which
    snGrad(sigma*kappa*alpha) == sigma*kappa*snGrad(alpha) identically -- they
    are a verification pair for operator co-location, and they pass
    (Courant max 1.2e-10, continuity errors ~1e-12 on the coupled stationary
    droplet). This class uses the VARYING computed kappa, so it answers the
    separate question: is the non-absorbable capillary flux measured for
    production (1.7e-5 of |phig| on the healthy plateau, at a pressure-solver
    residual of 1.5e-10, i.e. genuinely non-gradient) a property of kappa
    VARYING at all, or specifically of the sigma*kappa_f*snGrad(alpha)
    ASSEMBLY? If this arm drives that fraction to solver tolerance with the
    same kappa field, it is the assembly.

    Reads the same dictionary entries as reconstructedCurvature, so an existing
    case selects it by changing only `type`:
        levelSet { surfaceTensionForce { type curvaturePressurePotential;
                                         alpha alpha.water; kappa kappa; } }

SourceFiles
    curvaturePressurePotential.C

\*---------------------------------------------------------------------------*/

#ifndef curvaturePressurePotential_H
#define curvaturePressurePotential_H

#include "IOdictionary.H"
#include "surfaceTensionForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class curvaturePressurePotential Declaration
\*---------------------------------------------------------------------------*/

class curvaturePressurePotential
:
    public surfaceTensionForce
{
protected:

    // Protected data

        const fvSolution& fvSolutionDict_;
        const dictionary& levelSetDict_;
        const dictionary& surfTensionDict_;

        //- The COMPUTED curvature field (cell-centred).
        const volScalarField& curvature_;

        //- The phase indicator.
        const volScalarField& alpha_;

public:

    TypeName("curvaturePressurePotential");

    // Constructors

        explicit curvaturePressurePotential(const fvMesh& mesh);

    //- Destructor
    virtual ~curvaturePressurePotential() = default;

    // Member Functions

        //- snGrad(sigma*kappa*alpha)*|S_f|: an exact discrete gradient.
        virtual tmp<surfaceScalarField>
        calcFaceSurfaceTensionForceFlux() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
