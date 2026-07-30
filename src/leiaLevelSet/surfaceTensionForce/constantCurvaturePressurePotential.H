/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::constantCurvaturePressurePotential

Description
    Constant-curvature pressure-compatibility oracle.  The capillary pressure

        pSigma = sigma*kappa*alpha

    is constructed at cell centres and differentiated with the same snGrad
    operator used by the pressure equation.  The model returns the integrated
    oriented scalar face flux snGrad(pSigma)*magSf, as required by the strict
    surfaceTensionForce contract.

    Comparing this model with constantCurvatureSurfaceTension separates a
    pressure-compatible discrete potential from the nominally equivalent
    sigma*kappa*snGrad(alpha)*magSf assembly.

SourceFiles
    constantCurvaturePressurePotential.C

\*---------------------------------------------------------------------------*/

#ifndef constantCurvaturePressurePotential_H
#define constantCurvaturePressurePotential_H

#include "IOdictionary.H"
#include "surfaceTensionForce.H"

namespace Foam
{

class constantCurvaturePressurePotential
:
    public surfaceTensionForce
{
protected:

    const fvSolution& fvSolutionDict_;
    const dictionary& levelSetDict_;
    const dictionary& surfTensionDict_;

    const dimensionedScalar curvature_;
    const volScalarField& alpha_;

public:

    TypeName("constantCurvaturePressurePotential");

    constantCurvaturePressurePotential(const fvMesh& mesh);

    virtual ~constantCurvaturePressurePotential() = default;

    virtual tmp<surfaceScalarField>
        calcFaceSurfaceTensionForceFlux() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
