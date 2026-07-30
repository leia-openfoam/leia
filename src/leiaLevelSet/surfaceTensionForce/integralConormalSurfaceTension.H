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
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    Foam::integralConormalSurfaceTension

Description
    Curvature-free, conservative pseudo-2D prototype. A single vector traction
    is constructed on each internal mesh face cut by the reconstructed
    interface. Its finite-volume surface integral first supplies a conservative
    cell force density. The model then linearly interpolates that density to a
    face and returns its integrated normal force flux

        G_sigma,f = interpolate(f_sigma,c) . Sf.

    There is deliberately no direct cell-source delivery: the base-class API
    only permits the balanced face-flux path used by the pressure equation.

    This is deliberately a separate model from integralSurfaceTension: it is
    the direct surface-stress/conormal identity, not the structured staggered
    CST pressure-jump algorithm. Coupled/boundary cuts and 3D intersection
    polygons remain future work.

\*---------------------------------------------------------------------------*/

#ifndef integralConormalSurfaceTension_H
#define integralConormalSurfaceTension_H

#include "surfaceTensionForce.H"

namespace Foam
{

class integralConormalSurfaceTension
:
    public surfaceTensionForce
{
    const fvSolution& fvSolutionDict_;
    const dictionary& levelSetDict_;
    const dictionary& surfTensionDict_;
    const volScalarField& field_;
    const scalar isoValue_;
    const scalar sgn_;

    //- Conservative cell force produced by the conormal traction identity,
    //  independent of how that force is delivered to PIMPLE.
    tmp<volVectorField> conservativeCellSurfaceTensionForce() const;

public:

    TypeName("integralConormalSurfaceTension");

    integralConormalSurfaceTension(const fvMesh& mesh);
    virtual ~integralConormalSurfaceTension() = default;

    //- Integrated force flux linearInterpolate(fCell) & Sf [N/m].
    virtual tmp<surfaceScalarField> calcFaceSurfaceTensionForceFlux() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
