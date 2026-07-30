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
    Foam::pointValueScheme

Description
    Point-value semi-Lagrangian: psi^{n+1}(x_c) = psi^n(x_d). The Taylor
    backward-foot x_d is computed once, then the correction strategy assembles
    psi^{n+1} from the reconstruction evaluated at the feet. Non-conservative
    (a Lagrange interpolation of the departure value). The default scheme.

    The velocity inserted into the Taylor displacement is selectable as the
    supplied input field, the local normal projection (U.n)n, or the normal
    speed sampled at an on-demand closest point of the interpolated psi=0
    interface. Both time levels and grad(projected U) are formed consistently.

SourceFiles
    pointValueScheme.C

\*---------------------------------------------------------------------------*/

#ifndef pointValueScheme_H
#define pointValueScheme_H

#include "slScheme.H"

namespace Foam
{

class pointValueScheme
:
    public slScheme
{
    //- Velocity inserted into the Taylor departure-foot formula:
    //  input              -- use the fields passed by slAdvection unchanged;
    //  normalProjection   -- (U.n)n at each cell and time level;
    //  normalClosestPoint -- [U(piGamma).nGamma]nGamma, evaluated by an
    //                        on-demand closest-point walk before Taylor.
    word trajectoryVelocity_;

    //- Registry name of the physical velocity used by the two direct normal
    //  modes. This bypasses a separately constructed Eulerian Uext field.
    word trajectoryUName_;

    //- Smooth normal regularisation and closest-point guards.
    scalar normalEpsilon_;
    label projectionIterations_;
    scalar projectionTolerance_;
    scalar projectionMaxDistance_;

    //- Construct the normal trajectory velocity for one psi/U time level.
    autoPtr<volVectorField> normalVelocity
    (
        const word& name,
        const volScalarField& psiLevel,
        const volVectorField& ULevel,
        const bool closestPoint
    ) const;

public:

    TypeName("pointValue");

    explicit pointValueScheme(const fvMesh& mesh);

    virtual ~pointValueScheme() = default;

    virtual void advance
    (
        volScalarField& psi,
        const volVectorField& Unew,
        const volVectorField& Uold,
        slReconstruction& recon,
        slCorrector& corrector
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
