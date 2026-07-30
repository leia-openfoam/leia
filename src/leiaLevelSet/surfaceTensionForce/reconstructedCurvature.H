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
    Foam::reconstructedCurvature

Description
    CSF surface-tension force driven either by the solver's registered cell
    curvature or by a face-native curvature service.  In the latter path the
    curvature is attached to connected zero-surface elements, extended directly
    to faces, and is not interpolated again.  Both paths return the integrated,
    owner-oriented scalar flux required by the pressure equation:

        G_sigma,f = sigma * kappa_f * snGrad(weight) * magSf

    The curvature estimator and the force-localisation weight remain separate
    choices; the force model owns only the pressure-compatible face assembly.

    Selected with `surfaceTensionForce { type reconstructedCurvature; }`; the
    field names are dict-overridable (`kappa`, `alpha`).

SourceFiles
    reconstructedCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructedCurvature_H
#define reconstructedCurvature_H

#include "surfaceTensionForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class reconstructedCurvature Declaration
\*---------------------------------------------------------------------------*/

class reconstructedCurvature
:
    public surfaceTensionForce
{
protected:

    const fvSolution& fvSolutionDict_;
    const dictionary& levelSetDict_;
    const dictionary& surfTensionDict_;

    //- Mean curvature filled by the solver from the reconstruction (1/length).
    const volScalarField& kappa_;

    //- Phase fraction whose snGrad localises the force to the interface.
    const volScalarField& alpha_;

    //- Level set (needed by the interfaceWeighted face interpolation).
    const volScalarField& psi_;

    //- Face interpolation/delivery of kappa: "arithmetic" (fvc::interpolate),
    //  classical level-set CSF) or "interfaceWeighted" (Kang/GFM inverse-
    //  distance interpolation TO the psi=0 crossing,
    //  kappa_f = (kappa_P |psi_N| + kappa_N |psi_P|)/(|psi_P|+|psi_N|),
    //  which cancels the parallel-contour normal variation of the level-set
    //  curvature at the faces where the balance lives; Kang et al. 2000,
    //  Francois et al. 2006 (SSF), Abadie et al. 2015: exact balance).
    const word faceInterpolation_;

    //- Force localisation weight: "alpha" (snGrad of the geometric volume
    //  fraction, 2-3 face support) or "sharpHeaviside" (snGrad of the sign-
    //  based Heaviside of psi -- force ONLY on sign-change faces, the faithful
    //  LS-SSF/GFM pairing where the Kang interpolation is exact-in-structure).
    const word forceWeight_;

public:

    // Static Data Members
    TypeName("reconstructedCurvature");

    // Constructors
    reconstructedCurvature(const fvMesh& mesh);

    //- Destructor
    virtual ~reconstructedCurvature() = default;

    // Member Functions
    virtual tmp<surfaceScalarField> calcFaceSurfaceTensionForceFlux() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
