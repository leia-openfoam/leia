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
    Foam::surfaceTensionForce

Description
    Abstract base class for surface tension force modeling.

SourceFiles
    surfaceTensionForce.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTensionForce_H
#define surfaceTensionForce_H

#include "messageStream.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvSolution.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class surfaceTensionForce Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionForce
{
    
protected:
    
    const fvMesh& mesh_;
    
    const Time& runTime_;
    
    // Surface tension coefficient.
    IOdictionary transportProperties_;
    dimensionedScalar sigma_; 

    //- Model implementation of the integrated face force flux.
    virtual tmp<surfaceScalarField>
    calcFaceSurfaceTensionForceFlux() const = 0;

    //- Return the registered face-native curvature selected by a model, or
    //  nullptr when the model should use its own curvature estimator.
    //
    //  The common dictionary contract is
    //
    //      faceCurvatureSource model;       // model | registered
    //      faceCurvatureField  kappaInterfaceFace;
    //
    //  "connectedInterface" is accepted as a descriptive synonym for
    //  "registered".  Curvature reconstruction/regularisation is deliberately
    //  outside the force models: once kappa_f exists on the CSF faces, Kang or
    //  arithmetic cell-to-face interpolation must not be applied again.
    const surfaceScalarField* registeredFaceCurvature
    (
        const dictionary& modelDict,
        const bool legacyConnectedSelection = false
    ) const;

    //- Assemble the pressure-compatible, integrated scalar CSF flux directly
    //  from face curvature and a model-specific localisation field.
    tmp<surfaceScalarField> integratedCSFFlux
    (
        const surfaceScalarField& kappaFace,
        const volScalarField& forceWeight,
        const word& resultName = "GSigma"
    ) const;

public:

    // Static Data Members

    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceTensionForce, 
        Mesh, 
        (
            const fvMesh& mesh 
        ), 
        (mesh)
    );

    // Constructors
    
        //- Construct from components
        surfaceTensionForce(const fvMesh& mesh);

    // Selectors

        //- Select default constructed
        static autoPtr<surfaceTensionForce> New(const fvMesh& mesh);

    //- Destructor
    virtual ~surfaceTensionForce() = default;

    // Member Functions
    
    //- Integrated scalar capillary force flux G_sigma,f [N/m].
    //
    //  This is the only production delivery contract. It already includes the
    //  face-area integration (for a pointwise face force density f_sigma,f,
    //  G_sigma,f = f_sigma,f |Sf|). The pressure corrector therefore combines
    //  it directly with snGrad(p_rgh)|Sf|, and both enter the same
    //  fvc::reconstruct path. Models cannot bypass that balanced face-flux
    //  space with a cell-centred momentum source.
    tmp<surfaceScalarField> faceSurfaceTensionForceFlux() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
