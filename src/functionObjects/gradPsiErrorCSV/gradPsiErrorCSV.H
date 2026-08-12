/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Julian Reitzel, TU Darmstadt 
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
    Foam::functionObjects::gradPsiErrorCSV

Group
    grpFieldFunctionObjects

Description
    Looks up a signedDistance-error field (from gradPsiError functionObject for example).
    Calculates further properties and writes them to a CSV file "gradPsiError.csv"
    The following properties a written in every row.
        "TIME,"
        "MAX_DELTA_X,"
        "MEAN_DELTA_X,"
        "E_MAX_GRAD_PSI,"
        "E_MEAN_GRAD_PSI,"
        "E_NARROW_MAX_GRAD_PSI,"
        "E_NARROW_MEAN_GRAD_PSI, ..."

    Also writes volume-weighted L1 and L2 (RMS) norms of the signed-distance
    error e = ||grad psi| - 1|, globally and over the narrow band:
        L1 = sum(e V)/sum(V),  L2 = sqrt(sum(e^2 V)/sum(V))
    in columns E_L1_GRAD_PSI, E_L2_GRAD_PSI, E_NARROW_L1_GRAD_PSI,
    E_NARROW_L2_GRAD_PSI.

    When an SDPLS source is active it registers the interfacial normal strain
    a = n . grad(U) . n as the field \c R; its narrow-band min/mean/max are then
    written to NARROW_MIN_R, NARROW_MEAN_R and NARROW_MAX_R (NaN otherwise).
    These pair with NARROW_*_MAG_GRAD_PSI to test the sdplsBeta equilibrium
    g* = beta - a directly, per time step, rather than inferring it from a
    refinement study.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    gradPsiErrorCSV
    {
        type            gradPsiErrorCSV;
        libs            ( leiaFunctionObjects );
        writeControl    timeStep;
        field           gradPsiError;
    }

    \endverbatim



    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func gradPsiErrorCSV
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::gradPsiErrorCSV
    ...

SourceFiles
    gradPsiErrorCSV.C

\*---------------------------------------------------------------------------*/

#ifndef gradPsiErrorCSV_H
#define gradPsiErrorCSV_H

#include "fvCFD.H"
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class gradPsiErrorCSV Declaration
\*---------------------------------------------------------------------------*/

class gradPsiErrorCSV
:
    public fvMeshFunctionObject
{

    const volScalarField& field_;

    const volScalarField& psi_;

    // const volVectorField& gradPsi_;
    OFstream csvFile_;

public:

    //- Runtime type information
    TypeName("gradPsiErrorCSV");


    // Constructors

        //- Construct from Time and dictionary
        gradPsiErrorCSV
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        gradPsiErrorCSV(const gradPsiErrorCSV&) = delete;

        //- No copy assignment
        void operator=(const gradPsiErrorCSV&) = delete;

        virtual bool execute()
        {
            return true;
        }

        virtual bool write();


    //- Destructor
    virtual ~gradPsiErrorCSV() = default;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
