/*---------------------------------------------------------------------------*\
    leia :: davof :: davofAdvectionScheme
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    Derived from TwoPhaseFlow's advectionSchemes base class:
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::davofAdvectionScheme

Description
    The run-time selectable advection of a volume fraction with a geometric
    reconstruction, the base of the two schemes of leia's DAVOF library
    (docs/davof/PORT-DESIGN.md, 2026-11-06):

      davofSwept   the dual advection of the DAVOF state (the cell volume
                   fractions and the face area fractions) with one surface
                   per cell (davofSweptAdvection.H);
      isoAdvector  OpenFOAM's geometric VoF advection, Foam::isoAdvection of
                   libgeometricVoF (J. Roenby, H. Bredmose, H. Jasak, A
                   computational method for sharp interface advection, R.
                   Soc. Open Sci. 3 (2016) 160405, doi:10.1098/rsos.160405;
                   the bounding and the reconstruction schemes: H. Scheufler,
                   J. Roenby, Accurate and efficient surface reconstruction
                   from volume fraction data on general meshes, J. Comput.
                   Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009),
                   the comparison arm (isoAdvectorScheme.H).

    Selected by the key advectionScheme of the alpha solver dictionary
    (system/fvSolution, solvers/<alpha>), as TwoPhaseFlow's benchmark solver
    selected it. Registered in the object registry under the name
    davofAdvectionScheme, where the function objects (areaFractionError,
    volumeFractionError) find it and its reconstruction surf().

    The base holds the fields and the face volumetric transport alphaPhi_;
    the reconstruction is the scheme's own (surf() is virtual): the dual
    advection selects it through reconstructionSchemes::New, isoAdvection
    constructs its own.

SourceFiles
    davofAdvectionScheme.C
    davofAdvectionSchemeNew.C

\*---------------------------------------------------------------------------*/

#ifndef davofAdvectionScheme_H
#define davofAdvectionScheme_H

#include "typeInfo.H"
#include "fvc.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class davofAdvectionScheme Declaration
\*---------------------------------------------------------------------------*/

class davofAdvectionScheme
:
    public IOdictionary
{
    // Private data

        //- The alpha solver dictionary (system/fvSolution, solvers/<alpha>)
        dictionary advectionSchemesCoeffs_;


protected:

    // Protected data

        //- Reference to the VoF field
        volScalarField& alpha1_;

        //- Reference to the face fluxes
        const surfaceScalarField& phi_;

        //- Reference to the velocity field
        const volVectorField& U_;

        //- Face volumetric transport (the liquid volume flux per time)
        surfaceScalarField alphaPhi_;


    // Private Member Functions

        //- No copy construct
        davofAdvectionScheme(const davofAdvectionScheme&) = delete;

        //- No copy assignment
        void operator=(const davofAdvectionScheme&) = delete;


public:

    //- Runtime type information
    TypeName("davofAdvectionScheme");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        davofAdvectionScheme,
        components,
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        ),
        (alpha1, phi, U)
    );


    // Selectors

        //- Select by the key advectionScheme of the alpha solver dictionary
        static autoPtr<davofAdvectionScheme> New
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    // Constructors

        //- Construct from components
        davofAdvectionScheme
        (
            const word& type,
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~davofAdvectionScheme() = default;


    // Member Functions

        //- The alpha solver dictionary
        dictionary& modelDict();

        //- The alpha solver dictionary
        const dictionary& modelDict() const;

        //- Advection of the interface
        virtual void advect() = 0;

        //- Advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        ) = 0;

        //- The reconstruction of the scheme
        virtual reconstructionSchemes& surf() = 0;

        //- Reference to interface normals
        const volVectorField& normal()
        {
            return surf().normal();
        }

        //- Reference to interface centres
        const volVectorField& centre()
        {
            return surf().centre();
        }

        //- Return mass flux
        tmp<surfaceScalarField> getRhoPhi
        (
            const dimensionedScalar rho1,
            const dimensionedScalar rho2
        ) const
        {
            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    "rhoPhi",
                    (rho1 - rho2)*alphaPhi_ + rho2*phi_
                )
            );
        }

        //- Return mass flux
        tmp<surfaceScalarField> getRhoPhi
        (
            const volScalarField& rho1,
            const volScalarField& rho2
        )
        {
            return tmp<surfaceScalarField>
            (
               new surfaceScalarField
               (
                   "rhoPhi",
                   fvc::interpolate(rho1 - rho2)*alphaPhi_
                 + fvc::interpolate(rho2)*phi_
               )
            );
        }

        //- The face volumetric transport of the last step
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- The execution times (reconstruction, advection, unused) that
        //  TwoPhaseFlow's base class declared and no scheme ever filled
        //  (zero in every fork run); volumeFractionError divides them into
        //  its two per-step columns, which therefore stay zero here as in
        //  the fork's CSVs (the port reproduces those files bit for bit)
        vector runTime() const
        {
            return vector::zero;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
