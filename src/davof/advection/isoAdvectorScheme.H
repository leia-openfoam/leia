/*---------------------------------------------------------------------------*\
    leia :: davof :: isoAdvectorScheme
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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
    Foam::isoAdvectorScheme

Description
    The comparison arm of the DAVOF studies: OpenFOAM's own geometric VoF
    advection, Foam::isoAdvection of libgeometricVoF (OpenFOAM v2606), with
    its own reconstruction scheme (plicRDF, isoAlpha, gradAlpha, selected by
    reconstructionScheme of the alpha solver dictionary), behind the
    davofAdvectionScheme interface (advectionScheme isoAdvector).

    J. Roenby, H. Bredmose, H. Jasak, A computational method for sharp
    interface advection, R. Soc. Open Sci. 3 (2016) 160405,
    doi:10.1098/rsos.160405. H. Scheufler, J. Roenby, Accurate and efficient
    surface reconstruction from volume fraction data on general meshes,
    J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009.

    This is OpenFOAM proper's isoAdvection, not the TwoPhaseFlow fork's; the
    two differ (docs/davof/PORT-DESIGN.md section 1.3).

SourceFiles
    isoAdvectorScheme.C

\*---------------------------------------------------------------------------*/

#ifndef isoAdvectorScheme_H
#define isoAdvectorScheme_H

#include "davofAdvectionScheme.H"
#include "isoAdvection.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class isoAdvectorScheme Declaration
\*---------------------------------------------------------------------------*/

class isoAdvectorScheme
:
    public davofAdvectionScheme
{
    // Private data

        //- OpenFOAM's isoAdvection (it owns its reconstruction scheme)
        isoAdvection iso_;


public:

    //- Runtime type information
    TypeName("isoAdvector");


    // Constructors

        //- Construct given alpha, phi and the velocity field
        isoAdvectorScheme
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~isoAdvectorScheme() = default;


    // Member Functions

        //- Advection of the interface
        virtual void advect()
        {
            zeroField Sp;
            zeroField Su;
            iso_.advect<zeroField, zeroField>(Sp, Su);
            alphaPhi_ = iso_.alphaPhi();
        }

        //- Advection of the interface with source terms
        virtual void advect
        (
            const volScalarField::Internal& Sp,
            const volScalarField::Internal& Su
        )
        {
            iso_.advect<volScalarField::Internal, volScalarField::Internal>
            (
                Sp,
                Su
            );
            alphaPhi_ = iso_.alphaPhi();
        }

        //- The reconstruction of isoAdvection
        virtual reconstructionSchemes& surf()
        {
            return iso_.surf();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
