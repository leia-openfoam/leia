/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
    unofficial extension to OpenFOAM.

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
    Foam::functionObjects::reconstructionError

Description

Author
    Henning Scheufler, DLR, all rights reserved.
    Johan Roenby, DHI, all rights reserved.

SourceFiles
    reconstructionError.C

\*---------------------------------------------------------------------------*/

#ifndef reconstructionError_H
#define reconstructionError_H


#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

/*---------------------------------------------------------------------------*\
                       Class reconstructionError Declaration
\*---------------------------------------------------------------------------*/

class reconstructionError
:
    public functionObjects::writeFile
{
        const fvMesh& mesh_;

        scalar LNormalInf_;

        scalar LNormal1_;

        scalar LCentreInf_;

        scalar LCentre1_;

        scalar LCurvInf_;

        scalar LCurv1_;

        DynamicList<scalar> avgLNormal_;
        DynamicList<scalar> avgLCentre_;
        DynamicList<scalar> avgLCurv_;

        // Appended columns (davof, 2026-09-28): the normal error |n - n_exact|
        // (n and n_exact both unit, into alpha = 1) over the interface cells as
        // mean (1), root mean square (2) and maximum (Inf); the total
        // reconstructed area sum |normal|; the exact area (4 pi R^2 for a
        // sphere, 0 otherwise); the interface-cell count. Values of the last
        // calcError call, reduced over the ranks.
        scalar LNormalDiff1_;
        scalar LNormalDiff2_;
        scalar LNormalDiffInf_;
        scalar areaSum_;
        scalar areaExact_;
        label nInterface_;

        // Appended columns (davof, 2026-10-12): the position error, the
        // distance of centre() from the implicit surface (exact |R - |x - o||
        // for the sphere, the Newton-projected distance otherwise), over the
        // interface cells as mean (1), root mean square (2) and maximum (Inf).
        scalar LCentreDiff1_;
        scalar LCentreDiff2_;
        scalar LCentreDiffInf_;

        // Appended columns (davof, 2026-11-12): the angle between the unit
        // normals theta = 2 asin(|n - n_e|/2) and 1 - n.n_e = |n - n_e|^2/2,
        // over the interface cells as mean (1) and maximum (Inf). For a
        // normal whose angle converges at order p, 1 - n.n_e converges at
        // order 2p. LNormalOneMinusDotInf is the normal error of the plicRDF
        // benchmark (H. Scheufler, J. Roenby, Accurate and efficient surface
        // reconstruction from volume fraction data on general meshes,
        // J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009,
        // arXiv:1801.05382, Sec. 4.1, Eq. (21): L_inf^n = max_i(1 -
        // n_S,i . n_S,i^exact) over all interface cells), the original
        // LNormalInf_/LNormal1_ columns of this function object being the
        // same quantity accumulated over the calls.
        scalar LNormalAngle1_;
        scalar LNormalAngleInf_;
        scalar LNormalOneMinusDot1_;
        scalar LNormalOneMinusDotInf_;

    // Protected Member Functions

        //- File header information
        void writeFileHeader(Ostream& os) const;



private:

    // Private member functions

        //- Disallow default bitwise copy construct
        reconstructionError(const reconstructionError&);

        //- Disallow default bitwise assignment
        void operator=(const reconstructionError&);



public:


    // Constructors

        //- Construct from Time and dictionary
        reconstructionError
        (
            const fvMesh& mesh,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reconstructionError();


    // Member Functions


        //- Calculate the wall heat-flux
        void calcError
        (
            const dictionary& dict,
            const volVectorField& centre,
            const volVectorField& normal,
            const volScalarField& curv,
            bool randomMode
        );

        //- Write the wall heat-flux
        void write();

        //- Write the wall heat-flux
        void write
        (
            const scalar k1,
            const scalar k2
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
