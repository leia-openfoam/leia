/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCellCut
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
    Foam::davofCellCut

Description
    The cell cut of the davof reconstruction (the plane positioning and the
    explicit plane defects): OpenFOAM v2606's cutCellPLIC by default
    (grazingCut openfoam, Tomislav's decision of 2026-11-06), or the copy
    davofCutCellPLIC with the TwoPhaseFlow fork's grazing-cut threshold
    (grazingCut fork, 10*SMALL instead of ROOTVSMALL; the opt-in that proves
    the one arithmetic difference of the two libraries' cell cuts,
    docs/davof/PORT-DESIGN.md section 1.2 item 1). Both objects exist; one
    branch per call selects the one in use, so that the default path calls
    OpenFOAM's class exactly as before.

\*---------------------------------------------------------------------------*/

#ifndef davofCellCut_H
#define davofCellCut_H

#include "cutCellPLIC.H"
#include "davofCutCellPLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class davofCellCut Declaration
\*---------------------------------------------------------------------------*/

class davofCellCut
{
    // Private Data

        //- OpenFOAM's cut (the default)
        cutCellPLIC of_;

        //- The copy with the fork's threshold (grazingCut fork)
        davofCutCellPLIC fork_;

        //- Which one is in use
        const bool useFork_;


public:

    // Constructors

        //- Construct from the mesh and the rule (true: the fork's threshold)
        davofCellCut(const fvMesh& mesh, const bool useFork)
        :
            of_(mesh),
            fork_(mesh, 10*SMALL),
            useFork_(useFork)
        {}


    // Member Functions

        bool usesFork() const noexcept { return useFork_; }

        label calcSubCell(const label celli, const scalar cutValue, const vector& normal)
        {
            return useFork_ ? fork_.calcSubCell(celli, cutValue, normal)
                            : of_.calcSubCell(celli, cutValue, normal);
        }

        const point& subCellCentre() const { return useFork_ ? fork_.subCellCentre() : of_.subCellCentre(); }
        scalar subCellVolume() const { return useFork_ ? fork_.subCellVolume() : of_.subCellVolume(); }
        const DynamicList<point>& facePoints() { return useFork_ ? fork_.facePoints() : of_.facePoints(); }
        const point& faceCentre() const { return useFork_ ? fork_.faceCentre() : of_.faceCentre(); }
        const vector& faceArea() const { return useFork_ ? fork_.faceArea() : of_.faceArea(); }
        label cellStatus() const { return useFork_ ? fork_.cellStatus() : of_.cellStatus(); }
        scalar VolumeOfFluid() const { return useFork_ ? fork_.VolumeOfFluid() : of_.VolumeOfFluid(); }
        scalar cutValue() const { return useFork_ ? fork_.cutValue() : of_.cutValue(); }
        void clearStorage() { if (useFork_) fork_.clearStorage(); else of_.clearStorage(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
