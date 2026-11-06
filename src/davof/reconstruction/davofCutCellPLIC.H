/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCutCellPLIC
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    A copy of OpenFOAM v2606's cutCellPLIC (src/transportModels/geometricVoF):
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2018-2019 Johan Roenby
    Modified code Copyright (C) 2019-2020 DLR
    Modified code Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::davofCutCellPLIC

Description
    The cut of a cell by a plane (H. Scheufler, J. Roenby, Accurate and
    efficient surface reconstruction from volume fraction data on general
    meshes, J. Comput. Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009;
    the cutCell family of J. Roenby, H. Bredmose, H. Jasak, R. Soc. Open Sci.
    3 (2016) 160405, doi:10.1098/rsos.160405), as OpenFOAM v2606's
    cutCellPLIC computes it, with ONE parameter: the area below which a cut
    polygon counts as no cut and the cell as full or empty (the grazing cut,
    "a cell only touched at a point or a line"). OpenFOAM v2606 takes
    ROOTVSMALL (1e-150); the TwoPhaseFlow fork the DAVOF method was developed
    in takes 10*SMALL (1e-14). docs/davof/PORT-DESIGN.md section 1.2, item 1:
    the one arithmetic difference between the two libraries' cell cuts that
    the DAVOF code reaches, measured on the port's identity set (REPORT-swept.md,
    "The port into leia"): one grazing cut at one step of the rotation N 64
    changes a diagnostic count and anchorSkipConverged's decision for that
    cell, not the transported state.

    The opt-in `grazingCut fork` of the davof reconstruction dictionary
    (Tomislav, 2026-11-06) selects this class with the fork's threshold through
    davofCellCut (the adapter at the call sites) to prove that the threshold is
    the whole difference; the default stays OpenFOAM's own class.

SourceFiles
    davofCutCellPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef davofCutCellPLIC_H
#define davofCutCellPLIC_H

#include "cutCell.H"
#include "cutFacePLIC.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class davofCutCellPLIC Declaration
\*---------------------------------------------------------------------------*/

class davofCutCellPLIC
:
    public cutCell
{
    // Private Data

        //- Mesh whose cells and faces to cut
        const fvMesh& mesh_;

        //- The area below which a cut polygon is no cut (the grazing cut)
        const scalar areaTol_;

        //- Cell to cut
        label cellI_;

        //- Normal of the cutting plane
        vector normal_;

        //- Cutvalue used to cut cell
        scalar cutValue_;

        //- The face cutting
        cutFacePLIC cutFace_;

        //- List of face centres for CutFaces
        DynamicList<point> cutFaceCentres_;

        //- List of face area vectors for PLICCutFaces
        DynamicList<vector> cutFaceAreas_;

        //- Storage for subFace edges belonging to PLICFace
        DynamicList<DynamicList<point>> plicFaceEdges_;

        //- Points constituting the cell-PLICsurface intersection (isoface)
        DynamicList<point> facePoints_;

        //- Face centre of the cutFace
        point faceCentre_;

        //- Face normal of the PLICface (from high to low values)
        vector faceArea_;

        //- Cell centre of the submerged subcell
        point subCellCentre_;

        //- Volume of the submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- -1: fully below the plane, 0: cut, +1: fully above
        label cellStatus_;


public:

    // Constructors

        //- Construct from fvMesh and the grazing-cut area tolerance
        davofCutCellPLIC(const fvMesh& mesh, const scalar areaTol);


    // Member Functions

        //- Sets internal values and returns the cell status
        label calcSubCell
        (
            const label celli,
            const scalar cutValue,
            const vector& normal
        );

        const point& subCellCentre() const noexcept { return subCellCentre_; }
        scalar subCellVolume() const noexcept { return subCellVolume_; }
        const DynamicList<point>& facePoints();
        const point& faceCentre() const noexcept { return faceCentre_; }
        const vector& faceArea() const noexcept { return faceArea_; }
        label cellStatus() const noexcept { return cellStatus_; }
        scalar VolumeOfFluid() const noexcept { return VOF_; }
        scalar cutValue() const noexcept { return cutValue_; }
        scalar areaTol() const noexcept { return areaTol_; }

        //- Resets internal values
        void clearStorage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
