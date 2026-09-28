/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: plicInterfaceSurface

    The piecewise-planar interface of a mesh as one polygon per interface
    cell: the plane with unit normal n_c (pointing OUT of the liquid) through
    the point x_c is cut through cell c (OpenFOAM's cutCellPLIC, libgeometricVoF)
    and its polygon appended to a meshedSurface; meshCells() maps every polygon
    back to its cell, so cell fields can be attached to the polygons.

    Ported from TwoPhaseFlow's src/postProcessing/interface/interface.{H,C}
    (Henning Scheufler, DLR, GPL-3) with the reconstruction decoupled from the
    reconstructionSchemes registry: the caller passes the normal and the point
    on the plane, so the DAVOF plane (m_c/|m_c|, xPlane_c), a geometricVoF
    scheme (-normal_/|normal_|, centre_) or anything else can be drawn.

    writeLegacyVTK() writes the polygons as ASCII legacy VTK polydata with
    per-polygon CELL_DATA taken from cell fields (the error fields of
    leiaTestDavofNormal, the cell index, ...), one file per rank.

    Author: leia (2026-09-29), after Henning Scheufler's interface class.
\*---------------------------------------------------------------------------*/

#ifndef plicInterfaceSurface_H
#define plicInterfaceSurface_H

#include "volFields.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"
#include "boolList.H"
#include "fileName.H"

namespace Foam
{

class fvMesh;

class plicInterfaceSurface
:
    public meshedSurface
{
    typedef meshedSurface MeshStorage;

    const fvMesh& mesh_;

    //- For every polygon the cell it lies in
    labelList meshCells_;

public:

    TypeName("plicInterfaceSurface");

    //- One polygon per cell with |normalOut[c]| > 0 (and interfaceCell[c],
    //  when given): the plane through pointOnPlane[c] with the unit normal
    //  normalOut[c]/|normalOut[c]| pointing out of the liquid.
    plicInterfaceSurface
    (
        const fvMesh& mesh,
        const vectorField& normalOut,
        const vectorField& pointOnPlane,
        const boolList* interfaceCell = nullptr
    );

    const labelList& meshCells() const { return meshCells_; }

    //- Point interpolation is not provided (the polygons are disjoint);
    //  returns a zero field of the right size.
    template<class Type>
    tmp<Field<Type>> interpolate
    (
        const GeometricField<Type, fvPatchField, volMesh>& cCoords,
        const Field<Type>& pCoords
    ) const;

    //- ASCII legacy VTK polydata: POINTS, POLYGONS and, per name, a
    //  CELL_DATA SCALARS array with cellFields[i][meshCells_[polygon]].
    void writeLegacyVTK
    (
        const fileName& file,
        const wordList& names,
        const List<const scalarField*>& cellFields
    ) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "plicInterfaceSurfaceTemplates.C"
#endif

#endif

// ************************************************************************* //
