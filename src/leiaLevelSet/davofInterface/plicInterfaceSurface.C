/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: plicInterfaceSurface (see the header)
\*---------------------------------------------------------------------------*/

#include "plicInterfaceSurface.H"
#include "fvMesh.H"
#include "cutCellPLIC.H"
#include "OFstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plicInterfaceSurface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicInterfaceSurface::plicInterfaceSurface
(
    const fvMesh& mesh,
    const vectorField& normalOut,
    const vectorField& pointOnPlane,
    const boolList* interfaceCell
)
:
    MeshStorage(),
    mesh_(mesh)
{
    cutCellPLIC cellCut(mesh_);
    const vectorField& C = mesh_.C().primitiveField();

    DynamicList<List<point>> facePts(label(0.1*mesh_.nCells()) + 16);
    DynamicList<label> cellAddressing(label(0.1*mesh_.nCells()) + 16);

    forAll(normalOut, cellI)
    {
        if (interfaceCell && !(*interfaceCell)[cellI]) continue;
        const scalar nMag = mag(normalOut[cellI]);
        if (nMag <= 0) continue;
        const vector n = normalOut[cellI]/nMag;
        const scalar cutVal = (pointOnPlane[cellI] - C[cellI]) & n;
        cellCut.calcSubCell(cellI, cutVal, n);
        const auto& fPoints = cellCut.facePoints();
        if (fPoints.size() >= 3)
        {
            facePts.append(List<point>(fPoints));
            cellAddressing.append(cellI);
        }
    }

    meshCells_.transfer(cellAddressing);

    // Transfer to the mesh storage: disjoint polygons, no shared points.
    faceList faces(facePts.size());
    label nPoints = 0;
    forAll(facePts, i)
    {
        face f(facePts[i].size());
        forAll(f, fi)
        {
            f[fi] = nPoints + fi;
        }
        faces[i] = f;
        nPoints += facePts[i].size();
    }
    pointField points(nPoints);
    nPoints = 0;
    forAll(facePts, i)
    {
        forAll(facePts[i], fi)
        {
            points[nPoints++] = facePts[i][fi];
        }
    }
    MeshStorage updated(std::move(points), std::move(faces), surfZoneList());
    this->MeshStorage::transfer(updated);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plicInterfaceSurface::writeLegacyVTK
(
    const fileName& file,
    const wordList& names,
    const List<const scalarField*>& cellFields
) const
{
    if (names.size() != cellFields.size())
    {
        FatalErrorInFunction
            << "names (" << names.size() << ") and cellFields ("
            << cellFields.size() << ") differ in size" << exit(FatalError);
    }
    mkDir(file.path());
    OFstream os(file);
    os.precision(12);
    const pointField& pts = this->points();
    const faceList& fcs = this->surfFaces();

    os  << "# vtk DataFile Version 2.0" << nl
        << "leia plicInterfaceSurface: one polygon per interface cell" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << pts.size() << " double" << nl;
    forAll(pts, i)
    {
        os << pts[i].x() << ' ' << pts[i].y() << ' ' << pts[i].z() << nl;
    }
    label listSize = 0;
    forAll(fcs, i) listSize += 1 + fcs[i].size();
    os << "POLYGONS " << fcs.size() << ' ' << listSize << nl;
    forAll(fcs, i)
    {
        os << fcs[i].size();
        forAll(fcs[i], j) os << ' ' << fcs[i][j];
        os << nl;
    }
    os << "CELL_DATA " << fcs.size() << nl;
    os << "SCALARS cellId int 1" << nl << "LOOKUP_TABLE default" << nl;
    forAll(meshCells_, i) os << meshCells_[i] << nl;
    forAll(names, k)
    {
        const scalarField& fld = *cellFields[k];
        os << "SCALARS " << names[k] << " double 1" << nl
           << "LOOKUP_TABLE default" << nl;
        forAll(meshCells_, i) os << fld[meshCells_[i]] << nl;
    }
}


// ************************************************************************* //
