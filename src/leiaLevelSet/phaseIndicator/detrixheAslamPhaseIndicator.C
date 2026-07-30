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

\*---------------------------------------------------------------------------*/

#include "detrixheAslamPhaseIndicator.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "levelSetPlaneReconstruction.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(detrixheAslamPhaseIndicator, false);
addToRunTimeSelectionTable(phaseIndicator, detrixheAslamPhaseIndicator, Mesh);

// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace
{

//- Fraction of a tetrahedron occupied by Omega^- = {phi < 0}, given the four
//  signed-distance vertex values. Detrixhe & Aslam (2016), Eq. (7).
//
//  The two-in-two-out (middle) numerator below is the Omega^- form; the
//  numerator printed in the paper's Eq. (7) middle case evaluates to the
//  complement (Omega^+) fraction (verified against the divided-difference
//  formula and direct integration), so it has been corrected here.
//
//  In every non-trivial case the denominator is a product of differences of
//  vertex values that lie on opposite sides of the interface, so it is bounded
//  away from zero by the cross-interface gap -- there are no clipping
//  tolerances. stabilise() only guards the measure-zero case of a vertex lying
//  exactly on the interface; the result is clamped to [0,1] by the caller.
inline scalar tetNegativeFraction(scalar a, scalar b, scalar c, scalar d)
{
    scalar p[4] = {a, b, c, d};
    std::sort(p, p + 4);
    const scalar p1 = p[0], p2 = p[1], p3 = p[2], p4 = p[3];

    if (p4 <= 0) return 1;          // all four inside  (Omega^-)
    if (p1 >  0) return 0;          // all four outside (Omega^+)

    if (p3 <= 0)                    // three inside, one outside (corner at p4)
    {
        const scalar den = (p4 - p1)*(p4 - p2)*(p4 - p3);
        return 1 - (p4*p4*p4)/stabilise(den, VSMALL);
    }
    if (p2 <= 0)                    // two inside, two outside (middle)
    {
        const scalar num =
            p3*p4*(p1*p1 + p1*p2 + p2*p2)
          - p1*p2*(p1 + p2)*(p3 + p4)
          + p1*p1*p2*p2;
        const scalar den = (p3 - p1)*(p4 - p1)*(p3 - p2)*(p4 - p2);
        return num/stabilise(den, VSMALL);
    }
    // one inside, three outside (corner at p1)
    const scalar den = (p2 - p1)*(p3 - p1)*(p4 - p1);
    return -(p1*p1*p1)/stabilise(den, VSMALL);
}

} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

detrixheAslamPhaseIndicator::detrixheAslamPhaseIndicator(const fvMesh& mesh)
:
    phaseIndicator(mesh),
    narrowBand_(mesh.lookupObject<volScalarField>("NarrowBand")),
    geometrySource_
    (
        phaseIndDict_.getOrDefault<word>("geometrySource", "levelSetField")
    ),
    analyticSurface_()
{
    if
    (
        geometrySource_ != "levelSetField"
     && geometrySource_ != "analyticImplicitSurface"
    )
    {
        FatalIOErrorInFunction(phaseIndDict_)
            << "Unknown Detrixhe-Aslam geometrySource '" << geometrySource_
            << "'. Valid: levelSetField, analyticImplicitSurface."
            << exit(FatalIOError);
    }
    if (geometrySource_ == "analyticImplicitSurface")
    {
        const dictionary& surfaceDict = levelSetDict_.subDict("implicitSurface");
        analyticSurface_ = implicitSurface::New
        (
            surfaceDict.get<word>("type"),
            surfaceDict
        );
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void detrixheAslamPhaseIndicator::calcPhaseIndicator
(
    volScalarField& alpha,
    const volScalarField& psi
)
{
    const fvMesh& mesh = alpha.mesh();

    // Bulk: sign-based indicator.
    const bool analyticGeometry =
        geometrySource_ == "analyticImplicitSurface";
    forAll(alpha, cellID)
    {
        const scalar signValue = analyticGeometry
          ? analyticSurface_->value(mesh.cellCentres()[cellID])
          : psi[cellID];
        alpha[cellID] = (signValue < 0) ? 1 : 0;
    }

    const volScalarField& narrowBand = narrowBand_;

    const pointField&  points      = mesh.points();
    const faceList&    faces       = mesh.faces();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& cellCentres = mesh.cellCentres();

    // Coupled-patch neighbour psi / cell centres: evaluated once (matched MPI
    // call on every rank) so the per-cell reconstruction does no communication.
    const coupledFaceNeighbours coupledNei(mesh, psi);

    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] != 1) continue;

        // Production: the same linear least-squares reconstruction as
        // geometricPhaseIndicator.  Oracle: the exact signed-distance tangent
        // plane at the cell centre, phi(x)=n.(x-xc)+surface.value(xc).
        vector nc(Zero);
        scalar d = 0;
        if (analyticGeometry)
        {
            const point& xc = cellCentres[cellI];
            nc = analyticSurface_->grad(xc);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                nc /= nmag;
                d = analyticSurface_->value(xc) - (nc & xc);
            }
        }
        else
        {
            const scalarList pc =
                leastSquaresPlaneCoeffs(mesh, psi, cellI, coupledNei);
            nc = vector(pc[0], pc[1], pc[2]);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                d = pc[3]/nmag;
            }
        }

        const scalar nmag = mag(nc);
        if (nmag < SMALL)
        {
            // Degenerate (flat) reconstruction: fall back to the sign.
            const scalar signValue = analyticGeometry
              ? analyticSurface_->value(cellCentres[cellI])
              : psi[cellI];
            alpha[cellI] = (signValue < 0) ? 1 : 0;
            continue;
        }

        // Signed-distance plane: phi(x) = n & x + d, |n| = 1.
        const vector n = nc/nmag;

        const point& xc   = cellCentres[cellI];
        const scalar phic = (n & xc) + d;

        scalar liquidVol = 0;
        scalar totalVol  = 0;

        // Decompose the cell into tetrahedra (cell centroid, face centroid,
        // and the two end-points of each face edge) and fill each tet
        // analytically.
        const labelList& cFaces = mesh.cells()[cellI];
        forAll(cFaces, cf)
        {
            const label  faceI = cFaces[cf];
            const face&  f     = faces[faceI];
            const point& xf    = faceCentres[faceI];
            const scalar phif  = (n & xf) + d;

            forAll(f, ip)
            {
                const point& p0 = points[f[ip]];
                const point& p1 = points[f.nextLabel(ip)];
                const scalar phi0 = (n & p0) + d;
                const scalar phi1 = (n & p1) + d;

                // Tetrahedron volume (orientation independent).
                const scalar vol =
                    mag(((xf - xc) ^ (p0 - xc)) & (p1 - xc))/6.0;
                if (vol <= VSMALL) continue;

                liquidVol += tetNegativeFraction(phic, phif, phi0, phi1)*vol;
                totalVol  += vol;
            }
        }

        if (totalVol > VSMALL)
        {
            alpha[cellI] =
                min(max(liquidVol/totalVol, scalar(0)), scalar(1));
        }
    }

    alpha.correctBoundaryConditions();
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
