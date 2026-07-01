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

#include "IOobject.H"
#include "dimensionedScalarFwd.H"
#include "geometricPhaseIndicator.H"
#include "addToRunTimeSelectionTable.H"
#include "processorFvPatch.H"
#include "pTraits.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "levelSetImplicitSurfaces.H"
#include "foamGeometry.H"
#include "simpleMatrix.H"
#include "processorFvPatch.H"
#include "volFieldsFwd.H"
#include "levelSetPlaneReconstruction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(geometricPhaseIndicator, false);
addToRunTimeSelectionTable(phaseIndicator, geometricPhaseIndicator, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geometricPhaseIndicator::geometricPhaseIndicator(const fvMesh& mesh)
    :
        phaseIndicator(mesh),
        narrowBand_(mesh.lookupObject<volScalarField>("NarrowBand")),
        ncTmp_(new volVectorField
            (
                IOobject
                (
                    "nc",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE // FIXME(TM): readable IOOptions.
                ),
                mesh,
                dimensionedVector("nc", dimless, vector(0,0,0))
            )
        ),
        dcTmp_(new volScalarField
            (
                IOobject
                (
                    "dc",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE // FIXME(TM): readable IOOptions.
                ),
                mesh,
                dimensionedScalar("dc", dimless, 0)
            )
        )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void geometricPhaseIndicator::calcPhaseIndicator
(
    volScalarField& alpha,
    const volScalarField& psi
)
{
    alpha == dimensionedScalar("0", dimless, 0);
    const fvMesh& mesh = alpha.mesh();

    // Set phase-indicator to 1 and 0 in the bulk.
    forAll(alpha, cellID)
    {
        if (psi[cellID] < 0)
            alpha[cellID] = 1;
        else
            alpha[cellID] = 0;
    }

    // Approach: Lookup existing narrowBand
    const volScalarField& narrowBand = narrowBand_;

    // Linear Least-Squares Approximation of \Psi
    // \psi(x,y,z)= nc_x x  + nc_y y + nc_z z + dc
    volVectorField& nc_ = ncTmp_.ref();
    nc_ = dimensionedVector("nc", nc_.dimensions(), vector(0,0,0));
    volScalarField& dc_ = dcTmp_.ref();
    dc_ = dimensionedScalar("dc", dc_.dimensions(), 0.);

    // Coupled-patch neighbour psi / cell centres: evaluated once (matched MPI
    // call on every rank) so the per-cell reconstruction does no communication.
    const coupledFaceNeighbours coupledNei(mesh, psi);

    // Compute the phase indicator in the narrow band.
    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] == 1)
        {
            // Linear least-squares plane reconstruction, shared with
            // detrixheAslamPhaseIndicator: psi^l(x) = nc & x + dc.
            const scalarList planeCoeffs =
                leastSquaresPlaneCoeffs(mesh, psi, cellI, coupledNei);

            // Debugging
            nc_[cellI] = vector(planeCoeffs[0], planeCoeffs[1], planeCoeffs[2]);
            dc_[cellI] = planeCoeffs[3];

            // Ensure unit-normals
            hesseNormalPlane cutPlane(
                vector (planeCoeffs[0], planeCoeffs[1], planeCoeffs[2]),
                planeCoeffs[3]
            );

            // TODO: Use the cutPlane to test for intersection and do not
            //       intersect cells that do not require intersection. TM.
            auto cellIntersection = intersectCell(cellI, mesh, cutPlane);

            // TODO: Check intersection when the plane overlaps with the face. TM.
            //       - Bounding cell intersection volume is a quick fix.
            alpha[cellI] = max(cellIntersection.volume(), 0.) / mesh.V()[cellI];
        }
    }

    nc_.correctBoundaryConditions();
    dc_.correctBoundaryConditions();
    alpha.correctBoundaryConditions();
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
