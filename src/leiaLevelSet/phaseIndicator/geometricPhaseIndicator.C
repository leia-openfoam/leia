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

    const auto& psiBdryField = psi.boundaryField();
    const auto& patches = mesh.boundary();
    const auto& cellCells = mesh.cellCells(); 
    const auto& cellCenters = mesh.C();
    // Compute the phase indicator in the narrow band.
    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] == 1)
        {
            // Assemble the LLSQ linear system
            // TODO: Extend for 2D simulations in OpenFOAM. 
            //SquareMatrix<scalar> LLSQ(4,0);
            simpleMatrix LLSQ
            (
                4  /* size 4x4 */, 
                0. /* init coeff value*/, 
                0. /* init source value*/ 
            );
            scalarList planeCoeffs(4,0);
            scalarList LLSQsource(4,0);
        
            // NOTE: Nc doesn't contain cellI in OpenFOAM: cellI contributions
            // are handled additionallly. 
            const auto& Nc = cellCells[cellI];
                
            // equations 0,1,2
            // \frac{\partial e^{lsq}_c} {\partial n_{c,cmpt}} = 0 
            for (char row = 0; row < 3; ++row)
            {
                // Contributions from cellI : not in Nc
                // - nc coefficient contrib from cellI
                for(char col = 0; col < 3; ++col)
                    LLSQ(row,col) += 
                        cellCenters[cellI][col]*cellCenters[cellI][row];
                // - dc coefficient contrib from cellI 
                LLSQ(row, 3) += cellCenters[cellI][row];
                // - source contrib from cellI
                LLSQsource[row] += psi[cellI]*cellCenters[cellI][row];

                // Contributions from Nc 
                forAll(Nc, cellK)
                {
                    // - nc coefficient contrib from cellK
                    for(char col = 0; col < 3; ++col)
                        LLSQ(row,col) += cellCenters[Nc[cellK]][col]*cellCenters[Nc[cellK]][row];
                    // - dc coefficient contrib from cellK
                    LLSQ(row, 3) += cellCenters[Nc[cellK]][row];
                    // - source contrib from cellK
                    LLSQsource[row] += psi[Nc[cellK]]*cellCenters[Nc[cellK]][row];
                }
            }
            // equation 3
            // \frac{\partial e^{lsq}_c}{\partial d_c} = 0 
            // - nc coefficient contrib from cellI
            for(char col = 0; col < 3; ++col)
                LLSQ(3,col) += cellCenters[cellI][col];
            // - source contrib for cellI
            LLSQsource[3] += psi[cellI];

            // Contributions from Nc 
            forAll(Nc, cellK)
            {
                // - nc coefficient contrib from cellK
                for(char col = 0; col < 3; ++col)
                    LLSQ(3,col) += cellCenters[Nc[cellK]][col];
                // source term
                LLSQsource[3] += psi[Nc[cellK]];
            }
            // - dc coefficient is |Nc|: in OpenFOAM it is |Nc| + 1 (for cellI)
            LLSQ(3,3) = Nc.size() + 1;
            
            // Contributions from face-adjacent cells from coupled boundaries 
            const auto& nBandCell = mesh.cells()[cellI];
            const auto& cellCentersBdryField = cellCenters.boundaryField();
            forAll(nBandCell, faceI)
            {
                const label faceJ = nBandCell[faceI];
                // If the face is not internal. 
                if (! mesh.isInternalFace(faceJ))
                {
                    // Find the patch faceJ belongs to.
                    // - patch label.
                    label patchL = -1; 
                    // - face index local to the patch
                    label faceP = -1;

                    forAll(patches, patchI)
                    {
                        const fvPatch& patch = patches[patchI];
                        const polyPatch& ppatch = patch.patch();
                        faceP = ppatch.whichFace(faceJ);
                        // Patch is a coupled patch and the face faceP belongs to it.
                        if (isA<coupledFvPatch>(patch) && (faceP >= 0) && (faceP < patch.size())) 
                            patchL = patchI;
                    }
                    // If faceJ (faceP) belongs to the patchL coupled patch.
                    if ((patchL != -1) && (faceP >= 0)) 
                    {
                        // The face belongs to the coupled boundary. 
                        // Fetch the cell center of the cell.  
                        const auto& cellCentersPatchField = cellCentersBdryField[patchL]; 
                        auto cellCentersPatchNeiFieldTmp = cellCentersPatchField.patchNeighbourField();
                        const auto& cellCentersPatchNeiField = cellCentersPatchNeiFieldTmp();
                        const auto& cellCenter = cellCentersPatchNeiField[faceP];

                        // Fetch the psi value of the cell.  
                        const auto& psiPatchField = psiBdryField[patchL]; 
                        auto psiPatchNeiFieldTmp = psiPatchField.patchNeighbourField();
                        const auto& psiPatchNeiField = psiPatchNeiFieldTmp();
                        const auto& psiValue = psiPatchNeiField[faceP];
                        
                        // equations 0,1,2 contrib from coupled faceJ-adjacent cell
                        // \frac{\partial e^{lsq}_c} {\partial n_{c,cmpt}} = 0 
                        for (char row = 0; row < 3; ++row)
                        {
                            // - nc coefficient contrib from coupled faceJ-adjacent cell 
                            for(char col = 0; col < 3; ++col)
                                LLSQ(row,col) += cellCenter[col]*cellCenter[row];

                            // - dc coefficient contrib from coupled faceJ-adjacent cell 
                            LLSQ(row, 3) += cellCenter[row];

                            // - source contrib from coupled faceJ-adjacent cell 
                            LLSQsource[row] += psiValue*cellCenter[row];
                        }

                        // equation 3 contrib from coupled faceJ-adjacent cell
                        // \frac{\partial e^{lsq}_c}{\partial d_c} = 0 
                        
                        // - nc coefficient contrib 
                        for(char col = 0; col < 3; ++col)
                              LLSQ(3,col) += cellCenter[col];

                        // - dc coeff contrib for faceJ (faceP) coupled neighbor-cell
                        LLSQ(3,3) += 1;
                        
                        // - source term contrib
                        LLSQsource[3] += psiValue;
                    }
                }
            }

            //LUscalarMatrix LU (LLSQ); 
            // Crashes in parallel, uses PStream::parRun() even if the linear  
            // system solution is local to the coupled process. TM. 
            //LU.solve(planeCoeffs, LLSQsource);
            //Foam::LUsolve(LLSQ, LLSQsource);
            scalarField& source = LLSQ.source(); 
            source = LLSQsource;
            planeCoeffs = LLSQ.solve(); // TODO: Improve this. Gauss substitution. TM.
            
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
