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

Application
    
    leiaTestLevelSet

Description
    TODO

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "phaseIndicator.H"
#include "redistancer.H"
#include "narrowBand.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void redistancedPlaneIsUnchanged
(
    autoPtr<redistancer>& redist,
    autoPtr<narrowBand>& nBand,
    const fvMesh& mesh,
    volScalarField& psi
)
{
    // TEST: Redistancing a signed distance from a plane has no
    //       effect if the redistancing algorithm is second-order
    //       accurate. Boundary values are calculated exactly.
    //
    //       Exactness domain by model: the plane-anchored geometric models
    //       (planeFootWave, anchoredEikonal) are machine-exact on their
    //       ANCHOR set (narrow band + first face-neighbour ring; the LLS
    //       plane of planar data is the plane), while their far field is a
    //       foot-point-cloud distance / iterated Eikonal fill -- accurate to
    //       O(h^2/d) / solver tolerance BY DESIGN, not machine-exact. The
    //       other models (e.g. PDE) are asserted globally as before.

    // Set psi to a signed distance from a plane.
    const Time& runTime = mesh.time();
    const volVectorField& xc = mesh.C();

    dimensionedVector pTest = Foam::average(xc);  

    // TODO: extend the test with a parameterized normal orientation 

    // Zero the plane-normal component orthogonal to the 2D solution plane.
    // In the case of a 3D solution, all directions are 1.
    vector nVec (1,2.1,3);
    const auto& solutionDimensions = mesh.solutionD();
    for (int i = 0; i < 3; ++i)
    {
        if (solutionDimensions[i] < 0)
            nVec[i] = 0; 
    }
    // Normalize to obtain signed-distance
    nVec /= Foam::mag(nVec);
    dimensionedVector nTest ("nTest", dimless, nVec); 
    
    // Initialize the CASE psi as the signed distance to the plane. The
    // GeometricField operator== assigns internal AND patch values, so the
    // domain-boundary signed distances are prescribed exactly regardless of
    // the boundary-condition type.
    psi == (nTest & (xc - pTest));

    // Save the non-redistanced field for error evaluation
    volScalarField psi0("psiPlane0", psi);

    // Refresh the registered "NarrowBand" for the plane field: the
    // plane-anchored redistancer models read it to select anchor cells.
    nBand->calc();

    // |\nabla \psi0|
    volScalarField magGradPsi0("magGradPsi0", Foam::mag(fvc::grad(psi0)));
    Info << "mean(|\\nabla \\psi0|) : " << 
        Foam::average(magGradPsi0).value() << endl;

    // ||\nabla \psi0| - 1|
    volScalarField eGradPsi0 = Foam::mag(
        Foam::mag(magGradPsi0) - 
        dimensionedScalar("1", dimless, 1)
    );
    eGradPsi0.rename("eGradPsi0");

    // Redistance the signed distance to a plane (nothing should happen) 
    redist->redistance(psi);

    // |\nabla \psi| 
    volScalarField magGradPsi("magGradPsi", Foam::mag(fvc::grad(psi)));
    Info << "mean(|\\nabla \\psi|) : " << 
        Foam::average(magGradPsi).value() << endl;

    volScalarField ePsi = Foam::mag(psi - psi0);
    ePsi.rename("ePsi");

    // Valid errors available only at cell centers (internal field)
    volScalarField eGradPsi = Foam::mag(
        magGradPsi -
        dimensionedScalar("1", dimless, 1)
    );
    eGradPsi.rename("eGradPsi");
    eGradPsi.write();

    const fvSolution& fvSolutionDict(mesh);
    const word redistType =
        fvSolutionDict.subDict("levelSet").subDict("redistancer")
        .getOrDefault<word>("type", "noRedistancing");
    const bool planeAnchored =
        (redistType == "planeFootWave" || redistType == "anchoredEikonal");

    // The gradient checks are INTERIOR-cell checks: boundary-cell gradients
    // depend on the grad(psi) scheme's boundary treatment (e.g. noBvGrad
    // deliberately ignores patch values), not on the redistancer.
    boolList interiorCell(mesh.nCells(), true);
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<emptyPolyPatch>(pp) || pp.coupled()) continue;
        for (const label celli : pp.faceCells())
        {
            interiorCell[celli] = false;
        }
    }

    scalar ePsiLinf = Foam::max(ePsi).value();
    scalar eGradPsiLinf = 0;
    forAll(interiorCell, celli)
    {
        if (interiorCell[celli])
        {
            eGradPsiLinf = Foam::max(eGradPsiLinf, eGradPsi[celli]);
        }
    }
    reduce(eGradPsiLinf, maxOp<scalar>());
    Info << "eGradPsiLinf (interior): " << eGradPsiLinf << endl;

    if (planeAnchored)
    {
        // Machine exactness on the anchor set (band + first ring)...
        const volScalarField& nb =
            mesh.lookupObject<volScalarField>("NarrowBand");
        const labelListList& cellCells = mesh.cellCells();
        boolList anchorSet(mesh.nCells(), false);
        forAll(nb, celli)
        {
            if (nb[celli] == 1)
            {
                anchorSet[celli] = true;
                for (const label nbrI : cellCells[celli])
                {
                    anchorSet[nbrI] = true;
                }
            }
        }
        scalar ePsiLinfAnchors = 0;
        forAll(anchorSet, celli)
        {
            if (anchorSet[celli])
            {
                ePsiLinfAnchors = Foam::max(ePsiLinfAnchors, ePsi[celli]);
            }
        }
        reduce(ePsiLinfAnchors, maxOp<scalar>());
        Info << "ePsiLinf (band + ring): " << ePsiLinfAnchors << endl;
        if (ePsiLinfAnchors > 1e-12)
        {
            FatalErrorInFunction
                << redistType << ": planar signed distance not reproduced "
                << "on the anchor set (band + first ring): Linf = "
                << ePsiLinfAnchors
                << abort(FatalError);
        }

        // ... and sanity bounds over the TRIGGER BAND {|psi| <= 6h} -- the
        // region the phase indicator and the gradPsiThreshold criterion
        // read. The far field is each model's fill (donor-plane wave /
        // first-order Eikonal solve), approximate BY CONTRACT; it is
        // reported, not asserted.
        const scalar h =
            1.0/gMax(mesh.deltaCoeffs().primitiveField());
        scalar ePsiLinfBand6h = 0;
        scalar eGradPsiLinfBand6h = 0;
        forAll(psi0, celli)
        {
            if (mag(psi0[celli]) > 6.0*h) continue;
            ePsiLinfBand6h = Foam::max(ePsiLinfBand6h, ePsi[celli]);
            if (interiorCell[celli])
            {
                eGradPsiLinfBand6h =
                    Foam::max(eGradPsiLinfBand6h, eGradPsi[celli]);
            }
        }
        reduce(ePsiLinfBand6h, maxOp<scalar>());
        reduce(eGradPsiLinfBand6h, maxOp<scalar>());
        Info << "ePsiLinf (6h band): " << ePsiLinfBand6h
             << ", (global, reported only): " << ePsiLinf << endl;
        Info << "eGradPsiLinf (6h band, interior): " << eGradPsiLinfBand6h
             << endl;
        if (ePsiLinfBand6h > 0.5*h)
        {
            FatalErrorInFunction
                << redistType << ": planar 6h-band distance off by more "
                << "than h/2: Linf = " << ePsiLinfBand6h << " (h = " << h
                << ")"
                << abort(FatalError);
        }
        if (eGradPsiLinfBand6h > 0.15)
        {
            FatalErrorInFunction
                << redistType << ": planar 6h-band |grad psi| deviates by "
                << "more than 15%: Linf = " << eGradPsiLinfBand6h
                << abort(FatalError);
        }
    }
    else if (redistType == "PDE")
    {
        // The explicit reinitialization preserves a planar signed distance
        // exactly where its gradient is exact -- the INTERIOR. Boundary
        // cells drift by niterations*deltaT*(boundary-cell gradient error
        // of the grad(psi) scheme), bounded here by one cell size.
        const scalar h =
            1.0/gMax(mesh.deltaCoeffs().primitiveField());
        Info << "ePsiLinf (global): " << ePsiLinf << endl;
        if (eGradPsiLinf > 1e-11)
        {
            FatalErrorInFunction
                << "PDE: planar |grad psi| not preserved in the interior: "
                << eGradPsiLinf
                << abort(FatalError);
        }
        if (ePsiLinf > h)
        {
            FatalErrorInFunction
                << "PDE: planar signed distance drifted by more than h: "
                << ePsiLinf << " (h = " << h << ")"
                << abort(FatalError);
        }
    }
    else
    {
        if (ePsiLinf > 1e-13)
        {
            FatalErrorInFunction
                << " redistancing a planar signed distance fails: \n"
                << "Foam::max(Foam::mag(psi - psi0)) > 16*SMALL :"
                << ePsiLinf
                << abort(FatalError);
        }

        if (eGradPsiLinf > 1e-11)
        {
            FatalErrorInFunction
                << " redistancing a planar signed distance fails: \n"
                << "Foam::max(Foam::mag(Foam::mag(fvc::grad(psi)) - 1) > 16*SMALL) :"
                << eGradPsiLinf
                << abort(FatalError);
        }
    }
    
    // TODO: Only if test fails.
    Info << "Writing fields" << endl;

    psi0.write();
    magGradPsi0.write();

    psi.write();
    magGradPsi.write();

    ePsi.write();
    eGradPsi0.write();
    eGradPsi.write();

    Info << "Done." << endl;
};

void phaseIndicatorIsBounded(autoPtr<phaseIndicator>& phaseInd, const fvMesh& mesh)
{
    const Time& runTime = mesh.time();

    // TEST: The phase-indicator is bounded [0,1] 
    volScalarField psi 
    (
        IOobject
        (
            "psi", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimLength, 0),
        "calculated"
    );

    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0),
        "calculated"
    );

    phaseInd->calcPhaseIndicator(alpha, psi);
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Level set equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    redistancedPlaneIsUnchanged(redist, nBand, mesh, psi);

    phaseIndicatorIsBounded(phaseInd, mesh);
    
    return 0;
}


// ************************************************************************* //
