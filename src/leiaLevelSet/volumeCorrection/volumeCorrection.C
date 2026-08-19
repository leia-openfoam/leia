/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
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

#include "volumeCorrection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(volumeCorrection, false);
defineRunTimeSelectionTable(volumeCorrection, Mesh);
addToRunTimeSelectionTable(volumeCorrection, volumeCorrection, Mesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volumeCorrection::volumeCorrection(const fvMesh& mesh)
    :
        fvSolution_(mesh),
        levelSetDict_(fvSolution_.subDict("levelSet")),
        // subOrEmptyDict, NOT subDict: no case in this repository has a
        // volumeCorrection sub-dictionary, and the mandatory lookup would
        // abort every one of them at start-up.
        volCorrDict_(levelSetDict_.subOrEmptyDict("volumeCorrection")),
        targetVolume_(-1),
        lastShift_(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<Foam::volumeCorrection> volumeCorrection::New(const fvMesh& mesh)
{
    const fvSolution& fvSolution (mesh);
    const dictionary& levelSetDict = fvSolution.subDict("levelSet");

    // MUST be subOrEmptyDict (a COPY, held only for the duration of this
    // function): the absent sub-dictionary is the DEFAULT configuration of
    // this hierarchy, not an error. Same construction as narrowBand.C:48.
    const dictionary volCorrDict =
        levelSetDict.subOrEmptyDict("volumeCorrection");

    const word& modelType =
        volCorrDict.getOrDefault<word>("type", "noVolumeCorrection");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = MeshConstructorTable(modelType);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            fvSolution,
            "volumeCorrection",
            modelType,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Info<< "Selecting volumeCorrection type: " << modelType << nl << endl;

    // Construct the model and return the autoPtr to the object.
    return autoPtr<volumeCorrection>(ctorPtr(mesh));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void volumeCorrection::setTargetVolume(const volScalarField& alpha)
{
    // Idempotent: the target is captured exactly once per run. See the
    // declaration for why a second capture on a restart is a silent
    // re-targeting onto the already accumulated error.
    if (targetVolume_ >= 0)
    {
        return;
    }

    if (volCorrDict_.found("targetVolume"))
    {
        targetVolume_ = volCorrDict_.get<scalar>("targetVolume");

        Info<< "volumeCorrection: target volume V_target = " << targetVolume_
            << " (explicit levelSet/volumeCorrection/targetVolume entry)"
            << nl << endl;

        return;
    }

    // gSum over primitiveField()*mesh.V().field(): the SAME expression as
    // advectionErrors.H:45,57 and leiaLevelSetFoam.C:167, so the target and
    // the reported E_VOL_ALPHA_SIGNED baseline are the same number and not
    // merely two spellings that agree to round-off.
    targetVolume_ =
        gSum(alpha.primitiveField()*alpha.mesh().V().field());

    Info<< "volumeCorrection: target volume V_target = " << targetVolume_
        << " (captured from " << alpha.name() << " at t = "
        << alpha.mesh().time().timeOutputValue() << ")" << nl << endl;

    if (targetVolume_ <= 0)
    {
        // Phase 1 is empty (or alpha is not an indicator at all). Every
        // relative convergence test in this hierarchy divides by
        // targetVolume_, so state it here rather than let a corrector
        // produce inf/NaN inside a Newton loop.
        WarningInFunction
            << "captured target volume is " << targetVolume_
            << " <= 0: there is no phase-1 volume to conserve, and every "
            << "volumeCorrection model will refuse to act." << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
