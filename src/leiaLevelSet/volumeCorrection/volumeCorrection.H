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

Class
    Foam::volumeCorrection

Description
    Abstract base class opening the model hierarchy for corrections that
    restore the phase volume of the level set field \f$ \psi \f$ after
    advection, and at the same time the DEFAULT model \c noVolumeCorrection,
    whose correct() is a no-op.

    The conserved functional is the one the verification metrics report, and
    there is only one of it in this library
    (leiaLevelSetFoam/advectionErrors.H:45,57):

    \f{equation}{
        V(\psi) \;=\; \sum_c \alpha_c(\psi)\,|\Omega_c|,
        \qquad
        \alpha_c(\psi) \;=\;
        \frac{|\Omega_c \cap \{\psi < 0\}|}{|\Omega_c|},
    \f}

    with \f$ \alpha \f$ produced by the RUN-TIME SELECTED Foam::phaseIndicator
    (\c geometric by default), not by any indicator a corrector rolls itself.
    That is a hard requirement, not a convenience: the reported volume error
    \c E_VOL_ALPHA_SIGNED is \f$ V(\psi^n) - V(\psi^0) \f$ evaluated through
    the selected indicator, so a corrector that drives a DIFFERENT functional
    to zero -- an inline smoothed Heaviside, a hand-rolled cut-cell integral --
    leaves a residual in the reported column that it believes it removed.
    Hence Foam::volumeCorrection::correct() is handed the live
    Foam::phaseIndicator and evaluates \f$ V \f$ through it.

    Orientation, verified on every indicator in this library
    (geometricPhaseIndicator.C:101-107, detrixheAslamPhaseIndicator.C:152-159,
    heavisidePhaseIndicator.C:62-76): \f$ \alpha = 1 \f$ in
    \f$ \Omega^- = \{\psi < 0\} \f$. No sign flip is needed anywhere in this
    hierarchy.

    Selected at run time from the \c levelSet dictionary of \c system/fvSolution:

    \verbatim
        levelSet
        {
            volumeCorrection
            {
                type    newtonShift;
            }
        }
    \endverbatim

    The sub-dictionary is looked up with dictionary::subOrEmptyDict, so a
    \c levelSet dictionary WITHOUT a \c volumeCorrection entry -- which is
    every case in \c cases/ and every study the workflow generates -- selects
    \c noVolumeCorrection, whose correct() returns false and never touches
    \f$ \psi \f$. Existing studies run bit-identically.

    NOTE on the stored dictionary: dictionary::subOrEmptyDict returns a
    dictionary BY VALUE in OpenFOAM-v2512 (dictionary.H:836, "Find and return
    a sub-dictionary as a copy"). \c volCorrDict_ is therefore a VALUE member.
    A \c const \c dictionary& member -- the spelling used by
    Foam::redistancer and Foam::phaseIndicator, which are legal only because
    both use the mandatory subDict() -- would bind to a temporary and dangle
    at the first getOrDefault().

SourceFiles
    volumeCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef volumeCorrection_H
#define volumeCorrection_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvSolution.H"
#include "volFields.H"
#include "fvMesh.H"
#include "phaseIndicator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class volumeCorrection Declaration
\*---------------------------------------------------------------------------*/

class volumeCorrection
{
protected:

    const fvSolution& fvSolution_;
    const dictionary& levelSetDict_;

    //- BY VALUE. dictionary::subOrEmptyDict returns a COPY in v2512
    //  (dictionary.H:836); a reference member would dangle. See the class
    //  Description.
    const dictionary volCorrDict_;

    //- \f$ V_{target} = \sum_c \alpha_c |\Omega_c| \f$, the volume every
    //  correct() call restores. Negative means "not yet captured": it is set
    //  once, by setTargetVolume(), from the alpha field the solver has just
    //  built at start-up, so that it is bit-identical to the \c alphaV0
    //  baseline of the error metric (advectionErrors.H:45).
    scalar targetVolume_;

    //- The shift \f$ \varepsilon \f$ applied by the last correct() call
    //  (identically 0 for the no-op base). Kept so that a caller -- a CSV
    //  column, a diagnostic -- can report the correction the run actually
    //  needed, which is the only honest way to read a corrected volume error:
    //  a volume error of 0 with \f$ \varepsilon \f$ growing every step is a
    //  failing scheme with its symptom removed, not a converging one.
    scalar lastShift_;


public:

    // Static Data Members

        TypeName("noVolumeCorrection");


    declareRunTimeSelectionTable
    (
        autoPtr,
        volumeCorrection,
        Mesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from components
        volumeCorrection(const fvMesh& mesh);


    // Selectors

        //- Select from fvSolution levelSet.volumeCorrection.type,
        //  defaulting to the no-op noVolumeCorrection
        static autoPtr<volumeCorrection> New(const fvMesh& mesh);


    //- Destructor
    virtual ~volumeCorrection() = default;


    // Member Functions

        //- Capture \f$ V_{target} = \sum_c \alpha_c |\Omega_c| \f$ from the
        //  supplied indicator field, using gSum over the primitive field
        //  times mesh.V() -- the SAME expression as advectionErrors.H:45,57,
        //  so the target and the reported error are the same number.
        //
        //  Honours an explicit \c targetVolume dictionary entry, which is
        //  needed whenever the error metric's baseline is NOT the t = 0
        //  volume: errorCalculation.H:20-33 reads an optional \c alphaEnd
        //  field as \c alpha0, and a case shipping one has a baseline this
        //  routine cannot see.
        //
        //  IDEMPOTENT: a second call returns immediately. A restart from
        //  t > 0 would otherwise silently re-target onto the volume the run
        //  had already drifted to, and the corrector would then defend the
        //  accumulated error instead of removing it. The captured value is
        //  Info-logged unconditionally -- including by the no-op base, where
        //  it costs one global reduction at start-up -- so that every run
        //  states in its log which volume it enforced.
        virtual void setTargetVolume(const volScalarField& alpha);

        //- The captured target volume (< 0 before setTargetVolume())
        scalar targetVolume() const
        {
            return targetVolume_;
        }

        //- The shift applied by the last correct() call
        scalar lastShift() const
        {
            return lastShift_;
        }

        //- Apply a UNIFORM scalar shift to psi restoring targetVolume().
        //
        //  Returns true if and only if psi was modified. Callers use the
        //  return value exactly as they use redistancer::redistance() at
        //  leiaLevelSetTwoPhaseFoam/alphaEqn.H:50 -- to refresh the narrow
        //  band, since the zero crossing may have moved across a cell
        //  centre.
        //
        //  Base implementation: no-op returning false. This is what every
        //  case without a \c volumeCorrection sub-dictionary gets, and it is
        //  why adding this hierarchy changes no existing result.
        virtual bool correct
        (
            volScalarField& psi,
            phaseIndicator& phaseInd
        )
        {
            lastShift_ = 0;
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
