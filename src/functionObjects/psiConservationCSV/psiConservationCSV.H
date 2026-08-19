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
    Foam::functionObjects::psiConservationCSV

Group
    grpFieldFunctionObjects

Description
    Per-time-step conservation diagnostic of the level-set field psi. Writes one
    row per call into "psiConservation.csv" (opened in APPEND mode, rank 0
    writes) holding the global integral int(psi dV), its drift, a term-by-term
    budget of the discrete psi transport equation, and -- in the same row -- the
    phase volume vol{psi < 0}.

    WHY BOTH INTEGRALS ARE EMITTED, AND THE RETRACTION THAT FORCED IT.
    int(psi dV) and vol{psi < 0} = int(H(-psi) dV) are DIFFERENT functionals of
    the same field, and a scheme can conserve either one while destroying the
    other. A pure flux form
        d/dt int(psi dV) = -int_dOmega psi (v.n) dS
    conserves int(psi dV) to machine precision on a closed domain (see the
    telescoping argument below) and says NOTHING about the zero level set, which
    is where the phase volume lives: psi may drift by an O(1) constant, or tilt,
    with int(psi dV) unchanged, and the isosurface psi = 0 then sweeps through
    cells. The divergence-form SDPLS source sdplsRdiv (src/leiaLevelSet/
    sdplsSource/sdplsRdiv.C, commit c54efc8) was introduced with exactly that
    conflation in its rationale -- "the failure is a CONSERVATION failure, so the
    fix is a form whose fluxes cannot drain", supported by one N=32 run where the
    relative volume error dipped to -1.6% and recovered to -0.6%. That claim is
    RETRACTED: at N=64 the same form drained the droplet, and at N=128 it
    DIVERGED, relative volume error +16.9 at t = 0.031 s. Conservation of
    int(psi dV) was never evidence about vol{psi < 0}, and no single-column CSV
    could have shown the difference. This object therefore emits both, every
    step, in the same row: PSI_INT (= SUM_PSI_V) next to VOL_PSI_NEG, each with
    its own drift column.

    THE DISCRETE BUDGET THAT THE BUDGET_* COLUMNS DECOMPOSE.
    The Eulerian solver assembles (applications/solvers/leiaLevelSetTwoPhaseFoam/
    alphaEqn.H:20-29) the ADVECTIVE (Hamilton-Jacobi) form

        fvm::ddt(psi) + fvm::div(phi, psi) - fvm::Sp(divPhi, psi) == S,
        divPhi = fvc::div(phi).

    Multiplying the cell equation by V_c and summing over all cells c of the
    global mesh, with Euler ddt and every operator taken at the new time level:

        sum_c (psi_c^{n+1} - psi_c^n) V_c
          = -dt sum_c sum_{f in c} phi_f psi_f            (flux)
            +dt sum_c (div phi)_c psi_c^{n+1} V_c         (Sp volume term)
            +dt sum_c S_c V_c                             (source)

    Term 1 TELESCOPES. Every interior face f is owned by one cell and neighboured
    by another; fvm::div contributes +phi_f psi_f to the owner row and
    -phi_f psi_f to the neighbour row (the outward face-area vector flips sign),
    so the pair cancels exactly in the sum. A processor face is an interior face
    of the GLOBAL mesh -- its two halves live on different ranks and cancel after
    the parallel reduction -- which is why BUDGET_FLUX skips every
    patch().coupled() patch and sums only the non-coupled boundary. What remains,
    -dt sum_{boundary f} phi_f psi_f, is the ENTIRE non-conservation of the flux
    term, at ANY time order and for ANY div scheme (upwind, linearUpwind,
    limited): the cancellation is algebraic, not an accuracy statement. On the
    coupled droplet case (cases/stationaryDroplet2D) the walls are noSlip and
    frontAndBack is empty, so phi_f = 0 on every boundary face with faces at all
    and BUDGET_FLUX is identically zero -- any measured drift there is the Sp
    term, the source, or the residual, never the flux.

    Term 2 IS A VOLUME SOURCE, NOT A FLUX. -fvm::Sp(divPhi, psi) contributes
    +dt sum_c (div phi)_c psi_c V_c, which does not telescope and is bounded by

        |BUDGET_DIVPHI| <= dt * DIVPHI_LINF * int|psi| dV.

    It vanishes only for a DISCRETELY solenoidal flux, sum_{f in c} phi_f = 0.
    The intermediate PIMPLE fluxes are not solenoidal until the pressure
    corrector has converged, and that non-solenoidality is precisely why
    alphaEqn.H uses the advective form in the first place (the conservative form
    alone equals it only when div phi = 0). DIVPHI_LINF and DIVPHI_L1/L2_DIV_PHI
    measure that violation directly, in 1/s, and SUM_SP_DIVPHI_PSI reports the
    INSTANTANEOUS RATE sum_c (div phi)_c psi_c V_c (BUDGET_DIVPHI is the same
    quantity multiplied by the measured interval DELTA_T).

    Term 3, the SDPLS source. sdplsSource registers the coefficient f_nl(psi^n)
    as the volScalarField "SDPLS_nonLinearPart" (src/leiaLevelSet/sdplsSource/
    sdplsSource.C:117-127, dimensions 1/s, AUTO_WRITE, constructed even for
    noSource). The assembled per-cell contribution is V_c (Sc_c + Sp_c psi_c)
    (src/leiaLevelSet/sdplsSource/discretization/discretization.C:104-105 with
    the sign convention documented at :83-103), so
        simpleLinearImplicit:              Sc = 0, Sp = f_nl  -> V f_nl psi^{n+1}
        strictNegativeSpLinearImplicit:    V (max(f,0) psi^{iter} + min(f,0) psi^{n+1})
        explicit:                          V f_nl psi^{iter}
    and BUDGET_SOURCE = dt * gSum(f_nl psi V) is the exact source integral for
    the first and equals it up to the psi iterate for the other two.

    THE TWO NaN CASES, AND WHY THEY MUST BE NaN AND NOT A NUMBER.
      (a) MOLLIFIER. sdplsSource.C:139-148 multiplies nonLinearPart_ by
          mollifier_->field(psi) INSIDE fvmsdplsSource; the field left on the
          registry is the UNmollified f_nl. gSum(f_nl psi V) would then overstate
          the source by the cut-off it never applied.
      (b) THE DIVERGENCE-FORM FAMILY: Rdiv and RdivStrictSp. sdplsRdiv
          overrides fvmsdplsSource and assembles
          fvm::div(phiW, psi) - fvm::div(phi, psi) - fvm::Sp(fvc::div(phiW) - a, psi)
          with w = (n.v) n and phiW = fvc::flux(w); sdplsRdivStrictSp assembles
          the same two div terms with the coefficient f = a - div(w) sign-split,
          fvm::Su(max(f,0) psi^n) + fvm::Sp(min(f,0)). Both still fill
          SDPLS_nonLinearPart (= a, inherited from sdplsR), so the algebraic
          formula would report a number that is NOT the source that ran, and
          phiW is a local that is never registered, so this object cannot rebuild
          it without linking libleiaLevelSet.
    In both regimes BUDGET_SOURCE, BUDGET_SOURCE_BAND and (because the budget is
    then unattributable) BUDGET_RESIDUAL are quiet_NaN, and SOURCE_ALGEBRAIC
    records which regime the row is in: 1 = gSum(f_nl psi V) reproduces the
    assembled source, 0 = it does not. The two flags are read ONCE from
    fvSolution.levelSet.sdplsSource in the constructor, by the same dictionary
    route sdplsSource::New uses (sdplsSource.C:47-51), so no leia header and no
    extra link library are needed.

    WHAT BUDGET_RESIDUAL CONTAINS.
        BUDGET_RESIDUAL = D_PSI_INT - (BUDGET_FLUX + BUDGET_DIVPHI + BUDGET_SOURCE)
    is zero to the linear-solver tolerance ONLY when every assumption above
    holds. It absorbs, in order of expected size on the coupled droplet:
      - the flux mismatch. alphaEqn.H is included inside pimple.loop() and psi is
        re-solved on EVERY outer corrector with that corrector's phi, while pEqn
        updates phi afterwards -- so the registered phi read here is the LAST
        pressure-corrected flux of the step, not the flux that assembled the psi
        equation. With nOuterCorrectors > 1 the row also closes over several
        re-solves of the same step.
      - the time scheme. The identity above is exact for ddt(psi) Euler. It is
        NOT exact for CrankNicolson 1 (token PSI_DDT, cases/stationaryDroplet2D/
        system/aliasDict), where the ddt operator spans two old-time levels and
        the one-interval difference D_PSI_INT no longer matches the assembled
        temporal operator.
      - the linear-solver tolerance on psiEqn.solve().
    A residual much larger than DELTA_T * DIVPHI_LINF * PSI_ABS_INT is therefore
    a statement about the flux/time-scheme confound, not about the source.

    WHY PSI_ABS_INT IS THE NORMALISER. int(psi dV) is NOT a physical invariant:
    psi is a SIGNED DISTANCE, so the integral is dominated by the far field
    (|psi| ~ L away from the interface, |psi| ~ h next to it), can pass through
    zero for a symmetric configuration, and changes discontinuously whenever the
    redistancer fires. Relative columns are therefore formed with
    max(PSI_ABS_INT, VSMALL) = max(int|psi| dV, VSMALL), which cannot pass
    through zero, and the interface-local statement is made by the *_BAND
    columns, which restrict the same integrals to the registered "NarrowBand"
    mask (nbMask_c = NarrowBand_c > 0.5). SUM_PSI_V_REL is kept on the
    literally-requested denominator |int psi dV|(t=0) and is the column to
    distrust when that reference is near zero.

    UNDER THE SEMI-LAGRANGIAN SOLVER there is no fvm::div(phi, psi) and no
    linear solve for psi at all: psi^{n+1} is reconstructed at departure feet.
    D_PSI_INT, PSI_INT, VOL_PSI_NEG and the DIVPHI_* columns remain honest
    measurements (the first three of the field, the last of the flux the momentum
    equation carries), but BUDGET_FLUX/BUDGET_DIVPHI/BUDGET_SOURCE/
    BUDGET_RESIDUAL then decompose an equation that solver never assembles, and
    must be read as such.

    ABSENT FIELDS NEVER ABORT. phi, the SDPLS coefficient SDPLS_nonLinearPart,
    the NarrowBand mask and the phase indicator are all reached through the TYPED
    registry find mesh.cfindObject<Type>(name), which returns a null pointer
    instead of aborting -- typed and not the untyped objectRegistry::found for the
    reason recorded at gradPsiErrorCSV.C:197-201 (a name may be registered under
    another type, and the typed lookup that follows an untyped hit aborts).
    Every column a missing field feeds is quiet_NaN -- never 0.0, which a
    downstream "> 0" filter would silently keep and a plot would read as a
    perfect budget (the same reasoning as gradPsiErrorCSV.C:162-171).

    COLUMNS (psiConservation.csv, in order)
    \verbatim
      TIME                 runTime.timeOutputValue() of the completed step
      DELTA_T              t - t(previous row); MEASURED, not deltaTValue()
      PSI_INT              gSum(psi V)                       [== SUM_PSI_V]
      PSI_ABS_INT          gSum(|psi| V), the sign-free normaliser
      D_PSI_INT            PSI_INT - PSI_INT(previous row)
      D_PSI_INT_REL        D_PSI_INT / max(PSI_ABS_INT, VSMALL)
      BUDGET_FLUX          -DELTA_T sum_{non-coupled boundary f} phi_f psi_f
      BUDGET_DIVPHI        +DELTA_T gSum(divPhi psi V)
      BUDGET_SOURCE        +DELTA_T gSum(f_nl psi V), NaN if not algebraic
      BUDGET_RESIDUAL      D_PSI_INT - (FLUX + DIVPHI + SOURCE)
      BUDGET_RESIDUAL_REL  BUDGET_RESIDUAL / max(PSI_ABS_INT, VSMALL)
      DIVPHI_LINF          gMax(|divPhi|) [1/s]              [== MAX_DIV_PHI]
      DIVPHI_L1            gSum(|divPhi| V) / gSum(V) [1/s]
      PSI_INT_BAND         gSum(nbMask psi V)
      BUDGET_DIVPHI_BAND   +DELTA_T gSum(nbMask divPhi psi V)
      BUDGET_SOURCE_BAND   +DELTA_T gSum(nbMask f_nl psi V)
      SOURCE_ALGEBRAIC     1 = gSum(f_nl psi V) is the assembled source, else 0
      SUM_PSI_V_DRIFT      PSI_INT - PSI_INT(first row of this run)
      SUM_PSI_V_REL        SUM_PSI_V_DRIFT / max(|PSI_INT(first row)|, VSMALL)
      SUM_SP_DIVPHI_PSI    gSum(divPhi psi V) [1/s], the RATE of the Sp term
      L2_DIV_PHI           sqrt(gSum(divPhi^2 V) / gSum(V)) [1/s]
      VOL_PSI_NEG          gSum(alpha V), the phase volume the studies use
      VOL_PSI_NEG_REL      (VOL_PSI_NEG - VOL_PSI_NEG(first row)) / VOL_PSI_NEG(first row)
      VOL_PSI_NEG_SHARP    gSum([psi_c < 0] V), the cell-centre-sign volume
    \endverbatim
    SUM_PSI_V and MAX_DIV_PHI are NOT emitted as separate columns: they are
    PSI_INT and DIVPHI_LINF respectively, and a second copy of the same reduction
    in every row would be one more thing to keep consistent, not one more
    measurement.

    VOL_PSI_NEG uses the SAME measure the studies already report -- gSum(alpha V)
    on the registered phase indicator (applications/solvers/
    leiaLevelSetTwoPhaseFoam/advectionErrors.H:45,57 computes the volume error
    from exactly this), field name from the "alpha" entry, default alpha.water.
    That indicator is the SMEARED heaviside for PHASE_INDICATOR heaviside
    (src/leiaLevelSet/phaseIndicator/heavisidePhaseIndicator.C: epsilon =
    (nCells/2) max(1/deltaCoeffs), nCells = 2 by default), so it is NOT the sharp
    volume of {psi < 0}; VOL_PSI_NEG_SHARP is emitted beside it as the
    indicator-free counterpart, the cell-centre sign count sum_c [psi_c < 0] V_c,
    first-order in h and independent of PHASE_INDICATOR. Reading them together
    separates a real interface displacement (both move) from a change in the
    indicator's smearing width (only VOL_PSI_NEG moves).

    The t=0 reference for SUM_PSI_V_DRIFT / VOL_PSI_NEG_REL is the FIRST row this
    object writes, which the constructor emits. On a restart (startFrom
    latestTime) that row is the restart time, not 0 -- the drift columns are then
    drift since the restart, and the workflow deletes psiConservation.csv in
    rule preprocess so a re-run never appends to a stale reference.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    psiConservationCSV
    {
        type            psiConservationCSV;
        libs            ( leiaFunctionObjects );
        writeControl    timeStep;
        writeInterval   1;

        // Optional, shown with their defaults:
        // psi     psi;          // level-set field
        // phi     phi;          // mass/volumetric flux; NaN columns if absent
        // alpha   alpha.water;  // phase indicator for VOL_PSI_NEG
        // enabled true;         // false = never constructed, no CSV at all
    }
    \endverbatim

    Only scheme-free operators are used -- fvc::div(surfaceScalarField) is
    fvc::surfaceIntegrate and consults NO fvSchemes entry -- so the object runs
    unchanged on cases whose fvSchemes says \c divSchemes { default none; }
    (cases/stationaryDroplet2D/system/fvSchemes.template).

    Usage by the \c postProcess utility is not available: every interval column
    is defined against the PREVIOUS row this object wrote, which requires the
    object to have run over the time loop.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::gradPsiErrorCSV

SourceFiles
    psiConservationCSV.C

\*---------------------------------------------------------------------------*/

#ifndef psiConservationCSV_H
#define psiConservationCSV_H

#include "fvCFD.H"
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class psiConservationCSV Declaration
\*---------------------------------------------------------------------------*/

class psiConservationCSV
:
    public fvMeshFunctionObject
{
    //- The level set. A reference: without psi this object has nothing to
    //  measure, so a missing psi is a configuration error and lookupObject's
    //  abort is the correct report.
    const volScalarField& psi_;

    //- The flux, as a POINTER and not a reference. gradPsiErrorCSV can bind psi
    //  by reference because it only ever runs on level-set cases, but this
    //  object is meant to be switched on wherever int(psi dV) is in question,
    //  including runs that register no flux at all. A missing phi NaNs the
    //  BUDGET_FLUX / BUDGET_DIVPHI / DIVPHI_* / SUM_SP_DIVPHI_PSI columns and
    //  leaves PSI_INT, D_PSI_INT and the volume columns intact.
    const surfaceScalarField* phiPtr_;

    //- Phase-indicator field name for VOL_PSI_NEG (default alpha.water; the 3D
    //  cases name it alpha, cf. workflow/Snakefile SETFIELDS_ARGS).
    const word alphaName_;

    OFstream csvFile_;

    //- Previous-row state, so the drift is measured over the interval THIS row
    //  closes and never depends on GeometricField oldTime bookkeeping. psi
    //  .oldTime() would be wrong here for two independent reasons: the object
    //  may be written on a writeControl coarser than timeStep (the interval is
    //  then several steps), and the ddt rollover for the next step has not
    //  happened yet when function objects execute at the top of Time::run().
    scalar prevPsiInt_;
    scalar prevTime_;
    bool   havePrev_;

    //- First-row reference for the SUM_PSI_V_DRIFT / VOL_PSI_NEG_REL columns.
    //  Set on the first write() (the one the constructor issues), i.e. at
    //  startTime, which on a restart is the restart time and not 0.
    scalar psiInt0_;
    scalar volNeg0_;
    bool   haveRef_;

    //- True when gSum(f_nl psi V) reproduces the ASSEMBLED source (algebraic
    //  discretization, no mollifier); read once from fvSolution in the ctor,
    //  because fvSolution.levelSet.sdplsSource cannot change during a run
    //  (runTimeModifiable no in every study controlDict).
    bool   algebraicSource_;

    //- The CSV column schema, without the trailing newline. Single definition:
    //  written for a new file, and compared against the first line of an
    //  existing one so a file written by an older build cannot be appended to
    //  with a different number of fields.
    static const std::string& csvHeader();

public:

    //- Runtime type information
    TypeName("psiConservationCSV");


    // Constructors

        //- Construct from Time and dictionary
        psiConservationCSV
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        psiConservationCSV(const psiConservationCSV&) = delete;

        //- No copy assignment
        void operator=(const psiConservationCSV&) = delete;

        virtual bool execute()
        {
            return true;
        }

        virtual bool write();


    //- Destructor
    virtual ~psiConservationCSV() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
