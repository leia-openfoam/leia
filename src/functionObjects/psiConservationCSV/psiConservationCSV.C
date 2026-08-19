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


#include "psiConservationCSV.H"
#include "addToRunTimeSelectionTable.H"
#include "fvSolution.H"

#include <limits>
#include <fstream>
#include <string>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(psiConservationCSV, 0);
    addToRunTimeSelectionTable(functionObject, psiConservationCSV, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::psiConservationCSV::psiConservationCSV
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject (name, runTime, dict),
    psi_(mesh_.lookupObject<volScalarField>(dict.getOrDefault<word>("psi", "psi"))),
    // POINTER, not lookupObject: a kinematic or semi-Lagrangian run need not
    // register a flux, and the contract for this object is that a missing input
    // NaNs its columns instead of aborting the solver. cfindObject is the typed
    // find -- an untyped hit followed by a typed lookup aborts on a name reused
    // for another type (gradPsiErrorCSV.C:197-201 records that failure for `R`,
    // which OpenFOAM's turbulence models register as a symmTensor).
    phiPtr_(mesh_.cfindObject<surfaceScalarField>(dict.getOrDefault<word>("phi", "phi"))),
    alphaName_(dict.getOrDefault<word>("alpha", "alpha.water")),
    csvFile_(fileName("psiConservation.csv"), IOstreamOption(), IOstreamOption::appendType::APPEND),
    prevPsiInt_(0),
    prevTime_(0),
    havePrev_(false),
    psiInt0_(0),
    volNeg0_(0),
    haveRef_(false),
    algebraicSource_(true)
{
    // 16 significant digits, NOT the controlDict writePrecision (8 in the
    // droplet cases, 6 by OpenFOAM default). Every quantity in this row except
    // PSI_INT and VOL_PSI_NEG is a DIFFERENCE of numbers of order PSI_INT: on
    // cases/stationaryDroplet2D, int|psi| dV = O(1e-9) m^4 while a converged
    // one-step budget residual is many orders below it. Printing the integrals
    // at 6-8 significant digits would make the drift unrecoverable from the file
    // by any consumer that recomputes differences between rows, and would quantise
    // BUDGET_RESIDUAL to the last printed digit of PSI_INT.
    //
    // This is OSstream::precision(int) -- PER STREAM. It is deliberately NOT
    // IOstream::defaultPrecision(unsigned), which is static and would re-set the
    // precision of every subsequent field write in the run: that would change the
    // bytes of every written time directory the moment this diagnostic is switched
    // on, and a diagnostic that alters the data it observes is worthless.
    csvFile_.precision(16);

    if (Pstream::myProcNo() == 0)
    {
        const fileName csvName("psiConservation.csv");
        if (!fileSize(csvName))
        {
            csvFile_ << csvHeader().c_str() << "\n";
        }
        else
        {
            // The file is opened in APPEND mode and the header defines the
            // column schema, so a file written by an older build silently
            // receives rows with MORE fields than its header names -- every
            // column read by name downstream then refers to the wrong data.
            // This is reachable without a full re-materialisation: `rule solve`
            // removes only <solver>.csv, while `rm -f psiConservation.csv` lives
            // in `rule preprocess`, whose marker may still be up to date. Fail
            // loudly; the file is cheap to delete and impossible to repair.
            std::ifstream existing(csvName);
            std::string firstLine;
            std::getline(existing, firstLine);
            while (!firstLine.empty()
                && (firstLine.back() == '\r' || firstLine.back() == '\n'))
            {
                firstLine.pop_back();
            }
            if (firstLine != csvHeader())
            {
                FatalErrorInFunction
                    << "Existing " << csvName << " was written with a different"
                    << " column schema and is being appended to." << nl
                    << "  found  : " << firstLine.c_str() << nl
                    << "  expect : " << csvHeader().c_str() << nl
                    << "Appending would misalign every column read by name."
                    << " Delete the file (or re-run the case from a clean"
                    << " directory) and start again."
                    << exit(FatalError);
            }
        }
    }

    // Which SDPLS regime is this run in? Read ONCE, by the same dictionary route
    // sdplsSource::New uses (src/leiaLevelSet/sdplsSource/sdplsSource.C:47-51),
    // so this library needs no leia header and no extra link library -- the
    // alternative, holding the sdplsSource object, is impossible anyway: it is a
    // plain autoPtr local in the solver's createFields.H:153 and does not derive
    // from regIOobject, so it is not on the registry.
    //
    // subOrEmptyDict and not subDict for `levelSet`: this object must also be
    // usable on a case whose fvSolution has no levelSet dictionary at all (an
    // interFoam/interFlow comparison arm), where the answer is simply "no SDPLS
    // source", not a FatalIOError from a diagnostic.
    //
    // gSum(f_nl psi V) is the assembled source ONLY when the source is algebraic
    // and unmollified:
    //   - a mollifier multiplies nonLinearPart_ INSIDE fvmsdplsSource
    //     (sdplsSource.C:139-148) and the registered SDPLS_nonLinearPart stays
    //     UNmollified, so the algebraic integral overstates it by the cut-off;
    //   - the DIVERGENCE-FORM family assembles fluxes, not an algebraic
    //     coefficient, and still fills SDPLS_nonLinearPart with a, so the
    //     algebraic integral would be a number that is NOT the source that ran;
    //     phiW = fvc::flux((nHat & U)*nHat) is a local and is never registered,
    //     so it cannot be rebuilt here. Two types are in that family:
    //       `Rdiv`         (sdplsRdiv.C)          div(phiW,psi) - div(phi,psi)
    //                                             - Sp(div(phiW) - a, psi);
    //       `RdivStrictSp` (sdplsRdivStrictSp.C)  the same two div terms with
    //                                             the coefficient f = a -
    //                                             div(w) sign-split, Su for
    //                                             max(f,0) and Sp for min(f,0).
    //     BOTH must be excluded. Testing only for "Rdiv" would silently mark an
    //     RdivStrictSp run SOURCE_ALGEBRAIC 1 and publish a BUDGET_SOURCE
    //     column computed from a coefficient that never entered the matrix.
    // Both cases NaN the source and residual columns and set SOURCE_ALGEBRAIC 0.
    const fvSolution& fvSol(mesh_);
    const dictionary levelSetDict(fvSol.subOrEmptyDict("levelSet"));
    const dictionary sourceDict(levelSetDict.subOrEmptyDict("sdplsSource"));
    const word sourceType(sourceDict.getOrDefault<word>("type", "noSource"));
    const word mollifierType
    (
        sourceDict.subOrEmptyDict("mollifier").getOrDefault<word>("type", "none")
    );
    const bool divergenceFormSource =
        (sourceType == "Rdiv") || (sourceType == "RdivStrictSp");
    algebraicSource_ = !divergenceFormSource && (mollifierType == "none");

    // The t=0 row. Function objects are constructed from
    // functionObjectList::read(), which Time::run() reaches at the top of the
    // first loop iteration -- after createFields.H, which is what makes the psi
    // lookup above legal, and before the first ++runTime, so this row carries
    // the INITIAL condition and is the reference for SUM_PSI_V_DRIFT.
    write();
}


const std::string&
Foam::functionObjects::psiConservationCSV::csvHeader()
{
    // ONE definition of the column schema: it is both written and, for an
    // existing file, checked against. Keep in step with write() below.
    //
    // THE ROW DECOMPOSES ONE DISCRETE EQUATION. The Eulerian solver assembles
    // (alphaEqn.H:20-29)
    //
    //     fvm::ddt(psi) + fvm::div(phi, psi) - fvm::Sp(divPhi, psi) == S,
    //     divPhi = fvc::div(phi),
    //
    // and multiplying the cell equation by V_c and summing over all cells gives,
    // for Euler ddt with every operator at the new time level,
    //
    //     D_PSI_INT = BUDGET_FLUX + BUDGET_DIVPHI + BUDGET_SOURCE + BUDGET_RESIDUAL
    //
    // with BUDGET_RESIDUAL = 0 up to the linear-solver tolerance.
    //
    // BUDGET_FLUX IS THE WHOLE NON-CONSERVATION OF THE FLUX TERM. In
    // sum_c sum_{f in c} phi_f psi_f every interior face appears twice, +phi_f
    // psi_f in the owner row and -phi_f psi_f in the neighbour row, and cancels
    // pairwise; a processor face is an interior face of the GLOBAL mesh whose two
    // halves cancel across the reduction. Only non-coupled boundary faces survive,
    // so a pure flux form conserves int(psi dV) at ANY time order and for ANY div
    // scheme -- the cancellation is algebraic, not an accuracy statement. On
    // cases/stationaryDroplet2D (noSlip walls, empty front/back) BUDGET_FLUX is
    // identically zero and any drift there is Sp, source or residual.
    //
    // BUDGET_DIVPHI IS A VOLUME SOURCE, NOT A FLUX: +dt sum_c (div phi)_c psi_c
    // V_c does not telescope, is bounded by DELTA_T * DIVPHI_LINF * PSI_ABS_INT,
    // and vanishes only for a DISCRETELY solenoidal flux (sum_{f in c} phi_f = 0).
    // The intermediate PIMPLE fluxes are not solenoidal -- which is exactly why
    // alphaEqn.H:8-19 transports psi in the advective form instead of the
    // conservative one. SUM_SP_DIVPHI_PSI is the same quantity as an
    // instantaneous RATE, gSum(divPhi psi V) [1/s], defined on every row
    // including the first (where DELTA_T, and therefore BUDGET_DIVPHI, is not).
    //
    // BUDGET_SOURCE = dt gSum(f_nl psi V) from the registered
    // SDPLS_nonLinearPart. The assembled contribution is V (Sc + Sp psi)
    // (discretization.C:104-105), i.e. exactly V f_nl psi^{n+1} for
    // simpleLinearImplicit (Sc = 0, Sp = f_nl) and the same up to the psi iterate
    // for strictNegativeSpLinearImplicit and explicit. It is quiet_NaN, together
    // with BUDGET_SOURCE_BAND and BUDGET_RESIDUAL, whenever a mollifier is active
    // (sdplsSource.C:139-148 mollifies inside the assembly and leaves the
    // registered field UNmollified) or the source is one of the divergence-form
    // family, Rdiv (sdplsRdiv.C) or RdivStrictSp (sdplsRdivStrictSp.C), which
    // assemble fluxes whose phiW is never registered. SOURCE_ALGEBRAIC records
    // which regime the row is in.
    //
    // int(psi dV) IS NOT A PHYSICAL INVARIANT: psi is a signed distance, so the
    // far field (|psi| ~ L) dominates the integral, it can pass through zero for
    // a symmetric configuration, and it jumps whenever the redistancer fires.
    // PSI_ABS_INT = int|psi| dV is therefore the normaliser for the *_REL
    // columns, and the *_BAND columns restrict the same integrals to the
    // registered NarrowBand mask to make the interface-local statement.
    //
    // VOL_PSI_NEG IS A DIFFERENT FUNCTIONAL AND IS EMITTED BECAUSE CONFLATING THE
    // TWO PRODUCED A RETRACTED CLAIM. A flux form conserves int(psi dV) exactly
    // while the zero level set -- and with it the phase volume -- moves freely:
    // adding a constant to psi leaves neither invariant equal. The divergence-form
    // source sdplsRdiv was introduced (commit c54efc8) on the argument that its
    // fluxes "cannot drain" the droplet, evidenced by one N=32 run whose relative
    // volume error dipped to -1.6% and recovered to -0.6%. RETRACTED: at N=64 the
    // same form drained the droplet and at N=128 it DIVERGED, relative volume
    // error +16.9 at t = 0.031 s. PSI_INT and VOL_PSI_NEG therefore stand in the
    // same row, each with its own drift column, so that no future claim about one
    // can be made from a measurement of the other.
    //
    // NAME MAPPING for the columns as they were specified: SUM_PSI_V is PSI_INT
    // and MAX_DIV_PHI is DIVPHI_LINF. They are not duplicated -- a second copy of
    // the same reduction in every row is one more thing to keep consistent, not
    // one more measurement.
    static const std::string header =
        "TIME,"
        "DELTA_T,"
        "PSI_INT,"
        "PSI_ABS_INT,"
        "D_PSI_INT,"
        "D_PSI_INT_REL,"
        "BUDGET_FLUX,"
        "BUDGET_DIVPHI,"
        "BUDGET_SOURCE,"
        "BUDGET_RESIDUAL,"
        "BUDGET_RESIDUAL_REL,"
        "DIVPHI_LINF,"
        "DIVPHI_L1,"
        "PSI_INT_BAND,"
        "BUDGET_DIVPHI_BAND,"
        "BUDGET_SOURCE_BAND,"
        "SOURCE_ALGEBRAIC,"
        "SUM_PSI_V_DRIFT,"
        "SUM_PSI_V_REL,"
        "SUM_SP_DIVPHI_PSI,"
        "L2_DIV_PHI,"
        "VOL_PSI_NEG,"
        "VOL_PSI_NEG_REL,"
        "VOL_PSI_NEG_SHARP";
    return header;
}


bool Foam::functionObjects::psiConservationCSV::write()
{
    const fvMesh& mesh = this->mesh_;
    const Time& runTime = mesh.time();
    const scalar t = runTime.timeOutputValue();

    // Missing-input sentinel. NEVER 0.0: a zero budget term reads as a CLOSED
    // budget, which is the one conclusion this diagnostic exists to refuse, and
    // the ">0" filters in the plotting/order-fitting scripts drop NaN silently
    // while keeping a fabricated 0 (gradPsiErrorCSV.C:162-171 records the same
    // failure for the empty narrow band). NaN also propagates through every sum
    // below, so BUDGET_RESIDUAL is automatically NaN whenever any of its terms is
    // unavailable -- an unattributable budget must not be reported as a number.
    const scalar undefined = std::numeric_limits<scalar>::quiet_NaN();

    const scalarField& cellV = mesh.V().field();
    const scalarField& psiIn = psi_.primitiveField();

    // gSum(field*mesh.V().field()) on PRIMITIVE fields is the house form for a
    // volume integral (gradPsiErrorCSV.C:176-182, advectionErrors.H:45,57): one
    // MPI_Allreduce per integral, boundary faces excluded by construction (the
    // boundary contributes through BUDGET_FLUX and nowhere else).
    const scalar sumV = gSum(cellV);
    const scalar psiInt = gSum(psiIn*cellV);
    const scalar psiAbsInt = gSum(mag(psiIn)*cellV);
    const scalar psiNorm = Foam::max(psiAbsInt, VSMALL);

    // The interval THIS row closes, measured from the previous row rather than
    // taken from runTime.deltaTValue(). The two differ whenever writeControl is
    // coarser than timeStep, whenever setDeltaT changed the step, and on the
    // first row (no interval exists yet -> NaN, which propagates into every
    // BUDGET_* column of that row and leaves the state columns intact).
    const scalar deltaT = havePrev_ ? (t - prevTime_) : undefined;
    const scalar dPsiInt = havePrev_ ? (psiInt - prevPsiInt_) : undefined;
    const scalar dPsiIntRel = dPsiInt/psiNorm;

    // The narrow-band mask, built once and reused by every *_BAND integral.
    // Typed find for the same reason as `R`: an untyped registry hit on a name
    // reused for another type would abort in the lookup that follows it.
    const volScalarField* nbPtr = mesh.cfindObject<volScalarField>("NarrowBand");
    scalarField nbMask(psiIn.size(), Zero);
    bool haveBand = false;
    if (nbPtr)
    {
        const scalarField& nbField = nbPtr->primitiveField();
        forAll(nbMask, c)
        {
            nbMask[c] = (nbField[c] > 0.5) ? 1.0 : 0.0;
        }
        // The band is EMPTY once the phase is annihilated (no psi sign change
        // anywhere). Band integrals of an empty band are 0 by arithmetic, which
        // would read as "no drift in the band" for an interface that no longer
        // exists -- keep the NaN sentinel instead.
        haveBand = (gSum(nbMask*cellV) > VSMALL);
    }
    const scalar psiIntBand = haveBand ? gSum(nbMask*psiIn*cellV) : undefined;

    // divPhi = fvc::div(phi) is fvc::surfaceIntegrate for a surfaceScalarField
    // argument: it consults NO fvSchemes divSchemes entry, so this object runs on
    // cases whose fvSchemes says `divSchemes { default none; }`
    // (cases/stationaryDroplet2D/system/fvSchemes.template). It reproduces the
    // `divPhi` of alphaEqn.H:20 exactly. Held as a tmp and never named: the
    // solver registers its own volScalarField "divPhi" inside the PIMPLE loop,
    // and a second registration under that name is a registry clash.
    tmp<volScalarField> tdivPhi;
    scalar divPhiLinf = undefined;
    scalar divPhiL1 = undefined;
    scalar divPhiL2 = undefined;
    scalar spRate = undefined;        // gSum(divPhi psi V)          [1/s]
    scalar spRateBand = undefined;    // the same over the band      [1/s]
    scalar fluxOut = undefined;       // sum_{non-coupled bnd f} phi_f psi_f

    if (phiPtr_)
    {
        tdivPhi = fvc::div(*phiPtr_);
        const scalarField& divPhiIn = tdivPhi().primitiveField();

        divPhiLinf = gMax(mag(divPhiIn));
        divPhiL1 = gSum(mag(divPhiIn)*cellV)/Foam::max(sumV, VSMALL);
        divPhiL2 =
            Foam::sqrt(gSum(Foam::sqr(divPhiIn)*cellV)/Foam::max(sumV, VSMALL));
        spRate = gSum(divPhiIn*psiIn*cellV);
        if (haveBand)
        {
            spRateBand = gSum(nbMask*divPhiIn*psiIn*cellV);
        }

        // The surviving remainder of the telescoping flux sum: only NON-COUPLED
        // patches. A processor patch is an interior face of the global mesh --
        // its owner-side and neighbour-side contributions live on different ranks
        // and cancel in the reduce below -- so counting it here would add a term
        // that the discrete equation does not contain, of the order of the
        // through-partition flux (i.e. arbitrarily large, and decomposition
        // dependent, which would make the budget depend on `np`). `empty` and
        // `wedge` patches carry zero faces and contribute nothing either way.
        scalar bndFlux = 0;
        const surfaceScalarField::Boundary& phiBf = phiPtr_->boundaryField();
        const volScalarField::Boundary& psiBf = psi_.boundaryField();
        forAll(phiBf, patchi)
        {
            if (phiBf[patchi].patch().coupled())
            {
                continue;
            }
            const fvsPatchScalarField& pphi = phiBf[patchi];
            const fvPatchScalarField& ppsi = psiBf[patchi];
            forAll(pphi, facei)
            {
                bndFlux += pphi[facei]*ppsi[facei];
            }
        }
        reduce(bndFlux, sumOp<scalar>());
        fluxOut = bndFlux;
    }

    // The SDPLS source coefficient f_nl(psi^n), registered as a volScalarField by
    // the sdplsSource base class (sdplsSource.C:117-127) even for `noSource`,
    // where it is identically 0 and gSum(f_nl psi V) is a GENUINE zero -- not the
    // sentinel -- so the budget still closes on a sourceless run. Absent field =
    // no sdplsSource object was constructed at all (a kinematic or SL run), in
    // which case the equation this budget decomposes was never assembled and the
    // source column must not claim a value.
    const volScalarField* fnlPtr =
        mesh.cfindObject<volScalarField>("SDPLS_nonLinearPart");
    scalar sourceRate = undefined;      // gSum(f_nl psi V)     [1/s]
    scalar sourceRateBand = undefined;
    if (algebraicSource_ && fnlPtr)
    {
        const scalarField& fnlIn = fnlPtr->primitiveField();
        sourceRate = gSum(fnlIn*psiIn*cellV);
        if (haveBand)
        {
            sourceRateBand = gSum(nbMask*fnlIn*psiIn*cellV);
        }
    }

    // The SECOND conserved quantity, by the measure the studies already report:
    // gSum(alpha V) on the registered phase indicator, exactly as
    // advectionErrors.H:45,57 forms the volume error. Both indicators in the
    // library set alpha = 1 where psi < 0 (heavisidePhaseIndicator.C:62-77 ends
    // in `alpha == 1 - alpha`; geometricPhaseIndicator.C:103-107), so this IS the
    // {psi < 0} phase volume -- SMEARED over epsilon = (nCells/2) max(1/deltaCoeffs)
    // for PHASE_INDICATOR heaviside, geometric (plane-reconstructed) otherwise.
    const volScalarField* alphaPtr = mesh.cfindObject<volScalarField>(alphaName_);
    const scalar volNeg =
        alphaPtr ? gSum(alphaPtr->primitiveField()*cellV) : undefined;

    // The indicator-free counterpart: the cell-centre sign count
    // sum_c [psi_c < 0] V_c. First order in h and it cannot resolve a subcell
    // interface, but it depends on NOTHING except psi -- so a divergence between
    // it and VOL_PSI_NEG is a statement about the phase indicator (its smearing
    // width, or its reconstruction), while the two moving together is a real
    // displacement of the zero level set.
    scalarField negMask(psiIn.size(), Zero);
    forAll(negMask, c)
    {
        negMask[c] = (psiIn[c] < 0) ? 1.0 : 0.0;
    }
    const scalar volNegSharp = gSum(negMask*cellV);

    // Reference for the run-length drift columns, taken from the FIRST row this
    // object writes (the one the constructor issues at startTime).
    if (!haveRef_)
    {
        psiInt0_ = psiInt;
        volNeg0_ = volNeg;
        haveRef_ = true;
    }
    const scalar psiIntDrift = psiInt - psiInt0_;
    // Denominator |int psi dV|(t=0) as specified -- NOT PSI_ABS_INT. It is the
    // literal "relative to its initial value", and it is the column to distrust
    // when the initial configuration is nearly antisymmetric in psi: the signed
    // integral can then be orders of magnitude below int|psi| dV and this ratio
    // inflates without anything having drifted. D_PSI_INT_REL and
    // BUDGET_RESIDUAL_REL, normalised by PSI_ABS_INT, do not have that mode.
    const scalar psiIntDriftRel =
        psiIntDrift/Foam::max(Foam::mag(psiInt0_), VSMALL);
    const scalar volNegRel =
        (volNeg - volNeg0_)/Foam::max(volNeg0_, VSMALL);

    // The budget, all four terms in the same units as PSI_INT. NaN propagates:
    // no DELTA_T (first row), no phi, or a non-algebraic source each make
    // BUDGET_RESIDUAL NaN rather than a partial sum masquerading as a closure.
    const scalar budgetFlux = -deltaT*fluxOut;
    const scalar budgetDivPhi = deltaT*spRate;
    const scalar budgetSource = deltaT*sourceRate;
    const scalar budgetDivPhiBand = deltaT*spRateBand;
    const scalar budgetSourceBand = deltaT*sourceRateBand;
    const scalar budgetResidual =
        dPsiInt - (budgetFlux + budgetDivPhi + budgetSource);
    const scalar budgetResidualRel = budgetResidual/psiNorm;

    // Rank 0 writes; every reduction above has already been executed on ALL
    // ranks (a gSum inside this guard would deadlock in MPI_Allreduce).
    if (Pstream::myProcNo() == 0)
    {
        csvFile_ << t << ","
            << deltaT << ","
            << psiInt << ","
            << psiAbsInt << ","
            << dPsiInt << ","
            << dPsiIntRel << ","
            << budgetFlux << ","
            << budgetDivPhi << ","
            << budgetSource << ","
            << budgetResidual << ","
            << budgetResidualRel << ","
            << divPhiLinf << ","
            << divPhiL1 << ","
            << psiIntBand << ","
            << budgetDivPhiBand << ","
            << budgetSourceBand << ","
            << (algebraicSource_ ? 1 : 0) << ","
            << psiIntDrift << ","
            << psiIntDriftRel << ","
            << spRate << ","
            << divPhiL2 << ","
            << volNeg << ","
            << volNegRel << ","
            << volNegSharp << "\n";
    }

    // State for the NEXT row, on ALL ranks: prevPsiInt_ feeds a rank-local
    // subtraction whose result must be identical everywhere, and every rank
    // evaluates the havePrev_ branch above.
    prevPsiInt_ = psiInt;
    prevTime_ = t;
    havePrev_ = true;

    return true;
}


// ************************************************************************* //
