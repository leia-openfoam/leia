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

#include "newtonShiftVolumeCorrection.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
// mesh.deltaCoeffs() is a surfaceScalarField: without this include surfaceMesh
// is only forward-declared and gMax(...primitiveField()) fails to instantiate.
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(newtonShiftVolumeCorrection, false);
addToRunTimeSelectionTable(volumeCorrection, newtonShiftVolumeCorrection, Mesh);

// Alias. The method is named after its ALGORITHM here (a safeguarded
// Newton-Raphson root find for the shift eps), while the level-set literature
// and the design note for this work name it after its EFFECT ("global shift"
// volume correction). Registering both names costs one constructor-table
// entry and removes a silent-misconfiguration mode: a dictionary written as
// `type globalShift;` would otherwise abort with an unknown-type
// FatalIOError. Both spellings select this class; there is one implementation.
addNamedToRunTimeSelectionTable
(
    volumeCorrection,
    newtonShiftVolumeCorrection,
    Mesh,
    globalShift
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newtonShiftVolumeCorrection::newtonShiftVolumeCorrection(const fvMesh& mesh)
    :
        volumeCorrection(mesh),
        maxIter_(volCorrDict_.getOrDefault<label>("maxIter", 10)),
        relTol_(volCorrDict_.getOrDefault<scalar>("relTol", 1e-10)),
        maxShiftCells_(volCorrDict_.getOrDefault<scalar>("maxShiftCells", 0.5)),
        derivStepCells_
        (
            volCorrDict_.getOrDefault<scalar>("derivativeStepCells", 1e-3)
        ),
        psiBcChecked_(false)
{
    // The corrector solves f(eps) = V(eps) - V_target = 0 by Newton, so
    // V(eps) must be a differentiable function of eps. It is not for every
    // phase indicator in this library, and the failure is silent: a step
    // function has zero derivative almost everywhere, the Newton step
    // degenerates to bisection on a staircase, and the run reports a
    // "converged" shift that landed on a plateau. Refuse instead.
    const dictionary phaseIndDict =
        levelSetDict_.subOrEmptyDict("phaseIndicator");

    // Same default as phaseIndicator::New (phaseIndicator.C:31).
    const word phaseIndType =
        phaseIndDict.getOrDefault<word>("type", "geometric");

    if (phaseIndType == "sharpJump" || phaseIndType == "none")
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "volumeCorrection type " << type() << " requires a phase "
            << "indicator whose volume V(eps) = sum_c alpha_c(psi + eps) "
            << "|Omega_c| is a differentiable function of the uniform shift "
            << "eps." << nl
            << "The selected phaseIndicator type is '" << phaseIndType
            << "': 'sharpJump' makes V(eps) a piecewise-constant staircase "
            << "with jumps of one full cell volume |Omega_c|, and 'none' "
            << "leaves alpha untouched altogether "
            << "(phaseIndicator.H:97-104), so V does not depend on eps at "
            << "all. No root finder is well posed on either." << nl
            << "Use phaseIndicator type geometric, detrixheAslam or "
            << "heaviside, or remove the volumeCorrection sub-dictionary."
            << exit(FatalIOError);
    }

    // detrixheAslam with geometrySource analyticImplicitSurface computes
    // alpha from the ANALYTIC surface and ignores psi entirely
    // (detrixheAslamPhaseIndicator.C:152-159). V(eps) is then constant,
    // f'(eps) = 0 for every eps, and the correction can never converge -- it
    // would bisect the full bracket and apply a shift of 0.5 h chosen by
    // nothing.
    if
    (
        phaseIndType == "detrixheAslam"
     && phaseIndDict.getOrDefault<word>("geometrySource", "levelSetField")
        == "analyticImplicitSurface"
    )
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "volumeCorrection type " << type() << " cannot be used with "
            << "phaseIndicator detrixheAslam and geometrySource "
            << "analyticImplicitSurface: alpha is then evaluated from the "
            << "analytic implicit surface and is INDEPENDENT of psi, so "
            << "V(eps) is constant and f'(eps) = 0 identically." << nl
            << "Set geometrySource levelSetField, or remove the "
            << "volumeCorrection sub-dictionary." << exit(FatalIOError);
    }

    if (maxShiftCells_ <= 0)
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "maxShiftCells must be > 0 (it is the half-width of the "
            << "search bracket in cells: |eps| <= maxShiftCells*h), got "
            << maxShiftCells_ << "." << nl
            << "A non-positive value gives an empty bracket, which is not "
            << "the way to disable the correction -- remove the "
            << "volumeCorrection sub-dictionary, or set type "
            << "noVolumeCorrection." << exit(FatalIOError);
    }

    if (derivStepCells_ <= 0)
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "derivativeStepCells must be > 0 (it is the central-difference "
            << "step in cells: delta = derivativeStepCells*h), got "
            << derivStepCells_ << "." << exit(FatalIOError);
    }

    if (maxIter_ < 1)
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "maxIter must be >= 1, got " << maxIter_ << "."
            << exit(FatalIOError);
    }

    if (relTol_ <= 0)
    {
        FatalIOErrorInFunction(volCorrDict_)
            << "relTol must be > 0 (it is tested against the RELATIVE volume "
            << "residual |V(eps) - V_target|/V_target), got " << relTol_
            << "." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar newtonShiftVolumeCorrection::shiftedVolume
(
    const volScalarField& psi,
    volScalarField& alphaTmp,
    volScalarField& psiShift,
    phaseIndicator& phaseInd,
    const scalar eps
) const
{
    psiShift = psi + dimensionedScalar(psi.dimensions(), eps);

    // The coupled-patch values must follow the interior BEFORE the indicator
    // reads them: the geometric indicators fit their least-squares plane from
    // patchNeighbourField() (levelSetPlaneReconstruction.C:34-56), so a
    // processor-boundary cell fitted against an unshifted far side would
    // produce a decomposition-dependent alpha, and hence a
    // decomposition-dependent V(eps).
    psiShift.correctBoundaryConditions();

    phaseInd.calcPhaseIndicator(alphaTmp, psiShift);

    // gSum, not sum: identical on every rank, which is what keeps the Newton
    // iterate, the convergence test and the ITERATION COUNT identical on
    // every rank. calcPhaseIndicator issues matched MPI, so a rank-divergent
    // loop deadlocks rather than merely disagreeing.
    return gSum(alphaTmp.primitiveField()*psi.mesh().V().field());
}


bool newtonShiftVolumeCorrection::correct
(
    volScalarField& psi,
    phaseIndicator& phaseInd
)
{
    lastShift_ = 0;

    const fvMesh& mesh = psi.mesh();

    if (targetVolume_ <= 0)
    {
        WarningInFunction
            << "no usable target volume (V_target = " << targetVolume_
            << "): setTargetVolume() was not called, or phase 1 is empty. "
            << "No correction applied." << endl;
        return false;
    }

    // Checked ONCE per run, not once per step. A uniform shift leaves the
    // discrete gradient invariant only if the boundary values move with the
    // interior. zeroGradient (and the constraint types empty/symmetry/
    // processor) do that on correctBoundaryConditions(); a patch that FIXES
    // psi does not, and the near-boundary |grad psi| then changes by
    // eps/delta_n -- the one way this correction can damage the signed
    // distance property it is built to preserve.
    if (!psiBcChecked_)
    {
        psiBcChecked_ = true;

        forAll(psi.boundaryField(), patchi)
        {
            if (psi.boundaryField()[patchi].fixesValue())
            {
                WarningInFunction
                    << "psi patch '" << mesh.boundary()[patchi].name()
                    << "' is of type "
                    << psi.boundaryField()[patchi].type()
                    << ", which FIXES the boundary value. The uniform shift "
                    << "eps is applied to the internal field and the patch "
                    << "value does not follow, so |grad psi| changes by "
                    << "eps/delta_n in the boundary cells and the "
                    << "distance-preservation argument of this correction "
                    << "does not hold there." << endl;
            }
        }
    }

    // h and the bracket are GLOBAL reductions (redistancer.C:104 spelling),
    // so the bracket does not depend on the decomposition either.
    const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());
    const scalar epsMax = maxShiftCells_*h;
    const scalar delta = derivStepCells_*h;

    // Scratch fields. psiShift carries psi's boundary condition TYPES (copy
    // construction with a reset IOobject), so the shifted field is evaluated
    // exactly as psi would be. Neither field is registered: the object
    // registry must not gain a second "psi", and neither field may reach
    // runTime.write().
    volScalarField psiShift
    (
        IOobject
        (
            "psiVolumeCorrectionShift",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        psi
    );

    volScalarField alphaTmp
    (
        IOobject
        (
            "alphaVolumeCorrection",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        word("zeroGradient")
    );

    // f(eps) = V(eps) - V_target, evaluated through the SELECTED indicator.
    const scalar V0 =
        shiftedVolume(psi, alphaTmp, psiShift, phaseInd, 0);
    const scalar f0 = V0 - targetVolume_;

    if (mag(f0)/targetVolume_ <= relTol_)
    {
        // Already within tolerance: psi is not touched at all, so the caller
        // does not recalculate the narrow band and the step is bit-identical
        // to an uncorrected one.
        Info<< "volumeCorrection: V = " << V0 << ", target " << targetVolume_
            << ", dV/V " << f0/targetVolume_
            << " within relTol " << relTol_ << ", eps = 0 (no correction)"
            << endl;
        return false;
    }

    // Bracket the root on [-epsMax, +epsMax]. V(eps) is monotone decreasing
    // in eps wherever it is smooth (alpha = 1 in {psi < 0}, so raising psi
    // can only shrink Omega^-), so the two endpoint residuals normally
    // straddle zero. Evaluating them up front costs two indicator passes and
    // buys the bisection fallback used whenever the Newton step leaves the
    // bracket -- necessary because V(eps) is only PIECEWISE smooth: a bulk
    // cell centre crossing zero adds a jump of |Omega_c|.
    scalar lo = -epsMax;
    scalar hi =  epsMax;
    scalar fLo = shiftedVolume(psi, alphaTmp, psiShift, phaseInd, lo)
               - targetVolume_;
    scalar fHi = shiftedVolume(psi, alphaTmp, psiShift, phaseInd, hi)
               - targetVolume_;

    if (fLo*fHi > 0)
    {
        // No root within +-maxShiftCells*h: the volume error of this step is
        // larger than a uniform half-cell displacement can repair. Apply the
        // endpoint that reduces |f| the most rather than nothing: the
        // correction stays monotone toward the target, the per-step
        // distortion stays bounded by epsMax, and the residual is named in
        // the log. Refusing outright would let the deficit accumulate and
        // require an even larger shift on the next step, silently.
        const scalar epsSat = (mag(fLo) < mag(fHi)) ? lo : hi;
        const scalar fSat   = (mag(fLo) < mag(fHi)) ? fLo : fHi;

        WarningInFunction
            << "no root of V(eps) - V_target in the bracket [" << lo << ", "
            << hi << "] = +-" << maxShiftCells_ << " cells: f(-epsMax) = "
            << fLo << ", f(+epsMax) = " << fHi << "." << nl
            << "The volume error " << f0/targetVolume_
            << " (relative) exceeds what a uniform shift of "
            << maxShiftCells_ << " cells can restore. Applying the saturated "
            << "shift eps = " << epsSat << " and leaving a relative residual "
            << fSat/targetVolume_ << ". A volume error of this size is not a "
            << "bookkeeping error -- check the interface, not the corrector."
            << endl;

        psi += dimensionedScalar(psi.dimensions(), epsSat);
        psi.correctBoundaryConditions();
        lastShift_ = epsSat;

        Info<< "volumeCorrection: V = " << V0 << " -> " << fSat + targetVolume_
            << ", target " << targetVolume_
            << ", dV/V " << fSat/targetVolume_
            << ", eps = " << epsSat << " (" << epsSat/h
            << " cells, SATURATED), 0 Newton iterations" << endl;

        return (epsSat != 0);
    }

    // Insert the starting iterate eps = 0 into the bracket. The bracket
    // invariant maintained below is lo <= eps <= hi with fLo*fHi <= 0.
    scalar eps  = 0;
    scalar fEps = f0;

    if (fLo*fEps <= 0)
    {
        hi = eps;
        fHi = fEps;
    }
    else
    {
        lo = eps;
        fLo = fEps;
    }

    label iter = 0;
    bool converged = false;

    for (iter = 1; iter <= maxIter_; ++iter)
    {
        // Central difference through the SAME code path that produces the
        // residual: the derivative can never disagree with the function
        // being zeroed. Two indicator passes per iteration.
        const scalar dVdEps =
        (
            shiftedVolume(psi, alphaTmp, psiShift, phaseInd, eps + delta)
          - shiftedVolume(psi, alphaTmp, psiShift, phaseInd, eps - delta)
        )/(2*delta);

        scalar epsNew = 0.5*(lo + hi);

        if (mag(dVdEps) > VSMALL)
        {
            epsNew = eps - fEps/dVdEps;
        }

        // Written as a POSITIVE test so that a NaN or an infinity produced by
        // a degenerate derivative also falls through to bisection: the
        // comparison is false for NaN, and !false takes the bisection value.
        if (!(epsNew >= lo && epsNew <= hi))
        {
            epsNew = 0.5*(lo + hi);
        }

        eps = epsNew;
        fEps = shiftedVolume(psi, alphaTmp, psiShift, phaseInd, eps)
             - targetVolume_;

        if (fLo*fEps <= 0)
        {
            hi = eps;
            fHi = fEps;
        }
        else
        {
            lo = eps;
            fLo = fEps;
        }

        // The test is written on the REDUCED residual, so it evaluates
        // identically on every rank and the loop is left by every rank in
        // the same iteration.
        if (mag(fEps)/targetVolume_ <= relTol_)
        {
            converged = true;
            break;
        }
    }

    if (!converged)
    {
        WarningInFunction
            << "Newton-Raphson did not reach relTol = " << relTol_
            << " in " << maxIter_ << " iterations: remaining relative volume "
            << "residual |V(eps) - V_target|/V_target = "
            << mag(fEps)/targetVolume_ << " at eps = " << eps << " ("
            << eps/h << " cells)." << nl
            << "The shift IS applied -- it still reduces the volume error -- "
            << "but a correction that silently fails to converge is worse "
            << "than none, so the residual is named here. On the geometric "
            << "indicators the attainable floor is set by the polyhedron "
            << "clipping tolerance, not by this loop; raise relTol before "
            << "raising maxIter." << endl;
    }

    psi += dimensionedScalar(psi.dimensions(), eps);
    psi.correctBoundaryConditions();
    lastShift_ = eps;

    // iter is maxIter_ + 1 when the loop ran to exhaustion (the for-loop
    // increment fires before the test), so report the count of iterations
    // actually performed.
    Info<< "volumeCorrection: V = " << V0 << " -> " << fEps + targetVolume_
        << ", target " << targetVolume_
        << ", dV/V " << fEps/targetVolume_
        << ", eps = " << eps << " (" << eps/h << " cells), "
        << (converged ? iter : maxIter_) << " Newton iterations" << endl;

    return (eps != 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
