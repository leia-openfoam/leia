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

#include "pointValueScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvSolution.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointValueScheme, 0);
    addToRunTimeSelectionTable(slScheme, pointValueScheme, Mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointValueScheme::pointValueScheme(const fvMesh& mesh)
:
    slScheme(mesh),
    trajectoryVelocity_("input"),
    trajectoryUName_("U"),
    normalEpsilon_(1e-12),
    projectionIterations_(4),
    projectionTolerance_(0.1),
    projectionMaxDistance_(6.0)
{
    const fvSolution& fvSol(mesh);
    const dictionary& slDict =
        fvSol.subDict("levelSet").subOrEmptyDict("semiLagrangian");

    trajectoryVelocity_ =
        slDict.getOrDefault<word>("trajectoryVelocity", "input");
    trajectoryUName_ = slDict.getOrDefault<word>("trajectoryU", "U");
    normalEpsilon_ =
        slDict.getOrDefault<scalar>("trajectoryNormalEpsilon", 1e-12);
    projectionIterations_ =
        slDict.getOrDefault<label>("trajectoryProjectionIterations", 4);
    projectionTolerance_ =
        slDict.getOrDefault<scalar>("trajectoryProjectionTolerance", 0.1);
    projectionMaxDistance_ =
        slDict.getOrDefault<scalar>("trajectoryProjectionMaxDistance", 6.0);

    if
    (
        trajectoryVelocity_ != "input"
     && trajectoryVelocity_ != "normalProjection"
     && trajectoryVelocity_ != "normalClosestPoint"
    )
    {
        FatalIOErrorInFunction(slDict)
            << "Unknown semiLagrangian.trajectoryVelocity '"
            << trajectoryVelocity_ << "'. Valid values are input, "
            << "normalProjection and normalClosestPoint."
            << exit(FatalIOError);
    }
    if
    (
        trajectoryVelocity_ != "input"
     &&
        (
            normalEpsilon_ <= 0
         || projectionIterations_ < 1
         || projectionTolerance_ <= 0
         || projectionMaxDistance_ <= 0
        )
    )
    {
        FatalIOErrorInFunction(slDict)
            << "Direct normal trajectory controls must be positive: "
            << "trajectoryNormalEpsilon=" << normalEpsilon_ << ", "
            << "trajectoryProjectionIterations=" << projectionIterations_
            << ", trajectoryProjectionTolerance=" << projectionTolerance_
            << ", trajectoryProjectionMaxDistance="
            << projectionMaxDistance_ << exit(FatalIOError);
    }

    Info<< "pointValue trajectory velocity: " << trajectoryVelocity_;
    if (trajectoryVelocity_ != "input")
    {
        Info<< " (physical field " << trajectoryUName_ << ')';
    }
    Info<< endl;
}


Foam::autoPtr<Foam::volVectorField> Foam::pointValueScheme::normalVelocity
(
    const word& name,
    const volScalarField& psiLevel,
    const volVectorField& ULevel,
    const bool closestPoint
) const
{
    // Smoothly regularised unit normal. The projection (U.n)n is invariant to
    // a reversal of the level-set sign convention.
    const volVectorField gradPsi(fvc::grad(psiLevel));
    const dimensionedScalar eps
    (
        "trajectoryNormalEpsilon",
        gradPsi.dimensions(),
        normalEpsilon_
    );
    const volVectorField nHat
    (
        IOobject
        (
            name + "NHat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        gradPsi/sqrt(magSqr(gradPsi) + sqr(eps))
    );

    autoPtr<volVectorField> WPtr
    (
        new volVectorField
        (
            IOobject
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (ULevel & nHat)*nHat
        )
    );

    if (!closestPoint)
    {
        return WPtr;
    }

    // Characteristic size = minimum face-neighbour centre distance. This is
    // well scaled on a one-cell-thick 2D mesh, unlike cbrt(cell volume).
    scalarField cellSize(mesh_.nCells(), GREAT);
    const surfaceScalarField& dc = mesh_.deltaCoeffs();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, faceI)
    {
        const scalar d = 1.0/dc[faceI];
        cellSize[own[faceI]] = Foam::min(cellSize[own[faceI]], d);
        cellSize[nei[faceI]] = Foam::min(cellSize[nei[faceI]], d);
    }
    forAll(mesh_.boundary(), patchI)
    {
        if (!mesh_.boundary()[patchI].coupled())
        {
            continue;
        }
        const labelUList& fc = mesh_.boundary()[patchI].faceCells();
        const fvsPatchScalarField& pdc = dc.boundaryField()[patchI];
        forAll(fc, i)
        {
            cellSize[fc[i]] = Foam::min(cellSize[fc[i]], 1.0/pdc[i]);
        }
    }

    interpolationCellPoint<scalar> psiInterp(psiLevel);
    interpolationCellPoint<vector> gradInterp(gradPsi);
    interpolationCellPoint<vector> UInterp(ULevel);

    const volScalarField* bandPtr =
        mesh_.findObject<volScalarField>("NarrowBand");
    const pointField& C = mesh_.C();
    vectorField& W = WPtr->primitiveFieldRef();
    label nSuccess = 0;
    label nFallback = 0;

    forAll(C, c)
    {
        // Outside the phase-indicator/force band the local normal projection
        // already stored in W is sufficient. Avoid a long, ambiguous walk to a
        // remote interface which cannot affect the current geometric update.
        if (bandPtr && (*bandPtr)[c] < 0.5)
        {
            continue;
        }

        const scalar hc = cellSize[c];
        if (!std::isfinite(hc) || hc <= VSMALL)
        {
            ++nFallback;
            continue;
        }

        point x = C[c];
        label cx = c;
        scalar residual = Foam::mag(psiLevel[c]);
        bool ok = true;

        for (label iter = 0; iter < projectionIterations_; ++iter)
        {
            const scalar psiX = psiInterp.interpolate(x, cx);
            const vector gX = gradInterp.interpolate(x, cx);
            const scalar g2 = magSqr(gX);
            if
            (
                !std::isfinite(psiX)
             || !std::isfinite(g2)
             || g2 <= sqr(normalEpsilon_)
            )
            {
                ok = false;
                break;
            }

            const point xNew =
                x - psiX*gX/(g2 + sqr(normalEpsilon_));
            if (Foam::mag(xNew - C[c]) > projectionMaxDistance_*hc)
            {
                ok = false;
                break;
            }

            const label cNew = mesh_.findCell(xNew);
            if (cNew < 0)
            {
                // Processor-crossing projection needs a distributed point
                // sampler. Until that is added, retain the local projection.
                ok = false;
                break;
            }

            const scalar residualNew =
                Foam::mag(psiInterp.interpolate(xNew, cNew));
            if (!std::isfinite(residualNew) || residualNew > residual + SMALL*hc)
            {
                ok = false;
                break;
            }

            x = xNew;
            cx = cNew;
            residual = residualNew;
        }

        if (ok && residual <= projectionTolerance_*hc)
        {
            const vector gGamma = gradInterp.interpolate(x, cx);
            const scalar gm = Foam::sqrt
            (
                magSqr(gGamma) + sqr(normalEpsilon_)
            );
            if (std::isfinite(gm) && gm > normalEpsilon_)
            {
                const vector nGamma = gGamma/gm;
                const vector UGamma = UInterp.interpolate(x, cx);
                W[c] = (UGamma & nGamma)*nGamma;
                ++nSuccess;
                continue;
            }
        }

        ++nFallback;
    }

    reduce(nSuccess, sumOp<label>());
    reduce(nFallback, sumOp<label>());
    Info<< "normalClosestPoint trajectory sampling (" << psiLevel.name()
        << "): success/fallback = " << nSuccess << '/' << nFallback << endl;

    WPtr->correctBoundaryConditions();
    return WPtr;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pointValueScheme::advance
(
    volScalarField& psi,
    const volVectorField& Unew,
    const volVectorField& Uold,
    slReconstruction& recon,
    slCorrector& corrector
)
{
    const scalar dt = mesh_.time().deltaTValue()*dtScale_;
    const volVectorField& C = mesh_.C();

    const volVectorField* trajectoryNew = &Unew;
    const volVectorField* trajectoryOld = &Uold;
    autoPtr<volVectorField> WnewPtr;
    autoPtr<volVectorField> WoldPtr;

    if (trajectoryVelocity_ != "input")
    {
        // Direct normal modes deliberately sample the physical velocity from
        // the registry. This avoids accidentally projecting a precomputed
        // Eulerian Uext when a two-phase solver has supplied one to the legacy
        // input path.
        const volVectorField& physicalU =
            mesh_.lookupObject<volVectorField>(trajectoryUName_);
        const volVectorField& physicalUOld = physicalU.oldTime();
        const volScalarField& psiOld = psi.oldTime();
        const bool useClosest =
            (trajectoryVelocity_ == "normalClosestPoint");

        WnewPtr = normalVelocity
        (
            "trajectoryUNew",
            psi,
            physicalU,
            useClosest
        );
        WoldPtr = normalVelocity
        (
            "trajectoryUOld",
            psiOld,
            physicalUOld,
            useClosest
        );
        trajectoryNew = &WnewPtr();
        trajectoryOld = &WoldPtr();
    }

    // Gradient of the ACTUAL velocity inserted into Taylor. For the direct
    // normal modes this includes spatial variation of both the normal speed and
    // the projection tensor n*n; using grad(physical U) would lose second order.
    const volTensorField gradU(fvc::grad(*trajectoryNew, "gradU"));

    // ------------------------------------------------------------------ //
    // Departure (foot) points: computed ONCE. The Taylor backward-foot
    // integrator is UNCHANGED -- it is only hoisted out of the reconstruction
    // loop below and cached, so the deferred-correction passes reuse the same
    // characteristic feet. (KEEP the 1/2 on the dt^2 term.)
    //   x_d = x_c - u^{n+1} dt + 1/2 [ du/dt + (u^{n+1}.grad)u^{n+1} ] dt^2
    // ------------------------------------------------------------------ //
    pointField feet(mesh_.nCells());
    forAll(C, c)
    {
        const vector& uNew = (*trajectoryNew)[c];
        const vector& uOld = (*trajectoryOld)[c];
        const vector accel = (uNew - uOld)/dt + (uNew & gradU[c]);
        feet[c] = C[c] - uNew*dt + 0.5*accel*dt*dt;
    }

    // ------------------------------------------------------------------ //
    // The selected correction strategy assembles psi^{n+1} in place from psi^n,
    // evaluating the reconstruction at the fixed feet and (for deferredCorrection)
    // rebuilding it from the current iterate. It owns the value cap / non-finite
    // reset / foot-radius guard. The backtracking above is not its concern.
    // ------------------------------------------------------------------ //
    corrector.correct(psi, feet, recon);

    // Optional post-advection fix-up (band model: re-extend psi outside the band
    // as a clean signed distance so freshly-entered band cells get a good value).
    recon.postAdvect(psi);
}

// ************************************************************************* //
