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

#include "semiImplicitCapillaryForce.H"
#include "fvMatrices.H"
#include "fvcSnGrad.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(semiImplicitCapillaryForce, 0);
    addToRunTimeSelectionTable
    (
        option,
        semiImplicitCapillaryForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::fv::semiImplicitCapillaryForce::muFace
(
    const volScalarField& alpha
) const
{
    // coeff * sigma * dt * |snGrad(alpha)|_f : viscosity units, supported on
    // exactly the faces the balanced CSF flux GSigma differences.
    return tmp<surfaceScalarField>::New
    (
        "muSigmaFace",
        coeff_*sigma_*mesh_.time().deltaT()*mag(fvc::snGrad(alpha))
    );
}


void Foam::fv::semiImplicitCapillaryForce::addTerm(fvMatrix<vector>& eqn)
{
    const volVectorField& U = eqn.psi();

    const auto* alphaPtr = mesh_.cfindObject<volScalarField>(alphaName_);
    if (!alphaPtr)
    {
        FatalErrorInFunction
            << "semiImplicitCapillaryForce: field '" << alphaName_
            << "' (alphaName) is not registered." << exit(FatalError);
    }
    const volScalarField& alpha = *alphaPtr;

    const surfaceScalarField muf(muFace(alpha));

    // Diagnostics bookkeeping: one row per addSup call = per outer iteration.
    if (mesh_.time().timeIndex() != lastTimeIndex_)
    {
        lastTimeIndex_ = mesh_.time().timeIndex();
        outerCall_ = 0;
    }
    ++outerCall_;

    if (form_ == "increment")
    {
        // Deferred Jacobian: mu*Lap(u^{n+1} - u^(k)). U holds the current
        // outer iterate at assembly time, so no prevIter storage is needed.
        // The Laplace-Beltrami projector corrections cancel identically
        // between the implicit and the lagged part (both would be evaluated
        // at u^(k)), so this form needs neither kappa nor nHat.
        eqn += fvm::laplacian(muf, U);
        eqn -= fvc::laplacian(muf, U);

        if (diagnostics_)
        {
            // The lagged residual of the term: mu*Lap(u^(k)) in L2/Linf over
            // the domain. At outer-loop convergence the ADDED term is the
            // difference of two of these at consecutive iterates -> the
            // per-iterate CHANGE of this norm is the vanishing check.
            //
            // EVERY RANK must reach the reductions: gAverage/gMax are
            // COLLECTIVE. Guarding them behind Pstream::master() deadlocks the
            // run -- rank 0 blocks in the reduction and the others never enter
            // it. Measured 2026-08-28 on Lichtenberg: four arms sat alive but
            // silent for 76 minutes inside the first momentum assembly, which
            // is invisible in serial and is why local testing never saw it.
            // Reduce on all ranks; write on master only.
            const volVectorField lag(fvc::laplacian(muf, U));
            const scalarField magLag(mag(lag)());
            const scalar termL2 = Foam::sqrt(gAverage(magSqr(magLag)));
            const scalar termLinf = gMax(magLag);
            if (csvPtr_ && Pstream::master())
            {
                csvPtr_() << mesh_.time().timeOutputValue() << ","
                    << outerCall_ << ",increment,"
                    << termL2 << "," << termLinf << "," << 0.0 << "\n";
            }
        }
        return;
    }

    // ---- form == "value" (SAAMPLE / Raessi) -------------------------------

    eqn += fvm::laplacian(muf, U);

    if (laplaceBeltrami_)
    {
        const auto* kappaPtr = mesh_.cfindObject<volScalarField>(kappaName_);
        const auto* nHatPtr = mesh_.cfindObject<volVectorField>(nHatName_);
        if (!kappaPtr || !nHatPtr)
        {
            FatalErrorInFunction
                << "semiImplicitCapillaryForce (form value, laplaceBeltrami):"
                << " registered fields '" << kappaName_ << "' and '"
                << nHatName_ << "' are required. Present: kappa="
                << (kappaPtr != nullptr) << " nHat=" << (nHatPtr != nullptr)
                << ". The SL solvers register both; other solvers must name"
                << " their own or set laplaceBeltrami false."
                << exit(FatalError);
        }
        const volScalarField& kappa = *kappaPtr;
        const volVectorField& nHat = *nHatPtr;

        // Cell coefficient for the explicit corrections: the cell-native CSF
        // delta |grad(alpha)|. The face/cell delta mismatch sits inside a
        // deferred-corrected explicit term (one outer iteration lag), so it
        // perturbs the iteration path, not the converged operator; the
        // pairing-critical support is carried by the IMPLICIT face term.
        const volScalarField muC
        (
            "muSigmaCell",
            coeff_*sigma_*mesh_.time().deltaT()*mag(fvc::grad(alpha))
        );

        // SAAMPLE Eq. (23): Lap_Gamma(u) = Lap(u) - div((n.gradU) n)
        //                   + kappa*((gradU - (n.gradU) n) . n)
        // The implicit part above is Lap; subtract the normal parts
        // explicitly, evaluated at the current outer iterate. The tensor is
        // NAMED so its divergence scheme is the fixed dictionary entry
        // div(semiImplicitNormalFlux) (Gauss linear in the case templates)
        // instead of an unreadable expression-derived lookup key.
        tmp<volTensorField> tgradU = fvc::grad(U);
        const volTensorField& gradU = tgradU();

        const volTensorField nGradUn
        (
            "semiImplicitNormalFlux",
            (nHat & gradU)*nHat
        );

        const volVectorField normalPart
        (
            fvc::div(nGradUn)
          - kappa*((gradU - nGradUn) & nHat)
        );

        eqn -= muC*normalPart;

        if (diagnostics_)
        {
            // Damping power of the value form: integral mu_c |grad U|^2 dV --
            // the overdamping instrument (must be ~0 at equilibrium). As above:
            // gSum/gAverage/gMax are COLLECTIVE, so all ranks reduce and only
            // master writes.
            const scalar power =
                gSum((mesh_.V()*muC.primitiveField()
                     *magSqr(gradU.primitiveField()))());
            const scalarField magNP(mag(muC*normalPart)());
            const scalar termL2 = Foam::sqrt(gAverage(magSqr(magNP)));
            const scalar termLinf = gMax(magNP);
            if (csvPtr_ && Pstream::master())
            {
                csvPtr_() << mesh_.time().timeOutputValue() << ","
                    << outerCall_ << ",value,"
                    << termL2 << "," << termLinf << "," << power << "\n";
            }
        }
    }
    else if (diagnostics_ && csvPtr_ && Pstream::master())
    {
        csvPtr_() << mesh_.time().timeOutputValue() << ","
            << outerCall_ << ",value-bulk,0,0,0\n";
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::semiImplicitCapillaryForce::semiImplicitCapillaryForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    form_(coeffs_.getOrDefault<word>("form", "off")),
    coeff_(coeffs_.getOrDefault<scalar>("coeff", 1.0)),
    laplaceBeltrami_(coeffs_.getOrDefault<bool>("laplaceBeltrami", true)),
    alphaName_(coeffs_.getOrDefault<word>("alphaName", "alpha.water")),
    kappaName_(coeffs_.getOrDefault<word>("kappaName", "kappa")),
    nHatName_(coeffs_.getOrDefault<word>("nHatName", "nHatFit")),
    sigma_
    (
        "sigma",
        dimMass/sqr(dimTime),
        mesh_.lookupObject<IOdictionary>("transportProperties")
    ),
    diagnostics_(coeffs_.getOrDefault<bool>("diagnostics", true)),
    csvPtr_(nullptr),
    lastTimeIndex_(-1),
    outerCall_(0)
{
    if (form_ != "off" && form_ != "value" && form_ != "increment")
    {
        FatalIOErrorInFunction(coeffs_)
            << "Unknown form '" << form_
            << "'. Valid: off, value, increment." << exit(FatalIOError);
    }

    coeffs_.readEntry("fields", fieldNames_);
    if (fieldNames_.size() != 1)
    {
        FatalIOErrorInFunction(coeffs_)
            << "settings are: " << fieldNames_
            << " -- exactly one field (U) expected." << exit(FatalIOError);
    }
    fv::option::resetApplied();

    if (form_ != "off")
    {
        Info<< "semiImplicitCapillaryForce: form " << form_
            << ", coeff " << coeff_
            << ", laplaceBeltrami " << laplaceBeltrami_
            << ", sigma " << sigma_.value() << nl << endl;

        if (diagnostics_ && Pstream::master())
        {
            csvPtr_.reset
            (
                new OFstream("semiImplicitCapillary.csv")
            );
            csvPtr_() << "TIME,outerCall,form,termL2,termLinf,dampingPower\n";
            csvPtr_().precision(10);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::semiImplicitCapillaryForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (form_ == "off")
    {
        return;
    }
    addTerm(eqn);
}


void Foam::fv::semiImplicitCapillaryForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // The one-field momentum equation is assembled in conservative form with
    // rho folded into ddt/div; the capillary linearisation is a FORCE term
    // and carries no rho, exactly like GSigma itself.
    addSup(eqn, fieldi);
}


bool Foam::fv::semiImplicitCapillaryForce::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        coeffs_.readIfPresent("form", form_);
        coeffs_.readIfPresent("coeff", coeff_);
        coeffs_.readIfPresent("laplaceBeltrami", laplaceBeltrami_);
        if (form_ != "off" && form_ != "value" && form_ != "increment")
        {
            FatalIOErrorInFunction(coeffs_)
                << "Unknown form '" << form_
                << "'. Valid: off, value, increment." << exit(FatalIOError);
        }
        return true;
    }
    return false;
}


// ************************************************************************* //
