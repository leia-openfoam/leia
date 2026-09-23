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
    Foam::fv::semiImplicitCapillaryForce

Group
    grpFvOptionsSources

Description
    Semi-implicit (Hysing / Raessi / SAAMPLE sec. 3.4) capillary term for the
    momentum predictor, as an fvOption so ANY solver whose UEqn folds in
    fvOptions(rho, U) can select it from system/fvOptions with no solver code.

    The explicit CSF force evaluates sigma*kappa*n at the interface position of
    t^n. Linearising the force about the position the interface will have at
    t^{n+1} = x^n + dt*u^{n+1} adds
        dt * sigma * Lap_Gamma(u^{n+1}) * delta_Sigma,
    an interface-concentrated Laplace-Beltrami diffusion of the NEW velocity.
    Following SAAMPLE Eq. (23) the operator is split into one IMPLICIT bulk
    Laplacian plus EXPLICIT projector corrections (lagged by one outer
    iteration, i.e. a deferred correction whose converged limit is the true
    Lap_Gamma) -- legitimate only inside an iterative pressure-velocity loop.

    Two forms, selected by the `form` keyword:

    value (SAAMPLE / Raessi):
        eqn += fvm::laplacian(mu_f, U) - mu_c*normalCorrections(U)
      The term REMAINS in the converged momentum equation: it stabilises and
      it dissipates (the overdamping recorded in SAAMPLE sec. 5). Consistent
      pairing: interface FROZEN across the outer loop (PIMPLE.psiOuterCorrectors
      no) -- the term substitutes for the missing force update.

    increment (deferred Jacobian):
        eqn += fvm::laplacian(mu_f, U) - fvc::laplacian(mu_f, U)
      U on the fvc side is the CURRENT outer iterate u^(k), so the term acts on
      the outer-loop increment u^{n+1} - u^(k) and VANISHES to round-off at
      outer-loop convergence: zero artificial dissipation in the converged
      momentum equation. The projector corrections cancel identically between
      the implicit and lagged parts, so this form needs neither kappa nor the
      interface normal. On the first outer pass u^(0) = u^n, so the full
      sigma*dt damping acts as a continuation step and then self-ramps to zero
      -- no tuned decay factor. Consistent pairing: interface IN the outer loop
      (PIMPLE.psiOuterCorrectors yes) with residual-controlled outer
      iterations, where it is the Jacobian dt*dF/dx of the capillary coupling
      placed in the iteration matrix -- the preconditioner that makes the
      Picard iteration contract at capillary-scale time steps.

    The face coefficient is
        mu_f = coeff * sigma * dt * |snGrad(alpha)|_f
    -- the SAME discrete object the balanced CSF flux GSigma differences, so
    the damping support coincides with the force support exactly. No smoothing,
    no fvc::reconstruct, no artificial coefficient: sigma is read from
    transportProperties (fail-loudly), dt is the run's time step. The explicit
    projector corrections (value form only) take kappa and nHat from REGISTERED
    fields named in the dictionary -- for the SL solvers these are the quadratic
    fit's normal-extended curvature (`kappa`) and fit normal (`nHatFit`),
    consistent with the tubular-neighbourhood derivation (both constant along
    interface normals under the parallel-foliation hypothesis the SL chain
    already assumes).

    This is a linearisation of the capillary force, not a filter: the term is
    proportional to (an increment of) U and vanishes identically at the
    stationary Young-Laplace equilibrium, so it can neither shift the balance
    nor inject at grid scale where nothing moves.

Usage
    system/fvOptions:
    \verbatim
    semiImplicitCapillary
    {
        type            semiImplicitCapillaryForce;
        fields          (U);
        form            off;        // off | value | increment
        coeff           1.0;        // multiplies sigma*dt (2/3 = BDF2 Jacobian)
        laplaceBeltrami true;       // value form: subtract projector terms
        alphaName       alpha.water;
        kappaName       kappa;      // value form + laplaceBeltrami only
        nHatName        nHatFit;    // value form + laplaceBeltrami only
        diagnostics     true;       // per-outer-call CSV (master only)
    }
    \endverbatim

SourceFiles
    semiImplicitCapillaryForce.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitCapillaryForce_H
#define semiImplicitCapillaryForce_H

#include "fvOption.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class semiImplicitCapillaryForce Declaration
\*---------------------------------------------------------------------------*/

class semiImplicitCapillaryForce
:
    public fv::option
{
    // Private Data

        //- off | value | increment
        word form_;

        //- Multiplies sigma*dt in the coefficient (1 = backward Euler
        //  linearisation; 2/3 = the exact BDF2 Jacobian scale).
        scalar coeff_;

        //- value form: subtract the explicit projector corrections so the
        //  deferred-corrected operator is Lap_Gamma, not the bulk Laplacian.
        bool laplaceBeltrami_;

        word alphaName_;
        word kappaName_;
        word nHatName_;

        //- Surface tension coefficient, read once from the registered
        //  transportProperties (fail-loudly: no silent default).
        dimensionedScalar sigma_;

        //- Per-outer-call diagnostics CSV (master only)
        bool diagnostics_;
        autoPtr<OFstream> csvPtr_;
        label lastTimeIndex_;
        label outerCall_;


    // Private Member Functions

        //- The face coefficient mu_f = coeff*sigma*dt*|snGrad(alpha)|
        tmp<surfaceScalarField> muFace(const volScalarField& alpha) const;

        //- Assemble and add the selected form
        void addTerm(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("semiImplicitCapillaryForce");


    // Constructors

        //- Construct from explicit source name and mesh
        semiImplicitCapillaryForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        semiImplicitCapillaryForce(const semiImplicitCapillaryForce&) = delete;

        //- No copy assignment
        void operator=(const semiImplicitCapillaryForce&) = delete;


    //- Destructor
    virtual ~semiImplicitCapillaryForce() = default;


    // Member Functions

        //- Add implicit/explicit contributions to the momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Compressible/one-field variant (what UEqn's fvOptions(rho,U) calls)
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Read source dictionary (runtime-modifiable, unlike the
        //  sepMultiphaseFoam ancestor which aborts here)
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
