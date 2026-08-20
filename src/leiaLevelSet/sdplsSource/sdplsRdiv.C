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

#include "sdplsRdiv.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sdplsRdiv, false);
addToRunTimeSelectionTable(sdplsSource, sdplsRdiv, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sdplsRdiv::sdplsRdiv(const dictionary& dict, const fvMesh& mesh)
:
    sdplsR(dict, mesh)
{
    // The conservative assembly cannot compose with a cut-off: a mollifier
    // multiplying a FLUX destroys the conservation this form exists to
    // provide. Refuse loudly rather than run something the dictionary does
    // not describe.
    const word mollType =
        dict.subOrEmptyDict("mollifier").getOrDefault<word>("type", "none");
    if (mollType != "none")
    {
        FatalErrorInFunction
            << "sdplsSource type `Rdiv` is the CONSERVATIVE (divergence-form) "
            << "source; a mollifier (`" << mollType << "`) multiplying its "
            << "fluxes would break exactly the conservation it provides." << nl
            << "Use mollifier type `none` with Rdiv, or type `R` if a cut-off "
            << "is wanted."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvScalarMatrix>
sdplsRdiv::fvmsdplsSource(const volScalarField& psi, const volVectorField& U)
{
    // Refresh R_ (the strain a = n.grad(v).n) and the lagged normal, exactly
    // as the algebraic source does -- same Picard structure, measured
    // converged at nDefCorr = 3.
    update(psi, U);

    const fvMesh& mesh = psi.mesh();

    // The lagged interface normal, from the same dedicated gradient the
    // strain uses (gradPsi model + gradPsiSdpls scheme).
    volVectorField const gradPsi(grad(psi));
    dimensioned<scalar> const eps =
        dimensioned<scalar>("eps", gradPsi.dimensions(), SMALL);
    volVectorField const nHat(gradPsi/(mag(gradPsi) + eps));

    // Normal velocity w = (n.v) n and its face flux. The flux is NAMED so
    // fvSchemes can give div(phiW,psi) its own scheme; the solver's own mass
    // flux phi is looked up rather than rebuilt from U, so the same flux that
    // advects psi in the base equation appears in the correction terms.
    volVectorField const w(("wSdpls"), (nHat & U)*nHat);
    surfaceScalarField phiW("phiW", fvc::flux(w));

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>("phi");

    // S_matrix = div(phiW, psi) - div(phi, psi) - Sp(div(phiW) - a, psi).
    //
    // Derivation (see the header): substituting S = -div_s(psi v) into the
    // base equation and moving everything the base already assembles back to
    // the left leaves exactly these three terms on the right. At the continuum
    // they collapse to psi (a - div v): the algebraic source plus the same
    // compressibility correction the base equation carries. Discretely the
    // two div terms are CONSERVATIVE -- their global sum telescopes to
    // boundary fluxes -- which is the property the algebraic source lacks and
    // the coupled droplet measurements show it needs.
    // ---- PRODUCT-RULE RESIDUAL DIAGNOSTIC (added 2026-08-20) ----------------
    //
    // WHY THIS IS MEASURED HERE AND NOT INFERRED FROM A CONVERGENCE TABLE.
    // Rdiv and the algebraic source R are CONTINUUM-IDENTICAL for incompressible
    // flow, yet on the 2D reversed vortex with a prescribed divergence-free
    // velocity -- no capillary force, no curvature, no feedback of any kind --
    // R converges at +1.112 in the band gradient error while Rdiv DIVERGES at
    // -3.800 (0.391, 0.256, 1.56, 8.00, 20.5, 103, 669 over N = 32..256, every
    // run COMPLETING, relative volume error reaching 136% at N=256). The only
    // structural difference between them is that Rdiv routes part of the
    // operator through fvm::div FLUXES. So the defect must lie in that routing.
    //
    // The sign-split variant RdivStrictSp tested and FALSIFIED the first
    // hypothesis, that the unsigned Sp coefficient f = a - div(w) destroys
    // diagonal dominance: it demonstrably removed the M-matrix violation (512 of
    // 1024 cells carried f > 0 at N=32) and the divergence barely moved (-3.554
    // against -3.800, 635 against 669 at the finest rung).
    //
    // WHAT THIS RESIDUAL TESTS. At the continuum the product rule gives
    //     div(w psi) - psi div(w) - w . grad(psi) = 0
    // identically. Rdiv RELIES on that identity: the psi*div(w) content carried
    // implicitly inside fvm::div(phiW, psi) is supposed to cancel against the
    // -Sp(div(phiW), psi) term, leaving only the convective part. Discretely
    // there is NO guarantee that it does -- fvm::div(phiW, psi) interpolates psi
    // to faces with the div(phiW,psi) scheme while fvc::div(phiW) is a plain
    // flux divergence, and the two need not be consistent. If this residual
    // fails to vanish, and worse GROWS with refinement, then the cancellation
    // Rdiv is built on does not happen discretely, and that is the defect.
    //
    // The reference norm is printed alongside deliberately: the two terms nearly
    // cancel, so a residual is only meaningful relative to the magnitude of what
    // is being differenced.
    {
        volScalarField const prDiv("prDiv", fvc::div(phiW, psi));
        volScalarField const prPsiDivW("prPsiDivW", psi*fvc::div(phiW));
        volScalarField const prConv("prConv", w & gradPsi);
        volScalarField const prResid("prResid", prDiv - prPsiDivW - prConv);

        const scalarField& r = prResid.primitiveField();
        const scalarField& d = prDiv.primitiveField();
        const scalarField& V = mesh.V();
        const scalar Vtot = gSum(V);

        Info<< "sdplsRdiv product-rule residual"
            << " |div(w psi) - psi div w - w.grad psi|:"
            << " Linf = " << gMax(mag(r))
            << ", L2 = " << Foam::sqrt(gSum(sqr(r)*V)/Vtot)
            << " | reference |div(w psi)|: Linf = " << gMax(mag(d))
            << ", L2 = " << Foam::sqrt(gSum(sqr(d)*V)/Vtot)
            << endl;
    }
    // -------------------------------------------------------------------------

    return
    (
        fvm::div(phiW, psi)
      - fvm::div(phi, psi)
      - fvm::Sp(fvc::div(phiW) - strain(), psi)
    );
}

} // End namespace Foam

// ************************************************************************* //
