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

#include "sdplsRdivStrictSp.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(sdplsRdivStrictSp, false);
addToRunTimeSelectionTable(sdplsSource, sdplsRdivStrictSp, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sdplsRdivStrictSp::sdplsRdivStrictSp(const dictionary& dict, const fvMesh& mesh)
:
    sdplsRdiv(dict, mesh)
{
    // NO mollifier check here ON PURPOSE. sdplsRdiv's constructor, which this
    // one delegates to, raises a FatalError for any mollifier type other than
    // `none`: a cut-off multiplying a FLUX breaks the conservation the
    // divergence form exists to provide, and that is as true of the sign-split
    // variant as of the original. Inheriting the refusal keeps ONE copy of the
    // rule; replicating it here would create a second copy that can drift.
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvScalarMatrix>
sdplsRdivStrictSp::fvmsdplsSource
(
    const volScalarField& psi,
    const volVectorField& U
)
{
    // ------------------------------------------------------------------
    // Everything down to and including phiW is IDENTICAL to sdplsRdiv::
    // fvmsdplsSource. It is repeated rather than factored out so that the two
    // arms can be diffed line by line and the single-variable claim in the
    // header can be checked by reading, not by trusting a shared helper.
    // ------------------------------------------------------------------

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

    // ------------------------------------------------------------------
    // THE ONE DIFFERENCE. sdplsRdiv delivers the algebraic coefficient
    //
    //     f = a - div(w)          [1/s]
    //
    // fully implicitly, as  - fvm::Sp(div(phiW) - a, psi) = + f psi^{n+1}.
    //
    // fvm::Sp(sp, psi) does diag += V*sp, and this matrix sits on the RIGHT of
    // the solver's equation, which fvMatrix::operator== assembles as
    // LHS - RHS. The system diagonal therefore receives -V*f: a POSITIVE f
    // REDUCES the diagonal, and once f exceeds 1/dt (per unit volume) the
    // diagonal changes sign. The linear solve still converges -- it is a
    // well-posed system -- and returns an amplifying, sign-reversing update.
    // That is the failure mode measured on the 2D reversed vortex: every run
    // COMPLETED (solverFailed blank) while the band gradient error at t = T/2
    // grew from 0.391 at N = 32 to 668.9 at N = 256 (fitted order -3.800),
    // with the shape order at t = T at -1.306 and the volume-conservation
    // order at -1.544. See the class header for the full table and protocol.
    //
    // Split f in the sense of Patankar (1980), exactly as
    // strictNegativeSpLinearImplicit splits f_nl:
    //
    //     f = max(f,0) + min(f,0),
    //     max(f,0) psi^n     delivered EXPLICITLY through fvm::Su,
    //     min(f,0) psi^{n+1} delivered IMPLICITLY through fvm::Sp.
    //
    // The implicit branch contributes diag += V*min(f,0) <= 0 to THIS matrix,
    // hence -V*min(f,0) >= 0 to the system diagonal: the source can now only
    // ADD to the diagonal, never subtract from it, on every cell, at every
    // step, on every mesh. Nothing is dropped: the positive branch is still
    // delivered in full, only at time level n.
    //
    // psi.oldTime() and NOT the current psi: the advection loop re-assembles
    // nDefCorr times per step with psi updated in place, so an explicit branch
    // evaluated at the current iterate would converge back onto the fully
    // implicit sdplsRdiv equation and the experiment would be a null by
    // construction. With psi^n the split survives the loop. With nDefCorr = 1
    // the two choices are identical. Argued at length in the header.
    // ------------------------------------------------------------------

    volScalarField const f
    (
        "fSdplsRdivStrictSp",
        strain() - fvc::div(phiW)
    );

    dimensionedScalar const fZero("fZero", f.dimensions(), 0.0);

    // >= 0 everywhere: the branch that would reduce the system diagonal.
    volScalarField const fPos("fPosSdplsRdivStrictSp", max(f, fZero));

    // <= 0 everywhere: the branch that is safe to keep implicit.
    volScalarField const fNeg("fNegSdplsRdivStrictSp", min(f, fZero));

    // WITHOUT THIS LINE THE EXPERIMENT CANNOT BE INTERPRETED. If f <= 0 in
    // (almost) every cell then min(f,0) = f, max(f,0) = 0, and this class
    // assembles EXACTLY what sdplsRdiv assembles. An identical result would
    // then mean "the coefficient was already sign-safe", NOT "the diagonal was
    // not the defect" -- and those two readings imply opposite next steps. The
    // decisive number is f*dt, because the one-cell diagonal is 1/dt - f per
    // unit volume: max(f)*dt >= 1 marks a cell whose fully implicit diagonal
    // has been driven through zero, max(f)*dt << 1 a cell where the implicit
    // treatment was never at risk. Two global reductions per assembly, printed
    // once per defect-correction pass.
    const scalar dt = mesh.time().deltaTValue();
    const scalarField& fInternal = f.primitiveField();

    label nPositive = 0;
    forAll(fInternal, cellI)
    {
        if (fInternal[cellI] > 0)
        {
            ++nPositive;
        }
    }
    reduce(nPositive, sumOp<label>());

    Info<< "sdplsRdivStrictSp: f = a - div(w); cells with f > 0 (explicit"
        << " branch active): " << nPositive << "/"
        << returnReduce(mesh.nCells(), sumOp<label>())
        << ", max(f) dt = " << gMax(fInternal)*dt
        << ", min(f) dt = " << gMin(fInternal)*dt << endl;

    // S_matrix = div(phiW, psi) - div(phi, psi)
    //            + max(f,0) psi^n + min(f,0) psi^{n+1}.
    //
    // The two div terms are byte-identical to sdplsRdiv's, on the same fluxes
    // and the same fvSchemes entries: this variant isolates the treatment of
    // the algebraic coefficient and nothing else.
    //
    // Sign check on the two source calls, in the convention of this matrix
    // (its represented source is (diag*psi - source)/V, see
    // discretization.C): fvm::Sp(sp,psi) gives diag = V*sp, hence + sp*psi;
    // fvm::Su(su,psi) gives source = -V*su, hence + su. Both therefore enter
    // with a PLUS, unlike sdplsRdiv's `- fvm::Sp(div(phiW) - a, psi)`, whose
    // minus sign is absorbed by writing the coefficient as f = a - div(w)
    // instead of c = div(w) - a.
    return
    (
        fvm::div(phiW, psi)
      - fvm::div(phi, psi)
      + fvm::Su(fPos*psi.oldTime(), psi)
      + fvm::Sp(fNeg, psi)
    );
}

} // End namespace Foam

// ************************************************************************* //
