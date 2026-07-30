/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Julian Reitzel, TU Darmstadt
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

Application
    leiaTestRedistance

Description
    Apply the runtime-selected levelSet.redistancer ONCE to the read psi and
    write the result at time 1.

    If fvSolution has a levelSet.implicitSurface entry, the ANALYTIC signed
    distance psiExact = value(x)/|grad(value)(x)| (exact for planes at any
    'gradient' scaling and for spheres) is the error reference, and a
    single-row leiaTestRedistance.csv is written with pre/post redistancing
    error norms over the near band B = {|psiExact| <= 6 h}:

        H, E_LINF_BAND_PSI_PRE, E_LINF_BAND_PSI, E_L1_BAND_PSI,
        E_LINF_ALL_PSI, E_LINF_BAND_GRAD_PSI_PRE, E_LINF_BAND_GRAD_PSI,
        E_MEAN_BAND_GRAD_PSI

    (aggregated by the snakemake test-suite, e.g. the static
    config/redistanceCircle2D.yaml convergence study).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "redistancer.H"
#include "narrowBand.H"
#include "phaseIndicator.H"
#include "levelSetImplicitSurfaces.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct bandErrors
{
    scalar psiLinfBand;
    scalar psiL1Band;
    scalar psiLinfAll;
    scalar gradLinfBand;
    scalar gradMeanBand;
};


bandErrors computeBandErrors
(
    const volScalarField& psi,
    const scalarField& psiExact,
    const boolList& inBand
)
{
    const fvMesh& mesh = psi.mesh();
    const volVectorField gradPsi(fvc::grad(psi));

    bandErrors e{0, 0, 0, 0, 0};
    scalar l1Sum = 0, vSum = 0, gSumDev = 0;

    forAll(psiExact, celli)
    {
        const scalar ePsi = mag(psi[celli] - psiExact[celli]);
        e.psiLinfAll = max(e.psiLinfAll, ePsi);

        if (!inBand[celli]) continue;

        const scalar V = mesh.V()[celli];
        const scalar eGrad = mag(mag(gradPsi[celli]) - 1.0);

        e.psiLinfBand = max(e.psiLinfBand, ePsi);
        e.gradLinfBand = max(e.gradLinfBand, eGrad);
        l1Sum += V*ePsi;
        gSumDev += V*eGrad;
        vSum += V;
    }

    reduce(e.psiLinfBand, maxOp<scalar>());
    reduce(e.psiLinfAll, maxOp<scalar>());
    reduce(e.gradLinfBand, maxOp<scalar>());
    reduce(l1Sum, sumOp<scalar>());
    reduce(gSumDev, sumOp<scalar>());
    reduce(vSum, sumOp<scalar>());

    e.psiL1Band = (vSum > VSMALL) ? l1Sum/vSum : 0;
    e.gradMeanBand = (vSum > VSMALL) ? gSumDev/vSum : 0;

    return e;
}

} // End namespace Foam


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Redistance application."
    );



    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Analytic signed-distance reference from levelSet.implicitSurface (the
    // same entry leiaSetFields initialized psi from; value/|grad value| is
    // the exact SDF for planes and spheres regardless of profile scaling).
    const fvSolution& fvSolutionDict(mesh);
    const dictionary& levelSetDict = fvSolutionDict.subDict("levelSet");
    const bool haveExact = levelSetDict.found("implicitSurface");

    scalarField psiExact(mesh.nCells(), 0.0);
    boolList inBand(mesh.nCells(), false);
    Foam::bandErrors pre{0, 0, 0, 0, 0};

    if (haveExact)
    {
        const dictionary& surfDict = levelSetDict.subDict("implicitSurface");
        autoPtr<implicitSurface> surf =
            implicitSurface::New(surfDict.get<word>("type"), surfDict);

        const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());
        const vectorField& C = mesh.cellCentres();
        forAll(psiExact, celli)
        {
            const scalar v = surf->value(C[celli]);
            const scalar g = mag(surf->grad(C[celli]));
            psiExact[celli] = v/max(g, SMALL);
            inBand[celli] = (mag(psiExact[celli]) <= 6.0*h);
        }

        pre = Foam::computeBandErrors(psi, psiExact, inBand);
    }

    // Pre-event volume fraction (input reference for the volume/shape
    // metrics of ONE redistance event: on a signed-distance input a perfect
    // redistancer is IDEMPOTENT -- zero volume shift, zero shape change).
    phaseInd->calcPhaseIndicator(alpha, psi);
    const volScalarField alpha0("alpha0", alpha);
    const scalar Valpha0 =
        gSum(alpha.primitiveField()*mesh.V().field());

    redist->redistance(psi);

    // Post-event: refresh the band (cut cells may realign sub-h) + alpha.
    nBand->calc();
    phaseInd->calcPhaseIndicator(alpha, psi);
    const scalar Valpha1 =
        gSum(alpha.primitiveField()*mesh.V().field());
    const scalar eVol = mag(Valpha1 - Valpha0);
    const scalar eVolRel = eVol/max(Valpha0, VSMALL);
    const scalar eGeom =
        gSum(mag(alpha.primitiveField() - alpha0.primitiveField())
             *mesh.V().field());

    runTime.setTime(1,1);

    psi.write();
    alpha.write();

    if (haveExact)
    {
        const Foam::bandErrors post =
            Foam::computeBandErrors(psi, psiExact, inBand);
        const scalar h = 1.0/gMax(mesh.deltaCoeffs().primitiveField());

        Info<< "leiaTestRedistance: band Linf(psi - psiExact) "
            << pre.psiLinfBand << " -> " << post.psiLinfBand
            << ", band mean(||grad psi|-1|) "
            << pre.gradMeanBand << " -> " << post.gradMeanBand
            << ", volume " << Valpha0 << " -> " << Valpha1
            << " (E_VOL " << eVol << ", E_GEOM " << eGeom << ")" << endl;

        if (Pstream::master())
        {
            OFstream csv("leiaTestRedistance.csv");
            csv << "H,"
                << "E_LINF_BAND_PSI_PRE,"
                << "E_LINF_BAND_PSI,"
                << "E_L1_BAND_PSI,"
                << "E_LINF_ALL_PSI,"
                << "E_LINF_BAND_GRAD_PSI_PRE,"
                << "E_LINF_BAND_GRAD_PSI,"
                << "E_MEAN_BAND_GRAD_PSI,"
                << "E_VOL_ALPHA,"
                << "E_VOL_ALPHA_REL,"
                << "E_GEOM_ALPHA\n";
            csv << h << ","
                << pre.psiLinfBand << ","
                << post.psiLinfBand << ","
                << post.psiL1Band << ","
                << post.psiLinfAll << ","
                << pre.gradLinfBand << ","
                << post.gradLinfBand << ","
                << post.gradMeanBand << ","
                << eVol << ","
                << eVolRel << ","
                << eGeom << "\n";
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
