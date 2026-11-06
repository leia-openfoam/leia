/*---------------------------------------------------------------------------*\
    leia :: davofAdvectorFoam :: createFields.H: the fields, the mixture, the advection scheme, setFlow and the abort guard
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
    From the isoAdvector benchmark solver advectorVoF (Johan Roenby, DHI,
    2016; TwoPhaseFlow apps/benchmark/advectorVoF): J. Roenby, H. Bredmose,
    H. Jasak, A computational method for sharp interface advection, R. Soc.
    Open Sci. 3 (2016) 160405, doi:10.1098/rsos.160405
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

scalar V0 = fvc::domainIntegrate(alpha1).value();

Info << "sum(alpha*V) = " << V0
    << ", 1-max(alpha1) = " << 1 - gMax(alpha1.internalField())
    << "\t min(alpha1) = " << gMin(alpha1.internalField()) << endl;

//IsoAdvector related stuff
//isoAdvection advectorIso(alphaIso,phi,U);
// Construct incompressible turbulence model
autoPtr<davofAdvectionScheme> advector
(
    davofAdvectionScheme::New(alpha1,phi,U)
);

IOdictionary fvSolutionDict
(
    IOobject
    (
        "fvSolution",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


autoPtr<functionObjects::setFlow> flow;

if (fvSolutionDict.found("setFlow"))
{
    dictionary setFlowDict = fvSolutionDict.subDict("setFlow");
    flow.reset(new functionObjects::setFlow("setFlow",runTime,setFlowDict));
}

const dictionary& isoAdvectorDict = mesh.solutionDict().subDict("isoAdvector");
const scalar period = isoAdvectorDict.lookupOrDefault<scalar>("period", 0.0);
scalar reverseTime = isoAdvectorDict.lookupOrDefault<scalar>("reverseTime", 0.0);
scalar spirallingFlow = isoAdvectorDict.lookupOrDefault<scalar>("spirallingFlow", 0.0);

// The abort guard of the time loop (keys of the same isoAdvector dictionary;
// off by default, so that the testsuite is unchanged): the run stops with a
// normal "End" after an ABORT line when alpha leaves [abortMinAlpha,
// abortMaxAlpha] or the interface-cell count of the reconstruction scheme
// exceeds abortInterfaceCellFactor times its count at the first step
// (factor 0 disables that criterion).
const bool abortEnabled =
    isoAdvectorDict.lookupOrDefault<Switch>("abortEnabled", false);
const scalar abortMinAlpha =
    isoAdvectorDict.lookupOrDefault<scalar>("abortMinAlpha", -0.5);
const scalar abortMaxAlpha =
    isoAdvectorDict.lookupOrDefault<scalar>("abortMaxAlpha", 1.5);
const scalar abortInterfaceCellFactor =
    isoAdvectorDict.lookupOrDefault<scalar>("abortInterfaceCellFactor", 10);
label nInterface0 = -1;

//Copies of initial U and phi for use when flow is periodic
surfaceScalarField phi0 = phi;
volVectorField U0 = U;
