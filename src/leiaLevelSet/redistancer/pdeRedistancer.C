/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

#include "pdeRedistancer.H"
#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
// #include "fvCFD.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //
 
namespace Foam
{

static tmp<volScalarField> mysign(const volScalarField& field)
{
    tmp<volScalarField> tfield
    (
        new volScalarField
        (
            IOobject
            (
                word(),
                fileName(),
                field.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            field.mesh(),
            dimensionedScalar(dimless, 0.0)
        )
    );

    forAll(field, ID)
    {
        tfield.ref()[ID] = Foam::sign(field[ID]);
    }

    return tfield;
} 


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdeRedistancer, false);
addToRunTimeSelectionTable(redistancer, pdeRedistancer, Mesh);

} // End namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pdeRedistancer::pdeRedistancer(const fvMesh& mesh)
    :
        redistancer(mesh)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pdeRedistancer::redistance(volScalarField& psi) 
{
    dictionary pseudoTime_controlDict;

    pseudoTime_controlDict.add("startTime", 0);
    pseudoTime_controlDict.add("deltaT", 1);
    pseudoTime_controlDict.add("endTime", 4);
    pseudoTime_controlDict.add("writeControl", "timeStep");
    pseudoTime_controlDict.add("writeInterval", 1000);
    pseudoTime_controlDict.add("log", 2);

    const Time& runTime = psi.mesh().time(); 
    Time pseudoTime
        (
            pseudoTime_controlDict,
            runTime.rootPath(),
            runTime.caseName(),
            false,
            true
        );
    
    // Foam::autoPtr<Foam::Time> pseudoTime_ptr = Foam::Time::New();
    // Time& pseudoTime = pseudoTime_ptr.ref();

    Foam::fvMesh mesh_p
        (
            Foam::IOobject
            (
                Foam::polyMesh::defaultRegion,
                pseudoTime.timeName(),
                pseudoTime,
                Foam::IOobject::MUST_READ
            ),
            false
        );
    mesh_p.init(true);

    volScalarField psi_p = volScalarField
    (
        Foam::IOobject
        (
        "psi_p",
        pseudoTime.timeName(),
        pseudoTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh_p,
        psi.dimensions()
    );

    forAll(psi, ID)
    {
        psi_p[ID] = psi[ID];
    }

    fvScalarMatrix redistanceEqn
    (
        fvm::ddt(psi_p)
    ==
        sign(psi)*
        (1 - mag(fvc::grad(psi_p)))
        *dimensioned<scalar>(dimLength/dimTime, 1.0)
    );

    dictionary const& solverControls = fvSolution_.subDict("solvers").subDict(psi.name());

    redistanceEqn.solve(solverControls);

    // psi == psi_p; // Won't work: Different meshes
    forAll(psi, ID)
    {
        psi[ID] = psi_p[ID];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
