/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slRk2FootIntegrator.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
    defineTypeNameAndDebug(slRk2FootIntegrator, 0);
    addToRunTimeSelectionTable(slFootIntegrator, slRk2FootIntegrator, Mesh);
}


Foam::slRk2FootIntegrator::slRk2FootIntegrator(const fvMesh& mesh)
:
    slFootIntegrator(mesh)
{}


void Foam::slRk2FootIntegrator::feet
(
    const volVectorField& uNew,
    const volVectorField& uOld,
    const scalar dt,
    pointField& feet
) const
{
    // The rk2 branch of pointValueScheme, moved here unchanged.
    const volVectorField& C = mesh_.C();
    const volVectorField uHalf
    (
        "uHalfTrace",
        0.5*(uNew + uOld)
    );
    interpolationCellPoint<vector> uHalfInterp(uHalf);
    const vectorField& uH = uHalf.primitiveField();

    label nOffProc = 0;
    forAll(C, c)
    {
        const point xm = C[c] - 0.5*dt*uH[c];
        // The midpoint lies at most half a cell away at CFL <= 1, so the
        // arrival cell is nearly always the right one; findCell returns -1
        // only when it has crossed a processor boundary, where evaluating
        // in cell c is a local extrapolation of under half a cell rather
        // than a wrong-rank lookup. Counted, never silent.
        label cm = mesh_.findCell(xm);
        if (cm < 0) { cm = c; ++nOffProc; }
        const vector um = uHalfInterp.interpolate(xm, cm);
        feet[c] = C[c] - dt*(std::isfinite(magSqr(um)) ? um : uH[c]);
    }
    reduce(nOffProc, sumOp<label>());
    if (nOffProc > 0 && mesh_.time().writeTime())
    {
        Info<< "    rk2 foot: " << nOffProc
            << " midpoints outside the local decomposition, evaluated in"
            << " the arrival cell" << endl;
    }
}

// ************************************************************************* //
