/*---------------------------------------------------------------------------*\
License
    This file is part of the leia OpenFOAM module.
\*---------------------------------------------------------------------------*/

#include "slTaylorCellGradientFootIntegrator.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(slTaylorCellGradientFootIntegrator, 0);
    addToRunTimeSelectionTable
    (
        slFootIntegrator,
        slTaylorCellGradientFootIntegrator,
        Mesh
    );
}


Foam::slTaylorCellGradientFootIntegrator::slTaylorCellGradientFootIntegrator
(
    const fvMesh& mesh
)
:
    slFootIntegrator(mesh)
{}


Foam::tmp<Foam::tensorField>
Foam::slTaylorCellGradientFootIntegrator::cellGradient(const volVectorField& u)
{
    const fvMesh& mesh = u.mesh();
    const vectorField& C = mesh.C().primitiveField();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const label n = mesh.nCells();

    symmTensorField M(n, Zero);
    tensorField R(n, Zero);       // sum w d (u_n - u_c), d as the row index

    forAll(nei, f)
    {
        const label o = own[f], m = nei[f];
        const vector d = C[m] - C[o];
        const scalar w = 1.0/magSqr(d);
        const vector du = u[m] - u[o];
        M[o] += w*sqr(d); R[o] += w*(d*du);
        M[m] += w*sqr(d); R[m] += w*(d*du);    // (-d)(-du)
    }
    forAll(u.boundaryField(), patchi)
    {
        const fvPatchVectorField& pu = u.boundaryField()[patchi];
        if (!pu.coupled()) { continue; }
        const fvPatch& p = pu.patch();
        const labelUList& fc = p.faceCells();
        const vectorField uN(pu.patchNeighbourField());
        const vectorField delta(p.delta());          // owner centre to neighbour centre
        forAll(fc, i)
        {
            const label o = fc[i];
            const vector& d = delta[i];
            const scalar w = 1.0/magSqr(d);
            M[o] += w*sqr(d); R[o] += w*(d*(uN[i] - u[o]));
        }
    }

    // an empty direction (2D): unit diagonal, so that M is invertible and the
    // gradient has no component there (R is zero in that direction)
    const Vector<label>& sD = mesh.solutionD();
    auto tG = tmp<tensorField>::New(n, Zero);
    tensorField& G = tG.ref();
    forAll(G, c)
    {
        if (sD.x() < 0) { M[c].xx() += 1; }
        if (sD.y() < 0) { M[c].yy() += 1; }
        if (sD.z() < 0) { M[c].zz() += 1; }
        // G = inv(M) R: rows are d/dx_a of u
        G[c] = inv(M[c]) & R[c];
    }
    return tG;
}


void Foam::slTaylorCellGradientFootIntegrator::feet
(
    const volVectorField& uNew,
    const volVectorField& uOld,
    const scalar dt,
    pointField& feet
) const
{
    const volVectorField& C = mesh_.C();
    volVectorField uN("taylorCellGradientU", uNew);
    uN.correctBoundaryConditions();
    const tensorField G(cellGradient(uN));
    forAll(C, c)
    {
        const vector& un = uNew[c];
        const vector accel = (un - uOld[c])/dt + (un & G[c]);
        feet[c] = C[c] - un*dt + 0.5*accel*dt*dt;
    }
}

// ************************************************************************* //
