/*---------------------------------------------------------------------------*\
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, GPL-3.0-or-later.

Application
    davofInitContactLine

Description
    Writes the velocity U and the volumetric flux phi of the linear,
    divergence-free velocity field of the 2D contact line advection tests of
        M. Fricke, T. Maric, D. Bothe,
        Contact line advection using the geometrical Volume-of-Fluid method,
        Journal of Computational Physics 407 (2020) 109221,
        doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785
    at t = 0 in the x-z plane of a one-cell-thick mesh (x1 = x, x2 = z),
        v = (v0 + c1 x1 + c2 x2, -c1 x2)            (Eq. 10; Eqs. 36, 37)
    with v0, c1, c2 from system/contactLineDict (the time factor cos(pi t/tau)
    of Eq. 37 is applied by davofAdvectorFoam, spirallingFlow tau). The
    internal U is the field at the cell centres, every non-empty boundary
    patch carries the field at its face centres (forced assignment, also on
    the fixedValue wall), phi is the field at the face centres times the face
    area vector, which is exact for a linear field on planar faces and
    divergence-free to round-off. Reads U from the time directory (its
    boundary conditions). Prints the reference speed s = max over the cells
    of 0.5 sum_f |phi_f| / V_c times V_c^(1/3) and the largest
    |sum_f phi_f| / V_c.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write U and phi of the linear velocity of the contact line tests"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "contactLineDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const scalar v0 = dict.get<scalar>("v0");
    const scalar c1 = dict.get<scalar>("c1");
    const scalar c2 = dict.get<scalar>("c2");

    auto vel = [&](const point& p) -> vector
    {
        return vector(v0 + c1*p.x() + c2*p.z(), 0, -c1*p.z());
    };

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

    forAll(U, celli)
    {
        U[celli] = vel(mesh.C()[celli]);
    }
    volVectorField::Boundary& Ub = U.boundaryFieldRef();
    forAll(Ub, patchi)
    {
        if (isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi])) continue;
        const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
        vectorField Uf(Cf.size());
        forAll(Cf, i)
        {
            Uf[i] = vel(Cf[i]);
        }
        Ub[patchi] == Uf;
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVol/dimTime, Zero)
    );
    forAll(phi.primitiveFieldRef(), facei)
    {
        phi.primitiveFieldRef()[facei] =
            vel(mesh.Cf()[facei]) & mesh.Sf()[facei];
    }
    surfaceScalarField::Boundary& phib = phi.boundaryFieldRef();
    forAll(phib, patchi)
    {
        if (isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi])) continue;
        const label start = mesh.boundaryMesh()[patchi].start();
        forAll(phib[patchi], i)
        {
            phib[patchi][i] =
                vel(mesh.Cf()[start + i]) & mesh.Sf()[start + i];
        }
    }

    scalar s = 0;
    scalar divMax = 0;
    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    scalarField sumAbs(mesh.nCells(), Zero), sumPhi(mesh.nCells(), Zero);
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const scalar p = phi.primitiveField()[facei];
        sumAbs[own[facei]] += mag(p); sumAbs[nei[facei]] += mag(p);
        sumPhi[own[facei]] += p; sumPhi[nei[facei]] -= p;
    }
    forAll(phib, patchi)
    {
        const label start = mesh.boundaryMesh()[patchi].start();
        forAll(phib[patchi], i)
        {
            const scalar p = phib[patchi][i];
            sumAbs[own[start + i]] += mag(p); sumPhi[own[start + i]] += p;
        }
    }
    forAll(sumAbs, celli)
    {
        const scalar V = mesh.V()[celli];
        s = max(s, 0.5*sumAbs[celli]/V*Foam::cbrt(V));
        divMax = max(divMax, mag(sumPhi[celli])/V);
    }
    reduce(s, maxOp<scalar>());
    reduce(divMax, maxOp<scalar>());

    U.write();
    phi.write();

    // The wall velocity written (the check of the boundary values)
    forAll(Ub, patchi)
    {
        if (Ub[patchi].size() == 0) continue;
        Info<< "davofInitContactLine: patch " << mesh.boundaryMesh()[patchi].name()
            << " type " << Ub[patchi].type() << " fixesValue "
            << Ub[patchi].fixesValue() << endl;
    }

    Info<< "davofInitContactLine: v0 " << v0 << " c1 " << c1 << " c2 " << c2
        << ", reference speed (max over the cells of 0.5 sum_f |phi_f| / V"
        << " times V^(1/3)) " << s
        << ", max |sum_f phi_f| / V " << divMax << endl;

    return 0;
}
