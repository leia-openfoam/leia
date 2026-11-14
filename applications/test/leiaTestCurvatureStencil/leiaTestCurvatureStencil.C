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

Application
    leiaTestCurvatureStencil

Description
    How many cell layers of psi does the UFVM curvature
        kappa = div(grad(psi)/|grad(psi)|)
    read at the interface? (plan docs/plan-velocity-extension-verification.md,
    step A0.4; supervisor 2026-11-14: the band width is set by the stencil of
    the curvature operator, not by a free parameter.)

    Input: psi at the start time (the exact signed distance written by
    leiaSetFields) and levelSet.implicitSurface (implicitSphere: centre).
    The layer index l is topological: l = 0 on the cells with a face neighbour
    of the other sign of psi, l = k on the k-th face dilation. No distance is
    assumed.

    For each gradient scheme of levelSet.curvatureStencil.gradSchemes (fvSchemes
    keywords) and each perturbation, psi_k = psi + eps h r(c) on every cell
    with l >= k (or outside the dilation), with
        random: r in [-1, 1], reproducible per cell index (splitmix64),
        step:   r = 1 (an O(h) offset with a jump at the layer k edge).
    kappa = surfaceIntegrate(linearInterpolate(n) & Sf), n = g/|g|, g the
    scheme's gradient (Gauss linear divergence, no divSchemes entry needed).
    Exact: kappa* = (D - 1)/|x - x0|, the curvature of the level set of the
    signed distance through x, D = mesh.nSolutionD().

    CSV leiaTestCurvatureStencil.csv, one row per (scheme, type, k), k = 0 is
    the unperturbed baseline:
      E_KAPPA_CELL_L1/L2  over the cells with l = 0 (volume weighted)
      E_KAPPA_FACE_L1/L2  over the faces with an l = 0 owner or neighbour, of
                          linearInterpolate(kappa) (area weighted); these are
                          the faces where snGrad(alpha) is nonzero
      D_KAPPA_CELL_L1     L1 of kappa_k - kappa_0 over the l = 0 cells
      N_CHANGED_CELL/FACE the number of l = 0 cells / those faces whose kappa
                          differs from the baseline at all (exact test: 0)
    SERIAL only (internal faces).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

namespace
{
    // splitmix64 of (cell, seed) -> [-1, 1]
    double unitRandom(const uint64_t c, const uint64_t seed)
    {
        uint64_t z = c*0x9E3779B97F4A7C15ULL + seed + 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30))*0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27))*0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        return 2.0*(double(z >> 11)*(1.0/9007199254740992.0)) - 1.0;
    }
}

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField psi
    (
        IOobject("psi", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    const fvSolution& fvSol(mesh);
    const dictionary& lsDict = fvSol.subDict("levelSet");
    const dictionary& sDict = lsDict.subDict("implicitSurface");
    const point x0 = sDict.get<point>("center");
    const dictionary& tDict = lsDict.subOrEmptyDict("curvatureStencil");
    const wordList schemes =
        tDict.getOrDefault<wordList>("gradSchemes", wordList({"grad(psi)"}));
    const scalar eps = tDict.getOrDefault<scalar>("eps", 1.0);
    const label kMax = tDict.getOrDefault<label>("kMax", 6);
    const label seed = tDict.getOrDefault<label>("seed", 1);
    const scalar D = mesh.nSolutionD();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.Cf();
    const scalarField& V = mesh.V();
    const scalarField& magSf = mesh.magSf();

    // cell size: minimum face-neighbour centre distance
    scalarField h(mesh.nCells(), GREAT);
    forAll(own, f)
    {
        const scalar d = mag(C[nei[f]] - C[own[f]]);
        h[own[f]] = min(h[own[f]], d);
        h[nei[f]] = min(h[nei[f]], d);
    }

    // topological layer index from the sign change of psi
    labelList layer(mesh.nCells(), -1);
    forAll(own, f)
    {
        if (psi[own[f]]*psi[nei[f]] <= 0)
        {
            layer[own[f]] = 0;
            layer[nei[f]] = 0;
        }
    }
    for (label pass = 1; pass <= kMax + 2; ++pass)
    {
        labelList next(layer);
        forAll(own, f)
        {
            const label o = own[f], n = nei[f];
            if (layer[o] >= 0 && next[n] < 0) next[n] = pass;
            if (layer[n] >= 0 && next[o] < 0) next[o] = pass;
        }
        layer.transfer(next);
    }

    auto kExact = [&](const point& x)
    {
        vector r = x - x0;
        if (D < 3)
        {
            // the empty direction carries no curvature
            const Vector<label>& sol = mesh.solutionD();
            for (direction i = 0; i < vector::nComponents; ++i)
            {
                if (sol[i] < 0) r[i] = 0;
            }
        }
        return (D - 1)/max(mag(r), SMALL);
    };

    label nCell0 = 0;
    forAll(layer, c) if (layer[c] == 0) ++nCell0;
    label nFace0 = 0;
    forAll(own, f) if (layer[own[f]] == 0 || layer[nei[f]] == 0) ++nFace0;
    const scalar hRef = gAverage(h);

    OFstream csv("leiaTestCurvatureStencil.csv");
    csv << "DELTA_X,GRAD_SCHEME,TYPE,K,E_KAPPA_CELL_L1,E_KAPPA_CELL_L2,"
        << "E_KAPPA_FACE_L1,E_KAPPA_FACE_L2,D_KAPPA_CELL_L1,N_CHANGED_CELL,"
        << "N_CHANGED_FACE,N_CELL,N_FACE\n";

    for (const word& scheme : schemes)
    {
        scalarField k0;        // the baseline curvature (cells)
        scalarField kf0;       // the baseline face curvature

        for (const word& type : wordList({"random", "step"}))
        {
            for (label k = 0; k <= kMax; ++k)
            {
                if (k == 0 && type == "step" && k0.size()) continue;

                volScalarField psiK("psi", psi);
                if (k > 0)
                {
                    forAll(psiK, c)
                    {
                        if (layer[c] >= k || layer[c] < 0)
                        {
                            const scalar r = (type == "random")
                                ? unitRandom(uint64_t(c), uint64_t(seed)) : 1.0;
                            psiK[c] += eps*h[c]*r;
                        }
                    }
                    psiK.correctBoundaryConditions();
                }

                const volVectorField g(fvc::grad(psiK, scheme));
                const volVectorField n
                (
                    g/(mag(g) + dimensionedScalar(g.dimensions(), SMALL))
                );
                const volScalarField kappa
                (
                    fvc::surfaceIntegrate(linearInterpolate(n) & mesh.Sf())
                );
                const surfaceScalarField kf(linearInterpolate(kappa));

                scalar e1 = 0, e2 = 0, vs = 0, d1 = 0;
                label nCh = 0;
                forAll(layer, c)
                {
                    if (layer[c] != 0) continue;
                    const scalar e = mag(kappa[c] - kExact(C[c]));
                    e1 += e*V[c]; e2 += e*e*V[c]; vs += V[c];
                    if (k0.size())
                    {
                        const scalar dk = mag(kappa[c] - k0[c]);
                        d1 += dk*V[c];
                        if (dk != 0) ++nCh;
                    }
                }
                scalar f1 = 0, f2 = 0, as = 0;
                label nChF = 0;
                forAll(own, f)
                {
                    if (layer[own[f]] != 0 && layer[nei[f]] != 0) continue;
                    const scalar e = mag(kf[f] - kExact(Cf[f]));
                    f1 += e*magSf[f]; f2 += e*e*magSf[f]; as += magSf[f];
                    if (kf0.size() && kf[f] != kf0[f]) ++nChF;
                }
                if (k == 0)
                {
                    k0 = kappa.primitiveField();
                    kf0 = kf.primitiveField();
                }

                csv << hRef << "," << scheme << "," << (k == 0 ? word("none") : type)
                    << "," << k << "," << e1/max(vs, SMALL) << ","
                    << Foam::sqrt(e2/max(vs, SMALL)) << "," << f1/max(as, SMALL) << ","
                    << Foam::sqrt(f2/max(as, SMALL)) << "," << d1/max(vs, SMALL) << ","
                    << nCh << "," << nChF << "," << nCell0 << "," << nFace0 << "\n";

                Info<< scheme << " " << (k == 0 ? word("none") : type) << " k " << k
                    << ": E_kappa cell L1 " << e1/max(vs, SMALL)
                    << ", face L1 " << f1/max(as, SMALL)
                    << ", changed cells " << nCh << "/" << nCell0
                    << ", faces " << nChF << "/" << nFace0 << endl;
            }
        }
    }

    Info<< "End" << endl;
    return 0;
}

// ************************************************************************* //
