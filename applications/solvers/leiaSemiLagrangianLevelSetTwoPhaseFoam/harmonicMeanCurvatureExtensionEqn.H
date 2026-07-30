// Optional HARMONIC (Laplace) smoothing / extension of the mean curvature.
//
// Motivation: the cell-centred symbolic curvature (meanCurvatureNoExtension) is
// clean on the initialised signed-distance field, but during the coupled run it
// develops a grid-anisotropy-seeded, cos(4*theta) oscillation in the interface
// band that feeds the parasitic-current instability. The earlier remedy -- the
// foot-point Newton normal extension -- amplified the spikes (~40-85x) via its
// projection and was itself destabilising.
//
// This is a solver-friendly alternative with NO per-cell Newton projection: the
// curvature at the INTERFACE cells (those the psi=0 contour passes through) is
// pinned as a Dirichlet condition (fvMatrix::setValues), and laplacian(kappa)==0
// is solved once. The result is a smooth HARMONIC field that reproduces the
// trusted interface-cell curvature and extends it smoothly into the band,
// replacing the oscillating off-interface values that the CSF force interpolates.
//
// Enabled by the fvSolution entry:
//   levelSet { curvatureExtension { type harmonicLaplace; } }   // "none" disables
// (read into the bool harmonicCurvatureExtension in createSLFields.H). Requires a
// linear-solver entry for "kappa" in fvSolution/solvers.
if (harmonicCurvatureExtension)
{
    // --- interface cells: an internal- or coupled-face neighbour has psi of the
    //     opposite sign, i.e. the psi=0 contour passes through the cell pair. This
    //     is the one-cell-thick interface layer, needs no length scale, and is the
    //     set where the cell-centred curvature is most trustworthy (~1/R).
    const scalarField& psiIn = psi.primitiveField();
    boolList isInterfaceCell(mesh.nCells(), false);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    forAll(own, facei)
    {
        if (psiIn[own[facei]]*psiIn[nei[facei]] < 0)
        {
            isInterfaceCell[own[facei]] = true;
            isInterfaceCell[nei[facei]] = true;
        }
    }
    // coupled (processor / cyclic) boundaries
    forAll(psi.boundaryField(), patchi)
    {
        const fvPatchScalarField& psip = psi.boundaryField()[patchi];
        if (psip.coupled())
        {
            const scalarField psiNbr(psip.patchNeighbourField());
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, i)
            {
                if (psiIn[fc[i]]*psiNbr[i] < 0)
                {
                    isInterfaceCell[fc[i]] = true;
                }
            }
        }
    }

    DynamicList<label> ifaceCellsDyn(mesh.nCells()/8 + 1);
    forAll(isInterfaceCell, celli)
    {
        if (isInterfaceCell[celli]) { ifaceCellsDyn.append(celli); }
    }
    const labelList interfaceCells(ifaceCellsDyn);

    // pinned Dirichlet values = the current cell-centred curvature there
    scalarField pinnedKappa(interfaceCells.size());
    forAll(interfaceCells, i) { pinnedKappa[i] = kappa[interfaceCells[i]]; }

    // --- harmonic extension: laplacian(kappa) == 0 with kappa fixed on the
    //     interface cells. Unit (constant) diffusivity -> a pure Laplace/harmonic
    //     problem; the constant value is irrelevant to the harmonic solution. The
    //     field's own (zeroGradient) boundary conditions close the domain.
    fvScalarMatrix kappaExtEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("Dkappa", dimless, scalar(1)),
            kappa
        )
    );
    kappaExtEqn.setValues(interfaceCells, pinnedKappa);
    kappaExtEqn.solve();

    Info<< "Harmonic curvature extension: pinned "
        << returnReduce(interfaceCells.size(), sumOp<label>())
        << " interface cells; kappa in ["
        << gMin(kappa.primitiveField()) << ", "
        << gMax(kappa.primitiveField()) << "] 1/m" << endl;
}
