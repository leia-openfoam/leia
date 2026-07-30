// Optional grid-scale psi filter (see createSLFields.H for the dict and the
// rationale). Applied once per time step, AFTER the SL advection and the
// narrow-band update, BEFORE redistancing/indicator/curvature -- everything
// downstream reads the filtered field.
//
// psi_c <- psi_c - theta * L(L(psi))_c,   L(psi) = psi - avg(interp(psi)),
//
// restricted to the narrow band + one face-neighbour ring. L(L(.)) is a
// compact biharmonic-type high-pass: zero on linear AND (up to a constant
// absorbed by the second L) quadratic fields on uniform meshes, so a smooth
// curved signed distance is untouched to leading order -- only the grid-scale
// (2h) profile modes, the measured carrier of the parasitic feedback, are
// damped. This is the explicit delta-psi dissipation term the reinit-free
// semi-Lagrangian transport otherwise lacks (MULES-type algebraic advection
// gets an equivalent term for free).
if (psiFilterOn)
{
    // Band + 1 ring mask (face-neighbour dilation of the registered band).
    const volScalarField& nbMask =
        mesh.lookupObject<volScalarField>("NarrowBand");

    boolList inFilterBand(mesh.nCells(), false);
    forAll(nbMask, celli)
    {
        if (nbMask[celli] == 1) { inFilterBand[celli] = true; }
    }
    {
        boolList dilated(inFilterBand);
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        forAll(own, facei)
        {
            if (inFilterBand[own[facei]]) { dilated[nei[facei]] = true; }
            if (inFilterBand[nei[facei]]) { dilated[own[facei]] = true; }
        }
        inFilterBand.transfer(dilated);
    }

    const volScalarField Lpsi
    (
        psi - fvc::average(fvc::interpolate(psi))
    );
    const volScalarField L2psi
    (
        Lpsi - fvc::average(fvc::interpolate(Lpsi))
    );

    scalarField& psiIn = psi.primitiveFieldRef();
    const scalarField& L2In = L2psi.primitiveField();
    forAll(psiIn, celli)
    {
        if (inFilterBand[celli])
        {
            psiIn[celli] -= psiFilterTheta*L2In[celli];
        }
    }
    psi.correctBoundaryConditions();
}
