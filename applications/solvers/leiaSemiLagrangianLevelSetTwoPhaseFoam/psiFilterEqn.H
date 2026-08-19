#include "syncTools.H"

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
        // COUPLED PATCHES. The loop above walks INTERNAL faces only, so a cell
        // whose only band neighbour lives on another rank was never dilated into
        // and the filtered CELL SET itself depended on the decomposition.
        // Exchange the one-sided "my patch-internal cell is in the band" flag and
        // OR it in, so the band is the same set at any rank count.
        {
            scalarField nbrInBand(mesh.nBoundaryFaces(), 0.0);
            forAll(mesh.boundary(), patchI)
            {
                const fvPatch& patch = mesh.boundary()[patchI];
                if (!patch.coupled()) { continue; }
                const labelUList& fc = patch.faceCells();
                const label bStart = patch.start() - mesh.nInternalFaces();
                forAll(fc, i)
                {
                    nbrInBand[bStart + i] = inFilterBand[fc[i]] ? 1.0 : 0.0;
                }
            }
            syncTools::swapBoundaryFaceList(mesh, nbrInBand);
            forAll(mesh.boundary(), patchI)
            {
                const fvPatch& patch = mesh.boundary()[patchI];
                if (!patch.coupled()) { continue; }
                const labelUList& fc = patch.faceCells();
                const label bStart = patch.start() - mesh.nInternalFaces();
                forAll(fc, i)
                {
                    if (nbrInBand[bStart + i] > 0.5) { dilated[fc[i]] = true; }
                }
            }
        }
        inFilterBand.transfer(dilated);
    }

    // L(psi) = psi - average(interpolate(psi)), applied TWICE.
    //
    // THE HALO IS LOAD-BEARING AND WAS WRONG. Writing
    //     const volScalarField Lpsi(psi - fvc::average(fvc::interpolate(psi)));
    // makes Lpsi inherit its patch types from the SECOND operand, and
    // fvc::average is constructed with 'calculated' patch fields whose boundary
    // values are then overwritten with FACE values. So on a processor patch Lpsi
    // was UNCOUPLED, and the second application of L interpolated that stored
    // face value instead of the neighbour cell's Lpsi -- injecting a
    // seam-localised psi perturbation every step, in the very mechanism that
    // provides the damping.
    //
    // MEASURED, 3D N_L=60 on L=6R, 20 steps, serial against np=4:
    //                     filter OFF        filter ON (before this fix)
    //   max|U|            1.55e-06          2.43e-04     157x worse
    //   A2hL2Band         4.09e-09          8.38e-07     205x worse
    //   volumeRelError    3.64e-07          1.92e-05      53x worse
    // The rest of the solver is np-invariant at the solver floor; only the filter
    // was not. Constructing Lpsi from psi gives it psi's patch types (processor
    // on the seams) and correctBoundaryConditions() then performs a real
    // exchange, so L(L(psi)) sees the neighbour's L(psi) as it does internally.
    volScalarField Lpsi("Lpsi", psi);
    Lpsi = psi - fvc::average(fvc::interpolate(psi));
    Lpsi.correctBoundaryConditions();

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
