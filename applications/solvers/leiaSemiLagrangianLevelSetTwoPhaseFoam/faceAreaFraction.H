/*---------------------------------------------------------------------------*\
    Geometric FACE liquid-area fractions from the level-set reconstruction, for a
    consistent (rhoLENT-style) two-phase mass flux rhoPhi = rho_f * phi.

    For each internal face the fraction of the face area in Omega^- = {psi < 0}
    is computed TWICE -- once from the owner cell's linear least-squares plane and
    once from the neighbour cell's -- and the two are averaged:

        alpha_f = 1/2 ( a_owner + a_neighbour ),     a = area(psi<0)/area(face).

    The per-cell linear plane psi_c(x) = n_c . x + d_c is the SAME reconstruction
    used by detrixheAslamPhaseIndicator (leastSquaresPlaneCoeffs), so alpha_f is
    consistent with the cell Heaviside H. The face area fraction under a linear psi
    is evaluated by triangulating the face (centroid + edges) and filling each
    triangle analytically (triNegativeFraction, the 2D analogue of Detrixhe &
    Aslam's tetNegativeFraction). The linear density model then gives the face
    density rho_f = alpha_f rho1 + (1 - alpha_f) rho2.
\*---------------------------------------------------------------------------*/

#ifndef faceAreaFraction_H
#define faceAreaFraction_H

#include "fvCFD.H"
#include "levelSetPlaneReconstruction.H"
#include "levelSetImplicitSurfaces.H"
#include <algorithm>

namespace Foam
{

//- Fraction of a triangle in Omega^- = {phi < 0} from its three vertex phi
//  values (2D analogue of detrixheAslamPhaseIndicator::tetNegativeFraction).
inline scalar triNegativeFraction(scalar a, scalar b, scalar c)
{
    scalar p[3] = {a, b, c};
    std::sort(p, p + 3);
    const scalar p1 = p[0], p2 = p[1], p3 = p[2];

    if (p3 <= 0) return 1;                          // all three inside
    if (p1 >  0) return 0;                          // all three outside
    if (p2 <= 0)                                    // two inside, one out (corner p3)
    {
        return 1 - (p3*p3)/stabilise((p3 - p1)*(p3 - p2), VSMALL);
    }
    // one inside (p1), two outside
    return (p1*p1)/stabilise((p2 - p1)*(p3 - p1), VSMALL);
}


//- Liquid (Omega^-) area fraction of face f under the linear field
//  phi(x) = n & x + d, triangulating the face about its centroid xf.
inline scalar faceNegativeFraction
(
    const face& f,
    const pointField& points,
    const point& xf,
    const vector& n,
    const scalar d
)
{
    const scalar phif = (n & xf) + d;
    scalar liqA = 0, totA = 0;
    forAll(f, ip)
    {
        const point& p0 = points[f[ip]];
        const point& p1 = points[f.nextLabel(ip)];
        const scalar triA = 0.5*mag((p0 - xf) ^ (p1 - xf));
        if (triA <= VSMALL) continue;
        liqA += triNegativeFraction(phif, (n & p0) + d, (n & p1) + d)*triA;
        totA += triA;
    }
    return (totA > VSMALL) ? liqA/totA : (phif < 0 ? scalar(1) : scalar(0));
}


//- Face liquid-area-fraction field alpha_f from the level-set reconstruction.
//  upwind=true  : the UPWIND cell's plane intersected with the face (the phase
//                 composition of the fluid the flux phi carries through the face
//                 -- the consistent, stable donor choice; upwind = owner if
//                 phi>=0 else neighbour).
//  upwind=false : 1/2 (owner + neighbour) central average (comparison only).
//  Cells outside the narrow band use their sign (fully one phase). Boundary faces
//  use the owner-side fraction (exact for walls/empty away from the interface; a
//  serial-exact first cut on processor patches).
//- ADVECTED variant (alphaFSource donorPlaneAdvected): pass the OLD level set and a
//  face velocity, and the donor plane is displaced by the SAME semi-Lagrangian step
//  the CELL values undergo before it is cut against the face.
//
//  The cell update is psi_c^{n+1} = psi^n(x_c - u dt). Its face analogue is
//
//      alpha_f^{n+1} = area fraction of { x in f : psi^n(x - u_f dt) < 0 },
//
//  and because the reconstruction is a PLANE this needs no geometry at all:
//
//      psi^n(x - u_f dt) = n.x + d - (n . u_f) dt,
//
//  i.e. a constant shift of the plane offset per face, one dot product. No swept
//  volumes, no tet clipping -- that would be unsplit geometric VOF, a different and
//  far more expensive method.
//
//  WHY IT SHOULD HELP. Today the cell value comes from psi^n through the quadratic
//  fit at the traced-back point, while alpha_f comes from a LINEAR plane fit of the
//  ALREADY-UPDATED psi^{n+1} cut against the face: two different fields, two different
//  reconstructions, no shared displacement. Advecting the plane makes both "psi^n
//  sampled at the departure point", so they can then disagree only through u_f versus
//  u at the cell centre, instead of through the whole construction.
inline void computeFaceAreaFractions
(
    const fvMesh& mesh,
    const volScalarField& psi,
    const volScalarField& narrowBand,
    const surfaceScalarField& phi,
    surfaceScalarField& alphaf,
    const bool upwind = true,
    const surfaceVectorField* UfPtr = nullptr,   // non-null => advect the plane
    const scalar dt = 0
)
{
    const bool advectPlane = (UfPtr != nullptr && dt > 0);
    const pointField& points = mesh.points();
    const faceList&   faces  = mesh.faces();
    const vectorField& Cf    = mesh.faceCentres();
    const vectorField& C     = mesh.cellCentres();

    const dictionary& levelSetDict =
        static_cast<const fvSolution&>(mesh).subDict("levelSet");
    const dictionary& phaseIndicatorDict =
        levelSetDict.subDict("phaseIndicator");
    const word geometrySource = phaseIndicatorDict.getOrDefault<word>
    (
        "geometrySource",
        "levelSetField"
    );
    const bool analyticGeometry =
        geometrySource == "analyticImplicitSurface";
    if (geometrySource != "levelSetField" && !analyticGeometry)
    {
        FatalErrorInFunction
            << "Unknown face-density geometry source '" << geometrySource
            << "'. Valid: levelSetField, analyticImplicitSurface."
            << abort(FatalError);
    }
    autoPtr<implicitSurface> analyticSurface;
    if (analyticGeometry)
    {
        const dictionary& surfaceDict = levelSetDict.subDict("implicitSurface");
        analyticSurface = implicitSurface::New
        (
            surfaceDict.get<word>("type"),
            surfaceDict
        );
    }

    // Coupled-patch neighbour data (one matched MPI call), so the per-cell
    // reconstruction below does no communication.
    const coupledFaceNeighbours coupledNei(mesh, psi);

    // Per-cell linear plane n_c . x + d_c; planeSet=false -> use the cell sign.
    vectorField nCell(mesh.nCells(), Zero);
    scalarField dCell(mesh.nCells(), Zero);
    boolList    planeSet(mesh.nCells(), false);
    forAll(narrowBand, cellI)
    {
        if (narrowBand[cellI] != 1) continue;
        vector nc(Zero);
        scalar d = 0;
        if (analyticGeometry)
        {
            nc = analyticSurface->grad(C[cellI]);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                nc /= nmag;
                d = analyticSurface->value(C[cellI]) - (nc & C[cellI]);
            }
        }
        else
        {
            const scalarList pc =
                leastSquaresPlaneCoeffs(mesh, psi, cellI, coupledNei);
            nc = vector(pc[0], pc[1], pc[2]);
            const scalar nmag = mag(nc);
            if (nmag >= SMALL)
            {
                d = pc[3]/nmag;
            }
        }
        const scalar nmag = mag(nc);
        if (nmag < SMALL) continue;
        nCell[cellI] = nc/nmag;
        dCell[cellI] = d;
        planeSet[cellI] = true;
    }

    const scalarField& psiIn = psi.primitiveField();
    // dShift displaces the donor plane by the semi-Lagrangian step at THIS face:
    // psi^n(x - u_f dt) = n.x + d - (n . u_f) dt. Zero unless advecting, so the
    // default path is bit-identical.
    auto cellFrac = [&](const label c, const label faceI, const scalar dShift) -> scalar
    {
        return planeSet[c]
          ? faceNegativeFraction
            (
                faces[faceI], points, Cf[faceI], nCell[c], dCell[c] + dShift
            )
          :
            (
                (analyticGeometry ? analyticSurface->value(C[c]) : psiIn[c]) < 0
              ? scalar(1) : scalar(0)
            );
    };

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& phiIn = phi.primitiveField();
    scalarField& af = alphaf.primitiveFieldRef();
    forAll(own, faceI)
    {
        if (upwind)
        {
            // Donor cell: phi points owner->neighbour, so phi>=0 => owner upwind.
            const label up = (phiIn[faceI] >= 0) ? own[faceI] : nei[faceI];
            const scalar dShift = advectPlane && planeSet[up]
              ? -(nCell[up] & (*UfPtr)[faceI])*dt : scalar(0);
            af[faceI] = cellFrac(up, faceI, dShift);
        }
        else
        {
            const scalar dO = advectPlane && planeSet[own[faceI]]
              ? -(nCell[own[faceI]] & (*UfPtr)[faceI])*dt : scalar(0);
            const scalar dN = advectPlane && planeSet[nei[faceI]]
              ? -(nCell[nei[faceI]] & (*UfPtr)[faceI])*dt : scalar(0);
            af[faceI] =
                0.5*(cellFrac(own[faceI], faceI, dO) + cellFrac(nei[faceI], faceI, dN));
        }
    }

    surfaceScalarField::Boundary& afbf = alphaf.boundaryFieldRef();
    forAll(afbf, patchI)
    {
        const fvPatch& fvp = mesh.boundary()[patchI];
        const labelUList& fcs = fvp.faceCells();
        const label start = fvp.start();
        fvsPatchScalarField& afp = afbf[patchI];
        forAll(afp, i)
        {
            const label c = fcs[i];
            const scalar dShift = advectPlane && planeSet[c]
              ? -(nCell[c] & (*UfPtr).boundaryField()[patchI][i])*dt : scalar(0);
            afp[i] = cellFrac(c, start + i, dShift);
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
