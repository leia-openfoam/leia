// GEOMETRIC FACE DENSITY from the interface just built (alphaf, rhof, the face
// viscosity), and the geometric cell density for the mass-flux modes that solve
// no density equation. Shared by the semi-Lagrangian (slAlphaEqn.H) and the
// Eulerian (alphaEqn.H) two-phase solvers; the fields and switches come from
// createMassFluxFields.H. Moved here VERBATIM from slAlphaEqn.H on 2026-09-27.
        // Geometrical face density at the NEW interface location. In rhoLENT
        // this remains fixed during the outer iterations while phi^o changes;
        // the same rhof*phi^o then appears in the auxiliary density equation
        // and in div(rhoPhi,U). The older comparison modes also reuse rhof.
        if (!freezeRhoPhi && massFluxType != "interpolatedDensity")
        {
            const volScalarField& nb =
                mesh.lookupObject<volScalarField>("NarrowBand");
            if (alphaFAdvected)
            {
                // The face velocity that traces the plane is the SAME field the SL
                // scheme traces the cell feet with -- reconstruct(phi) for the
                // projectedFlux setting -- interpolated to faces. Using phi's own
                // reconstruction keeps the face and cell displacements from being two
                // different velocity fields, which is the whole point.
                const surfaceVectorField Uf
                (
                    "UfSLtrace", fvc::interpolate(fvc::reconstruct(phi))
                );
                computeFaceAreaFractions
                (
                    mesh, psi.oldTime(), nb, phi, alphaf, massFluxUpwind,
                    &Uf, runTime.deltaTValue()
                );
            }
            else
            {
                computeFaceAreaFractions
                (
                    mesh, psi, nb, phi, alphaf, massFluxUpwind
                );
            }
            if (alphaFTrapezoid)
            {
                // TIME-AVERAGED face density (author's proposal, 2026-09-02). The
                // density that crosses a face during the step is the average over the
                // SWEPT region, not the plane value at t^{n+1}; the trapezoid of the
                // two endpoint reconstructions is its second-order approximation.
                //
                // psi.oldTime() is the interface at t^n and is maintained by the
                // semi-Lagrangian scheme (pointValueScheme uses it directly). The
                // donor direction is taken from the CURRENT phi for both levels: phi
                // changes by O(dt) across the step, so using it for the old level is
                // a second-order-consistent choice and keeps the two reconstructions
                // on the same upwind side, which a mixed choice would not.
                computeFaceAreaFractions
                (
                    mesh, psi.oldTime(), nb, phi, alphafOld, massFluxUpwind
                );
                const surfaceScalarField alphafBar
                (
                    "alphafBar", 0.5*(alphaf + alphafOld)
                );
                rhof = alphafBar*rho1 + (scalar(1) - alphafBar)*rho2;
            }
            else
            {
                rhof = alphaf*rho1 + (scalar(1) - alphaf)*rho2;
            }

            // Face dynamic viscosity: alpha_f source and mixing model are two
            // independent choices, both handled in one place.
            #include "faceViscosity.H"
        }

        // The non-rhoLENT modes use the geometrical cell density immediately.
        // rhoLENT deliberately keeps rho.oldTime() as the previous geometrical
        // state and computes the current auxiliary rho below.
        if (!freezeRhoPhi && !useRhoLENT)
        {
            rho == alpha1*rho1 + alpha2*rho2;
            rho.correctBoundaryConditions();
        }
