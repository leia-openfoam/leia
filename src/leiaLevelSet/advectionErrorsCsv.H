/*---------------------------------------------------------------------------*\
    Single definition of the advection-error CSV schema.

Why this exists
---------------
Three solvers write the SAME error table -- leiaLevelSetFoam,
leiaRedistancedLevelSetFoam (which includes leiaLevelSetFoam's advectionErrors.H
via `-I../leiaLevelSetFoam`) and leiaSemiLagrangeLevelSetFoam -- and each one
carried its OWN copy of the header string. The comment in the semi-Lagrangian
copy read "same columns as leiaLevelSetFoam.csv for direct comparison", which is
a statement no compiler was checking.

It broke the moment the schema changed: adding the signed volume-error columns
to the shared advectionErrors.H made leiaRedistancedLevelSetFoam write 14 fields
under a 12-field header, silently shifting every column after E_VOL_ALPHA_REL.
A curated CSV built from that would have read the geometric error as the signed
volume error, and nothing would have complained.

So the header lives here, once. A solver that writes these rows includes this
file and cannot disagree with the writer about what the columns are.

Keep this in lockstep with the write block in advectionErrors.H.
\*---------------------------------------------------------------------------*/

#ifndef advectionErrorsCsv_H
#define advectionErrorsCsv_H

namespace Foam
{

//- Column header for the per-step advection-error CSV, newline included.
inline const char* advectionErrorsCsvHeader()
{
    return
        "TIME,"
        "DELTA_X,"
        "L_INF_E_PSI,"
        // Magnitude of the phase-volume error, and relative to the initial
        // volume. Kept for continuity: every published table reads these.
        "E_VOL_ALPHA,"
        "E_VOL_ALPHA_REL,"
        // ...and the SIGNED counterparts. A magnitude cannot separate a
        // systematic volume LOSS from an error that oscillates in sign, and
        // those have different causes. Measured on the 2D steady vortex, the
        // sourceless baseline is negative at all seven resolutions while
        // SDPLS:R changes sign twice.
        "E_VOL_ALPHA_SIGNED,"
        "E_VOL_ALPHA_SIGNED_REL,"
        "E_GEOM_ALPHA,"
        "E_GEOM_ALPHA_REL,"
        "E_BOUND_ALPHA,"
        "ELAPSED_CPU_TIME,"
        "ELAPSED_CLOCK_TIME,"
        "MAX_CFL,"
        "MAX_ABS_SUM_VOL_FLUXES\n";
}

} // End namespace Foam

#endif

// ************************************************************************* //
