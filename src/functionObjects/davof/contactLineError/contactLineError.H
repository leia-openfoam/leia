/*---------------------------------------------------------------------------*\
    leia :: functionObjects :: contactLineError
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, GPL-3.0-or-later.

Class
    Foam::functionObjects::contactLineError

Description
    The contact points and the contact angles of the 2D contact line
    advection tests (the planar wall z = 0 of a one-cell-thick mesh, the
    interface a circular cap), measured on the reconstruction of the step in
    the way of
        M. Fricke, T. Maric, D. Bothe,
        Contact line advection using the geometrical Volume-of-Fluid method,
        Journal of Computational Physics 407 (2020) 109221,
        doi:10.1016/j.jcp.2019.109221, arXiv:1907.01785.
    For every cell that owns a face of the wall patch and is an interface cell
    of the reconstruction (davofAdvectionScheme::surf()), the PLIC plane
    (normal n, centre c) is intersected with the wall; the contact point
    x_cp = c_x + (n_z c_z - n_y (y_f - c_y))/n_x (y_f the wall face centre)
    is counted when it lies in the wall face's x-range. The angle through
    the liquid follows cos(theta) = n_Sigma,z with n_Sigma = -n/|n| the unit
    normal of the interface pointing out of the liquid (the reconstruction's
    normal points into the liquid), theta in (0, pi). The left contact point
    has n_Sigma,x < 0, the right one n_Sigma,x > 0; when several cells of a
    side qualify the outermost one (the leftmost for the left point, the
    rightmost for the right point) is taken and the number of candidates is
    reported. The reconstruction of a step is that of alpha^n at the start of
    the advection step, so the row of the call at t^{n+1} is labelled with
    t^n = t^{n+1} - dt. Per row: the liquid volume, the L1 distance
    sum |alpha - alpha(0)| V of the current alpha from the initial one
    (Eq. 35 of the paper at the end of a period), the extremes of alpha.
    The reference values are not computed here (the post-processing script
    scores every row).

SourceFiles
    contactLineError.C
\*---------------------------------------------------------------------------*/

#ifndef functionObjects_contactLineError_H
#define functionObjects_contactLineError_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

namespace Foam
{
namespace functionObjects
{

class contactLineError
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of the volume fraction field
        word fieldName_;

        //- Name of the wall patch
        word wallPatch_;

        //- The initial alpha (internal field)
        scalarField alpha0_;

        //- Whether alpha0_ is stored
        bool haveAlpha0_;

        //- Steps without a detected contact cell (left, right)
        label nMissing_[2];

        //- Print-only key dumpCandidates: every candidate contact cell of
        //  every step, the chosen one marked, with alpha of the step's
        //  reconstruction, |m|/V^(2/3), the plane, the wall face fraction
        //  alphaf (-1 where the field does not exist) and the wall cells
        //  around the chosen one
        bool dump_;

        //- alpha^n (the field the step's reconstruction used)
        scalarField alphaPrev_;

        //- alphaf of the wall faces before the step (-1: no such field)
        scalarField afPrevWall_;

    virtual void writeFileHeader(Ostream& os) const;

    //- No copy construct
    contactLineError(const contactLineError&) = delete;

    //- No copy assignment
    void operator=(const contactLineError&) = delete;

public:

    //- Runtime type information
    TypeName("contactLineError");

    //- Construct from Time and dictionary
    contactLineError
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~contactLineError() = default;

    virtual bool read(const dictionary& dict);

    virtual bool execute();

    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif
