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

Class
    Foam::functionObjects::interfaceMarkers

Description
    The marker read-out of step B of the PDE level-set plan (section 6; the
    marker test of D. Bothe, T. Maric and K. Soga, manuscript of 2026-11-15,
    section 5). Verification instrument, not a method model.

    Markers Y_l are placed on the zero set of the initial surface
    (levelSet.implicitSurface; a circle in 2D, a Fibonacci sphere in 3D, from
    its `center` and `radius`) and moved with the PRESCRIBED velocity model
    (velocityModel, with its cos(pi t/tau) factor when it oscillates) by RK4,
    nSubSteps per solver step. On the interface the Newton point is the point
    itself, so the characteristics of the velocity-extension model are the
    material markers there, and WITHOUT feedback q = |grad psi| is constant
    along them for any flow:

        q_l(t) = q_l(0) = |grad psi_0(Y_l(0))|   (analytic, from the surface).

    At every execution the field is sampled at the markers by cell-point
    interpolation (unstructured, compact) of psi and of the cell gradient of
    psi (gradScheme, default the fvSchemes entry gradPsiMetric):

        E_Q    = q_h(Y_l(t)) - q_l(t)     (the target of step B)
        E_POS  = psi_h(Y_l(t))/q_l(t)     (the interface position error, to
                                           first order in the distance)

    reported as L1 and L2 means over the markers, never L_inf. One row per
    execution in interfaceMarkers.csv (master), the row at the start time
    included (the interpolation error of the initial field). Parallel: the
    markers are moved identically on every rank; the rank whose cell contains
    a marker samples it; the sums are reduced; N_FOUND counts the samples.

    Dictionary (system/controlDict functions):
        interfaceMarkers
        {
            type        interfaceMarkers;
            libs        (leiaCore);
            field       psi;
            nMarkers    720;             // 2D: on the circle; 3D: Fibonacci
            nSubSteps   4;               // RK4 sub-steps per solver step
            gradScheme  gradPsiMetric;   // fvSchemes gradSchemes key
        }

SourceFiles
    interfaceMarkers.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceMarkers_H
#define interfaceMarkers_H

#include "fvMeshFunctionObject.H"
#include "velocityModel.H"
#include "levelSetImplicitSurfaces.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class interfaceMarkers
:
    public fvMeshFunctionObject
{
    word fieldName_;
    label nMarkers_;
    label nSubSteps_;
    word gradScheme_;

    autoPtr<velocityModel> velModel_;
    autoPtr<implicitSurface> surface_;

    //- Marker positions and their exact q (no feedback: constant)
    pointField markers_;
    scalarField q0_;

    //- The time the markers belong to
    scalar tMarkers_;

    autoPtr<OFstream> csv_;

    vector velocityAt(const point& x, const scalar t) const;
    void advance(const scalar t1);
    void sampleAndWrite();

public:

    TypeName("interfaceMarkers");

    interfaceMarkers
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~interfaceMarkers() = default;

    virtual bool read(const dictionary& dict);
    virtual bool execute();
    virtual bool write();
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
