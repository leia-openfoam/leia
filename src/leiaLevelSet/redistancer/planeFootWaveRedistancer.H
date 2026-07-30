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
    Foam::planeFootWaveRedistancer

Description
    Plane-anchored geometric redistancer: the narrow-band cells' NORMALIZED
    linear least-squares planes (planeAnchors.H -- the same planes the
    Detrixhe-Aslam indicator and the geometric mass flux read) seed a
    MeshWave/FaceCellWave nearest-donor propagation. The wave orders donors
    by distance to the plane FOOT point; every cell then evaluates the
    CONTINUOUS signed distance |(x - foot) & nHat| to its donor's plane
    (geometricVoF/plicRDF's trick; the payload carries the donor normal).
    Evaluating the plane instead of the discrete foot-point cloud removes
    the scalloping floor point-cloud distances put on the band deviation
    (measured on the static gate: a redistance event must LOWER the band
    L1 of ||grad psi|-1|, never raise it). Anchor cells (band + first
    ring, planeAnchors.H) keep their signed plane distances exactly; the
    wave region is re-signed with sign(psi_old). No linear solves and no
    parameters -- monotone, causal, parallel-native. Requires the
    registered "NarrowBand" field (levelSet.narrowBand.type signChange).

SourceFiles
    planeFootWaveRedistancer.C

\*---------------------------------------------------------------------------*/

#ifndef planeFootWaveRedistancer_H
#define planeFootWaveRedistancer_H

#include "redistancer.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class planeFootWaveRedistancer Declaration
\*---------------------------------------------------------------------------*/

class planeFootWaveRedistancer
:
    public redistancer
{
private:

    const fvMesh& mesh_;

public:

    TypeName("planeFootWave");

    planeFootWaveRedistancer(const fvMesh& mesh);

    virtual ~planeFootWaveRedistancer() = default;

    virtual bool doRedistance(volScalarField& psi);
};

} // End namespace Foam

#endif

// ************************************************************************* //
