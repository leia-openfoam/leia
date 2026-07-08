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
    Foam::meshWaveExt

Description
    Velocity extension by OpenFOAM's native unstructured wave propagation:
    MeshWave<wallPointData<vector>> (the wallDist machinery carrying a vector
    payload) seeded at every psi-sign-change face with (origin = face foot
    point x_f - psi_f n_f, payload = interface velocity there). The wave
    delivers the payload of the NEAREST seed to every cell -- a discrete
    closest-point extension, ordered by distance like fast marching
    (Adalsteinsson & Sethian JCP 148 (1999)), parallel/cyclic/AMI-native via
    FaceCellWave. Accuracy is O(h) (nearest SEED, not the true foot point);
    the parallel-robust fallback of the extension family. Selected as
    `type meshWave;`.

SourceFiles
    meshWaveExt.C

\*---------------------------------------------------------------------------*/

#ifndef meshWaveExt_H
#define meshWaveExt_H

#include "interfaceExtension.H"

namespace Foam
{

class meshWaveExt
:
    public interfaceExtension
{
protected:

    virtual void extend() override;

public:

    TypeName("meshWave");

    meshWaveExt(const fvMesh& mesh);

    virtual ~meshWaveExt() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //
