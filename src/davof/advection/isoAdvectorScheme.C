/*---------------------------------------------------------------------------*\
    leia :: davof :: isoAdvectorScheme (see isoAdvectorScheme.H)
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tomislav Maric, TU Darmstadt
-------------------------------------------------------------------------------
License
    This file is part of leia, which is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    leia is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with leia.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "isoAdvectorScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoAdvectorScheme, 0);
    addToRunTimeSelectionTable
    (
        davofAdvectionScheme,
        isoAdvectorScheme,
        components
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoAdvectorScheme::isoAdvectorScheme
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    davofAdvectionScheme(typeName, alpha1, phi, U),
    iso_(alpha1, phi, U)
{
    Info<< "isoAdvector: OpenFOAM's isoAdvection (libgeometricVoF) with the"
        << " reconstruction scheme " << iso_.surf().type() << endl;
}


// ************************************************************************* //
