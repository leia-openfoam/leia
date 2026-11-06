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
    Foam::davofNormalModels::plicRDF

Description
    Reference: the helper is OpenFOAM's reconstruction::plicRDF of
    libgeometricVoF (H. Scheufler, J. Roenby, Accurate and efficient surface
    reconstruction from volume fraction data on general meshes, J. Comput.
    Phys. 383 (2019) 1-23, doi:10.1016/j.jcp.2019.01.009), OpenFOAM v2606's
    version since the port into leia (2026-11-06), not the TwoPhaseFlow fork's.

    The normal of a helper plicRDF reconstruction of the same alpha, owned
    by the model and reconstructed whenever the davof scheme is. The
    helper's interfaceNormal/interfaceCentre fields and its
    reconstructionScheme dictionary carry the names the davof scheme has
    registered: objectRegistry::checkIn refuses the duplicates and the
    registry keeps the davof scheme's objects. Default applyTo check.

    The helper marks only the cells with alpha strictly inside
    (surfCellTol, 1 - surfCellTol); the davof scheme marks its interface
    cells by their face fractions (|m_c| > wispTol V^(2/3)), so the cells at
    the edge of the band that the scheme marks and the helper does not have
    no model normal. What they get is the base class's noNormalFallback
    (keep: m_c, counted "kept m"; ringArea: the ring normal, counted "ring
    fallbacks").

SourceFiles
    plicRDFNormal.C

\*---------------------------------------------------------------------------*/

#ifndef davofNormalModels_plicRDF_H
#define davofNormalModels_plicRDF_H

#include "davofNormalModel.H"
#include "plicRDF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace davofNormalModels
{

/*---------------------------------------------------------------------------*\
                           Class plicRDF Declaration
\*---------------------------------------------------------------------------*/

class plicRDF
:
    public davofNormalModel
{
    // Private Data

        //- The helper reconstruction
        autoPtr<reconstruction::plicRDF> helperPtr_;


public:

    //- Runtime type information
    TypeName("plicRDF");


    //- Construct from the davof dictionary and the scheme
    plicRDF
    (
        const dictionary& dict,
        const reconstruction::davof& scheme
    );


    //- Destructor
    virtual ~plicRDF() = default;


    // Member Functions

        virtual void modelNormals(const scalar tState, vectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace davofNormalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
