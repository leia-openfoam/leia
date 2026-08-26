/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Tomislav Maric & Julian Reitzel, TU Darmstadt
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
    Foam::neighbourNarrowBand

Description
    Concrete class which calculates a field with {0,1} which mark cells near the zero level-set interface.

SourceFiles
    neighbourNarrowBand.C


\*---------------------------------------------------------------------------*/

#ifndef neighbourNarrowBand_H
#define neighbourNarrowBand_H

#include "signChangeNarrowBand.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class neighbourNarrowBand Declaration
\*---------------------------------------------------------------------------*/

class neighbourNarrowBand
    : public signChangeNarrowBand
{
    int ntimes_;
    
public:

    TypeName("neighbours");
    neighbourNarrowBand(const dictionary& dict, const volScalarField& psi);
    virtual ~neighbourNarrowBand() = default;

    //- The SEED that the dilation in calc() grows from. Factored out as a
    //  virtual so a derived class can change WHICH cells start the band
    //  without duplicating the dilation, which is subtle and
    //  parallel-critical: it exchanges once per layer because layer k+1 on
    //  this rank may only be reachable from a cell that layer k placed on the
    //  far side of a processor boundary. Duplicating that loop would let the
    //  two copies drift.
    //
    //  Default: the psi sign change, so neighbourNarrowBand itself is
    //  behaviour-identical to before this was factored out.
    virtual void seed() { signChangeNarrowBand::calc(); }

    virtual void calc();
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
