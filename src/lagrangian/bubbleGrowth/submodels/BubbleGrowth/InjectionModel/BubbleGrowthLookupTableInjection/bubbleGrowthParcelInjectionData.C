/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "bubbleGrowthParcelInjectionData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bubbleGrowthParcelInjectionData, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleGrowthParcelInjectionData::bubbleGrowthParcelInjectionData()
:
    x_(point::zero),
    U_(Zero),
    d_(0.0),
    rho_(0.0),
    mDot_(0.0)
{}


Foam::bubbleGrowthParcelInjectionData::bubbleGrowthParcelInjectionData
(
    const dictionary& dict
)
:
    x_(dict.lookup("x")),
    U_(dict.lookup("U")),
    d_(readScalar(dict.lookup("d"))),
    rho_(readScalar(dict.lookup("rho"))),
    mDot_(readScalar(dict.lookup("mDot")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubbleGrowthParcelInjectionData::~bubbleGrowthParcelInjectionData()
{}


// ************************************************************************* //
