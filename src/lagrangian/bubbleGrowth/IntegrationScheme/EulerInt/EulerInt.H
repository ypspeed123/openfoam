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

Class
    Foam::EulerInt

Description
    EulerInt-implicit integration

\*---------------------------------------------------------------------------*/

#ifndef EulerInt_H
#define EulerInt_H

#include "IntegrationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class EulerInt Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerInt
:
    public IntegrationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("EulerInt");


    // Constructors

        //- Construct from components
        EulerInt(const word& phiName, const dictionary& dict);

        //- Copy constructor
        EulerInt(const EulerInt& is);

        //- Construct and return clone
        virtual autoPtr<IntegrationScheme<Type>> clone() const
        {
            return autoPtr<IntegrationScheme<Type>>(new EulerInt<Type>(*this));
        }


    //- Destructor
    virtual ~EulerInt();


    // Member Functions

        //- Perform the integration
        virtual typename IntegrationScheme<Type>::integrationResult integrate
        (
            const Type& phi,
            const scalar dt,
            const Type& alphaBeta,
            const scalar beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EulerInt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
