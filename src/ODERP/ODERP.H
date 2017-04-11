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

Description

\*---------------------------------------------------------------------------*/
#ifndef ODERP_H
#define ODERP_H

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "IOdictionary.H"
#include "liquidProperties.H"
#include "dictionary.H"
#include <fstream>
namespace Foam
{

template<class Type>
class Cloud;

template<class ParcelType>
class BubbleGrowthParcel;

template<class CloudType>
class BubbleGrowthCloud;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class ParcelType>
class ODERP
:
    public ODESystem
{
//private members
	IOdictionary ODEProperties_;
	scalar R0_;
	scalar gamma_;
	scalar p0_;
	const word ODESolverType_;
	scalar ODESolverStep_;
    scalar ODESolverSubStep_;
    scalar ODESolverTol_;
	scalar pL_;
	scalar Tv_;
	scalar Tinf_;
	autoPtr<liquidProperties> fuel_;
	scalar writeControl_;
	
	// coefficients for vapor density
	scalar a_;
	scalar b_;
	scalar c_;
	scalar d_;
	
//public members	
public:
	friend class BubbleGrowthParcel<ParcelType>;
	friend class Cloud<BubbleGrowthParcel<ParcelType>>;
	friend class BubbleGrowthCloud<Cloud<BubbleGrowthParcel<ParcelType>>>;
	
	//constuctor	
	inline ODERP(const polyMesh&);
	inline ODERP(const ODERP<ParcelType>&);

	//destructor
    inline virtual ~ODERP();
    
//member functions

    label nEqns() const
    {
        return 2;
    }

    inline void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const;
    
    inline void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const;
   	
   	inline void operator=(const ODERP<ParcelType>&);
   	
   	inline void setEnvironment
   	(
   		const scalar& pL,
   		const scalar& Tv
   	);
   	inline const liquidProperties& fuel() const;
   	inline const word& ODESolverType() const;
   	inline scalar ODESolverStep() const;
   	inline scalar ODESolverSubStep() const;
   	inline scalar ODESolverTol() const;
   	inline scalar Tinf() const;
   	inline scalar rhov(scalar);
   	inline scalar gradT(scalar,scalar,scalar);
   	inline scalar delta
   	(
		scalar P,
		scalar Tv,
		scalar r,
		scalar rDot,
		scalar rDotPre,
		scalar delta,
		scalar Tinf
	);
   	inline void getTv
   	(
	scalar delta,
	scalar& Tv,
	scalar& Tref,
	scalar P,
	scalar r,
	scalar rDot
	);
	inline void output(scalar t,scalar R,scalar Rdot,scalar Tv,scalar delta,scalar pL,scalar Id)
	{
		std::ofstream outfile;
		outfile.open("bubble.txt",std::ios_base::out | std::ios_base::app);		
		outfile<<Id<<"\t"<<t<<"\t"<<R<<"\t"<<Rdot<<"\t"<<Tv<<"\t"<<delta<<"\t"<<pL<<"\n";
		outfile.close();
	}
	inline bool writeControl(scalar t)
	{
		if(writeControl_==0)
		writeControl_=ODESolverStep_;
		if (t>writeControl_)
		{
			writeControl_=writeControl_*pow(10,0.01);
			return true;
		}
		else
		return false;
	}
	inline scalar Tsat(scalar P)
	{
		if (P>10)
		{
		scalar a,b,c,d,e;
		a=73.649;
		b=-7258.2;
		c=-7.3037;
		d=4.1653e-6;
		e=2.0;
		scalar x=log10(P);
		scalar T=0.5197969075*pow(x,4)-7.4169394106*pow(x,3)+45.9989060649*pow(x,2)-99.612229501*x+323.1018292924;
		scalar pSat=fuel_->pv(P,T);
		while (abs(P-pSat)>1.0)
		{
			scalar temp=(P-pSat)/(pSat*(-2*b/T/T+c/T+d*e*pow(T,e-1)));
			T+=temp;
			pSat=fuel_->pv(P,T);
		}
		return max(T, 1.0);
		}
		else
		return 1.0;
	}
	inline scalar deri(scalar t,scalar A,scalar B)
	{
		return A*sqrt(A*A/B/B*t+1)-A*A/B*sqrt(t);
	}
};

}

#include "ODERPI.H"

#endif
// ************************************************************************* //
