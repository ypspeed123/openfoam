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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(dictionary::null),
    parcelTypeId_(dict_, -1),
    rhoMin_(dict_, 0.0),
    rho0_(dict_, 0.0),
    minParcelMass_(dict_, 0.0)
{}


template<class ParcelType>
inline Foam::BubbleGrowthParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParcelMass_(cp.minParcelMass_)
{}


template<class ParcelType>
inline Foam::BubbleGrowthParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(dict_, "parcelTypeId", -1),
    rhoMin_(dict_, "rhoMin", 1e-15),
    rho0_(dict_, "rho0"),
    minParcelMass_(dict_, "minParcelMass", 1e-15)
{}


template<class ParcelType>
inline Foam::BubbleGrowthParcel<ParcelType>::BubbleGrowthParcel
(

    const polyMesh& owner,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI
)
:
    ParcelType(owner, position, celli, tetFacei, tetPtI),
    oderp_(owner),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    rDot_(0.0),
    rDotPre_(0.0),
    delta_(0.0),
    Tv_(0.0),
    Tref_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    rhoc_(0.0),
    Uc_(Zero),
    muc_(0.0)
{}


template<class ParcelType>
inline Foam::BubbleGrowthParcel<ParcelType>::BubbleGrowthParcel
(

    const polyMesh& owner,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar rDot0,
    const scalar rDotPre0,
    const scalar delta0,
    const scalar Tv,
    const scalar Tref,
    const scalar dTarget0,
    const vector& U0,
    const constantProperties& constProps
)
:
    ParcelType(owner, position, celli, tetFacei, tetPtI),
    oderp_(owner),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    rDot_(rDot0),
    rDotPre_(rDotPre0),
    delta_(delta0),
    Tv_(Tv),
    Tref_(Tref),
    dTarget_(dTarget0),
    U_(U0),
    rho_(constProps.rho0()),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    rhoc_(0.0),
    Uc_(Zero),
    muc_(0.0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::BubbleGrowthParcel<ParcelType>::constantProperties::minParcelMass() const
{
    return minParcelMass_.value();
}


// * * * * * * * BubbleGrowthParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::BubbleGrowthParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::BubbleGrowthParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::rDot() const
{
    return rDot_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::rDotPre() const
{
    return rDotPre_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::delta() const
{
    return delta_;
}


template<class ParcelType>
inline void Foam::BubbleGrowthParcel<ParcelType>::update
(
	scalar d,
	scalar rDot,
	scalar rDotPre,
	scalar delta,
	scalar Tv,
	scalar Tref
)
{
	d_=d;
	rDotPre_=rDotPre;
	rDot_=rDot;
	delta_=delta;
	Tv_=Tv;
	Tref_=Tref;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::Tv() const
{
    return Tv_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::Tref() const
{
    return Tref_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::BubbleGrowthParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::age() const
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::BubbleGrowthParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::rhoc() const
{
    return rhoc_;
}


template<class ParcelType>
inline const Foam::vector& Foam::BubbleGrowthParcel<ParcelType>::Uc() const
{
    return Uc_;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::muc() const
{
    return muc_;
}


template<class ParcelType>
inline bool& Foam::BubbleGrowthParcel<ParcelType>::active()
{
    return active_;
}


template<class ParcelType>
inline Foam::label& Foam::BubbleGrowthParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::d()
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::rDot()
{
    return rDot_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::rDotPre()
{
    return rDotPre_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::delta()
{
    return delta_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::Tv()
{
    return Tv_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::Tref()
{
    return Tref_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::dTarget()
{
    return dTarget_;
}


template<class ParcelType>
inline Foam::vector& Foam::BubbleGrowthParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::rho()
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::age()
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar& Foam::BubbleGrowthParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template<class ParcelType>
inline Foam::vector& Foam::BubbleGrowthParcel<ParcelType>::UTurb()
{
    return UTurb_;
}


template<class ParcelType>
inline Foam::label Foam::BubbleGrowthParcel<ParcelType>::faceInterpolation() const
{
    // Use volume-based interpolation if dealing with external faces
    if (this->cloud().internalFace(this->face()))
    {
        return this->face();
    }
    else
    {
        return -1;
    }
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::massCell
(
    const label celli
) const
{
    return rhoc_*this->mesh().cellVolumes()[celli];
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::Re
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar muc
) const
{
    return rhoc*mag(U - Uc_)*d/(muc + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::We
(
    const vector& U,
    const scalar d,
    const scalar rhoc,
    const scalar sigma
) const
{
    return rhoc*magSqr(U - Uc_)*d/(sigma + ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::BubbleGrowthParcel<ParcelType>::Eo
(
    const vector& a,
    const scalar d,
    const scalar sigma
) const
{
    vector dir = U_/(mag(U_) + ROOTVSMALL);
    return mag(a & dir)*(rho_ - rhoc_)*sqr(d)/(sigma + ROOTVSMALL);
}


// ************************************************************************* //
