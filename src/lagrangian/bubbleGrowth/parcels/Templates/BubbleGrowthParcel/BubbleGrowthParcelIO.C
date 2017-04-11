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

#include "BubbleGrowthParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::BubbleGrowthParcel<ParcelType>::propertyList_ =
    Foam::BubbleGrowthParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::BubbleGrowthParcel<ParcelType>::sizeofFields_
(
    offsetof(BubbleGrowthParcel<ParcelType>, rhoc_)
  - offsetof(BubbleGrowthParcel<ParcelType>, active_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::BubbleGrowthParcel<ParcelType>::BubbleGrowthParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    oderp_(mesh),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
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
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            rDot_ = readScalar(is);
            rDotPre_ = readScalar(is);
            delta_ = readScalar(is);
            Tv_ = readScalar(is);
            Tref_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "BubbleGrowthParcel<ParcelType>::BubbleGrowthParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleGrowthParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<label> active(c.fieldIOobject("active", IOobject::MUST_READ));
    c.checkFieldIOobject(c, active);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<scalar>
        nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);
    
    IOField<scalar> rDot(c.fieldIOobject("rDot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rDot);
    
    IOField<scalar> rDotPre(c.fieldIOobject("rDotPre", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rDotPre);
    
    IOField<scalar> delta(c.fieldIOobject("delta", IOobject::MUST_READ));
    c.checkFieldIOobject(c, delta);

    IOField<scalar> Tv(c.fieldIOobject("Tv", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Tv);
    
    IOField<scalar> Tref(c.fieldIOobject("Tref", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Tref);

    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age(c.fieldIOobject("age", IOobject::MUST_READ));
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UTurb);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        BubbleGrowthParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.rDot_ = rDot[i];
        p.rDotPre_ = rDotPre[i];
        p.delta_ = delta[i];
        p.Tv_ = Tv[i];
        p.Tref_ = Tref[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::BubbleGrowthParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> rDot(c.fieldIOobject("rDot", IOobject::NO_READ), np);
    IOField<scalar> rDotPre(c.fieldIOobject("rDotPre", IOobject::NO_READ), np);
    IOField<scalar> delta(c.fieldIOobject("delta", IOobject::NO_READ), np);
    IOField<scalar> Tv(c.fieldIOobject("Tv", IOobject::NO_READ), np);
    IOField<scalar> Tref(c.fieldIOobject("Tref", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const BubbleGrowthParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        rDot[i] = p.rDot();
        rDotPre[i] = p.rDotPre();
        delta[i] = p.delta();
        Tv[i] = p.Tv_;
        Tref[i] = p.Tref_;
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();

        i++;
    }

    active.write();
    typeId.write();
    nParticle.write();
    d.write();
    rDot.write();
    rDotPre.write();
    delta.write();
    Tv.write();
    Tref.write();
    dTarget.write();
    U.write();
    rho.write();
    age.write();
    tTurb.write();
    UTurb.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BubbleGrowthParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.rDot()
            << token::SPACE << p.rDotPre()
            << token::SPACE << p.delta()
            << token::SPACE << p.Tv()
            << token::SPACE << p.Tref()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            BubbleGrowthParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const BubbleGrowthParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
