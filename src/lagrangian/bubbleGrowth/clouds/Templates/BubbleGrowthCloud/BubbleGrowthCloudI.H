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

#include "fvmSup.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::BubbleGrowthCloud<CloudType>&
Foam::BubbleGrowthCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::BubbleGrowthCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::BubbleGrowthCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::BubbleGrowthCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::BubbleGrowthCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::BubbleGrowthCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::BubbleGrowthCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::BubbleGrowthCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::BubbleGrowthCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::BubbleGrowthCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::BubbleGrowthCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::BubbleGrowthCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::BubbleGrowthCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::BubbleGrowthCloud<CloudType>::p() const
{
    return p_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::BubbleGrowthCloud<CloudType>::g() const
{
    return g_;
}





template<class CloudType>
inline Foam::scalar Foam::BubbleGrowthCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::BubbleGrowthCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::BubbleGrowthCloud<CloudType>::forceType&
Foam::BubbleGrowthCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::BubbleGrowthCloud<CloudType>::forceType&
Foam::BubbleGrowthCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::BubbleGrowthCloud<CloudType>::functionType&
Foam::BubbleGrowthCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::stochasticCollision() const
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::stochasticCollision()
{
    return stochasticCollisionModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::BubbleGrowthCloud<CloudType>>&
Foam::BubbleGrowthCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::BubbleGrowthCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::label Foam::BubbleGrowthCloud<CloudType>::nParcels() const
{
    return this->size();
}


template<class CloudType>
inline Foam::scalar Foam::BubbleGrowthCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::BubbleGrowthCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::BubbleGrowthCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::BubbleGrowthCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::BubbleGrowthCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


template<class CloudType>
inline Foam::cachedRandom& Foam::BubbleGrowthCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*>>&
Foam::BubbleGrowthCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline const Foam::scalarField&
Foam::BubbleGrowthCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::BubbleGrowthCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::BubbleGrowthCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::BubbleGrowthCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::BubbleGrowthCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::BubbleGrowthCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            const DimensionedField<scalar, volMesh>
                Vdt(mesh_.V()*this->db().time().deltaT());

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm.ref();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::BubbleGrowthCloud<CloudType>::vDotSweep() const
{
    tmp<volScalarField> tvDotSweep
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":vDotSweep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& vDotSweep = tvDotSweep.ref();
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        vDotSweep[celli] += p.nParticle()*p.areaP()*mag(p.U() - U_[celli]);
    }

    vDotSweep.primitiveFieldRef() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::BubbleGrowthCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& theta = ttheta.ref();
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        theta[celli] += p.nParticle()*p.volume();
    }

    theta.primitiveFieldRef() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::BubbleGrowthCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha.ref().primitiveFieldRef();
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        alpha[celli] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::BubbleGrowthCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + ":rhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff.ref().primitiveFieldRef();
    forAllConstIter(typename BubbleGrowthCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label celli = p.cell();

        rhoEff[celli] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


// ************************************************************************* //
