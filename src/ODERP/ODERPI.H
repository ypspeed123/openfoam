
template<class ParcelType>
inline Foam::ODERP<ParcelType>::ODERP
(
	const polyMesh& mesh
)
:
   	ODESystem(),
   	ODEProperties_
   	(
   		IOobject
   		(
   			"ODEProperties",
   			mesh.time().constant(),
   			mesh,
   			IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
   		)
   	),
   	R0_(dimensionedScalar(ODEProperties_.lookup("R0")).value()),
    gamma_(dimensionedScalar(ODEProperties_.lookup("gamma")).value()),
    p0_(dimensionedScalar(ODEProperties_.lookup("p0")).value()),
    ODESolverType_(ODEProperties_.lookup("ODESolverType")),
    ODESolverStep_(dimensionedScalar(ODEProperties_.lookup("ODESolverStep")).value()),
    ODESolverSubStep_(dimensionedScalar(ODEProperties_.lookup("ODESolverSubStep")).value()),
    ODESolverTol_(dimensionedScalar(ODEProperties_.lookup("ODESolverTol")).value()),
    pL_(0),
    Tv_(0),
    Tinf_(dimensionedScalar(ODEProperties_.lookup("Tinf")).value()),
    writeControl_(0),
    a_(dimensionedScalar(ODEProperties_.lookup("a")).value()),
    b_(dimensionedScalar(ODEProperties_.lookup("b")).value()),
    c_(dimensionedScalar(ODEProperties_.lookup("c")).value()),
    d_(dimensionedScalar(ODEProperties_.lookup("d")).value())
{
	const word liquid(ODEProperties_.lookup("fuelType"));
	fileName fuel(liquid);
	dictionary fuelType(fuel);
	fuelType.add("defaultCoeffs", true);
	fuel_=liquidProperties::New(fuelType);
}

template<class ParcelType>
inline Foam::ODERP<ParcelType>::ODERP(const ODERP<ParcelType>& bt)
:
    ODESystem(),
    ODEProperties_(bt.ODEProperties_),
    R0_(bt.R0_),
    gamma_(bt.gamma_),
    p0_(bt.p0_),
    ODESolverType_(bt.ODESolverType_),
    ODESolverStep_(bt.ODESolverStep_),
    ODESolverSubStep_(bt.ODESolverSubStep_),
    ODESolverTol_(bt.ODESolverTol_),
    pL_(bt.pL_),
	Tv_(bt.Tv_),
	Tinf_(bt.Tinf_),
    fuel_(bt.fuel_),
    writeControl_(bt.writeControl_),
    a_(bt.a_),
    b_(bt.b_),
    c_(bt.c_),
    d_(bt.d_)    
{}

template<class ParcelType>
inline Foam::ODERP<ParcelType>::~ODERP()
{}

//member functions
template<class ParcelType>
inline void Foam::ODERP<ParcelType>::setEnvironment
(
	const scalar& pL,
	const scalar& Tv
)
{
	Tv_=Tv;
	pL_=pL;
}

template<class ParcelType>
inline void Foam::ODERP<ParcelType>::derivatives
(
	const scalar x,
	const scalarField& y,
	scalarField& dydx
) const
{
	
	scalar kappa=gamma_;
	scalar sigma=fuel_->sigma(pL_,Tv_);
	scalar muL=fuel_->mu(pL_,Tv_);
	scalar pv=fuel_->pv(pL_,Tv_);
	scalar rhoL=fuel_->rho(pL_,Tv_);
/*	Info<<"sigma: "<<sigma<<endl;
	Info<<"muL: "<<muL<<endl;
	Info<<"pv: "<<pv<<endl;
	Info<<"rhoL: "<<rhoL<<endl;
	Info<<"R0: "<<R0_<<endl;
	Info<<"R: "<<y[0]<<endl;
	Info<<"Rdot: "<<y[1]<<endl;
	Info<<"pL: "<<pL_<<endl;
*/
	if(y[0]>=R0_)
	kappa=1.0;
//	scalar pB=(p0_+2.0*sigma/R0_-pv)*Foam::pow(R0_/y[0],3.0*kappa)+pv;
	scalar pB=pv;
//Rayleigh-Plesset Equation

	dydx[0] = y[1] ;
    dydx[1] = (-3.0/2.0*Foam::sqr(y[1])+(pB-pL_)/rhoL)/y[0]-(4.0*muL*y[1]+2.0*sigma)/(rhoL*Foam::sqr(y[0]));
//    Info<<"pB: "<<pB<<" R: "<<dydx[0]<<" "<<dydx[1]<<endl;
}

template<class ParcelType>
inline void Foam::ODERP<ParcelType>::jacobian
(
	const scalar x,
	const scalarField& y,
	scalarField& dfdx,
	scalarSquareMatrix& dfdy
) const
{

	scalar kappa=gamma_;
	scalar sigma=fuel_->sigma(pL_,Tv_);
	scalar muL=fuel_->mu(pL_,Tv_);
	scalar pv=fuel_->pv(pL_,Tv_);
	scalar rhoL=fuel_->rho(pL_,Tv_);
	if(y[0]>=R0_)
	kappa=1.0;
	scalar pB=(p0_+2.0*sigma/R0_-pv)*Foam::pow(R0_/y[0],3.0*kappa)+pv;
	dfdx[0] = 0.0;
    dfdx[1] = 0.0;

    dfdy[0][0] = 0.0;
    dfdy[0][1] = 1.0;


    dfdy[1][0] =  -(-3.0/2.0*sqr(y[1])-(pB-pL_)/rhoL)/sqr(y[0])+(8.0*muL*y[1]+4.0*sigma)/(rhoL*pow(y[0],3));
    dfdy[1][1] =  -3.0*y[1] / y[0] - 4.0*muL / (rhoL*sqr(y[0])) ;
}

template<class ParcelType>
inline const Foam::liquidProperties& Foam::ODERP<ParcelType>::fuel() const
{
	return fuel_();
}


template<class ParcelType>
inline const Foam::word& Foam::ODERP<ParcelType>::ODESolverType() const
{
	return ODESolverType_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::ODESolverStep() const
{
	return ODESolverStep_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::ODESolverSubStep() const
{
	return ODESolverSubStep_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::ODESolverTol() const
{
	return ODESolverTol_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::Tinf() const
{
	return Tinf_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::rhov(scalar T)
{
	return ((a_*T+b_)*T+c_)*T+d_;
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::gradT(scalar P,scalar T,scalar Rdot)
{
	return fuel_->hl(P,T)*rhov(T)*Rdot/fuel_->K(P,T);
}


template<class ParcelType>
inline Foam::scalar Foam::ODERP<ParcelType>::delta
(
	scalar P,
	scalar Tv,
	scalar r,
	scalar rDot,
	scalar rDotPre,
	scalar delta,
	scalar Tinf
)
{
	scalar alpha=fuel_->K(P,Tv)/fuel_->rho(P,Tv)/fuel_->Cp(P,Tv);
	scalar dTdr=gradT(P,Tv,rDot);
//	Info<<"alpha: "<<alpha<<" dTdr: "<<dTdr<<endl;
	if (rDot-rDotPre<0 && rDot>1e-4)
	{
		scalar a1=(Tv-Tinf_)*exp(dTdr*delta/(Tv-Tinf_));
		scalar a2=((Tv-Tinf_)+alpha*ODESolverStep()*dTdr*dTdr/(Tv-Tinf_)-2*alpha*ODESolverStep()*dTdr/(r+delta));
		return (Tv-Tinf_)/dTdr*log(a1/a2);
	}
	
	
	else 
	return (Tv-Tinf_)/dTdr*log(0.103);
}

template<class ParcelType>
inline void Foam::ODERP<ParcelType>::getTv
(
	scalar delta,
	scalar& Tv,
	scalar& Tref,
	scalar P,
	scalar r,
	scalar rDot
)
{
	
	scalar dTdr=gradT(P,Tv,rDot);
	scalar alpha=fuel_->K(P,Tv)/fuel_->rho(P,Tv)/fuel_->Cp(P,Tv);
	scalar d2Tdr2_ref,dTdr_ref;
	scalar Rref=1e-8;
	if(delta!=0)
	{
	label n=3;
	scalarSquareMatrix a(n);
	a[0][0]=delta*delta;
	a[0][1]=pow(delta,3);
	a[0][2]=pow(delta,4);
	a[1][0]=2*delta;
	a[1][1]=3*delta*delta;
	a[1][2]=4*pow(delta,3);
	a[2][0]=2;
	a[2][1]=6*delta;
	a[2][2]=12*delta*delta;

	scalarField b(n),c(n);
	b[0]=Tinf_-Tv-dTdr*delta;
	b[1]=-dTdr;
	b[2]=0;

	//Gaussian elimination
	for (int i=0;i<2;i++)
	{
		for( int j=i+1;j<3;j++)
		{
			scalar temp=a[j][i]/a[i][i];
			for(int k=0;k<3;k++)
			{
				a[j][k]-=temp*a[i][k];
			}
			b[j]-=temp*b[i];
		}
	}
	for (int i=2;i>=0;i--)
	{
		c[i]=b[i]/a[i][i];
		for (int k=2;k>i;k--)
		c[i]-=a[i][k]*c[k]/a[i][i];
	}
	
	//find the 2rd spatial derivatives at a distance form the bubble interface
		
	d2Tdr2_ref=c[2]*12*Rref*Rref+6*c[1]*Rref+2*c[0];
	dTdr_ref=4*c[2]*pow(Rref,3)+3*c[1]*Rref*Rref+2*c[0]*Rref+dTdr;
	}
	else
	{
	d2Tdr2_ref=0;
	dTdr_ref=dTdr;
	}
	scalar dTdt_ref=alpha*(d2Tdr2_ref+2/r*dTdr_ref);
	Tref+=dTdt_ref*ODESolverStep();
	Tv=Tref-dTdr*Rref;
	
}



template<class ParcelType>
inline void Foam::ODERP<ParcelType>::operator=(const ODERP<ParcelType>& rhs)
{
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "ODERP::operator=(const ODERP&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
    ODEProperties_=rhs.ODEProperties_;
    R0_=rhs.R0_;
    gamma_=rhs.gamma_;
    p0_=rhs.p0_;
    ODESolverType_=rhs.ODESolverType_;
    ODESolverStep_=rhs.ODESolverStep_;
    ODESolverSubStep_=rhs.ODESubSolverStep_;
    ODESolverTol_=rhs.ODESolverTol_;
    pL_ = rhs.pL_;
    Tv_=  rhs. Tv_;
    Tinf_ = rhs.Tinf_;
    fuel_ = rhs.fuel_;
}
