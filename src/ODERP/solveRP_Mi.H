
label cellB = mesh.findCell(p.position());
vector posB= p.position();
const scalar& pLD =  td.pInterp().interpolate(posB, cellB);
const scalar& rho = td.rhoInterp().interpolate(posB, cellB);
scalar pL=pLD*rho;


scalar R=this->d()/2;
scalar Rdot=this->rDot();
scalar Tv= this->Tv();
scalar rDotPre= this->rDotPre();
scalar delta=this->delta();
scalar Tref= this->Tref();

if (Rdot==0)
{
	R=2*oderp_.fuel().sigma(pL,oderp_.Tinf())/(oderp_.fuel().pv(pL,oderp_.Tinf())-pL);
}

scalar age=this->age();
scalar tStep=oderp_.ODESolverStep();
scalar Tsat=oderp_.Tsat(pL);
scalar hfg=oderp_.fuel().hl(pL,oderp_.Tinf());
scalar rhov=oderp_.rhov(oderp_.Tinf());
scalar rhol=oderp_.fuel().rho(pL,oderp_.Tinf());
scalar cpl=oderp_.fuel().Cp(pL,oderp_.Tinf());
scalar alpha=oderp_.fuel().K(pL,oderp_.Tinf())/oderp_.fuel().rho(pL,oderp_.Tinf())/oderp_.fuel().Cp(pL,oderp_.Tinf());
scalar A=max(1e-12,(2*hfg*rhov*(oderp_.Tinf()-Tsat))/(3*rhol*Tsat));
scalar B=max(1e-12,sqrt(12/3.14159*alpha)*(oderp_.Tinf()-Tsat)*cpl*rhol/hfg/rhov);
scalar xStart=tNow;

while ((xStart+oderp_.ODESolverStep())<=(tNow+dt))
{
	scalar k1=oderp_.deri(age,A,B);
	scalar k2=oderp_.deri(age+tStep/2.0,A,B);
	scalar k3=oderp_.deri(age+tStep,A,B);
	R=R+tStep/6*(k1+4*k2+k3);
	Rdot=k3;
	age=age+tStep;
	xStart+=tStep;
//	if (oderp_.writeControl(xStart))
//	oderp_.output(xStart,R,Rdot,Tv,delta,pL,this->origId());

}

/*
scalarField y(yStart);

while ((xStart+oderp_.ODESolverStep())<=(tNow+dt))
{
//	Info<<Tv<<" "<<rDotPre<<" "<<delta<<" "<<Tref<<" "<<pL<<endl;
	oderp_.setEnvironment( pL, Tv);
	scalar x = xStart;
    scalar xEnd = x + oderp_.ODESolverStep();
    if (oderp_.writeControl(x))
	oderp_.output(x,yStart[0],yStart[1],Tv,delta,pL,this->origId());

	//solve the RP equation and find the new vapor temp
    scalar hEst = oderp_.ODESolverSubStep();
    solveRP -> relTol()= oderp_.ODESolverTol();
    solveRP -> solve(x, xEnd, y,hEst);
//    Info<<"starting calculating delta with Tv="<<Tv<<endl;
    scalar deltaNew=oderp_.delta(pL,Tv,yStart[0],yStart[1],rDotPre,delta,oderp_.Tinf());

//    Info<<"R: "<<y[0]<<" Rdot: "<<y[1]<<" Tv: "<<Tv<<" delta: "<<deltaNew<<endl;
//    Info<<"Tref: "<<Tref<<" rDOtPre "<<rDotPre<<endl;
    //update for the next step
    oderp_.getTv(deltaNew,Tv,Tref,pL,yStart[0],yStart[1]);
    rDotPre=yStart[1];  
    yStart[0] = y[0];
    yStart[1] = y[1];
//    Info<<y[0]<<" "<<y[1]<<" "<<Tv<<endl;
    delta=deltaNew;
    xStart = xEnd;
}
*/
this->update(2*R,Rdot,rDotPre,delta,Tv,Tref);

tNow +=dt;
