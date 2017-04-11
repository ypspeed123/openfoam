
label cellB = mesh.findCell(p.position());
vector posB= p.position();
const scalar& pLD =  td.pInterp().interpolate(posB, cellB);
const scalar& rho = td.rhoInterp().interpolate(posB, cellB);
scalar pL=pLD*rho;
scalar Tv= this->Tv();
scalar rDotPre= this->rDotPre();
scalar delta=this->delta();
scalar Tref= this->Tref();

dictionary dict;
dict.add("solver",  oderp_.ODESolverType());

autoPtr<ODESolver> solveRP = ODESolver::New(oderp_,dict);
Info<<"Now caculating time t="<<tNow<<endl;
Info<<"lagrangian tracking time dt="<<dt<<endl;
scalar xStart=tNow;
scalarField yStart(oderp_.nEqns());
yStart[0] = this->d()/2;
yStart[1] = this->rDot();
if (yStart[1]==0)
{
	yStart[1]=1e-6;
	yStart[0]=2*oderp_.fuel().sigma(pL,oderp_.Tinf())/(oderp_.fuel().pv(pL,oderp_.Tinf())-pL);
	Tv=oderp_.Tinf();
	Tref=oderp_.Tinf();
//	Info<<yStart[0]<<" "<<yStart[1]<<endl;
}

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

this->update(2*y[0],y[1],rDotPre,delta,Tv,Tref);

tNow +=dt;
