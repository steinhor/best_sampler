#include "msu_sampler/hyper.h"

using namespace msu_sampler;

Chyper::Chyper(){
	muB=muI=muS=0.0;
	Rvisc_calculated=false;
}

void Chyper::Print(){
	printf("----- hyper info -------\n");
	printf("T=%g, sigma=%g, rhoB=%g, rhoI=%g\n",T,sigma,rhoB,rhoI);
	printf("muB=%g, muI=%g, muS=%g\n",muB,muI,muS);
	printf("epsilon=%g, P=%g, s=%g\n",epsilon,P,GetEntropyDensity());
	printf("dOmega=(%g,%g,%g,%g), udotdOmega=%g\n",dOmega[0],dOmega[1],dOmega[2],dOmega[3],udotdOmega);
	printf("u=(%g,%g,%g,%g).    tau=%g,eta=%g\n",u[0],u[1],u[2],u[3],tau,eta);
	for(int alpha=0;alpha<4;alpha++){
			printf("%10.3e %10.3e %10.3e %10.3e\n",
			pitilde[alpha][0],pitilde[alpha][1],pitilde[alpha][2],pitilde[alpha][3]);
	}
	printf("-----------------------\n");
}

double Chyper::GetEntropyDensity(){
	double s=((epsilon+P)/T)-muB*rhoB-muI*rhoI-muS*rhoS;
	return s;	
}
