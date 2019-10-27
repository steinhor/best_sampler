#include "pratt_sampler/hyper.h"

Chyper::Chyper(){
	muB=muI=muS=0.0;
}

void Chyper::Print(){
	printf("----- hyper info -------\n");
	printf("T=%g, sigma=%g, rhoB=%g, rhoI=%g\n",T,sigma,rhoB,rhoI);
	printf("dOmega=(%g,%g,%g,%g)\n",dOmega[0],dOmega[1],dOmega[2],dOmega[3]);
	printf("u=(%g,%g,%g,%g).    tau=%g,eta=%g\n",u[0],u[1],u[2],u[3],tau,eta);
	for(int alpha=0;alpha<4;alpha++){
			printf("%10.3e %10.3e %10.3e %10.3e\n",
			pitilde[alpha][0],pitilde[alpha][1],pitilde[alpha][2],pitilde[alpha][3]);
	}
	double T,sigma,rhoB,rhoI;
	FourVector dOmega; // hyper volume
	FourVector u;  // collective velocity
	FourVector r; // position
	double tau,eta; // Bjorken tau and spatial rapidity
	double pitilde[4][4]; // shear tensor
}
