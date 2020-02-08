#ifndef __INCLUDE_MISC_CC
#define __INCLUDE_MISC_CC
#include <iostream>
#include <cmath>
#include <complex>
#include "pratt_sampler/misc.h"
using namespace pratt_sampler;

void Misc::Boost(FourVector &u,FourVector &p,FourVector &pprime){
	int alpha;
	if (std::isnan(pprime[1]*pprime[2]*pprime[3])) printf("3 p= %lf %lf %lf %lf\n",pprime[0],pprime[1],pprime[2],pprime[3]);
	double pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
	double A=-(pdotu+p[0])/(1.0+u[0]);
	//printf("A=%lf\n",A);
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]=p[alpha];
		//printf("pprime[alpha]=%lf\tp[alpha]=%lf\n",pprime[alpha],p[alpha]);
	pprime[0]+=A;
	for(alpha=0;alpha<4;alpha++)
		pprime[alpha]+=(A+2*p[0])*u[alpha];
}

void Misc::BoostToCM(FourVector &u,FourVector &p,FourVector &ptilde){
	int alpha;
	double pdotu=p[0]*u[0]-p[1]*u[1]-p[2]*u[2]-p[3]*u[3];
	double A=-(pdotu+p[0])/(1.0+u[0]);
	for(alpha=0;alpha<4;alpha++)
		ptilde[alpha]=p[alpha];
	ptilde[0]+=A+2*pdotu;
	for(alpha=0;alpha<4;alpha++)
		ptilde[alpha]+=A*u[alpha];
}

void Misc::Pause(){
	double pausedummy;
	printf("PAUSED, enter anything to continue: ");
	scanf("%lf",&pausedummy);
}

void Misc::Pause(int seconds){
	//const double CLOCKS_PER_SEC=1.0;
	clock_t endwait;
	endwait = clock () + seconds * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}

#endif
